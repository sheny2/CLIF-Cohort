# Shared ICD-9-CM/ICD-10-CM helpers for the adult pneumonia cohort.
#
# The expanded encounter definition follows Lindenauer et al. (JAMA 2012):
#   1. a principal diagnosis of pneumonia; or
#   2. a principal diagnosis of sepsis or respiratory failure together with a
#      secondary diagnosis of pneumonia in the same component hospitalization.
#
# The original publication used ICD-9-CM. The ICD-10-CM pneumonia families here
# are the corresponding bacterial, other-organism, and organism-unspecified
# categories. Influenza pneumonia, viral pneumonia, aspiration pneumonia, and
# pneumonia in diseases classified elsewhere are intentionally excluded.

normalize_icd_format <- function(x) {
  value <- toupper(gsub("[^A-Z0-9]", "", as.character(x)))
  data.table::fcase(
    value %in% c("ICD9", "ICD9CM", "9", "09"), "ICD9CM",
    value %in% c("ICD10", "ICD10CM", "10"), "ICD10CM",
    grepl("ICD9", value), "ICD9CM",
    grepl("ICD10", value), "ICD10CM",
    default = NA_character_
  )
}

normalize_icd_code <- function(x) {
  toupper(gsub("[^A-Z0-9]", "", as.character(x)))
}

normalize_diagnosis_primary <- function(x) {
  value <- tolower(trimws(as.character(x)))
  as.integer(value %in% c("1", "true", "t", "yes", "y", "primary", "principal"))
}

classify_pneumonia_icd <- function(code, format) {
  code_clean <- normalize_icd_code(code)
  format_clean <- normalize_icd_format(format)

  icd9 <- format_clean == "ICD9CM"
  icd10 <- format_clean == "ICD10CM"

  pneumonia <-
    (icd9 & grepl("^(481|482|483|485|486)[0-9]*$", code_clean)) |
    (icd10 & grepl("^J(13|14|15|16|18)[A-Z0-9]*$", code_clean))

  sepsis <-
    (icd9 & grepl("^(038[0-9]*|99591|99592|78552)$", code_clean)) |
    (icd10 & grepl("^(A4[01][A-Z0-9]*|R652[01])$", code_clean))

  respiratory_failure <-
    (icd9 & grepl("^(51881|51882|51884|7991)$", code_clean)) |
    (icd10 & grepl("^(J80|J96(0|2|9)[A-Z0-9]*|R092)$", code_clean))

  pneumonia[is.na(pneumonia)] <- FALSE
  sepsis[is.na(sepsis)] <- FALSE
  respiratory_failure[is.na(respiratory_failure)] <- FALSE

  data.table::data.table(
    diagnosis_code_clean = code_clean,
    diagnosis_code_format_clean = format_clean,
    is_pneumonia_code = as.integer(pneumonia),
    is_sepsis_code = as.integer(sepsis),
    is_respiratory_failure_code = as.integer(respiratory_failure)
  )
}

derive_pneumonia_encounter_flags <- function(diagnosis) {
  required <- c(
    "hospitalization_id", "diagnosis_code", "diagnosis_code_format",
    "diagnosis_primary"
  )
  missing <- setdiff(required, names(diagnosis))
  if (length(missing)) {
    stop(
      "Diagnosis data are missing required column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  x <- data.table::as.data.table(data.table::copy(diagnosis))
  x[, hospitalization_id := as.character(hospitalization_id)]
  classification <- classify_pneumonia_icd(
    x$diagnosis_code,
    x$diagnosis_code_format
  )
  x[, `:=`(
    diagnosis_code_original = as.character(diagnosis_code),
    diagnosis_code_clean = classification$diagnosis_code_clean,
    diagnosis_code_format_clean = classification$diagnosis_code_format_clean,
    diagnosis_primary_binary = normalize_diagnosis_primary(diagnosis_primary),
    is_pneumonia_code = classification$is_pneumonia_code,
    is_sepsis_code = classification$is_sepsis_code,
    is_respiratory_failure_code = classification$is_respiratory_failure_code
  )]

  encounter_flags <- x[
    !is.na(hospitalization_id) & hospitalization_id != "",
    .(
      principal_pneumonia = as.integer(any(
        is_pneumonia_code == 1L & diagnosis_primary_binary == 1L
      )),
      any_pneumonia = as.integer(any(is_pneumonia_code == 1L)),
      secondary_pneumonia = as.integer(any(
        is_pneumonia_code == 1L & diagnosis_primary_binary == 0L
      )),
      principal_sepsis = as.integer(any(
        is_sepsis_code == 1L & diagnosis_primary_binary == 1L
      )),
      principal_respiratory_failure = as.integer(any(
        is_respiratory_failure_code == 1L & diagnosis_primary_binary == 1L
      ))
    ),
    by = hospitalization_id
  ]
  encounter_flags[, expanded_pneumonia_case := as.integer(
    principal_pneumonia == 1L |
      (
        secondary_pneumonia == 1L &
          (principal_sepsis == 1L | principal_respiratory_failure == 1L)
      )
  )]
  encounter_flags[, pneumonia_case_pathway := data.table::fcase(
    principal_pneumonia == 1L, "principal pneumonia",
    principal_sepsis == 1L & secondary_pneumonia == 1L,
      "principal sepsis + secondary pneumonia",
    principal_respiratory_failure == 1L & secondary_pneumonia == 1L,
      "principal respiratory failure + secondary pneumonia",
    any_pneumonia == 1L, "other/secondary pneumonia only",
    default = "no qualifying pneumonia code"
  )]

  list(diagnoses = x, encounters = encounter_flags)
}
