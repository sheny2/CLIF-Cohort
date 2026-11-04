# Shared ICD-10-CM helpers for the Hopkins severe-CAP cohort.
#
# The primary definition agreed for this target-trial emulation is an explicit
# pneumonia diagnosis from J12.*, J13-J16.*, or J18.* that is documented as
# present on admission. Diagnosis position is not restricted. Codes for
# influenza pneumonia (J09-J11), pneumonia in diseases classified elsewhere
# (J17), aspiration (J69), and ventilator-associated pneumonia (J95.851) do not
# qualify. Etiology or sepsis codes (for example U07.1, B95.3, B96.1, B97.29,
# and A40.3) are not sufficient without a qualifying pneumonia diagnosis.

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

normalize_poa_present <- function(x) {
  value <- tolower(trimws(as.character(x)))
  as.integer(value %in% c("1", "true", "t", "yes", "y"))
}

classify_pneumonia_icd <- function(code, format) {
  code_clean <- normalize_icd_code(code)
  format_clean <- normalize_icd_format(format)
  icd10 <- format_clean == "ICD10CM"

  qualifying_pneumonia <- icd10 & grepl(
    "^J(12|13|14|15|16|18)[A-Z0-9]*$", code_clean
  )
  influenza_pneumonia <- icd10 & grepl("^J(09|10|11)[A-Z0-9]*$", code_clean)
  excluded_pneumonia <- icd10 & (
    grepl("^J17[A-Z0-9]*$", code_clean) |
      grepl("^J69[A-Z0-9]*$", code_clean) |
      code_clean == "J95851"
  )

  qualifying_pneumonia[is.na(qualifying_pneumonia)] <- FALSE
  influenza_pneumonia[is.na(influenza_pneumonia)] <- FALSE
  excluded_pneumonia[is.na(excluded_pneumonia)] <- FALSE

  data.table::data.table(
    diagnosis_code_clean = code_clean,
    diagnosis_code_format_clean = format_clean,
    is_qualifying_pneumonia_code = as.integer(qualifying_pneumonia),
    is_influenza_pneumonia_code = as.integer(influenza_pneumonia),
    is_explicitly_excluded_pneumonia_code = as.integer(excluded_pneumonia)
  )
}

derive_pneumonia_encounter_flags <- function(diagnosis) {
  required <- c(
    "hospitalization_id", "diagnosis_code", "diagnosis_code_format",
    "diagnosis_primary", "poa_present"
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
    poa_present_binary = normalize_poa_present(poa_present),
    is_qualifying_pneumonia_code =
      classification$is_qualifying_pneumonia_code,
    is_influenza_pneumonia_code =
      classification$is_influenza_pneumonia_code,
    is_explicitly_excluded_pneumonia_code =
      classification$is_explicitly_excluded_pneumonia_code
  )]
  x[, qualifies_poa_pneumonia := as.integer(
    is_qualifying_pneumonia_code == 1L & poa_present_binary == 1L
  )]

  encounter_flags <- x[
    !is.na(hospitalization_id) & hospitalization_id != "",
    .(
      poa_pneumonia_case = as.integer(any(qualifies_poa_pneumonia == 1L)),
      poa_pneumonia_primary = as.integer(any(
        qualifies_poa_pneumonia == 1L & diagnosis_primary_binary == 1L
      )),
      poa_pneumonia_secondary = as.integer(any(
        qualifies_poa_pneumonia == 1L & diagnosis_primary_binary == 0L
      )),
      influenza_pneumonia_present = as.integer(any(
        is_influenza_pneumonia_code == 1L & poa_present_binary == 1L
      )),
      explicitly_excluded_pneumonia_present = as.integer(any(
        is_explicitly_excluded_pneumonia_code == 1L & poa_present_binary == 1L
      ))
    ),
    by = hospitalization_id
  ]

  list(diagnoses = x, encounters = encounter_flags)
}
