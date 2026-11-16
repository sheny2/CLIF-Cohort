# Shared preparation for this standalone dated analysis only.
prepare_analysis_inputs <- function(data_dir, cohort_dir, treatment_dir, sensitive_dir,
                               shareable_dir, grace_hours, cutoff_text, max_patients) {
  step_hours <- 6L; horizon_days <- 28L
  tau <- as.integer(grace_hours / step_hours)
  M <- as.integer(horizon_days * 24 / step_hours)
  sample_seed <- 20260814L
  data_cutoff_dttm <- as.POSIXct(paste0(as.Date(cutoff_text), " 23:59:59"), tz="UTC")
  if (is.na(data_cutoff_dttm)) stop("DATA_CUTOFF_DATE must be YYYY-MM-DD")
  cohort_path <- file.path(cohort_dir, "adult_pneumonia_cohort.parquet")
  crosswalk_path <- file.path(cohort_dir, "adult_pneumonia_hospitalization_crosswalk.parquet")
  steroid_path <- file.path(treatment_dir, "qualifying_steroid_administrations.parquet")
  for (f in c(cohort_path, crosswalk_path, steroid_path))
    if (!file.exists(f)) stop("Missing upstream output: ", f, "; run scripts 01--04 first.")
  write_shareable_csv <- function(x, filename) {
    forbidden <- intersect(c("patient_id", "hospitalization_id", "hospital_block_id"), names(x))
    if (length(forbidden)) stop("Identifier in aggregate preparation output")
    fwrite(as.data.table(x), file.path(shareable_dir, filename), na="")
  }
  write_sensitive_csv <- function(x, filename) fwrite(as.data.table(x), file.path(sensitive_dir, filename), na="")
required_file <- function(filename) {
  path <- file.path(data_dir, filename)
  if (!file.exists(path)) stop("Required CLIF file not found: ", path)
  path
}

read_selected <- function(filename, columns) {
  path <- required_file(filename)
  dataset <- open_dataset(path)
  schema_names <- dataset$schema$names
  missing <- setdiff(columns, schema_names)
  if (length(missing)) {
    stop(filename, " is missing: ", paste(missing, collapse = ", "))
  }
  query <- dataset |> select(all_of(columns))
  if ("hospitalization_id" %in% columns && exists("crosswalk", inherits = TRUE)) {
    ids <- unique(get("crosswalk", inherits = TRUE)$hospitalization_id)
    query <- query |> filter(.data$hospitalization_id %in% ids)
  }
  if ("patient_id" %in% columns && !"hospitalization_id" %in% columns &&
      exists("cohort", inherits = TRUE)) {
    ids <- unique(get("cohort", inherits = TRUE)$patient_id)
    query <- query |> filter(.data$patient_id %in% ids)
  }
  query |> collect() |> as.data.table()
}

clean_code <- function(x) toupper(gsub("[^A-Z0-9]", "", as.character(x)))

to_interval <- function(hours) {
  result <- rep.int(M + 1L, length(hours))
  observed <- is.finite(hours)
  result[observed] <- pmin(
    M + 1L,
    pmax(1L, as.integer(ceiling(hours[observed] / step_hours)))
  )
  result
}

# ---- Hospitalization, demographics, outcome, comorbidities ------------------

message("Building baseline hospital-block cohort and 28-day in-hospital death/hospice outcome...")
cohort <- as.data.table(read_parquet(cohort_path, as_data_frame = TRUE))
required_cohort_columns <- c(
  "hospital_block_id", "patient_id", "admission_dttm", "discharge_dttm",
  "age_at_admission", "respiratory_failure_start",
  "respiratory_support_category", "source_site_id"
)
missing_cohort_columns <- setdiff(required_cohort_columns, names(cohort))
if (length(missing_cohort_columns)) {
  stop(
    "The cohort is not the revised block-level output; missing: ",
    paste(missing_cohort_columns, collapse = ", "),
    "\nRerun 01_build_adult_pneumonia_cohort.R."
  )
}
cohort[, `:=`(
  hospital_block_id = as.character(hospital_block_id),
  patient_id = as.character(patient_id),
  admission_dttm = as.POSIXct(admission_dttm, tz = "UTC"),
  discharge_dttm = as.POSIXct(discharge_dttm, tz = "UTC"),
  respiratory_failure_start = as.POSIXct(respiratory_failure_start, tz = "UTC"),
  index_dttm = as.POSIXct(respiratory_failure_start, tz = "UTC")
)]
if (anyDuplicated(cohort$patient_id)) {
  stop("The revised cohort must contain exactly one block per patient.")
}

crosswalk <- as.data.table(read_parquet(crosswalk_path, as_data_frame = TRUE))
required_crosswalk_columns <- c(
  "patient_id", "hospital_block_id", "hospitalization_id"
)
missing_crosswalk_columns <- setdiff(required_crosswalk_columns, names(crosswalk))
if (length(missing_crosswalk_columns)) {
  stop(
    "Hospital-block crosswalk is missing: ",
    paste(missing_crosswalk_columns, collapse = ", ")
  )
}
crosswalk[, `:=`(
  patient_id = as.character(patient_id),
  hospital_block_id = as.character(hospital_block_id),
  hospitalization_id = as.character(hospitalization_id)
)]
crosswalk <- unique(crosswalk)
if (!all(crosswalk$hospital_block_id %in% cohort$hospital_block_id)) {
  stop("Crosswalk contains blocks absent from the selected cohort.")
}
if (crosswalk[, uniqueN(hospital_block_id), by = hospitalization_id][V1 > 1L, .N]) {
  stop("A hospitalization_id maps to more than one selected hospital block.")
}

# Restrict the original selected cohort using BASELINE information only.
sf_selection <- select_baseline_sf(data_dir,cohort,crosswalk,shareable_dir,sensitive_dir)
cohort <- sf_selection$cohort;crosswalk <- sf_selection$crosswalk
source_cohort_patients <- nrow(sf_selection$audit)
encoded_sf <- encode_baseline_sf(cohort)
cohort <- encoded_sf$cohort
write_shareable_csv(encoded_sf$processing,"baseline_sf_processing.csv")
rm(sf_selection,encoded_sf)

# Sampling occurs before auxiliary CLIF tables are loaded or the longitudinal
# grid is built. Because 01 already selected one block per patient, the sampled
# patient, block, and independent analysis-unit counts are identical.
full_patients <- uniqueN(cohort$patient_id)
full_blocks <- uniqueN(cohort$hospital_block_id)
full_component_hospitalizations <- uniqueN(crosswalk$hospitalization_id)
if (max_patients > 0L && full_patients > max_patients) {
  set.seed(sample_seed)
  sampled_patient_ids <- sample(
    sort(unique(cohort$patient_id)),
    size = max_patients,
    replace = FALSE
  )
  cohort <- cohort[patient_id %in% sampled_patient_ids]
  crosswalk <- crosswalk[patient_id %in% sampled_patient_ids]
  message(
    "Random sample retained ", uniqueN(cohort$patient_id),
    " patients, ", uniqueN(cohort$hospital_block_id), " hospital blocks, and ",
    uniqueN(crosswalk$hospitalization_id), " component hospitalization IDs."
  )
} else {
  message(
    "Full cohort retained because it contains only ", full_patients,
    " patients."
  )
}

sampling_summary <- data.table(
  measure = c(
    "sample_seed", "maximum_patients_requested", "full_cohort_patients",
    "full_cohort_hospital_blocks", "full_component_hospitalizations",
    "sampled_patients", "sampled_hospital_blocks",
    "sampled_component_hospitalizations"
  ),
  value = c(
    sample_seed, max_patients, full_patients, full_blocks,
    full_component_hospitalizations, uniqueN(cohort$patient_id),
    uniqueN(cohort$hospital_block_id), uniqueN(crosswalk$hospitalization_id)
  )
)
write_shareable_csv(sampling_summary, "analysis_sampling_summary.csv")

map_to_selected_blocks <- function(x) {
  x[, hospitalization_id := as.character(hospitalization_id)]
  x <- x[hospitalization_id %in% crosswalk$hospitalization_id]
  merge(
    x,
    crosswalk[, .(hospitalization_id, hospital_block_id)],
    by = "hospitalization_id",
    all = FALSE,
    sort = FALSE
  )
}

patient <- read_selected(
  "clif_patient.parquet",
  c(
    "patient_id", "sex_category", "race_category", "ethnicity_category",
    "death_dttm"
  )
)
patient[, patient_id := as.character(patient_id)]
patient <- unique(patient, by = "patient_id")

new_patient_columns <- setdiff(names(patient), names(cohort))
cohort <- merge(
  cohort,
  patient[, c("patient_id", new_patient_columns), with = FALSE],
  by = "patient_id",
  all.x = TRUE,
  sort = FALSE
)

diagnosis <- read_selected(
  "clif_hospital_diagnosis.parquet",
  c("hospitalization_id", "diagnosis_code", "diagnosis_code_format",
    "poa_present")
)
diagnosis[, `:=`(
  hospitalization_id = as.character(hospitalization_id),
  code = clean_code(diagnosis_code),
  format = toupper(gsub("[^A-Z0-9]", "", as.character(diagnosis_code_format)))
)]
diagnosis <- diagnosis[
  tolower(trimws(as.character(poa_present))) %chin%
    c("1", "true", "t", "yes", "y")
]
diagnosis <- map_to_selected_blocks(diagnosis)

# Claims-based proxies for the chronic conditions listed in the screenshot.
diagnosis[, chronic_lung := as.integer(
  (format == "ICD10CM" & grepl("^J(4[0-7]|6[0-7])", code)) |
    (format == "ICD9CM" & grepl("^(49[0-6]|50[0-5])", code))
)]
diagnosis[, uncontrolled_diabetes := as.integer(
  (format == "ICD10CM" & grepl("^E(08|09|10|11|13)(1|65)", code)) |
    (format == "ICD9CM" & grepl("^250[0-9][23]$", code))
)]
diagnosis[, immunosuppressed := as.integer(
  (format == "ICD10CM" & grepl("^(D8[0-4]|B20|Z94|Z796)", code)) |
    (format == "ICD9CM" & grepl("^(042|279|V42|V5865)", code))
)]
diagnosis[, metastatic_solid_tumor := as.integer(
  (format == "ICD10CM" & grepl("^C(7[7-9]|80)", code)) |
    (format == "ICD9CM" & grepl("^19[6-9]", code))
)]
diagnosis[, hematologic_malignancy := as.integer(
  (format == "ICD10CM" & grepl("^(C8[1-9]|C9[0-6]|D4[5-7])", code)) |
    (format == "ICD9CM" & grepl("^(20[0-8]|2387)", code))
)]
diagnosis[, heart_failure := as.integer(
  (format == "ICD10CM" & grepl("^I50", code)) |
    (format == "ICD9CM" & grepl("^428", code))
)]
diagnosis[, chronic_kidney_disease := as.integer(
  (format == "ICD10CM" & grepl("^N18", code)) |
    (format == "ICD9CM" & grepl("^585", code))
)]

comorbidities <- diagnosis[
  ,
  .(
    chronic_lung = max(chronic_lung),
    uncontrolled_diabetes = max(uncontrolled_diabetes),
    immunosuppressed = max(immunosuppressed),
    metastatic_solid_tumor = max(metastatic_solid_tumor),
    hematologic_malignancy = max(hematologic_malignancy),
    heart_failure = max(heart_failure),
    chronic_kidney_disease = max(chronic_kidney_disease)
  ),
  by = hospital_block_id
]
cohort <- merge(
  cohort,
  comorbidities,
  by = "hospital_block_id",
  all.x = TRUE,
  sort = FALSE
)
comorbidity_names <- c(
  "chronic_lung", "uncontrolled_diabetes", "immunosuppressed",
  "metastatic_solid_tumor", "hematologic_malignancy", "heart_failure",
  "chronic_kidney_disease"
)
for (variable in comorbidity_names) setnafill(cohort, "const", 0, cols = variable)

# The primary endpoint excludes deaths after the selected hospital block.
cohort <- derive_primary_outcome(cohort,"index_dttm")
write_outcome_definition(cohort,shareable_dir,data_cutoff_dttm)
cohort[, censor_hours := as.numeric(difftime(
  data_cutoff_dttm, index_dttm, units = "hours"
))]
if (any(cohort$censor_hours < 0, na.rm = TRUE)) {
  stop("At least one cohort time zero occurs after the data cutoff.")
}
cohort[, delta := as.integer(
  is.finite(event_hours) & event_hours <= horizon_days * 24 &
    event_hours <= censor_hours
)]
cohort[, T_event := fifelse(delta == 1L, to_interval(event_hours), M + 1L)]
cohort[, T_censor := fifelse(
  censor_hours < horizon_days * 24,
  pmin(M, pmax(1L, as.integer(floor(censor_hours / step_hours)) + 1L)),
  M + 1L
)]
cohort[, T_obs := pmin(T_event, T_censor, M)]

cohort[, `:=`(
  age_z = (as.numeric(age_at_admission) - 65) / 10,
  female = fifelse(tolower(sex_category) == "female", 1L, 0L, 0L),
  sex_missing = as.integer(
    is.na(sex_category) | !tolower(sex_category) %in% c("female", "male")
  ),
  baseline_support_imv = as.integer(respiratory_support_category == "imv"),
  baseline_support_nippv = as.integer(respiratory_support_category %chin% c("nippv", "cpap")),
  admission_year_z = (as.integer(format(index_dttm, "%Y")) - 2021) / 5,
  admission_season_sin = sin(2 * pi * as.integer(format(index_dttm, "%m")) / 12),
  admission_season_cos = cos(2 * pi * as.integer(format(index_dttm, "%m")) / 12)
)]

# Initial location proxies. Hospital identity is the federating/site-stratifying
# variable and therefore is not duplicated as a patient-level model covariate.
adt <- read_selected(
  "clif_adt.parquet",
  c("hospitalization_id", "in_dttm", "out_dttm", "location_category")
)
adt <- map_to_selected_blocks(adt)
adt <- merge(
  adt,
  cohort[, .(hospital_block_id, index_dttm)],
  by = "hospital_block_id",
  all.x = TRUE,
  sort = FALSE
)
adt[, `:=`(
  in_dttm = as.POSIXct(in_dttm, tz = "UTC"),
  out_dttm = as.POSIXct(out_dttm, tz = "UTC")
)]
adt[, `:=`(
  contains_index = as.integer(
    !is.na(in_dttm) & in_dttm <= index_dttm &
      (is.na(out_dttm) | out_dttm >= index_dttm)
  ),
  prior_to_index = as.integer(!is.na(in_dttm) & in_dttm <= index_dttm),
  distance = abs(as.numeric(difftime(in_dttm, index_dttm, units = "hours")))
)]
setorder(adt, hospital_block_id, -contains_index, -prior_to_index, distance)
initial_location <- adt[, .SD[1L], by = hospital_block_id]
initial_location[, location_text := tolower(as.character(location_category))]
initial_location[, `:=`(
  initial_location_icu = as.integer(grepl("icu|intensive", location_text)),
  initial_location_ed = as.integer(grepl("emerg|ed", location_text))
)]
cohort <- merge(
  cohort,
  initial_location[, .(
    hospital_block_id,
    initial_location_icu,
    initial_location_ed
  )],
  by = "hospital_block_id",
  all.x = TRUE,
  sort = FALSE
)
setnafill(
  cohort,
  "const",
  0,
  cols = c("initial_location_icu", "initial_location_ed")
)

# ---- Qualifying glucocorticoid initiation -----------------------------------

message("Reading qualifying systemic glucocorticoid administrations...")
steroid <- as.data.table(read_parquet(steroid_path, as_data_frame = TRUE))
required_steroid <- c("hospital_block_id", "admin_dttm")
missing_steroid <- setdiff(required_steroid, names(steroid))
if (length(missing_steroid)) {
  stop("Qualifying steroid file is missing: ",
       paste(missing_steroid, collapse = ", "))
}
steroid[, `:=`(
  hospital_block_id = as.character(hospital_block_id),
  admin_dttm = as.POSIXct(admin_dttm, tz = "UTC")
)]
steroid <- unique(steroid[
  hospital_block_id %in% cohort$hospital_block_id,
  .(hospital_block_id, admin_dttm)
])
steroid <- merge(
  steroid,
  cohort[, .(hospital_block_id, index_dttm, discharge_dttm, event_dttm)],
  by = "hospital_block_id", all.x = TRUE, sort = FALSE
)
steroid[, treatment_hours := as.numeric(difftime(
  admin_dttm, index_dttm, units = "hours"
))]
steroid <- steroid[
  is.finite(treatment_hours) & treatment_hours >= 0 &
    treatment_hours <= horizon_days * 24 & admin_dttm <= discharge_dttm &
    admin_dttm <= data_cutoff_dttm &
    (is.na(event_dttm) | admin_dttm <= event_dttm)
]
first_steroid <- steroid[order(admin_dttm), .SD[1L], by = hospital_block_id]
cohort <- merge(
  cohort,
  first_steroid[, .(hospital_block_id, treatment_dttm = admin_dttm,
                     treatment_hours)],
  by = "hospital_block_id", all.x = TRUE, sort = FALSE
)
cohort[, S := to_interval(treatment_hours)]
cohort[S > T_event | S > T_censor, S := M + 1L]
cohort[, A_tau := as.integer(S <= tau)]

# Raw outcome counts. These are descriptive cohort quantities, not adjusted
# causal estimates. In particular, treatment-group comparisons can be affected
# by treatment-confounding and by primary events that occur before the grace period ends.
events_within_grace <- cohort[
  delta == 1L & is.finite(event_hours) & event_hours <= grace_hours,
  .N
]
events_after_grace <- cohort[
  delta == 1L & is.finite(event_hours) & event_hours > grace_hours,
  .N
]
events_28d <- cohort[delta == 1L, .N]
cohort_outcome_summary <- data.table(
  measure = c(
    "final sampled cohort",
    paste0("primary event within ", grace_hours, " hours of time zero"),
    paste0("primary event after ", grace_hours, " hours through day ", horizon_days),
    paste0("all primary events within ", horizon_days, " days"),
    paste0("no observed primary event within ", horizon_days, " days")
  ),
  n = c(
    nrow(cohort),
    events_within_grace,
    events_after_grace,
    events_28d,
    nrow(cohort) - events_28d
  )
)
cohort_outcome_summary[, percent_of_cohort := 100 * n / nrow(cohort)]
write_shareable_csv(cohort_outcome_summary, "cohort_outcome_summary.csv")

outcome_by_observed_treatment <- cohort[
  ,
  .(
    hospital_blocks = .N,
    events_within_grace = sum(
      delta == 1L & is.finite(event_hours) & event_hours <= grace_hours
    ),
    events_28d = sum(delta == 1L)
  ),
  by = .(glucocorticoid_within_grace = A_tau)
]
outcome_by_observed_treatment[, `:=`(
  events_within_grace_percent = 100 * events_within_grace / hospital_blocks,
  events_28d_percent = 100 * events_28d / hospital_blocks,
  interpretation = paste(
    "Descriptive only; observed treatment groups are not baseline randomized",
    "groups and an early primary event can prevent treatment initiation."
  )
)]
write_shareable_csv(
  outcome_by_observed_treatment,
  "cohort_outcome_by_observed_grace_period_treatment.csv"
)

# ---- Actual Hopkins hospitals as sites --------------------------------------

cohort[, site := as.character(source_site_id)]
if (anyNA(cohort$site) || any(!nzchar(cohort$site))) {
  stop("Every analysis unit must have a nonmissing source_site_id.")
}
if (uniqueN(cohort$site) < 2L) stop("At least two hospitals are required.")
patient_sites <- unique(cohort[, .(patient_id, site)])

# ---- Person-interval grid ----------------------------------------------------

grid <- cohort[
  ,
  .(m = seq_len(T_obs)),
  by = .(
    hospital_block_id,
    patient_id,
    site,
    index_dttm,
    admission_dttm,
    discharge_dttm,
    event_dttm,
    S,
    A_tau,
    T_event,
    T_censor,
    T_obs,
    delta,
    age_z,
    baseline_sf,
    baseline_sf_model,
    female,
    sex_missing,
    admission_year_z,
    admission_season_sin,
    admission_season_cos,
    baseline_support_imv,
    baseline_support_nippv,
    initial_location_icu,
    initial_location_ed,
    chronic_lung,
    uncontrolled_diabetes,
    immunosuppressed,
    metastatic_solid_tumor,
    hematologic_malignancy,
    heart_failure,
    chronic_kidney_disease
  )
]
grid[, decision_dttm := index_dttm + (m - 1L) * step_hours * 3600]

rolling_last <- function(measurement, time_col, value_col, output_col,
                         lookback_hours = 24) {
  if (!nrow(measurement)) {
    answer <- unique(grid[, .(hospital_block_id, m)])
    answer[, (output_col) := NA_real_]
    return(answer)
  }
  measurement <- copy(measurement)
  block <- match(measurement$hospital_block_id,cohort$hospital_block_id)
  measured <- as.numeric(as.POSIXct(measurement[[time_col]],tz="UTC"))
  end <- pmin(as.numeric(cohort$discharge_dttm[block]),as.numeric(data_cutoff_dttm))
  death <- as.numeric(cohort$event_dttm[block]); end[is.finite(death)] <- pmin(end[is.finite(death)],death[is.finite(death)])
  keep <- is.finite(measured)&measured>=as.numeric(cohort$admission_dttm[block])&measured<=end
  measurement <- measurement[which(keep)]
  measurement[, (value_col) := screen_measurements(get(value_col), output_col)]
  x <- measurement[
    !is.na(get(time_col)) & is.finite(as.numeric(get(value_col))),
    .(
      hospital_block_id = as.character(hospital_block_id),
      recorded_dttm = as.POSIXct(get(time_col), tz = "UTC"),
      measurement_seconds = as.numeric(get(time_col)),
      value = as.numeric(get(value_col))
    )
  ]
  if(output_col=="support_code") {
    setorder(x,hospital_block_id,recorded_dttm,value)
    x <- unique(x,by=c("hospital_block_id","recorded_dttm"),fromLast=TRUE)
  }
  setkey(x, hospital_block_id, recorded_dttm)
  answer <- x[
    grid[, .(hospital_block_id, m, decision_dttm)],
    on = .(hospital_block_id, recorded_dttm <= decision_dttm),
    mult = "last",
    .(
      hospital_block_id = i.hospital_block_id,
      m = i.m,
      decision_dttm = i.decision_dttm,
      measurement_seconds = x.measurement_seconds,
      value = x.value
    )
  ]
  answer[, age_hours := as.numeric(decision_dttm) / 3600 - measurement_seconds / 3600]
  answer[!is.finite(age_hours) | age_hours < 0 | age_hours > lookback_hours,
         value := NA_real_]
  answer[, c("decision_dttm", "measurement_seconds", "age_hours") := NULL]
  setnames(answer, "value", output_col)
  answer
}

# ---- Vitals ------------------------------------------------------------------

message("Aligning temperature, SpO2, and MAP...")
vitals <- read_selected(
  "clif_vitals.parquet",
  c("hospitalization_id", "recorded_dttm", "vital_category", "vital_value")
)
vitals <- map_to_selected_blocks(vitals)
vitals[, vital_category := gsub(
  "[^a-z0-9]+", "_", tolower(trimws(as.character(vital_category)))
)]
vitals <- vitals[
  vital_category %chin% c(
    "temperature", "temperature_c", "temp_c", "spo2", "map",
    "mean_arterial_pressure"
  )
]
# Raw range audit precedes cleaning; invalid readings never hide older valid ones.
raw_screen_audit <- list()
audit_raw <- function(x, value_col, variable) {
  y <- merge(x, cohort[,.(hospital_block_id,site,index_dttm)], by="hospital_block_id")
  y[, hours := as.numeric(difftime(recorded_dttm,index_dttm,units="hours"))]
  y <- y[is.finite(hours) & hours >= -24 & hours <= 672]
  raw_screen_audit[[variable]] <<- y[, {
    z <- as.numeric(get(value_col)); cleaned <- screen_measurements(z,variable)
    .(records=.N, finite_records=sum(is.finite(z)), rejected=sum(is.finite(z)&!is.finite(cleaned)),
      raw_min=if(any(is.finite(z)))min(z[is.finite(z)]) else NA_real_,
      raw_max=if(any(is.finite(z)))max(z[is.finite(z)]) else NA_real_)
  }, by=site][,variable:=variable]
}
audit_raw(vitals[vital_category %chin% c("temperature","temperature_c","temp_c")],"vital_value","temperature")
audit_raw(vitals[vital_category=="spo2"],"vital_value","spo2")
audit_raw(vitals[vital_category %chin% c("map","mean_arterial_pressure")],"vital_value","map_value")
temperature <- rolling_last(
  vitals[vital_category %chin% c("temperature", "temperature_c", "temp_c")],
  "recorded_dttm", "vital_value", "temperature", 12
)
temperature <- fill_nearest_temperature(temperature, vitals, grid, cohort,
  data_cutoff_dttm, shareable_dir, sensitive_dir)
spo2 <- rolling_last(
  vitals[vital_category == "spo2"],
  "recorded_dttm", "vital_value", "spo2", 12
)
map_value <- rolling_last(
  vitals[vital_category %chin% c("map", "mean_arterial_pressure")],
  "recorded_dttm", "vital_value", "map_value", 12
)
rm(vitals); gc(verbose = FALSE)

# ---- Labs --------------------------------------------------------------------

message("Aligning laboratory covariates...")
lab_path <- required_file("clif_labs.parquet")
lab_dataset <- open_dataset(lab_path)
lab_schema <- lab_dataset$schema$names
lab_time_column <- if ("lab_collect_dttm" %in% lab_schema) {
  "lab_collect_dttm"
} else if ("lab_result_dttm" %in% lab_schema) {
  "lab_result_dttm"
} else stop("clif_labs.parquet lacks lab_collect_dttm/lab_result_dttm.")
lab_categories <- lab_dataset |>
  distinct(lab_category) |>
  collect() |>
  pull(lab_category) |>
  as.character()

normalize_lab_category <- function(x) {
  gsub("[^a-z0-9]+", "_", tolower(trimws(as.character(x))))
}
# Explicit blood/serum aliases: do not combine venous and arterial blood gases,
# bilirubin fractions, or urine and serum measurements through broad regexes.
lab_aliases <- list(
  wbc = c("wbc", "white_blood_cell_count", "white_blood_cells"),
  glucose = c("glucose_serum", "serum_glucose", "glucose", "blood_glucose"),
  albumin = c("albumin", "albumin_serum", "serum_albumin"),
  platelet = c("platelet_count", "platelets", "platelet"),
  bilirubin = c("bilirubin_total", "total_bilirubin", "bilirubin"),
  creatinine = c("creatinine", "creatinine_serum", "serum_creatinine")
)
lab_map <- lapply(lab_aliases, function(aliases) {
  lab_categories[normalize_lab_category(lab_categories) %in% aliases]
})
lab_mapping_audit <- rbindlist(lapply(names(lab_map), function(variable) {
  matched <- lab_map[[variable]]
  data.table(
    variable = variable,
    matched_lab_category = if (length(matched)) paste(matched, collapse = ";") else NA_character_,
    mapping_available = length(matched) > 0L,
    definition = if (variable %in% c("ph", "paco2")) "arterial only" else "blood/serum"
  )
}))
write_shareable_csv(lab_mapping_audit, "laboratory_category_mapping.csv")
selected_lab_categories <- unique(unlist(lab_map, use.names = FALSE))
if (length(selected_lab_categories)) {
  labs <- lab_dataset |>
    filter(
      .data$hospitalization_id %in% crosswalk$hospitalization_id,
      .data$lab_category %in% selected_lab_categories
    ) |>
    select(
      hospitalization_id,
      all_of(lab_time_column),
      lab_category,
      lab_value_numeric
    ) |>
    collect() |>
    as.data.table()
  labs[, hospitalization_id := as.character(hospitalization_id)]
  setnames(labs, lab_time_column, "lab_dttm")
  labs <- map_to_selected_blocks(labs)
} else {
  labs <- data.table(
    hospitalization_id = character(),
    hospital_block_id = character(),
    lab_dttm = as.POSIXct(character()),
    lab_category = character(),
    lab_value_numeric = numeric()
  )
}

lab_roll <- function(name, lookback = 24) {
  rolling_last(
    labs[lab_category %in% lab_map[[name]]],
    "lab_dttm",
    "lab_value_numeric",
    name,
    lookback
  )
}
for (v in names(lab_map)) {
 z <- copy(labs[lab_category %in% lab_map[[v]]]); setnames(z,"lab_dttm","recorded_dttm")
 audit_raw(z,"lab_value_numeric",v)
}
wbc <- lab_roll("wbc")
glucose <- lab_roll("glucose")
albumin <- lab_roll("albumin", 48)
platelet <- lab_roll("platelet", 48)
bilirubin <- lab_roll("bilirubin", 48)
creatinine <- lab_roll("creatinine", 48)
rm(labs); gc(verbose = FALSE)

# ---- Respiratory support and GCS --------------------------------------------

message("Aligning respiratory support and GCS...")
respiratory <- read_selected(
  "clif_respiratory_support.parquet",
  c("hospitalization_id", "recorded_dttm", "device_category", "fio2_set")
)
respiratory <- map_to_selected_blocks(respiratory)
respiratory[, device_category_clean := gsub(
  "[^a-z0-9]+", "_", tolower(trimws(as.character(device_category)))
)]
respiratory[, support_code := match(map_support(device_category_clean), SUPPORT_LEVELS)]
audit_raw(respiratory,"fio2_set","fio2")
fio2 <- rolling_last(respiratory, "recorded_dttm", "fio2_set", "fio2", 12)
support <- rolling_last(respiratory, "recorded_dttm", "support_code", "support_code", 12)
rm(respiratory); gc(verbose = FALSE)

assessment <- read_selected(
  "clif_patient_assessments.parquet",
  c("hospitalization_id", "recorded_dttm", "assessment_category", "numerical_value")
)
assessment <- map_to_selected_blocks(assessment)
assessment <- assessment[
  assessment_category == "gcs_total"
]
audit_raw(assessment,"numerical_value","gcs")
gcs <- rolling_last(
  assessment,
  "recorded_dttm",
  "numerical_value",
  "gcs",
  24
)
rm(assessment); gc(verbose = FALSE)

# ---- Merge and construct analysis covariates --------------------------------

message("Constructing longitudinal covariates and SOFA proxy...")
long <- copy(grid)
for (object_name in c(
  "temperature", "spo2", "map_value", "wbc", "glucose", "albumin",
  "platelet", "bilirubin", "creatinine", "fio2", "support", "gcs"
)) {
  long <- merge(
    long,
    get(object_name),
    by = c("hospital_block_id", "m"),
    all.x = TRUE,
    sort = FALSE
  )
}
setorder(long, site, hospital_block_id, m)

long[fio2 > 1.5, fio2 := fio2 / 100]
long[fio2 < 0.21 | fio2 > 1, fio2 := NA_real_]
long[spo2 < 50 | spo2 > 100, spo2 := NA_real_]
long[, sf_ratio := fifelse(is.finite(spo2) & spo2 <= 97, spo2 / fio2, NA_real_)]
long[, support_category := SUPPORT_LEVELS[as.integer(support_code)]]
long[is.na(support_category),support_category:="Unknown"]
for (category in SUPPORT_LEVELS[-1]) {
  column <- paste0("support_",category)
  long[, (column) := as.integer(support_category==category)]
}
long[,support_unknown:=as.integer(support_category=="Unknown")]
long[,support_code:=NULL]
write_shareable_csv(rbindlist(raw_screen_audit,fill=TRUE),"raw_measurement_screening_counts.csv")
write_shareable_csv(measurement_rules(),"screening_ranges.csv")
write_shareable_csv(data.table(variable=c("arterial_pH","arterial_PaCO2"),
 reason="Excluded a priori because of severe missingness; no values, missingness indicators or lags enter models"),
 "excluded_covariates.csv")

# Simplified non-respiratory SOFA: coagulation, liver, cardiovascular MAP,
# renal creatinine, and neurologic GCS. Vasopressor dose and urine-output
# criteria are not available in this first-pass construction.
long[, sofa_platelet := fcase(
  platelet < 20, 4, platelet < 50, 3, platelet < 100, 2,
  platelet < 150, 1, is.finite(platelet), 0, default = NA_real_
)]
long[, sofa_bilirubin := fcase(
  bilirubin >= 12, 4, bilirubin >= 6, 3, bilirubin >= 2, 2,
  bilirubin >= 1.2, 1, is.finite(bilirubin), 0, default = NA_real_
)]
long[, sofa_map := fcase(
  map_value < 70, 1, is.finite(map_value), 0, default = NA_real_
)]
long[, sofa_creatinine := fcase(
  creatinine >= 5, 4, creatinine >= 3.5, 3, creatinine >= 2, 2,
  creatinine >= 1.2, 1, is.finite(creatinine), 0, default = NA_real_
)]
long[, sofa_gcs := fcase(
  gcs < 6, 4, gcs < 10, 3, gcs < 13, 2, gcs < 15, 1,
  is.finite(gcs), 0, default = NA_real_
)]

sofa_parts <- c(
  "sofa_platelet", "sofa_bilirubin", "sofa_map",
  "sofa_creatinine", "sofa_gcs"
)
long[, sofa_components_observed := rowSums(is.finite(as.matrix(.SD))), .SDcols=sofa_parts]
long[, nonresp_sofa_proxy := rowSums(.SD), .SDcols = sofa_parts]

raw_tv <- c(
  "temperature", "wbc", "glucose", "albumin",
  "nonresp_sofa_proxy"
)

plausible_limits <- lapply(setNames(raw_tv,raw_tv),function(v) {
 r <- measurement_rules()[variable==v]; c(r$lower,r$upper)
})
# Estimate transformations/imputation from untreated decision rows during grace,
# so post-discharge follow-up does not dominate those choices.
processing_rows <- long$m <= tau & long$m <= long$S

covariate_qa <- rbindlist(lapply(raw_tv, function(variable) {
  limits <- plausible_limits[[variable]]
  x <- as.numeric(long[[variable]])
  x[x < limits[1L] | x > limits[2L]] <- NA_real_
  observed <- is.finite(x)
  observed_fraction <- mean(observed)
  if (any(observed & processing_rows)) {
    winsor <- quantile(x[observed & processing_rows], c(0.005, 0.995), na.rm = TRUE, names = FALSE)
    x[observed] <- pmin(pmax(x[observed], winsor[1L]), winsor[2L])
    fill <- median(x[observed & processing_rows], na.rm = TRUE)
  } else {
    winsor <- c(NA_real_, NA_real_)
    fill <- 0
  }
  missing_name <- paste0(variable, "_missing")
  long[, (missing_name) := as.integer(!observed)]
  x[!observed] <- fill
  if (variable %in% c("wbc", "glucose")) x <- log1p(x)
  scale_value <- sd(x[processing_rows])
  if (!is.finite(scale_value) || scale_value <= 0) scale_value <- 1
  center_value <- mean(x[processing_rows])
  model_name <- paste0(variable, "_model")
  if (variable == "nonresp_sofa_proxy") {
    long[, (model_name) := x]
  } else {
    long[, (model_name) := (x - center_value) / scale_value]
  }
  data.table(
    variable = variable,
    observed_fraction = observed_fraction,
    observed_fraction_untreated_during_grace = mean(observed[processing_rows]),
    imputation_value = fill,
    winsor_lower = winsor[1L],
    winsor_upper = winsor[2L],
    model_variable = model_name,
    usable = any(observed & processing_rows) && uniqueN(x[processing_rows]) > 1L
  )
}))
write_shareable_csv(covariate_qa, "covariate_availability_and_processing.csv")
write_shareable_csv(
  data.table(
    target = names(lab_map),
    matched_clif_categories = vapply(
      lab_map,
      function(x) paste(x, collapse = ";"),
      character(1)
    )
  ),
  "matched_lab_categories.csv"
)

baseline_covariates <- c(
  "baseline_sf_model", "age_z", "female", "sex_missing", "baseline_support_imv",
  "baseline_support_nippv", "initial_location_icu", "initial_location_ed",
  "admission_year_z", "admission_season_sin", "admission_season_cos",
  comorbidity_names
)
baseline_covariates <- baseline_covariates[
  vapply(long[, ..baseline_covariates], function(x) uniqueN(x) > 1L, logical(1))
]
usable_tv <- covariate_qa[usable == TRUE, variable]
time_varying_covariates <- c(
  covariate_qa[usable == TRUE, model_variable],
  paste0(usable_tv, "_missing"),
  paste0("support_",SUPPORT_LEVELS[-1]), "support_unknown"
)
time_varying_covariates <- time_varying_covariates[
  vapply(long[, ..time_varying_covariates], function(x) uniqueN(x) > 1L,
         logical(1))
]

if (!length(time_varying_covariates)) {
  stop("No requested time-varying covariates passed availability checks.")
}

# Record the exact, post-availability-check terms used to estimate the
# stabilized IPCW components. This is more informative than the candidate list
# because constant baseline variables and poorly observed longitudinal
# variables have already been removed.
ipcw_terms <- list(
  initiation_denominator = c(
    ".time", baseline_covariates, time_varying_covariates,
    paste0("lag_", time_varying_covariates)
  ),
  initiation_numerator = ".time",
  censoring_denominator = c(".time", baseline_covariates, ".trt"),
  censoring_numerator = ".time"
)
ipcw_model_specification <- rbindlist(lapply(
  names(ipcw_terms),
  function(component) {
    data.table(
      model_component = component,
      term = ipcw_terms[[component]],
      term_type = fcase(
        ipcw_terms[[component]] == ".time", "interval time",
        ipcw_terms[[component]] == ".trt", "current treatment",
        ipcw_terms[[component]] %chin% baseline_covariates, "baseline",
        grepl("^lag_", ipcw_terms[[component]]), "lagged time-varying",
        default = "current time-varying"
      ),
      note = if (grepl("^censoring", component)) {
        paste(
          "Administrative censoring occurs only when the 28-day horizon",
          "extends beyond the prespecified data cutoff."
        )
      } else {
        paste(
          "The analysis uses the centrally aggregated interval-specific",
          "initiation numerator broadcast to all sites; denominator models",
          "are hospital-specific for both federated and pooled",
          "site-stratified CCW."
        )
      }
    )
  }
))
write_shareable_csv(ipcw_model_specification, "ipcw_model_specification.csv")

long[, init_m := as.integer(S == m & S <= M)]
long[, event_m := as.integer(delta == 1L & T_event == m)]
long[, censor_m := as.integer(
  T_censor <= M & T_censor == m & (delta == 0L | T_event > T_censor)
)]

required_analysis <- c(baseline_covariates, time_varying_covariates)
if (anyNA(long[, ..required_analysis]) ||
    any(!is.finite(as.matrix(long[, ..required_analysis])))) {
  stop("Non-finite values remain in analysis covariates.")
}
stopifnot(
  all(long$event_m + long$censor_m <= 1L),
  all(long[, .N, by = hospital_block_id]$N <= M),
  uniqueN(cohort$patient_id) == uniqueN(patient_sites$patient_id)
)

# ---- Save prepared data and diagnostics -------------------------------------

design <- list(
  source = "Hopkins multisite EHR data in CLIF 2.1 format",
  population = BASELINE_SF_DEFINITION,
  sf_role = "baseline_only",
  source_cohort_patients = source_cohort_patients,
  analysis_unit = "one eligible six-hour hospital block per patient",
  time_zero = "respiratory-failure start",
  treatment = "first qualifying systemic glucocorticoid administration",
  strategy_1 = "initiate qualifying glucocorticoid within the grace period",
  strategy_0 = "do not initiate qualifying glucocorticoid within the grace period; later initiation allowed",
  outcome = PRIMARY_OUTCOME,
  alive_discharge = "Other alive discharge is event-free, not censoring; subsequent post-discharge death excluded",
  step_hours = step_hours,
  grace_hours = grace_hours,
  tau_intervals = tau,
  horizon_days = horizon_days,
  horizon_intervals = M,
  data_cutoff = format(data_cutoff_dttm, "%Y-%m-%d"),
  observed_hospitals = sort(unique(cohort$site)),
  maximum_sampled_patients = max_patients,
  sample_seed = sample_seed,
  full_cohort_patients = full_patients,
  full_cohort_hospital_blocks = full_blocks,
  full_component_hospitalizations = full_component_hospitalizations,
  sampled_patients = uniqueN(cohort$patient_id),
  sampled_hospital_blocks = uniqueN(cohort$hospital_block_id),
  sampled_component_hospitalizations = uniqueN(crosswalk$hospitalization_id),
  baseline_covariates = baseline_covariates,
  time_varying_covariates = time_varying_covariates,
  omitted_requested_features = c(
    "Arterial pH and PaCO2 excluded a priori for severe missingness",
    "full non-respiratory SOFA: vasopressor dose and urine-output components omitted"
  )
)

  stopifnot(!any(grepl("^(lag_)?sf_ratio(_model|_missing)?$",c(baseline_covariates,time_varying_covariates))),
            all(is.finite(long$baseline_sf)),all(cohort$analysis_baseline_spo2<=97))
  write_covariate_eda(long, cohort, shareable_dir, tau, step_hours)
  write_parquet(as.data.frame(cohort), file.path(sensitive_dir, "pneumonia_ccw_cohort.parquet"))
  write_parquet(as.data.frame(long), file.path(sensitive_dir, "pneumonia_ccw_long.parquet"))
  write_parquet(as.data.frame(crosswalk), file.path(sensitive_dir, "sampled_hospital_block_crosswalk.parquet"))
  saveRDS(design, file.path(sensitive_dir, "analysis_design.rds"))
  list(design=design, long=as.data.frame(long), cohort=cohort,
       baseline_covariates=baseline_covariates, time_varying_covariates=time_varying_covariates,
       tau=tau, M=M, step_hours=step_hours, horizon_days=horizon_days,
       cohort_outcome_summary=cohort_outcome_summary)
}
