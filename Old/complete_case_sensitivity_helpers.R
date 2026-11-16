# Sensitivity-local compatibility for a legitimately empty clinical covariate
# set. Base R paste0("lag_", character()) returns "lag_", not character().
# Do not create that phantom predictor; the source method library is unchanged.
sensitivity_lags <- function(x) if(length(x)) paste0("lag_", x) else character()
g6_prepare_data <- local({
  original <- g6_prepare_data
  function(...) {
    answer <- original(...)
    if(identical(answer$lag_names, "lag_") && !"lag_" %in% names(answer$data))
      answer$lag_names <- character()
    answer
  }
})

# Shared rules for 0701/0901. All outputs are separate from the primary analyses.
# Clinical mappings/lookbacks follow 09; missing values are never median-filled.
# Prior-interval model terms retain the existing initial-lag=0 convention.
select_complete_cases <- function(cohort,long,baseline,tv,output) {
  invalid <- function(x) !is.finite(as.numeric(x))
  first <- long[m == 1L]
  audit <- rbindlist(c(lapply(baseline,function(v) data.table(variable=v,type="baseline",
    denominator=nrow(first),missing=sum(invalid(first[[v]])))),
    lapply(tv,function(v) data.table(variable=v,type="time-varying",
      denominator=nrow(long),missing=sum(invalid(long[[v]]))))))
  audit[, missing_fraction := missing/denominator]
  audit[, retain_by_missingness := missing_fraction < .30]
  fwrite(audit,file.path(output,"sensitivity_variable_selection.csv"))
  b <- audit[type=="baseline" & retain_by_missingness==TRUE,variable]
  t <- audit[type=="time-varying" & retain_by_missingness==TRUE,variable]
  required <- c(b,t)
  bad <- if(length(required)) Reduce(`|`, lapply(required,function(v) invalid(long[[v]]))) else rep(FALSE,nrow(long))
  excluded <- unique(long$hospital_block_id[bad])
  block <- unique(long[,.(hospital_block_id,patient_id,site)])
  block[, retained := !hospital_block_id %in% excluded]
  fwrite(block[,.(before_patients=.N,excluded_patients=sum(!retained),
    retained_patients=sum(retained)),by=site],file.path(output,"sensitivity_attrition_by_site.csv"))
  fwrite(data.table(step=c("Before complete-case restriction","Excluded: missing retained confounder in any observed interval","Final complete-case cohort"),
    patients=c(nrow(block),sum(!block$retained),sum(block$retained))),file.path(output,"sensitivity_attrition.csv"))
  if(!any(block$retained)) stop("No complete-case blocks remain. See sensitivity_variable_selection.csv and sensitivity_attrition.csv; rule was not relaxed.")
  cohort <- cohort[!hospital_block_id %in% excluded]
  long <- long[!hospital_block_id %in% excluded]
  b <- b[vapply(b,function(v) uniqueN(long[[v]])>1L,logical(1))]
  t <- t[vapply(t,function(v) uniqueN(long[[v]])>1L,logical(1))]
  audit[, used_in_models := variable %in% c(b,t)]
  audit[, reason := fifelse(!retain_by_missingness,"omitted: missingness >=30%",
    fifelse(used_in_models,"retained; complete-case restriction applied","omitted: constant after restriction"))]
  fwrite(audit,file.path(output,"sensitivity_variable_selection.csv"))
  if(any(grepl("_missing$",c(b,t)))) stop("Missingness indicator entered sensitivity models")
  if(length(c(b,t)) && any(!is.finite(as.matrix(long[,c(b,t),with=FALSE])))) stop("Missing sensitivity model inputs remain")
  list(cohort=cohort,long=long,baseline=b,tv=t)
}

prepare_complete_case_inputs <- function(data_dir, cohort_dir, treatment_dir, sensitive_dir,
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
    fwrite(as.data.table(x), file.path(shareable_dir, paste0("preselection_", filename)), na="")
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

message("Building baseline hospital-block cohort and 28-day mortality outcome...")
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

# Use patient death_dttm for 28-day death after respiratory-failure start,
# including deaths after discharge. Fall back to an expired-discharge timestamp
# only when death_dttm is absent. Administrative censoring uses the explicit
# study data cutoff.
cohort[, event_dttm := as.POSIXct(death_dttm, tz = "UTC")]
cohort[
  is.na(event_dttm) & tolower(discharge_category) == "expired",
  event_dttm := as.POSIXct(discharge_dttm, tz = "UTC")
]
cohort[, event_hours := as.numeric(difftime(
  event_dttm,
  index_dttm,
  units = "hours"
))]
cohort[
  !is.finite(event_hours) | event_hours < 0,
  `:=`(event_dttm = as.POSIXct(NA, tz = "UTC"), event_hours = NA_real_)
]
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
  baseline_support_nippv = as.integer(respiratory_support_category == "nippv"),
  admission_year_z = (as.integer(format(index_dttm, "%Y")) - 2021) / 5,
  admission_season_sin = sin(2 * pi * as.integer(format(index_dttm, "%m")) / 12),
  admission_season_cos = cos(2 * pi * as.integer(format(index_dttm, "%m")) / 12)
)]

cohort[sex_missing == 1L, female := NA_real_]

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
initial_location[is.na(location_text) | !nzchar(trimws(location_text)),
                 c("initial_location_icu", "initial_location_ed") := NA_integer_]
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
# by treatment-confounding and by deaths that occur before the grace period ends.
deaths_within_grace <- cohort[
  delta == 1L & is.finite(event_hours) & event_hours <= grace_hours,
  .N
]
deaths_after_grace <- cohort[
  delta == 1L & is.finite(event_hours) & event_hours > grace_hours,
  .N
]
deaths_28d <- cohort[delta == 1L, .N]
cohort_mortality_summary <- data.table(
  measure = c(
    "final sampled cohort",
    paste0("death within ", grace_hours, " hours of time zero"),
    paste0("death after ", grace_hours, " hours through day ", horizon_days),
    paste0("all deaths within ", horizon_days, " days"),
    paste0("no observed death within ", horizon_days, " days")
  ),
  n = c(
    nrow(cohort),
    deaths_within_grace,
    deaths_after_grace,
    deaths_28d,
    nrow(cohort) - deaths_28d
  )
)
cohort_mortality_summary[, percent_of_cohort := 100 * n / nrow(cohort)]
write_shareable_csv(cohort_mortality_summary, "cohort_mortality_summary.csv")

mortality_by_observed_treatment <- cohort[
  ,
  .(
    hospital_blocks = .N,
    deaths_within_grace = sum(
      delta == 1L & is.finite(event_hours) & event_hours <= grace_hours
    ),
    deaths_28d = sum(delta == 1L)
  ),
  by = .(glucocorticoid_within_grace = A_tau)
]
mortality_by_observed_treatment[, `:=`(
  deaths_within_grace_percent = 100 * deaths_within_grace / hospital_blocks,
  deaths_28d_percent = 100 * deaths_28d / hospital_blocks,
  interpretation = paste(
    "Descriptive only; observed treatment groups are not baseline randomized",
    "groups and early death can prevent treatment initiation."
  )
)]
write_shareable_csv(
  mortality_by_observed_treatment,
  "cohort_mortality_by_observed_grace_period_treatment.csv"
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
    S,
    A_tau,
    T_event,
    T_censor,
    T_obs,
    delta,
    age_z,
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
  x <- measurement[
    !is.na(get(time_col)) & is.finite(as.numeric(get(value_col))),
    .(
      hospital_block_id = as.character(hospital_block_id),
      recorded_dttm = as.POSIXct(get(time_col), tz = "UTC"),
      measurement_seconds = as.numeric(get(time_col)),
      value = as.numeric(get(value_col))
    )
  ]
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
temperature <- rolling_last(
  vitals[vital_category %chin% c("temperature", "temperature_c", "temp_c")],
  "recorded_dttm", "vital_value", "temperature", 12
)
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
  ph = c("ph_arterial", "arterial_ph", "ph_blood_arterial"),
  paco2 = c("pco2_arterial", "arterial_pco2", "paco2",
            "carbon_dioxide_arterial", "arterial_carbon_dioxide"),
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
wbc <- lab_roll("wbc")
glucose <- lab_roll("glucose")
albumin <- lab_roll("albumin", 48)
ph <- lab_roll("ph")
paco2 <- lab_roll("paco2")
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
respiratory[, support_score := fcase(
  device_category_clean == "room_air", 0,
  device_category_clean %chin% c("nasal_cannula", "face_mask", "trach_collar"), 1,
  device_category_clean %chin% c("high_flow_nc", "nippv", "cpap"), 2,
  device_category_clean == "imv", 3,
  default = NA_real_
)]
fio2 <- rolling_last(respiratory, "recorded_dttm", "fio2_set", "fio2", 12)
support <- rolling_last(
  respiratory,
  "recorded_dttm",
  "support_score",
  "support_score",
  12
)
rm(respiratory); gc(verbose = FALSE)

assessment <- read_selected(
  "clif_patient_assessments.parquet",
  c("hospitalization_id", "recorded_dttm", "assessment_category", "numerical_value")
)
assessment <- map_to_selected_blocks(assessment)
assessment <- assessment[
  assessment_category == "gcs_total"
]
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
  "temperature", "spo2", "map_value", "wbc", "glucose", "albumin", "ph",
  "paco2", "platelet", "bilirubin", "creatinine", "fio2", "support", "gcs"
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
long[, sf_ratio := spo2 / fio2]

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
long[, nonresp_sofa_proxy := rowSums(.SD), .SDcols = sofa_parts]

raw_tv <- c(
  "temperature", "wbc", "glucose", "albumin", "ph", "paco2",
  "sf_ratio", "support_score", "nonresp_sofa_proxy"
)

plausible_limits <- list(
  temperature = c(25, 45),
  wbc = c(0.1, 200),
  glucose = c(20, 1000),
  albumin = c(0.5, 7),
  ph = c(6.5, 8.0),
  paco2 = c(5, 200),
  sf_ratio = c(20, 700),
  support_score = c(0, 3),
  nonresp_sofa_proxy = c(0, 17)
)

covariate_qa <- rbindlist(lapply(raw_tv, function(variable) {
  limits <- plausible_limits[[variable]]
  x <- as.numeric(long[[variable]])
  x[x < limits[1L] | x > limits[2L]] <- NA_real_
  observed <- is.finite(x)
  observed_fraction <- mean(observed)
  if (any(observed)) {
    winsor <- quantile(x[observed], c(0.005, 0.995), na.rm = TRUE, names = FALSE)
    x[observed] <- pmin(pmax(x[observed], winsor[1L]), winsor[2L])
    fill <- NA_real_
  } else {
    winsor <- c(NA_real_, NA_real_)
    fill <- NA_real_
  }
  x[!observed] <- NA_real_
  if (variable %in% c("wbc", "glucose")) x <- log1p(x)
  scale_value <- sd(x, na.rm = TRUE)
  if (!is.finite(scale_value) || scale_value <= 0) scale_value <- 1
  center_value <- if(any(is.finite(x))) mean(x, na.rm = TRUE) else 0
  model_name <- paste0(variable, "_model")
  if (variable %in% c("support_score", "nonresp_sofa_proxy")) {
    long[, (model_name) := x]
  } else {
    long[, (model_name) := (x - center_value) / scale_value]
  }
  data.table(
    variable = variable,
    observed_fraction = observed_fraction,
    imputation_value = fill,
    winsor_lower = winsor[1L],
    winsor_upper = winsor[2L],
    model_variable = model_name,
    usable = observed_fraction > 0.70 && uniqueN(x[is.finite(x)]) > 1L
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
  "age_z", "female", "baseline_support_imv",
  "baseline_support_nippv", "initial_location_icu", "initial_location_ed",
  "admission_year_z", "admission_season_sin", "admission_season_cos",
  comorbidity_names
)
# Screen baseline at patient level, clinical histories at observed person-interval level.
# Fractions are computed ONCE before deleting cases; exactly 30% is excluded.
time_varying_covariates <- paste0(raw_tv, "_model")
selected <- select_complete_cases(cohort, long, baseline_covariates,
                                  time_varying_covariates, shareable_dir)
cohort <- selected$cohort; long <- selected$long
baseline_covariates <- selected$baseline
 time_varying_covariates <- selected$tv
crosswalk <- crosswalk[hospital_block_id %in% cohort$hospital_block_id]
patient_sites <- unique(cohort[, .(patient_id, site)])
full_patients <- uniqueN(cohort$patient_id); full_blocks <- nrow(cohort)
full_component_hospitalizations <- nrow(crosswalk)

# Record the exact, post-availability-check terms used to estimate the
# stabilized IPCW components. This is more informative than the candidate list
# because constant baseline variables and poorly observed longitudinal
# variables have already been removed.
ipcw_terms <- list(
  initiation_denominator = c(
    ".time", baseline_covariates, time_varying_covariates,
    sensitivity_lags(time_varying_covariates)
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
long <- long[, unique(c("hospital_block_id", "patient_id", "site", "m", "init_m", "event_m", "censor_m", baseline_covariates, time_varying_covariates)), with=FALSE]

design <- list(
  source = "Hopkins multisite EHR data in CLIF 2.1 format",
  analysis_unit = "one eligible six-hour hospital block per patient",
  time_zero = "respiratory-failure start",
  treatment = "first qualifying systemic glucocorticoid administration",
  strategy_1 = "initiate qualifying glucocorticoid within the grace period",
  strategy_0 = "do not initiate qualifying glucocorticoid through 28 days",
  outcome = "death within 28 days of respiratory-failure start",
  alive_discharge = paste(
    "not censoring; patient death_dttm provides 28-day follow-up and",
    "discharge Expired is used only when death_dttm is absent"
  ),
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
    "hospital type: not available as a patient-level varying field",
    "full non-respiratory SOFA: vasopressor dose and urine-output components omitted"
  )
)

  write_parquet(as.data.frame(cohort), file.path(sensitive_dir, "pneumonia_ccw_cohort.parquet"))
  write_parquet(as.data.frame(long), file.path(sensitive_dir, "pneumonia_ccw_long.parquet"))
  write_parquet(as.data.frame(crosswalk), file.path(sensitive_dir, "sampled_hospital_block_crosswalk.parquet"))
  saveRDS(design, file.path(sensitive_dir, "analysis_design.rds"))
  design$sensitivity <- "strict complete-case; confounders with missingness >=30% omitted; no missingness indicators or covariate imputation"
  design$missingness_window <- "full observed follow-up, before complete-case deletion"
  saveRDS(design, file.path(sensitive_dir, "analysis_design.rds"))
  list(design=design, long=as.data.frame(long), cohort=cohort,
       crosswalk=crosswalk, baseline_covariates=baseline_covariates,
       time_varying_covariates=time_varying_covariates, step_hours=step_hours,
       grace_hours=grace_hours, tau=tau, M=M, horizon_days=horizon_days)
}

