# Baseline eligibility sensitivity: no selection using future saturations.
BASELINE_SF_DEFINITION <- "Baseline S/F complete: latest valid prior SpO2 50-97% and FiO2 0.21-1, each within 12 hours of original time zero"

select_baseline_sf <- function(data_dir,cohort,crosswalk,aggregate_dir,sensitive_dir) {
 c <- data.table::copy(data.table::as.data.table(cohort))
 cw <- data.table::copy(data.table::as.data.table(crosswalk))
 c[,hospital_block_id:=as.character(hospital_block_id)]
 cw[,`:=`(hospitalization_id=as.character(hospitalization_id),hospital_block_id=as.character(hospital_block_id))]
 cw<-unique(cw[hospital_block_id%in%c$hospital_block_id,.(hospitalization_id,hospital_block_id)])
 if(anyDuplicated(cw$hospitalization_id))stop("Encounter maps to multiple hospital blocks")
 ids<-cw$hospitalization_id
 grid<-c[,.(hospital_block_id,source_site_id,admission_dttm,discharge_dttm,
             decision_dttm=as.POSIXct(respiratory_failure_start,tz="UTC"))]
 norm<-function(x)gsub("[^a-z0-9]+","_",tolower(trimws(as.character(x))))
 ds<-arrow::open_dataset(file.path(data_dir,"clif_vitals.parquet"))
 cats<-ds |> dplyr::distinct(vital_category) |> dplyr::collect()
 categories<-cats$vital_category[which(norm(cats$vital_category)=="spo2")]
 sp<-ds |> dplyr::filter(.data$hospitalization_id%in%ids,.data$vital_category%in%categories) |>
   dplyr::select(hospitalization_id,recorded_dttm,vital_value) |> dplyr::collect() |> data.table::as.data.table()
 fi<-arrow::open_dataset(file.path(data_dir,"clif_respiratory_support.parquet")) |>
   dplyr::filter(.data$hospitalization_id%in%ids) |>
   dplyr::select(hospitalization_id,recorded_dttm,fio2_set) |> dplyr::collect() |> data.table::as.data.table()
 latest<-function(x,value_col,variable) {
  x[,hospitalization_id:=as.character(hospitalization_id)]
  x<-merge(x,cw,by="hospitalization_id",allow.cartesian=FALSE)
  x<-merge(x,grid,by="hospital_block_id",allow.cartesian=TRUE)
  x[,`:=`(seconds=as.numeric(as.POSIXct(recorded_dttm,tz="UTC")),
           value=screen_measurements(get(value_col),variable))]
  # Keep the latest usable physical measurement, including saturation 98-100.
  # The <=97 eligibility restriction is applied only AFTER choosing that reading.
  x<-x[is.finite(value)&is.finite(seconds)&seconds>=as.numeric(admission_dttm)&
        seconds<=as.numeric(discharge_dttm)&seconds<=as.numeric(decision_dttm)&
        seconds>=as.numeric(decision_dttm)-12*3600]
  x[,measured_seconds:=seconds]
  data.table::setkey(x,hospital_block_id,seconds)
  if(!nrow(x))return(data.table::data.table(hospital_block_id=grid$hospital_block_id,value=NA_real_,age_hours=NA_real_))
  query<-copy(grid);query[,decision_seconds:=as.numeric(decision_dttm)]
  x[query,on=.(hospital_block_id,seconds<=decision_seconds),mult="last",
    .(hospital_block_id=i.hospital_block_id,value=x.value,
      age_hours=(i.decision_seconds-x.measured_seconds)/3600)]
 }
 sp<-latest(sp,"vital_value","spo2");fi<-latest(fi,"fio2_set","fio2")
 audit<-c[,.(hospital_block_id,patient_id,source_site_id)]
 audit[sp,on="hospital_block_id",`:=`(analysis_baseline_spo2=i.value,baseline_spo2_age_hours=i.age_hours)]
 audit[fi,on="hospital_block_id",`:=`(analysis_baseline_fio2=i.value,baseline_fio2_age_hours=i.age_hours)]
 audit[,baseline_sf:=analysis_baseline_spo2/analysis_baseline_fio2]
 audit[,baseline_sf_exclusion:=data.table::fcase(
   !is.finite(analysis_baseline_spo2),"No usable prior SpO2 within 12h",
   analysis_baseline_spo2>97,"Latest valid baseline SpO2 >97%",
   !is.finite(analysis_baseline_fio2),"SpO2 eligible but no usable prior FiO2 within 12h",
   !is.finite(screen_measurements(baseline_sf,"sf_ratio")),"Baseline S/F outside screening range",
   default="Retained")]
 audit[,baseline_sf_eligible:=baseline_sf_exclusion=="Retained"]
 aggregate<-audit[,.(patients=.N),by=.(site=source_site_id,baseline_sf_exclusion)]
 aggregate[,source_site_patients:=sum(patients),by=site]
 aggregate[,percent_of_site:=100*patients/source_site_patients]
 data.table::fwrite(aggregate,file.path(aggregate_dir,"baseline_sf_selection_by_site.csv"),na="")
 total<-audit[,.(patients=.N),by=baseline_sf_exclusion]
 total[,`:=`(source_cohort_patients=nrow(audit),percent_of_source=100*patients/nrow(audit))]
 data.table::fwrite(total,file.path(aggregate_dir,"baseline_sf_selection_summary.csv"),na="")
 data.table::fwrite(data.table::data.table(parameter=c("definition","selection_time","subsequent_high_spo2","model_role","sf_missing_indicators"),
  value=c(BASELINE_SF_DEFINITION,"Original time zero, unchanged","Retain all subsequent follow-up; do not censor/exclude for later SpO2 >97%",
   "Baseline covariate only; no time-varying S/F or lagged S/F","None")),file.path(aggregate_dir,"baseline_sf_sensitivity_definition.csv"))
 arrow::write_parquet(audit,file.path(sensitive_dir,"baseline_sf_selection_rows.parquet"))
 selected<-audit[baseline_sf_eligible==TRUE]
 if(!nrow(selected))stop("No patients have eligible complete baseline S/F; inspect baseline_sf_selection_summary.csv")
 # Avoid duplicating fields when the function is called with a previously selected cohort.
 added<-setdiff(names(selected),names(c))
 c<-merge(c,selected[,c("hospital_block_id",added),with=FALSE],by="hospital_block_id",all=FALSE,sort=FALSE)
 cw<-crosswalk[hospital_block_id%in%c$hospital_block_id]
 message("Baseline S/F sensitivity: retained ",nrow(c)," / ",nrow(audit)," patients")
 list(cohort=c,crosswalk=cw,audit=audit)
}

encode_baseline_sf <- function(cohort) {
 c<-data.table::copy(data.table::as.data.table(cohort))
 if(any(!is.finite(c$baseline_sf)))stop("Baseline S/F must be complete; no imputation is permitted")
 limits<-as.numeric(quantile(c$baseline_sf,c(.005,.995),names=FALSE))
 value<-pmin(pmax(c$baseline_sf,limits[1]),limits[2])
 center<-mean(value);scale<-sd(value)
 if(!is.finite(scale)||scale<=0)scale<-1
 c[,baseline_sf_model:=(value-center)/scale]
 list(cohort=c,processing=data.table::data.table(variable="baseline_sf",patients=nrow(c),missing=0L,
   winsor_lower=limits[1],winsor_upper=limits[2],center=center,scale=scale,imputation="None"))
}
