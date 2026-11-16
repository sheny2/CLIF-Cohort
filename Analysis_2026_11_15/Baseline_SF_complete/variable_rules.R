# Single source of analytic screening and categorical respiratory support rules.
# Bounds are inclusive, intentionally broad, and are analytic exclusions rather
# than assertions that every excluded extreme is a verified source-data error.
measurement_rules <- function() data.table::data.table(
 variable=c("temperature","spo2","map_value","wbc","glucose","albumin",
            "platelet","bilirubin","creatinine","gcs","fio2","sf_ratio","nonresp_sofa_proxy"),
 lower=c(25,50,20,.1,20,.5,0,0,0,3,.21,20,0),
 upper=c(45,100,200,200,1000,7,2000,100,30,15,1,700,17),
 unit=c("C","percent","mmHg","10^3/uL","mg/dL","g/dL","10^3/uL","mg/dL","mg/dL",
        "3-15","fraction","percent / fraction","proxy points"),
 action="Outside range becomes missing BEFORE temporal selection; never delete a patient solely for one invalid record",
 additional_rule=c(rep("",10),"Convert 21-100 percent form to fraction before screening",
                   "Only construct when the selected SpO2 is <=97%; never cap saturation", "Requires all five proxy components"))

screen_measurements <- function(x, variable) {
 z <- suppressWarnings(as.numeric(x))
 if(variable=="fio2") z[is.finite(z)&z>=21&z<=100] <- z[is.finite(z)&z>=21&z<=100]/100
 rule <- measurement_rules(); v <- variable; rule <- rule[rule$variable==v]
 if(nrow(rule)) z[!is.finite(z)|z<rule$lower|z>rule$upper] <- NA_real_
 else z[!is.finite(z)] <- NA_real_
 z
}
SUPPORT_LEVELS <- c("room_air","conventional_oxygen","hfno","nippv_cpap","imv")
map_support <- function(x) {
 x <- gsub("[^a-z0-9]+","_",tolower(trimws(as.character(x))))
 data.table::fcase(x=="room_air","room_air",
  x %in% c("nasal_cannula","face_mask","trach_collar"),"conventional_oxygen",
  x=="high_flow_nc","hfno", x %in% c("nippv","cpap"),"nippv_cpap",
  x=="imv","imv",default=NA_character_)
}

# Keep available past-12-hour temperatures first. Fill remaining INPATIENT gaps
# by nearest valid same-block reading in either direction, with earlier tie-break.
# No maximum distance is imposed for the fallback; distances and future use are audited.
fill_nearest_temperature <- function(temperature,vitals,grid,cohort,cutoff,aggregate_dir,sensitive_dir) {
 x <- copy(vitals[vital_category %chin% c("temperature","temperature_c","temp_c")])
 x[,value:=screen_measurements(vital_value,"temperature")]
 x <- merge(x,cohort[,.(hospital_block_id,admission_dttm,discharge_dttm,event_dttm)],by="hospital_block_id")
 x <- x[is.finite(value)&!is.na(recorded_dttm)&recorded_dttm>=admission_dttm &
   recorded_dttm<=discharge_dttm & recorded_dttm<=cutoff &
   (is.na(event_dttm)|recorded_dttm<=event_dttm)]
 x <- x[,.(hospital_block_id,seconds=as.numeric(recorded_dttm),measured=as.numeric(recorded_dttm),value)]
 setorder(x,hospital_block_id,seconds,value)
 x <- unique(x,by=c("hospital_block_id","seconds"));setkey(x,hospital_block_id,seconds)
 a <- merge(grid[,.(hospital_block_id,m,site,decision_dttm,admission_dttm,discharge_dttm)],temperature,
            by=c("hospital_block_id","m"),sort=FALSE)
 a[,`:=`(before_fill_missing=!is.finite(temperature),fallback_used=FALSE,
         future_reading=FALSE,distance_hours=NA_real_)]
 target <- a[before_fill_missing & decision_dttm>=admission_dttm & decision_dttm<=discharge_dttm]
 if(nrow(x)&&nrow(target)) {
  target[,seconds:=as.numeric(decision_dttm)]
  prior <- x[target,on=.(hospital_block_id,seconds<=seconds),mult="last",
   .(hospital_block_id=i.hospital_block_id,m=i.m,prior_value=x.value,prior_time=x.measured)]
  nextval <- x[target,on=.(hospital_block_id,seconds>=seconds),mult="first",
   .(hospital_block_id=i.hospital_block_id,m=i.m,next_value=x.value,next_time=x.measured)]
  target <- merge(target,prior,by=c("hospital_block_id","m"))
  target <- merge(target,nextval,by=c("hospital_block_id","m"))
  target[, use_prior := is.finite(prior_time)&(!is.finite(next_time)|seconds-prior_time<=next_time-seconds)]
  target[,`:=`(fill=fifelse(use_prior,prior_value,next_value),
   measured=fifelse(use_prior,prior_time,next_time))]
  target[,`:=`(distance=abs(measured-seconds)/3600,future=is.finite(measured)&measured>seconds)]
  a[target,on=.(hospital_block_id,m),`:=`(temperature=i.fill,fallback_used=is.finite(i.fill),
    future_reading=i.future,distance_hours=i.distance)]
 }
 a[, decision_hours := (m-1)*6]
 a[, period := fifelse(decision_hours<as.numeric(Sys.getenv("TREATMENT_GRACE_HOURS","24")),"during_grace","after_grace")]
 summary <- a[,.(rows=.N,missing_before=sum(before_fill_missing),filled_nearest=sum(fallback_used),
  filled_future=sum(fallback_used&future_reading),remaining_missing=sum(!is.finite(temperature)),
  median_fill_distance_hours=if(any(fallback_used))median(distance_hours[fallback_used]) else NA_real_,
  maximum_fill_distance_hours=if(any(fallback_used))max(distance_hours[fallback_used]) else NA_real_,
  fills_over_24h=sum(fallback_used&distance_hours>24,na.rm=TRUE)),by=.(site,period)]
 fwrite(summary,file.path(aggregate_dir,"temperature_nearest_fill_audit.csv"),na="")
 write_parquet(a,file.path(sensitive_dir,"temperature_nearest_fill_rows.parquet"))
 a[,.(hospital_block_id,m,temperature)]
}
