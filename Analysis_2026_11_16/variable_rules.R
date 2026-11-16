# Single source of analytic screening and categorical respiratory support rules.
# Bounds are inclusive, intentionally broad, and are analytic exclusions rather
# than assertions that every excluded extreme is a verified source-data error.
measurement_rules <- function() data.table::data.table(
 variable=c("temperature","spo2","map_value","wbc","glucose","albumin",
            "platelet","bilirubin","creatinine","gcs","fio2","sf_ratio","nonresp_sofa_proxy"),
 lower=c(30,70,0,0,0,0,0,0,0,3,.21,70,0),
 upper=c(45,100,250,500,2000,15,2000,80,30,15,1,462,17),
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

# Fill remaining gaps AFTER bounded prior carry-forward. Future measurements are
# permitted only within the SAME selected block and before its primary event /
# discharge / administrative cutoff (source x is already filtered in preparation).
fill_adjacent_measurements <- function(answer,x,grid,variable,tau,aggregate_dir,sensitive_dir) {
 a<-merge(answer,grid[,.(hospital_block_id,m,site,discharge_dttm)],by=c("hospital_block_id","m"),sort=FALSE)
 a[,`:=`(missing_before=!is.finite(value),filled=FALSE,future=FALSE,distance_hours=NA_real_)]
 target<-a[missing_before==TRUE & decision_dttm<=discharge_dttm]
 z<-x[,.(hospital_block_id,seconds=measurement_seconds,measured=measurement_seconds,val=value)]
 setkey(z,hospital_block_id,seconds)
 if(nrow(target)&&nrow(z)) {
  target[,seconds:=as.numeric(decision_dttm)]
  prior<-z[target,on=.(hospital_block_id,seconds<=seconds),mult="last",
   .(hospital_block_id=i.hospital_block_id,m=i.m,pv=x.val,pt=x.measured)]
  nextval<-z[target,on=.(hospital_block_id,seconds>=seconds),mult="first",
   .(hospital_block_id=i.hospital_block_id,m=i.m,nv=x.val,nt=x.measured)]
  target<-merge(target,prior,by=c("hospital_block_id","m"))
  target<-merge(target,nextval,by=c("hospital_block_id","m"))
  target[,use_prior:=is.finite(pt)&(!is.finite(nt)|seconds-pt<=nt-seconds)]
  target[,`:=`(fill=fifelse(use_prior,pv,nv),measured=fifelse(use_prior,pt,nt))]
  target[,`:=`(gap=abs(measured-seconds)/3600,used_future=is.finite(measured)&measured>seconds)]
  a[target,on=.(hospital_block_id,m),`:=`(value=i.fill,filled=is.finite(i.fill),future=i.used_future,distance_hours=i.gap)]
 }
 a[,period:=fifelse(m<=tau,"during_grace","after_grace")]
 summary<-a[,.(rows=.N,missing_before=sum(missing_before),filled_adjacent=sum(filled),
  filled_future=sum(filled&future),remaining_missing=sum(!is.finite(value)),
  median_fill_distance_hours=if(any(filled))median(distance_hours[filled]) else NA_real_,
  maximum_fill_distance_hours=if(any(filled))max(distance_hours[filled]) else NA_real_,
  fills_over_24h=sum(filled&distance_hours>24,na.rm=TRUE)),by=.(site,period)]
 summary[,variable:=variable]
 fwrite(summary,file.path(aggregate_dir,paste0("adjacent_fill_audit_",variable,".csv")),na="")
 arrow::write_parquet(a[missing_before==TRUE],file.path(sensitive_dir,paste0("adjacent_fill_rows_",variable,".parquet")))
 a[,names(answer),with=FALSE]
}
