#!/usr/bin/env Rscript
# Raw measurement quality EDA for the selected script-01 cohort, before
# imputation, winsorization or plausibility exclusions. No analysis is changed.
# Rscript 0501_measurement_quality_eda.R CLIF_DATA_DIR [COHORT_DIR] [NEW_OUTPUT_DIR]
# Defaults: cohort from COHORT_OUTPUT_DIR; -24 to +672 hours around time zero.
# Optional EDA_PRE_HOURS=24, EDA_FOLLOWUP_HOURS=672.
# aggregate/: unsuppressed distributions/ranges/units/time trends and plots.
# sensitive/: flagged and extreme individual records; keep on the server.
# Flags describe review candidates, NOT verified errors or clinical normal ranges.
suppressPackageStartupMessages({library(arrow);library(data.table);library(dplyr);library(ggplot2)})
a<-commandArgs(TRUE)
sf<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
root<-dirname(normalizePath(gsub("~+~"," ",sf,fixed=TRUE)))
arg<-function(i,default)if(length(a)>=i && nzchar(a[i]))a[i] else default
data_dir<-arg(1,Sys.getenv("CLIF_DATA_DIR","/projects/emulate_icu/rclif_2_1"))
cohort_dir<-arg(2,Sys.getenv("COHORT_OUTPUT_DIR",file.path(root,"pipeline_output/sensitive/cohort")))
out<-arg(3,file.path(root,"pipeline_output",paste0("0501_measurement_eda_",format(Sys.time(),"%Y%m%d_%H%M%S"))))
pre<-as.numeric(Sys.getenv("EDA_PRE_HOURS","24"));horizon<-as.numeric(Sys.getenv("EDA_FOLLOWUP_HOURS","672"))
if(!is.finite(pre)||pre<0||!is.finite(horizon)||horizon<=0)stop("Invalid EDA time window")
if(normalizePath(out,mustWork=FALSE)%in%vapply(c(data_dir,cohort_dir),normalizePath,character(1),mustWork=FALSE))stop("Use a separate output directory")
agg<-file.path(out,"aggregate");sensitive<-file.path(out,"sensitive");plots<-file.path(agg,"plots")
dir.create(plots,recursive=TRUE,showWarnings=FALSE);dir.create(sensitive,recursive=TRUE,showWarnings=FALSE)
write_agg<-function(x,f){
 if(any(c("patient_id","hospitalization_id","hospital_block_id","recorded_dttm")%in%names(x)))stop("Identifier/timestamp in aggregate output")
 fwrite(x,file.path(agg,f),na="")
}
norm<-function(x)gsub("[^a-z0-9]+","_",tolower(trimws(as.character(x))))
cohort<-as.data.table(read_parquet(file.path(cohort_dir,"adult_pneumonia_cohort.parquet")))
cw<-as.data.table(read_parquet(file.path(cohort_dir,"adult_pneumonia_hospitalization_crosswalk.parquet")))
required<-c("hospital_block_id","patient_id","source_site_id","respiratory_failure_start","discharge_dttm")
if(!all(required%in%names(cohort)))stop("Rerun script 01: required block fields missing")
if(anyDuplicated(cohort$hospital_block_id)||anyDuplicated(cohort$patient_id))stop("Expected one block per patient")
cohort[,`:=`(hospital_block_id=as.character(hospital_block_id),site=as.character(source_site_id),
 time_zero=as.POSIXct(respiratory_failure_start,tz="UTC"),discharge=as.POSIXct(discharge_dttm,tz="UTC"))]
cw[,`:=`(hospitalization_id=as.character(hospitalization_id),hospital_block_id=as.character(hospital_block_id))]
cw<-unique(cw[hospital_block_id%in%cohort$hospital_block_id,.(hospitalization_id,hospital_block_id)])
if(anyDuplicated(cw$hospitalization_id))stop("Encounter maps to multiple blocks")
links<-merge(cw,cohort[,.(hospital_block_id,site,time_zero,discharge)],by="hospital_block_id")
ids<-unique(links$hospitalization_id)

# These screening bounds reproduce the main pipeline where available; MAP,
# GCS and auxiliary venous gases use broad review ranges. They are not exclusion rules.
rules<-data.table(variable=c("temperature","spo2","map","wbc","glucose","albumin","ph_arterial","pco2_arterial",
 "ph_venous","pco2_venous","platelet","bilirubin_total","creatinine","gcs","fio2_raw","support_score"),
 lower=c(25,50,20,.1,20,.5,6.5,5,6.5,5,0,0,0,3,.21,0),
 upper=c(45,100,200,200,1000,7,8,200,8,200,2000,100,30,15,1,3),
 assumed_convention=c("C","percent","mmHg","10^9/L or 10^3/uL","mg/dL","g/dL","unitless","mmHg",
 "unitless","mmHg","10^9/L or 10^3/uL","mg/dL","mg/dL","3-15","fraction (percent form separately flagged)","ordinal 0-3"))
write_agg(rules,"review_bounds_and_assumed_scales.csv")
aliases<-list(
 temperature=c("temp_c","temperature","temperature_c"),spo2=c("spo2"),map=c("map","mean_arterial_pressure"),
 wbc=c("wbc","white_blood_cell_count","white_blood_cells"),glucose=c("glucose_serum","serum_glucose","blood_glucose","glucose"),
 albumin=c("albumin","serum_albumin","albumin_serum"),ph_arterial=c("ph_arterial","arterial_ph","ph_blood_arterial"),
 pco2_arterial=c("pco2_arterial","arterial_pco2","paco2"),ph_venous=c("ph_venous","venous_ph"),pco2_venous=c("pco2_venous","venous_pco2"),
 platelet=c("platelet_count","platelets","platelet"),bilirubin_total=c("bilirubin_total","total_bilirubin","bilirubin"),
 creatinine=c("creatinine","serum_creatinine","creatinine_serum"),gcs=c("gcs_total"))
sources<-list(
 list(file="clif_vitals.parquet",category="vital_category",time="recorded_dttm",value="vital_value",variables=c("temperature","spo2","map")),
 list(file="clif_labs.parquet",category="lab_category",time="lab_collect_dttm",value="lab_value_numeric",variables=c("wbc","glucose","albumin","ph_arterial","pco2_arterial","ph_venous","pco2_venous","platelet","bilirubin_total","creatinine")),
 list(file="clif_patient_assessments.parquet",category="assessment_category",time="recorded_dttm",value="numerical_value",variables="gcs"))
q<-function(x,p){x<-x[is.finite(x)];if(length(x))as.numeric(quantile(x,p,names=FALSE)) else rep(NA_real_,length(p))}
stats<-function(x){
 finite<-is.finite(x);z<-x[finite];qs<-q(z,c(.001,.01,.25,.5,.75,.99,.999))
 list(records=length(x),finite_records=length(z),nonfinite_records=sum(!finite),
 minimum=if(length(z))min(z) else NA_real_,p001=qs[1],p01=qs[2],q1=qs[3],median=qs[4],q3=qs[5],p99=qs[6],p999=qs[7],
 maximum=if(length(z))max(z) else NA_real_,mean=if(length(z))mean(z) else NA_real_,sd=if(length(z)>1)sd(z) else NA_real_)
}
summaries<-time6<-time_day<-calendar<-flags<-bins<-mappings<-list();pair_inputs<-list()
theme_set(theme_minimal(base_size=12))
attach_records<-function(x){
 x[,hospitalization_id:=as.character(hospitalization_id)]
 x<-merge(x,links,by="hospitalization_id",allow.cartesian=FALSE)
 x[,recorded_dttm:=as.POSIXct(recorded_dttm,tz="UTC")]
 x[,hours:=as.numeric(difftime(recorded_dttm,time_zero,units="hours"))]
 x<-x[is.finite(hours)&hours>=-pre&hours<=horizon]
 x[,`:=`(source_unit=fifelse(is.na(source_unit)|!nzchar(trimws(source_unit)),"[not supplied]",trimws(source_unit)),
   relative_6h=floor(hours/6)*6,relative_day=floor(hours/24),calendar_month=format(recorded_dttm,"%Y-%m"),
   after_block_discharge=!is.na(discharge)&recorded_dttm>discharge)]
 x
}
review_variable<-function(x,v){
 if(!nrow(x))return(invisible(NULL))
 message("Summarizing raw measurements: ",v," (",nrow(x)," records)")
 x<-copy(x);x[,variable:=v]
 x[,duplicate_record:=duplicated(.SD),.SDcols=c("hospitalization_id","recorded_dttm","category","value","source_unit")]
 b<-rules[variable==v];lo<-b$lower[1];hi<-b$upper[1]
 x[,outside_review_range:=is.finite(value)&(value<lo|value>hi)]
 x[,scale_candidate:=fcase(
   v=="temperature" & value>=80 & value<=115,"possible Fahrenheit numeric scale",
   v=="spo2" & value>0 & value<=1,"possible proportion rather than percent",
   v=="fio2_raw" & value>=21 & value<=100,"possible percent rather than fraction",
   v=="glucose" & value>0 & value<20,"low value or possible mmol/L scale; verify unit",
   v=="wbc" & value>200 & value<=200000,"high value or possible cells/uL scale; verify unit",
   v=="albumin" & value>7 & value<=80,"high value or possible g/L scale; verify unit",
   default="")]
 keys<-c("variable","site","source_unit")
 summaries[[v]]<<-x[,c(stats(value),list(patients=uniqueN(hospital_block_id))),by=keys]
 flags[[v]]<<-x[,.(records=.N,nonfinite=sum(!is.finite(value)),outside_review_range=sum(outside_review_range),
   outside_review_percent=100*mean(outside_review_range),scale_review_records=sum(nzchar(scale_candidate)),
   exact_duplicate_records=sum(duplicate_record),
   spo2_above_97=if(v=="spo2")sum(is.finite(value)&value>97) else NA_integer_,
   after_discharge_records=sum(after_block_discharge)),by=keys]
 time6[[v]]<<-x[,stats(value),by=c(keys,"relative_6h")]
 time_day[[v]]<<-x[,stats(value),by=c(keys,"relative_day")]
 calendar[[v]]<<-x[,stats(value),by=c(keys,"calendar_month")]
 # Never delete suspect values. Save row-level examples separately on-server.
 bad<-x[!is.finite(value)|outside_review_range|nzchar(scale_candidate)|after_block_discharge]
 write_parquet(bad,file.path(sensitive,paste0(v,"_review_records.parquet")))
 extreme<-x[is.finite(value)][order(value),rbind(head(.SD,20),tail(.SD,20)),by=.(site,source_unit)]
 write_parquet(extreme,file.path(sensitive,paste0(v,"_extreme_examples.parquet")))
 finite<-x[is.finite(value)]
 if(!nrow(finite))return(invisible(NULL))
 # Aggregate bins use ALL records, not a sample. Signed-log keeps zeros and negatives visible.
 finite[,signed_log:=sign(value)*log10(1+abs(value))]
 bound<-range(finite$signed_log);width<-max(diff(bound)/100,.001)
 h<-finite[,.(records=.N),by=.(site,source_unit,bin=floor((signed_log-bound[1])/width))]
 h[,signed_log_mid:=bound[1]+(bin+.5)*width];h[,variable:=v]
 bins[[v]]<<-h
 panel_count<-uniqueN(paste(h$site,h$source_unit));height<-max(5,ceiling(panel_count/3)*2.4)
 p<-ggplot(h,aes(signed_log_mid,records))+geom_col(width=width*.95,fill="#2878A5")+
  scale_y_continuous(trans="log1p")+
  facet_wrap(~site+source_unit,ncol=3,scales="free_y")+
  labs(title=paste(v,": all raw values by site and reported unit"),
   subtitle="No clipping or imputation; all observations contribute",
   x="sign(value) × log10(1 + |value|)",y="Measurement records (log1p axis)")
 ggsave(file.path(plots,paste0(v,"_full_distribution.png")),p,width=12,height=height,dpi=160,limitsize=FALSE)
 # A separate central-range view is purely a display zoom; extreme values remain in summaries/full plot.
 raw_range<-q(finite$value,c(.01,.99));w<-max(diff(raw_range)/80,1e-6)
 z<-finite[value>=raw_range[1]&value<=raw_range[2],.(records=.N),by=.(site,source_unit,bin=floor((value-raw_range[1])/w))]
 z[,mid:=raw_range[1]+(bin+.5)*w]
 p<-ggplot(z,aes(mid,records))+geom_col(width=w*.95,fill="#4B9C78")+
  facet_wrap(~site+source_unit,ncol=3,scales="free_y")+
  labs(title=paste(v,": central raw-value distribution"),subtitle="Display limited to pooled 1st–99th percentiles; full ranges are preserved separately",x="Raw numeric value",y="Measurement records")
 ggsave(file.path(plots,paste0(v,"_central_distribution.png")),p,width=12,height=height,dpi=160,limitsize=FALSE)
 trend<-time_day[[v]]
 p<-ggplot(trend,aes(relative_day,median))+geom_ribbon(aes(ymin=q1,ymax=q3),fill="#2878A5",alpha=.25)+
  geom_line(color="#2878A5")+geom_line(aes(y=minimum),color="#BC554A",linetype=3)+
  geom_line(aes(y=maximum),color="#BC554A",linetype=3)+facet_wrap(~site+source_unit,ncol=3,scales="free_y")+
  labs(title=paste(v,": daily raw ranges and distributions"),subtitle="Blue: median and IQR; red dotted lines: observed minimum/maximum",x="Days relative to time zero",y="Raw numeric value")
 ggsave(file.path(plots,paste0(v,"_daily_min_max.png")),p,width=12,height=height,dpi=160,limitsize=FALSE)
 monthly<-copy(calendar[[v]])
 monthly[,month_start:=as.Date(paste0(calendar_month,"-01"))]
 p<-ggplot(monthly,aes(month_start,median))+geom_ribbon(aes(ymin=p01,ymax=p99),fill="#2878A5",alpha=.2)+
  geom_line()+geom_point(size=.6)+facet_wrap(~site+source_unit,ncol=3,scales="free_y")+
  labs(title=paste(v,": recording scale over calendar time"),
   subtitle="Monthly median and 1st–99th percentiles; exact monthly extrema are in the CSV",x="Calendar month",y="Raw numeric value")
 ggsave(file.path(plots,paste0(v,"_calendar_month.png")),p,width=12,height=height,dpi=160,limitsize=FALSE)
 early<-time6[[v]][relative_6h>=0&relative_6h<96]
 if(nrow(early)){
  p<-ggplot(early,aes(relative_6h,median))+geom_ribbon(aes(ymin=p01,ymax=p99),fill="#2878A5",alpha=.18)+
   geom_line()+geom_point(size=.7)+facet_wrap(~site+source_unit,ncol=3,scales="free_y")+
   geom_vline(xintercept=c(24,96),linetype=3)+labs(title=paste(v,": first 96 hours"),
    subtitle="Six-hour measurement bins; median and 1st–99th percentiles",x="Bin start, hours",y="Raw numeric value")
  ggsave(file.path(plots,paste0(v,"_first96h.png")),p,width=12,height=height,dpi=160,limitsize=FALSE)
 }
}
for(src in sources){
 path<-file.path(data_dir,src$file)
 if(!file.exists(path)){warning("Missing optional measurement source: ",src$file);next}
 ds<-open_dataset(path);schema<-ds$schema$names
 required<-c("hospitalization_id",src$category,src$time,src$value)
 if(!all(required%in%schema))stop("Missing required measurement columns in ",src$file)
 catcol<-src$category
 available<-ds |> select(all_of(catcol)) |> distinct() |> collect()
 available<-as.character(available[[catcol]])
 map<-lapply(aliases[src$variables],function(z)available[norm(available)%in%z])
 mappings[[src$file]]<-rbindlist(lapply(names(map),function(v)data.table(source=src$file,variable=v,
   matched_categories=paste(map[[v]],collapse=";"),mapping_available=length(map[[v]])>0)))
 wanted<-unique(unlist(map,use.names=FALSE));if(!length(wanted))next
 cols<-unique(c(required,intersect(c("reference_unit","lab_specimen_category"),schema)))
 raw<-ds |> filter(.data$hospitalization_id%in%ids,.data[[catcol]]%in%wanted) |>
  select(all_of(cols)) |> collect() |> as.data.table()
 setnames(raw,c(src$category,src$time,src$value),c("category","recorded_dttm","value"))
 raw[,value:=as.numeric(value)]
 if(!"reference_unit"%in%names(raw))raw[,reference_unit:=NA_character_]
 raw[,source_unit:=as.character(reference_unit)]
 raw<-attach_records(raw)
 for(v in names(map)){
  x<-raw[category%in%map[[v]]]
  review_variable(x,v)
  if(v=="spo2")pair_inputs$spo2<-copy(x)
 }
 rm(raw);invisible(gc())
}

path<-file.path(data_dir,"clif_respiratory_support.parquet")
if(file.exists(path)){
 raw<-open_dataset(path) |> filter(.data$hospitalization_id%in%ids) |>
  select(hospitalization_id,recorded_dttm,device_category,fio2_set) |> collect() |> as.data.table()
 raw[,`:=`(category=as.character(device_category),source_unit=NA_character_)]
 raw<-attach_records(raw);raw[,device_clean:=norm(device_category)]
 raw[,support_score:=fcase(device_clean=="room_air",0,
  device_clean%in%c("nasal_cannula","face_mask","trach_collar"),1,
  device_clean%in%c("high_flow_nc","nippv","cpap"),2,device_clean=="imv",3,default=NA_real_)]
 write_agg(raw[,.(records=.N,patients=uniqueN(hospital_block_id)),by=.(site,device_clean,support_score)],"respiratory_device_mapping_counts.csv")
 x<-copy(raw);x[,value:=as.numeric(fio2_set)];review_variable(x,"fio2_raw");pair_inputs$fio2<-copy(x)
 x<-copy(raw);x[,value:=support_score];review_variable(x,"support_score")
 rm(raw,x);invisible(gc())
}

# Audit the S/F construction on a six-hour grid using the main 12-hour last-value
# rule for BOTH inputs. Keep raw and normalized FiO2; do not convert SpO2 proportions.
if(all(c("spo2","fio2")%in%names(pair_inputs))){
 grid<-cohort[,.(m=seq_len(ceiling(horizon/6))),by=.(hospital_block_id,site,time_zero,discharge)]
 grid[,`:=`(row_id=.I,decision_seconds=as.numeric(time_zero)+(m-1)*6*3600)]
 last_input<-function(x){
  z<-x[is.finite(value),.(hospital_block_id,seconds=as.numeric(recorded_dttm),measured_seconds=as.numeric(recorded_dttm),val=value)]
  if(!nrow(z))return(data.table(row_id=grid$row_id,value=NA_real_,age_hours=NA_real_))
  setkey(z,hospital_block_id,seconds)
  ans<-z[grid,on=.(hospital_block_id,seconds<=decision_seconds),mult="last",
   .(row_id=i.row_id,value=x.val,age_hours=(i.decision_seconds-x.measured_seconds)/3600)]
  ans[!is.finite(age_hours)|age_hours<0|age_hours>12,value:=NA_real_]
  ans
 }
 sp<-last_input(pair_inputs$spo2);fi<-last_input(pair_inputs$fio2)
 grid[,`:=`(spo2_raw=sp$value[match(row_id,sp$row_id)],fio2_raw=fi$value[match(row_id,fi$row_id)])]
 grid[,fio2_fraction:=fifelse(fio2_raw>1.5,fio2_raw/100,fio2_raw)]
 grid[fio2_fraction<.21|fio2_fraction>1,fio2_fraction:=NA_real_]
 grid[,sf_ratio:=fifelse(spo2_raw>=50&spo2_raw<=100,spo2_raw/fio2_fraction,NA_real_)]
 grid[,`:=`(sf_eligible_le97=is.finite(sf_ratio)&spo2_raw<=97,
   relative_6h=(m-1)*6,after_discharge=decision_seconds>as.numeric(discharge))]
 write_agg(grid[,c(stats(sf_ratio),list(paired_rows=sum(is.finite(sf_ratio)),
  paired_rows_spo2_gt97=sum(is.finite(sf_ratio)&spo2_raw>97),eligible_le97_rows=sum(sf_eligible_le97))),
  by=.(site,relative_6h,after_discharge)],"sf_grid_pairing_audit.csv")
 write_parquet(grid,file.path(sensitive,"sf_grid_pairing_rows.parquet"))
}
for(item in list(c("summaries","raw_distributions_by_site_and_unit.csv"),c("time6","raw_six_hour_ranges.csv"),
 c("time_day","raw_daily_ranges.csv"),c("calendar","raw_calendar_month_ranges.csv"),c("flags","review_flag_counts.csv"),
 c("bins","full_distribution_plot_bins.csv"),c("mappings","category_mapping_audit.csv"))){
 obj<-get(item[1]);if(length(obj))write_agg(rbindlist(obj,fill=TRUE),item[2])
}
if(length(summaries)){
 u<-rbindlist(summaries)
 write_agg(u[,.(reported_unit_groups=.N,nonblank_unit_groups=sum(source_unit!="[not supplied]"),
   reported_units=paste(source_unit,collapse=" | ")),by=.(variable,site)],"reported_unit_inventory.csv")
}
write_agg(data.table(parameter=c("cohort_patients","pre_hours","followup_hours","data_dir","cohort_dir"),
 value=as.character(c(nrow(cohort),pre,horizon,data_dir,cohort_dir))),"eda_settings.csv")
writeLines(c(
 "Raw measurement EDA for the selected script-01 cohort; no treatment assignment or 07/09 outputs needed.",
 "Start with raw_distributions_by_site_and_unit.csv, review_flag_counts.csv, reported_unit_inventory.csv and plots/.",
 "Six-hour and daily summaries are relative to time zero; calendar-month summaries detect changes in recording scale over years.",
 "Min/max and quantiles include all finite recorded values, including suspect values. Raw records are not imputed, clipped or deleted for being implausible.",
 "Distribution counts are measurement-weighted, not patient-weighted; frequently monitored patients contribute more observations.",
 "Unit fields are reported as delivered. No lab-unit conversions are made. Vitals generally lack an explicit unit field; their scales are assumptions to investigate.",
 "Multiple unit strings may represent spelling variants rather than incompatible scales. Even the same unit label can conceal a mapping error.",
 "Review bounds and scale candidates are screening heuristics, not verified errors or exclusion criteria. SpO2 >97% is physiologically valid; its separate flag concerns use in S/F.",
 "Full-range plots use signed log10 to show extreme/negative/zero values. Central plots are display-only 1st–99th percentile zooms.",
 "Each site/unit has its own distribution panel. Daily raw extrema may compress central patterns; consult the central distributions and first-96h percentile plots too.",
 "S/F pairing uses independent last values within 12h at six-hour decision times; it is not necessarily simultaneous. FiO2 percent normalization follows the main pipeline; original SpO2/FIO2 values remain saved.",
 "S/F grid is a measurement audit, not a death-censored analysis risk set. Its post-discharge rows are explicitly labelled and must not be confused with inpatient measurement gaps.",
 "The current main support score mapping is reproduced only for audit; unknown device categories remain unmapped.",
 "Arterial and venous gases are summarized separately. No source measurement or analysis file has been modified.",
 "Aggregate summaries and plots are unsuppressed. Individual review records/extreme examples/exact timestamps are only in sensitive/ and must remain on server."
),file.path(out,"START_HERE.txt"))
writeLines(capture.output(sessionInfo()),file.path(out,"SESSION_INFO.txt"))
writeLines("Complete",file.path(out,"COMPLETE.txt"))
message("Measurement quality EDA complete: ",normalizePath(out))
