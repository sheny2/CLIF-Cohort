suppressPackageStartupMessages(library(data.table))
proposal_limits <- list(temperature=c(25,45),wbc=c(.1,200),glucose=c(20,1000),
 albumin=c(.5,7),ph=c(6.5,8),paco2=c(5,200),sf_ratio=c(20,700),support_score=c(0,4),
 nonresp_sofa_proxy=c(0,17),spo2=c(50,100),fio2=c(.21,1),map_value=c(20,200),
 platelet=c(0,2000),bilirubin=c(0,100),creatinine=c(.1,30),gcs=c(3,15))
proposal_normalize <- function(x) gsub('[^a-z0-9]','',tolower(trimws(as.character(x))))
proposal_location <- function(x) {
 x <- proposal_normalize(x)
 fcase(x %in% c('ed','emergency','emergencydepartment'),'ED',
       x %in% c('icu','intensivecare','intensivecareunit'),'ICU',
       x %in% c('imc','intermediate','intermediatecare','intermediatecareunit','stepdown','stepdownunit'),'IMC',
       x %in% c('ward','floor','medsurg','acute','acutecare','acutecareunit','generalward'),'Ward',
       default='Unknown')
}
proposal_latest <- function(grid,measurement,time_col,value_col,output_col,baseline=FALSE) {
 g <- copy(as.data.table(grid))
 x <- copy(as.data.table(measurement))
 x <- x[,.(hospital_block_id=as.character(hospital_block_id),
            source_time=as.POSIXct(get(time_col),tz='UTC'),
            value=suppressWarnings(as.numeric(get(value_col))))]
 if(output_col=='fio2') x[value>1 & value<=100,value:=value/100]
 if(output_col %in% names(proposal_limits)) {
   bounds <- proposal_limits[[output_col]]
   x[value<bounds[1]|value>bounds[2],value:=NA_real_]
 }
 x <- x[!is.na(source_time)&is.finite(value)]
 bounds <- unique(g[,.(hospital_block_id,admission_dttm)])
 x <- merge(x,bounds,by='hospital_block_id')
 x <- x[source_time>=admission_dttm]
 # Baseline means nearest valid observation AT OR BEFORE time zero; no backfill.
 q <- if(baseline) unique(g[,.(hospital_block_id,index_dttm)]) else
   g[,.(hospital_block_id,m,index_dttm=decision_dttm)]
 setkey(x,hospital_block_id,source_time)
 ans <- x[q,on=.(hospital_block_id,source_time<=index_dttm),mult='last',
          .(hospital_block_id=i.hospital_block_id,query_time=i.index_dttm,value=x.value)]
 if(baseline) ans <- merge(g[,.(hospital_block_id,m)],ans[,.(hospital_block_id,value)],by='hospital_block_id',all.x=TRUE)
 else ans[,m:=q$m]
 setnames(ans,'value',output_col)
 ans[,c('hospital_block_id','m',output_col),with=FALSE]
}
proposal_sf_series <- function(spo2,respiratory,max_pair_hours=4) {
 s <- copy(as.data.table(spo2))[,.(hospital_block_id,spo2_time=as.POSIXct(recorded_dttm,tz='UTC'),sp=as.numeric(vital_value))]
 f <- copy(as.data.table(respiratory))[,.(hospital_block_id,ftime=as.POSIXct(recorded_dttm,tz='UTC'),fio=as.numeric(fio2_set))]
 f[fio>1 & fio<=100,fio:=fio/100]
 f <- f[is.finite(fio)&fio>=.21&fio<=1&!is.na(ftime)]
 s <- s[is.finite(sp)&sp>=50&sp<=100&!is.na(spo2_time)]
 f[,source_seconds:=as.numeric(ftime)]
 setkey(f,hospital_block_id,ftime)
 paired <- f[s,on=.(hospital_block_id,ftime<=spo2_time),mult='last',
  .(hospital_block_id=i.hospital_block_id,recorded_dttm=i.spo2_time,
    spo2=i.sp,fio2=x.fio,age_seconds=as.numeric(i.spo2_time)-x.source_seconds)]
 paired <- paired[is.finite(age_seconds)&age_seconds>=0&age_seconds<=max_pair_hours*3600]
 paired[,sf_ratio:=spo2/fio2]
 paired[,.(hospital_block_id,recorded_dttm,sf_ratio)]
}
proposal_process <- function(long,grace_intervals,max_missing=.70) {
 stopifnot(is.finite(max_missing),max_missing>=0,max_missing<1)
 x <- copy(as.data.table(long))
 features <- c('glucose','albumin','temperature','wbc','ph','paco2','sf_ratio','support_score','nonresp_sofa_proxy')
 baseline <- c('glucose','albumin')
 # Before initiation and within grace: the risk set used for treatment decisions.
 decision <- x$m<=grace_intervals & x$m<=x$S
 qa <- rbindlist(lapply(features,function(v) {
  limits <- proposal_limits[[v]]
  value <- as.numeric(x[[v]])
  value[!is.finite(value)|value<limits[1]|value>limits[2]] <- NA_real_
  ref <- if(v %in% baseline) value[!duplicated(x$hospital_block_id)] else value[decision]
  observed <- is.finite(ref);miss <- if(length(ref)) mean(!observed) else 1
  eligible <- any(observed) && miss<=max_missing
  fill <- if(any(observed)) median(ref[observed]) else NA_real_
  winsor <- if(any(observed)) as.numeric(quantile(ref[observed],c(.005,.995))) else c(NA_real_,NA_real_)
  if(v=="support_score") winsor <- c(0,4)
  missing <- as.integer(!is.finite(value))
  # Residual missingness (no prior measurement) cannot be filled by LOCF.
  # Impute only retained features; indicator remains tied to original missingness.
  value[is.finite(value)] <- pmin(pmax(value[is.finite(value)],winsor[1]),winsor[2])
  value[missing==1L] <- fill
  if(v %in% c('wbc','glucose')) value <- log1p(value)
  center <- if(any(is.finite(value))) mean(value[is.finite(value)]) else 0
  scale <- sd(value,na.rm=TRUE);if(!is.finite(scale)||scale<=0) scale<-1
  if(!v %in% c('support_score','nonresp_sofa_proxy')) value <- (value-center)/scale
  x[,(paste0(v,'_model')):=value]
  x[,(paste0(v,'_missing')):=missing]
  varies <- length(unique(value[is.finite(value)]))>1
  data.table(variable=v,role=if(v %in% baseline)'baseline' else 'time_varying',
   assessment_population=if(v %in% baseline)'one baseline per patient' else 'untreated decision intervals within grace',
   missing_fraction_after_locf=miss,maximum_missing_fraction=max_missing,
   usable=eligible&&varies,
   reason=if(!any(observed))'no observed values' else if(miss>max_missing)'exceeds missingness threshold' else if(!varies)'constant after processing' else 'included',
   imputation_value=fill,winsor_lower=winsor[1],winsor_upper=winsor[2],model_variable=paste0(v,'_model'))
 }))
 base <- qa[role=='baseline' & usable==TRUE,variable]
 tv <- qa[role=='time_varying' & usable==TRUE,variable]
 support_terms <- character()
 if("support_score" %in% tv) {
   # Encode support type categorically; room air is the reference.
   for(k in 1:4) {
     nm <- paste0(c("support_low_flow","support_hfno","support_nippv","support_imv")[k],"_model")
     x[,(nm):=as.integer(support_score==k)]
     x[is.na(get(nm)),(nm):=0L]
     support_terms <- c(support_terms,nm)
   }
   support_terms <- c(support_terms,"support_score_missing")
   qa[variable=="support_score",model_variable:=paste(support_terms,collapse=";")]
   tv <- setdiff(tv,"support_score")
 }
 make_terms <- function(f) if(length(f)) as.vector(rbind(paste0(f,'_model'),paste0(f,'_missing'))) else character()
 list(data=x,qa=qa,baseline=make_terms(base),tv=c(make_terms(tv),support_terms))
}
