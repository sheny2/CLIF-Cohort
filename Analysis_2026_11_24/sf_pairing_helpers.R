# Strict chronological S/F pairing; no future/device backfill.
SF_DEFINITION<-"Prior same-device-episode FiO2 <=4h at SpO2 time; latest prior SpO2 <=4h at decision; categorical S/F"
SF_LEVELS<-c("Undefined",">300","250-300","200-250","120-200","<=120")
SF_DUMMIES<-c("sf_undefined","sf_250_300","sf_200_250","sf_120_200","sf_le120")
sf_categories<-function(x) fcase(!is.finite(x),"Undefined",x>300,">300",x>250,"250-300",
 x>200,"200-250",x>120,"120-200",default="<=120")
build_sf_pairs<-function(sp,resp) {
 sp<-copy(as.data.table(sp));resp<-copy(as.data.table(resp))
 sp[,`:=`(hospital_block_id=as.character(hospital_block_id),spo2_seconds=as.numeric(recorded_dttm),
  spo2=screen_measurements(vital_value,"spo2"))]
 sp<-sp[is.finite(spo2)&is.finite(spo2_seconds)]
 resp[,`:=`(hospital_block_id=as.character(hospital_block_id),seconds=as.numeric(recorded_dttm),
  device=gsub("[^a-z0-9]+","_",tolower(trimws(as.character(device_category)))),fi=screen_measurements(fio2_set,"fio2"))]
 resp[is.na(device)|device=="",device:=NA_character_]
 resp<-resp[is.finite(seconds)]
 resp[,priority:=match(map_support(device),SUPPORT_LEVELS)];resp[is.na(priority),priority:=0L]
 setorder(resp,hospital_block_id,seconds,priority)
 # Explicit device wins over blank; most intensive device wins timestamp ties.
 r<-resp[,{
  chosen<-if(any(!is.na(device)))tail(device[!is.na(device)],1) else NA_character_
  values<-fi[if(is.na(chosen))is.na(device) else !is.na(device)&device==chosen]
  list(device=chosen,fi=if(any(is.finite(values)))tail(values[is.finite(values)],1) else NA_real_)
 },by=.(hospital_block_id,seconds)]
 if(!nrow(r)) {
  sp[,`:=`(paired_fio2=NA_real_,fio2_age_hours=NA_real_,paired_device=NA_character_,episode=NA_integer_,measured_fio2=NA_real_,sf_reason="No prior device/FiO2",sf_ratio=NA_real_)]
  attr(sp,"device_timeline")<-data.table(hospital_block_id=character(),seconds=numeric(),episode=integer())
  return(sp)
 }
 setorder(r,hospital_block_id,seconds)
 r[,device:= {ii<-cummax(ifelse(!is.na(device),seq_len(.N),0L));ans<-rep(NA_character_,.N);ans[ii>0]<-device[ii[ii>0]];ans},by=hospital_block_id]
 r[,episode:=rleid(device),by=hospital_block_id]
 r[,fio2_seconds:=fifelse(is.finite(fi),seconds,NA_real_)]
 r[,`:=`(paired_fio2=nafill(fi,type="locf"),measured_fio2=nafill(fio2_seconds,type="locf")),by=.(hospital_block_id,episode)]
 setkey(r,hospital_block_id,seconds)
 a<-r[sp,on=.(hospital_block_id,seconds<=spo2_seconds),mult="last",
  .(hospital_block_id=i.hospital_block_id,spo2_seconds=i.spo2_seconds,spo2=i.spo2,
    paired_device=x.device,episode=x.episode,paired_fio2=x.paired_fio2,measured_fio2=x.measured_fio2)]
 a[,fio2_age_hours:=(spo2_seconds-measured_fio2)/3600]
 a[,sf_reason:=fcase(is.na(paired_device),"No prior device",!is.finite(paired_fio2),"No FiO2 in current device episode",
  !is.finite(fio2_age_hours)|fio2_age_hours<0|fio2_age_hours>4,"FiO2 older than 4h",spo2>97,"SpO2 >97%",default="Eligible")]
 a[,sf_ratio:=fifelse(sf_reason=="Eligible",screen_measurements(spo2/paired_fio2,"sf_ratio"),NA_real_)]
 attr(a,"device_timeline")<-r[,.(hospital_block_id,seconds,episode)]
 a
}
sf_pairs_to_grid<-function(pairs,grid) {
 p<-copy(pairs);setorder(p,hospital_block_id,spo2_seconds)
 p[,measured_spo2:=spo2_seconds];setkey(p,hospital_block_id,spo2_seconds)
 g<-copy(grid);g[,decision_seconds:=as.numeric(decision_dttm)]
 a<-p[g,on=.(hospital_block_id,spo2_seconds<=decision_seconds),mult="last",
  .(hospital_block_id=i.hospital_block_id,m=i.m,spo2=x.spo2,fio2=x.paired_fio2,sf_ratio=x.sf_ratio,
    sf_reason=x.sf_reason,sf_episode=x.episode,paired_device=x.paired_device,
    sf_spo2_seconds=x.measured_spo2,sf_fio2_seconds=x.measured_fio2,sf_spo2_age_hours=(i.decision_seconds-x.measured_spo2)/3600,
    sf_fio2_age_at_pair_hours=x.fio2_age_hours)]
 timeline<-attr(pairs,"device_timeline")
 if(!is.null(timeline)&&nrow(timeline)) {
  setkey(timeline,hospital_block_id,seconds)
  now<-timeline[g,on=.(hospital_block_id,seconds<=decision_seconds),mult="last",.(episode=x.episode)]
  a[!is.na(sf_episode)&!is.na(now$episode)&sf_episode!=now$episode,
    `:=`(sf_ratio=NA_real_,sf_reason="Device changed since SpO2")]
 }
 # Four-hour recency for the SpO2 observation itself; no future or best-value selection.
 a[!is.finite(sf_spo2_age_hours)|sf_spo2_age_hours<0|sf_spo2_age_hours>4,
   `:=`(sf_ratio=NA_real_,sf_reason="No usable SpO2 in preceding 4h")]
 a[,sf_category:=sf_categories(sf_ratio)]
 a
}
