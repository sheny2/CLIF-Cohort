# Patient-weighted six-hour support summaries and cleaned, pre-imputation EDA.
write_covariate_eda <- function(long,cohort,out,tau,step_hours) {
 dir.create(file.path(out,"plots"),recursive=TRUE,showWarnings=FALSE)
 write <- function(x,f) fwrite(x,file.path(out,f),na="")
 # Baseline S/F contributes one observation per patient, not repeated rows.
 sf <- cohort[,.(site,baseline_sf)]
 write(sf[,.(patients=.N,missing=sum(!is.finite(baseline_sf)),minimum=min(baseline_sf),
   p01=as.numeric(quantile(baseline_sf,.01)),median=median(baseline_sf),
   p99=as.numeric(quantile(baseline_sf,.99)),maximum=max(baseline_sf)),by=site],"baseline_sf_distribution.csv")
 p <- ggplot(sf,aes(baseline_sf))+geom_histogram(bins=35,fill="#2878A5")+facet_wrap(~site,ncol=3)+
  labs(title="Baseline S/F in the restricted cohort",subtitle="One baseline value per patient; no baseline S/F imputation",
   x="Baseline SpO2 / FiO2",y="Patients")+theme_minimal()
 ggsave(file.path(out,"plots","baseline_sf_distribution.png"),p,width=12,height=7,dpi=160)
 # Full baseline cohort denominator: explicit death/discharge/admin censoring states.
 g <- cohort[,.(m=seq_len(112L)),by=.(hospital_block_id,site,index_dttm,
             discharge_dttm,event_hours,event_type,delta,censor_hours)]
 g[,hours:=(m-1)*step_hours]
 g <- merge(g,long[,.(hospital_block_id,m,support_category)],by=c("hospital_block_id","m"),all.x=TRUE)
 g[, state := fcase(delta==1L & is.finite(event_hours)&event_hours<=hours & event_type=="in_hospital_death","In-hospital death",
   delta==1L & is.finite(event_hours)&event_hours<=hours & event_type=="hospice_discharge","Hospice discharge",
   is.finite(censor_hours)&censor_hours<=hours,"Administratively censored",
   as.numeric(index_dttm)+hours*3600>as.numeric(discharge_dttm),"Discharged alive",
   !is.na(support_category),support_category,default="Unknown")]
 # Other alive discharge is a resolved non-event for the selected block.
 # Death outside that block does not change this state; hospice is separate.
 g[,state:=factor(state,levels=c(SUPPORT_LEVELS,"Unknown","Discharged alive","In-hospital death","Hospice discharge","Administratively censored"))]
 summarize_support <- function(z,group) {
  counts <- z[,.(patients=.N),by=c(group,"state")]
  denom <- z[,.(denominator=.N),by=group]
  a <- merge(counts,denom,by=group);a[,proportion:=patients/denominator];a
 }
 overall <- summarize_support(g,"hours");overall[,site:="All sites"]
 bysite <- summarize_support(g,c("site","hours"))
 all <- rbindlist(list(overall,bysite),use.names=TRUE)
 write(all,"respiratory_support_proportions_full_cohort.csv")
 inhospital <- g[!as.character(state)%chin%c("Discharged alive","In-hospital death","Hospice discharge","Administratively censored")]
 if(nrow(inhospital)) {
  b <- summarize_support(inhospital,c("site","hours"))
  pooled <- summarize_support(inhospital,"hours");pooled[,site:="All sites"]
  b <- rbindlist(list(pooled,b),use.names=TRUE)
  write(b,"respiratory_support_proportions_alive_in_hospital.csv")
 }
 # Complete combinations prevent ggplot from bridging an absent category.
 for (kind in c("full_cohort","alive_in_hospital")) {
  d <- if(kind=="full_cohort")all else if(exists("b"))b else NULL
  if(is.null(d)||!nrow(d))next
  complete <- CJ(site=unique(d$site),hours=sort(unique(d$hours)),state=levels(g$state))
  d[,state:=as.character(state)]
  complete <- merge(complete,d,by=c("site","hours","state"),all.x=TRUE)
  complete[is.na(proportion),proportion:=0]
  complete[,state:=factor(state,levels=levels(g$state))]
  for(window in c(96,672)) {
   p <- ggplot(complete[hours<window],aes(hours/24,proportion,fill=state))+
    geom_area(position="stack")+
    scale_fill_manual(values=c(room_air="#66C2A5",conventional_oxygen="#3288BD",hfno="#FEE08B",
      nippv_cpap="#FDAE61",imv="#D53E4F",Unknown="#999999",`Discharged alive`="#ABD9E9",
      `In-hospital death`="#333333",`Hospice discharge`="#AA4499",`Administratively censored`="#DDDDDD"),drop=FALSE,
      labels=c(room_air="Room air",conventional_oxygen="Conventional oxygen",hfno="HFNO",
       nippv_cpap="NIPPV/CPAP",imv="IMV",Unknown="Unknown",`Discharged alive`="Discharged alive",
       `In-hospital death`="In-hospital death",`Hospice discharge`="Hospice discharge",`Administratively censored`="Administratively censored"))+
    guides(fill=guide_legend(nrow=3,byrow=TRUE))+facet_wrap(~site,ncol=2)+
    geom_vline(xintercept=c(24,96)/24,linetype=3)+
    scale_y_continuous(labels=function(x)paste0(round(100*x),"%"))+coord_cartesian(ylim=c(0,1))+
    labs(title="Respiratory support over time",subtitle=if(kind=="full_cohort")
      "One state per original patient per time point; full baseline cohort denominator" else
      "One state per patient still alive and within the selected hospital block",
      x="Days since time zero",y="Patient proportion",fill="Support / observation state")+
    theme_minimal(base_size=11)+theme(legend.position="bottom")
   ggsave(file.path(out,"plots",paste0("respiratory_support_",kind,"_",window,"h.png")),p,width=12,height=10,dpi=160)
  }
 }
 # Consecutive 6-hour states, rather than transitions between every raw chart row.
 setorder(g,site,hospital_block_id,m)
 g[,from_state:=shift(as.character(state)),by=hospital_block_id]
 write(g[!is.na(from_state),.(transitions=.N),by=.(site,from_state,to_state=as.character(state))],
       "respiratory_support_six_hour_transitions.csv")
 numeric_vars <- c("temperature","wbc","glucose","albumin","nonresp_sofa_proxy")
 z <- melt(long,id.vars=c("hospital_block_id","site","m","S","decision_dttm","discharge_dttm"),
           measure.vars=numeric_vars,variable.name="variable",value.name="value")
 z[,`:=`(hours=(m-1)*step_hours,in_hospital=decision_dttm<=discharge_dttm,
          decision_risk_set=m<=tau&m<=S)]
 q <- function(x,p)if(any(is.finite(x)))as.numeric(quantile(x[is.finite(x)],p,names=FALSE)) else NA_real_
 stats <- z[,.(rows=.N,observed=sum(is.finite(value)),missing=sum(!is.finite(value)),
   missing_percent=100*mean(!is.finite(value)),minimum=q(value,0),p01=q(value,.01),
   median=q(value,.5),p99=q(value,.99),maximum=q(value,1)),by=.(site,variable,hours,in_hospital)]
 write(stats,"cleaned_covariate_ranges_and_missingness_by_time.csv")
 miss <- rbindlist(list(
  z[m==1,.(scope="baseline",rows=.N,missing=sum(!is.finite(value))),by=.(site,variable)],
  z[decision_risk_set==TRUE,.(scope="untreated_decisions_during_grace",rows=.N,missing=sum(!is.finite(value))),by=.(site,variable)],
  z[in_hospital==TRUE,.(scope="in_hospital_followup",rows=.N,missing=sum(!is.finite(value))),by=.(site,variable)],
  z[in_hospital==FALSE,.(scope="after_discharge_followup",rows=.N,missing=sum(!is.finite(value))),by=.(site,variable)]))
 miss[,missing_percent:=100*missing/rows];write(miss,"model_relevant_missingness.csv")
 for(v in numeric_vars) {
  p <- ggplot(stats[variable==v&in_hospital],aes(hours/24,median))+
   geom_ribbon(aes(ymin=p01,ymax=p99),fill="#2878A5",alpha=.2)+geom_line()+
   geom_line(aes(y=minimum),linetype=3,color="#BC554A")+
   geom_line(aes(y=maximum),linetype=3,color="#BC554A")+facet_wrap(~site,ncol=3,scales="free_y")+
   labs(title=paste(v,": screened values over time"),
    subtitle="Before median imputation/winsorization; temperature includes audited nearest fills; inpatient rows only",
    x="Days since time zero",y="Value: median, 1st–99th percentiles, min/max")+theme_minimal()
  ggsave(file.path(out,"plots",paste0(v,"_screened_time_ranges.png")),p,width=12,height=7,dpi=160)
 }
}
