#!/usr/bin/env Rscript
# Independent local CCW analyses followed by inverse-variance meta-analysis.
# Run AFTER 07: Rscript 08_local_ccw_meta_analysis.R RUN_ROOT
# Only reads prepared 07 outputs; writes new aggregate/local_ccw_meta and
# sensitive/local_ccw_meta outputs. Existing five-method files are unchanged.
suppressPackageStartupMessages({library(arrow);library(data.table);library(ggplot2)})
if(!requireNamespace("metafor",quietly=TRUE))stop("08 requires metafor: install.packages('metafor')")
data.table::setDTthreads(1L)
a<-commandArgs(TRUE)
f<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
root<-dirname(normalizePath(gsub("~+~"," ",f,fixed=TRUE)))
for(helper in c("06_generalized_five_methods.R","variable_rules.R","outcome_helpers.R","local_meta_helpers.R"))source(file.path(root,helper))
run_root<-if(length(a)&&nzchar(a[1]))a[1] else Sys.getenv("RUN_ROOT","")
if(!nzchar(run_root)||!dir.exists(run_root))stop("Supply the existing RUN_ROOT from a completed 07 run")
input<-file.path(run_root,"sensitive","analysis");original<-file.path(run_root,"aggregate","analysis")
required<-c(file.path(input,c("pneumonia_ccw_long.parquet","analysis_design.rds")),
 file.path(original,c("five_method_results.csv","COMPLETE.txt")))
if(any(!file.exists(required)))stop("07 must finish first. Missing: ",paste(required[!file.exists(required)],collapse=", "))
design<-readRDS(file.path(input,"analysis_design.rds"))
if(!identical(design$outcome,PRIMARY_OUTCOME))stop("07 outcome differs from the current composite endpoint; rerun 07 first")
out<-file.path(run_root,"aggregate","local_ccw_meta");sensitive<-file.path(run_root,"sensitive","local_ccw_meta")
dir.create(out,recursive=TRUE,showWarnings=FALSE);dir.create(sensitive,recursive=TRUE,showWarnings=FALSE)
unlink(file.path(out,"COMPLETE.txt"))
write<-function(x,f) {
 if(any(c("patient_id","hospital_block_id","hospitalization_id",".uid")%in%names(x)))stop("Identifier in aggregate output")
 fwrite(x,file.path(out,f),na="")
}
long<-as.data.table(read_parquet(file.path(input,"pneumonia_ccw_long.parquet")))
if(any(long[,uniqueN(site),by=hospital_block_id]$V1!=1L))stop("Hospital block assigned to more than one site")
results<-fits<-term_lists<-processing<-diagnostics<-curves<-weight_checks<-list();all_warnings<-character()
for(s in sort(unique(long$site))) {
 message("Fitting local CCW at ",s)
 site_long<-long[site==s];units<-unique(site_long[,.(hospital_block_id,S,delta)])
 early<-sum(units$S<=design$tau_intervals);other<-nrow(units)-early
 row<-data.table(site=s,patients=nrow(units),early_initiators=early,later_or_never=other,
   primary_events_28d=sum(units$delta),outcome=PRIMARY_OUTCOME,site_valid=FALSE,site_status="",
   outcome_risk_treatment=NA_real_,outcome_risk_no_treatment=NA_real_,risk_difference=NA_real_,risk_difference_se=NA_real_,
   risk_difference_lower=NA_real_,risk_difference_upper=NA_real_,risk_ratio=NA_real_,log_risk_ratio_se=NA_real_,
   risk_ratio_lower=NA_real_,risk_ratio_upper=NA_real_,odds_ratio=NA_real_,log_odds_ratio_se=NA_real_,
   odds_ratio_lower=NA_real_,odds_ratio_upper=NA_real_,rmst_treatment_days=NA_real_,rmst_no_treatment_days=NA_real_,
   rmst_difference_days=NA_real_,rmst_difference_se_days=NA_real_,rmst_difference_lower_days=NA_real_,rmst_difference_upper_days=NA_real_)
 if(early==0L||other==0L) {
  row[,site_status:="No observed early or later/never group; cannot support local contrast"]
  results[[s]]<-row;next
 }
 local_warnings<-character()
 attempt<-tryCatch(withCallingHandlers({
  p<-local_ccw_prepare(site_long,design)
  # Numerator uses THIS SITE's empirical interval-specific initiation rates.
  num<-g6_common_init_numerator(p$data,p$cfg$horizon)
  component<-g6_local_components(p$data,p$cfg,common_init_num=num)
  fit<-g6_central_ccw(list(g6_local_stats(component)))
  list(prepared=p,component=component,fit=fit,numerator=num)
 },warning=function(w){local_warnings<<-c(local_warnings,conditionMessage(w));invokeRestart("muffleWarning")}),
 error=function(err)err)
 if(length(local_warnings))all_warnings<-c(all_warnings,paste(s,local_warnings))
 if(inherits(attempt,"error")) {
  row[,site_status:=paste("Local fit failed:",conditionMessage(attempt))];results[[s]]<-row;next
 }
 p<-attempt$prepared;c<-attempt$component;fit<-attempt$fit
 fits[[s]]<-attempt
 converged<-function(fit)!inherits(fit$model,"glm")||isTRUE(fit$model$converged)
 status<-if(!converged(c$initiation_fit)||!converged(c$censoring_fit))"Nuisance model did not converge" else "ok"
 # A flat curve after loss of every at-risk clone with nonzero survival is not
 # sufficient support for a 28-day contrast. Distinguish this from all dying.
 if(any(colSums(c$W1*c$R1)<=0 & head(fit$S1,-1)>1e-8)||
    any(colSums(c$W0*c$R0)<=0 & head(fit$S0,-1)>1e-8))status<-"Empty clone risk set before horizon with nonzero event-free probability"
 row[,`:=`(site_valid=status=="ok",site_status=status,
  outcome_risk_treatment=fit$psi1,outcome_risk_no_treatment=fit$psi0,
  risk_difference=fit$RD,risk_difference_se=fit$SE_RD,risk_difference_lower=fit$CI_RD[1],risk_difference_upper=fit$CI_RD[2],
  risk_ratio=fit$RR,log_risk_ratio_se=fit$SE_logRR,risk_ratio_lower=fit$CI_RR[1],risk_ratio_upper=fit$CI_RR[2],
  odds_ratio=fit$OR,log_odds_ratio_se=fit$SE_logOR,odds_ratio_lower=fit$CI_OR[1],odds_ratio_upper=fit$CI_OR[2],
  rmst_treatment_days=fit$RMST1*design$step_hours/24,rmst_no_treatment_days=fit$RMST0*design$step_hours/24,
  rmst_difference_days=fit$RMST_diff*design$step_hours/24,rmst_difference_se_days=fit$SE_RMSTdiff*design$step_hours/24,
  rmst_difference_lower_days=fit$CI_RMSTdiff[1]*design$step_hours/24,rmst_difference_upper_days=fit$CI_RMSTdiff[2]*design$step_hours/24)]
 results[[s]]<-row
 processing[[s]]<-copy(p$processing)[,site:=s]
 term_lists[[s]]<-rbindlist(lapply(names(p$cfg$terms),function(nm)data.table(site=s,component=nm,term=p$cfg$terms[[nm]],
   type=fcase(p$cfg$terms[[nm]]%in%p$cfg$baseline,"baseline",grepl("^lag_",p$cfg$terms[[nm]]),"lagged time-varying",default="current/time/treatment"),
   note=if(nm=="initiation_numerator")"Actual numerator is the site's empirical interval-specific initiation rate" else "")))
 for(nm in c("initiation_fit","censoring_fit")) {
  model<-c[[nm]]$model
  diagnostics[[paste(s,nm)]]<-data.table(site=s,component=nm,
   converged=converged(c[[nm]]),rank=if(inherits(model,"glm"))model$rank else 0L,
   coefficients=if(inherits(model,"glm"))length(coef(model)) else 0L,
   aliased=if(inherits(model,"glm"))sum(is.na(coef(model))) else 0L)
 }
 curves[[s]]<-data.table(site=s,day=rep(0:c$M*design$step_hours/24,2),
   strategy=rep(c("Early initiation","Later or never"),each=c$M+1L),event_free_probability=c(fit$S1,fit$S0))
 for(m in unique(c(design$tau_intervals,design$horizon_intervals)))for(arm in c("early","later_or_never"))for(version in c("raw","clipped")) {
  risk<-if(arm=="early")c$R1[,m] else c$R0[,m]
  mat<-if(version=="raw"){if(arm=="early")c$raw_W1 else c$raw_W0} else {if(arm=="early")c$W1 else c$W0}
  w<-mat[risk,m]
  weight_checks[[length(weight_checks)+1L]]<-data.table(site=s,interval=m,arm=arm,version=version,n=length(w),
   maximum=if(length(w))max(w) else NA_real_,ess=if(length(w)&&sum(w^2)>0)sum(w)^2/sum(w^2) else NA_real_)
 }
}
site_results<-rbindlist(results,fill=TRUE);write(site_results,"local_ccw_site_results.csv")
write(rbindlist(processing,fill=TRUE),"local_covariate_processing.csv")
write(rbindlist(term_lists,fill=TRUE),"local_ipcw_model_specification.csv")
write(rbindlist(diagnostics,fill=TRUE),"local_nuisance_fit_diagnostics.csv")
write(rbindlist(weight_checks,fill=TRUE),"local_weight_checks.csv")
if(length(curves))write(rbindlist(curves),"local_event_free_curves.csv")
meta<-pool_local_contrasts(meta_contrasts(site_results))
write(meta$summary,"meta_analysis_results.csv");write(meta$weights,"meta_analysis_site_weights_and_exclusions.csv")
all_warnings<-unique(c(all_warnings,meta$warnings));writeLines(all_warnings,file.path(out,"warnings.txt"))
saveRDS(list(design=design,sites=fits,meta_models=meta$models),file.path(sensitive,"local_ccw_meta_fits.rds"))
# Append contrast estimates, not separately pooled arm risks: each estimand has
# its own inverse-variance weights, so arm risks would not define these contrasts.
base<-fread(file.path(original,"five_method_results.csv"))
if(nrow(base)!=5L||any(base$outcome!=PRIMARY_OUTCOME))stop("Expected five matching-outcome methods from 07")
extra<-copy(base[1]);for(nm in names(extra))set(extra,j=nm,value=if(is.numeric(extra[[nm]]))NA_real_ else NA_character_)
extra[,`:=`(method_key="local_ccw_meta",method="Local CCW + random-effects meta-analysis",outcome=PRIMARY_OUTCOME,
 variance_type="Site influence-function variances (nuisance fixed); REML heterogeneity; Hartung-Knapp t intervals")]
map<-list(risk_difference=c("risk_difference","risk_difference_se","risk_difference_lower","risk_difference_upper"),
 log_risk_ratio=c("risk_ratio","log_risk_ratio_se","risk_ratio_lower","risk_ratio_upper"),
 log_odds_ratio=c("odds_ratio","log_odds_ratio_se","odds_ratio_lower","odds_ratio_upper"),
 rmst_difference_days=c("rmst_difference_days","rmst_difference_se_days","rmst_difference_lower_days","rmst_difference_upper_days"))
for(e in names(map)) {
 z<-meta$summary[estimand==e&meta_method=="REML_HK"]
 set(extra,j=map[[e]][1],value=z$estimate);set(extra,j=map[[e]][2],value=z$analysis_scale_se)
 set(extra,j=map[[e]][3],value=z$lower);set(extra,j=map[[e]][4],value=z$upper)
}
six<-rbindlist(list(base,extra),fill=TRUE);write(six,"six_method_results.csv")
plot_data<-function(x,labelcol)rbindlist(list(
 x[,.(method=get(labelcol),estimand="28-day primary-outcome risk difference",estimate=risk_difference,lower=risk_difference_lower,upper=risk_difference_upper)],
 x[,.(method=get(labelcol),estimand="28-day restricted mean event-free time difference (days)",estimate=rmst_difference_days,lower=rmst_difference_lower_days,upper=rmst_difference_upper_days)]))
forest<-function(z,levels,title,filename) {
 z[,method:=factor(method,levels=rev(levels))]
 z[,meta_result:=grepl("meta-analysis",as.character(method))]
 p<-ggplot(z,aes(estimate,method,color=meta_result,shape=meta_result))+
  geom_vline(xintercept=0,linetype=3,color="grey50")+
  geom_errorbar(aes(xmin=lower,xmax=upper),orientation="y",width=.18,na.rm=TRUE)+geom_point(size=2.7,na.rm=TRUE)+
  scale_color_manual(values=c("FALSE"="#2878A5","TRUE"="#9B3675"),guide="none")+
  scale_shape_manual(values=c("FALSE"=16,"TRUE"=18),guide="none")+
  facet_wrap(~estimand,scales="free_x",ncol=1)+
  labs(title=title,subtitle=PRIMARY_OUTCOME,x="Estimate and 95% CI",y=NULL,
    caption="Positive RD: more adverse events; positive event-free time difference: benefit. Meta-analysis: REML + Hartung–Knapp.")+theme_minimal(base_size=11)
 ggsave(file.path(out,filename),p,width=12,height=9,dpi=200)
}
forest(plot_data(six,"method"),six$method,"Six-method treatment-effect comparison","six_method_effect_comparison.png")
# Site forest includes both pooling choices, with eligibility/exclusions in CSV.
site_plot<-copy(site_results)[,method:=ifelse(site_valid,site,paste0(site," (excluded from pooling)"))]
for(kind in c("REML_HK","fixed_effect")) {
 row<-copy(extra);row[,method:=if(kind=="REML_HK")"Pooled random-effects meta-analysis" else "Pooled fixed-effect meta-analysis"]
 for(e in c("risk_difference","rmst_difference_days")) {
  z<-meta$summary[estimand==e&meta_method==kind]
  for(j in c(1,3,4))set(row,j=map[[e]][j],value=z[[c("estimate","lower","upper")[match(j,c(1,3,4))]]])
 }
 site_plot<-rbindlist(list(site_plot,row),fill=TRUE)
}
forest(plot_data(site_plot,"method"),site_plot$method,"Local CCW estimates and meta-analysis","local_ccw_site_forest.png")
rr<-six[,.(method,risk_ratio,risk_ratio_lower,risk_ratio_upper)]
rr[,method:=factor(method,levels=rev(six$method))]
p<-ggplot(rr,aes(risk_ratio,method))+geom_vline(xintercept=1,linetype=3)+
 geom_errorbar(aes(xmin=risk_ratio_lower,xmax=risk_ratio_upper),orientation="y",width=.18,na.rm=TRUE)+
 geom_point(size=2.7,na.rm=TRUE)+scale_x_log10()+labs(title="Six-method primary-outcome risk ratios",subtitle=PRIMARY_OUTCOME,
 x="Risk ratio (log scale), 95% CI",y=NULL)+theme_minimal()
ggsave(file.path(out,"six_method_risk_ratio_forest.png"),p,width=12,height=6,dpi=200)
writeLines(c("08 fits weights AND weighted event-free outcome curves within each site independently.",
 "Screened raw values/nearest temperature and selected variable definitions come from 07; imputation, winsorization and scaling are recalculated locally.",
 "Initiation numerator rates, denominator fits and clipping are site-local; no common cross-site numerator is reused.",
 "Site effects/SEs are pooled: RD and event-free time on original scales, RR/OR on log scales.",
 "Primary pooling: REML random effects with Hartung-Knapp 95% t intervals; fixed-effect inverse-variance pooling is supplementary.",
 "No continuity correction is added for zero-risk ratio contrasts; invalid/nonpositive variances are excluded per estimand with reasons.",
 "Nonconverged local nuisance fits and unsupported clone risk sets are excluded from pooling; their point estimates remain visible for review.",
 "At least two valid sites are required; no fallback pooled result is fabricated if REML fails.",
 "The meta-analysis targets an average site effect and need not equal the participant-pooled CCW effect.",
 "With five sites, heterogeneity and prediction intervals are imprecise. Causal interpretation still requires exchangeability and positivity.",
 "Influence-function site variances treat estimated nuisance parameters as fixed, matching 07; no bootstrap is performed.",
 "Arm-specific risks/RMSTs are not separately pooled in the sixth row because contrast weights differ by estimand.",
 "07 outputs are read-only; all new figures/results live under aggregate/local_ccw_meta."),file.path(out,"METHOD_NOTES.txt"))
writeLines(capture.output(sessionInfo()),file.path(out,"SESSION_INFO.txt"))
writeLines("Complete",file.path(out,"COMPLETE.txt"))
print(meta$summary)
cat("\n08 aggregate results:",normalizePath(out),"\n")
