# Diagnostics use the exact raw and clipped matrices produced by script 06.
write_weight_diagnostics <- function(analysis,out,sensitive,step_hours) {
 write <- function(x,f) fwrite(x,file.path(out,f),na="")
 summarize <- function(w) {
  w <- w[is.finite(w)]
  if(!length(w))return(list(n=0,mean=NA_real_,sd=NA_real_,minimum=NA_real_,median=NA_real_,
    p95=NA_real_,p99=NA_real_,maximum=NA_real_,ess=NA_real_,ess_fraction=NA_real_,
    over_10=0,over_20=0,over_50=0,over_100=0,top_one_percent_weight_share=NA_real_))
  n <- length(w);ess <- sum(w)^2/sum(w^2)
  list(n=n,mean=mean(w),sd=sd(w),minimum=min(w),median=median(w),
    p95=as.numeric(quantile(w,.95)),p99=as.numeric(quantile(w,.99)),maximum=max(w),
    ess=ess,ess_fraction=ess/n,over_10=sum(w>10),over_20=sum(w>20),over_50=sum(w>50),over_100=sum(w>100),
    top_one_percent_weight_share=sum(sort(w,decreasing=TRUE)[seq_len(max(1,ceiling(n*.01)))])/sum(w))
 }
 trajectory <- probabilities <- coefficients <- fit_checks <- profiles <- clipping <- risk_counts <- list()
 fit_audit <- function(fit,site,component) {
  if(is.null(fit))return(invisible(NULL))
  model <- fit$model
  if(inherits(model,"glm")) {
   beta <- coef(model);se <- rep(NA_real_,length(beta));names(se)<-names(beta)
   tab <- coef(summary(model));se[rownames(tab)]<-tab[,2]
   coefficients[[paste(site,component)]] <<- data.table(site=site,component=component,
     term=names(beta),coefficient=as.numeric(beta),standard_error=as.numeric(se),aliased=is.na(beta))
   fit_checks[[paste(site,component)]] <<- data.table(site=site,component=component,
     n=nobs(model),parameters=length(beta),rank=model$rank,aliased=sum(is.na(beta)),
     converged=isTRUE(model$converged),iterations=model$iter)
  } else fit_checks[[paste(site,component)]] <<- data.table(site=site,component=component,
     n=NA_integer_,parameters=0L,rank=0L,aliased=0L,converged=TRUE,iterations=0L)
 }
 d <- as.data.table(analysis$data)
 for(site_name in names(analysis$site_components)) {
  c <- analysis$site_components[[site_name]]
  fit_audit(c$initiation_fit,site_name,"initiation_denominator")
  fit_audit(c$censoring_fit,site_name,"natural_censoring_denominator")
  for(m in seq_len(c$M))for(arm in c("early","later_or_never")) {
   r <- if(arm=="early")c$R1[,m] else c$R0[,m]
   raw <- if(arm=="early")c$raw_W1[,m] else c$raw_W0[,m]
   clipped <- if(arm=="early")c$W1[,m] else c$W0[,m]
   for(version in c("raw","clipped")) {
    w <- if(version=="raw")raw[r] else clipped[r]
    trajectory[[length(trajectory)+1L]] <- as.data.table(c(list(site=site_name,interval=m,
       hours=(m-1)*step_hours,arm=arm,version=version),summarize(w)))
   }
  }
  for(arm in c("early","later_or_never")) {
   raw <- if(arm=="early")c$raw_W1 else c$raw_W0
   clipped <- if(arm=="early")c$W1 else c$W0
   risk <- if(arm=="early")c$R1 else c$R0
   limits <- quantile(raw[c$obs & is.finite(raw)],analysis$config$trunc,names=FALSE)
   clipping[[paste(site_name,arm)]] <- data.table(site=site_name,arm=arm,
    lower=limits[1],upper=limits[2],risk_set_rows=sum(risk),
    clipped_rows=sum(risk&abs(raw-clipped)>1e-12),
    clipped_percent=100*sum(risk&abs(raw-clipped)>1e-12)/sum(risk))
  }
  site <- d[.site==site_name]
  decision <- site[.time<=analysis$config$tau & .time<=.S]
  decision[,in_hospital:=decision_dttm<=discharge_dttm]
  risk_counts[[site_name]] <- decision[,.(rows=.N,patients=uniqueN(.uid),initiations=sum(.init)),
    by=.(site=.site,interval=.time,in_hospital)]
  cell <- cbind(match(as.character(decision$.uid),c$uid),decision$.time)
  decision[,`:=`(probability=c$Hden[cell],early=as.integer(.S<=analysis$config$tau))]
  probabilities[[site_name]] <- decision[,.(rows=.N,initiations=sum(.init),
    minimum=min(probability),p01=as.numeric(quantile(probability,.01)),median=median(probability),
    p99=as.numeric(quantile(probability,.99)),maximum=max(probability),
    below_001=sum(probability<.01),above_099=sum(probability>.99),
    near_clamp=sum(probability<=1.00001e-6|probability>=1-1.00001e-6)),by=.(site=.site,interval=.time,early)]
  # Small fixed probability bins preserve aggregate-only plots.
  decision[,bin:=pmin(49L,floor(probability*50))]
  bins <- decision[,.(records=.N),by=.(site=.site,early,bin)]
  write(bins,paste0("initiation_probability_bins_",site_name,".csv"))
  vars <- unique(c(analysis$config$baseline,analysis$config$tv))
  for(arm in c("early","later_or_never")) {
   mat <- if(arm=="early")c$raw_W1 else c$raw_W0
   risk <- if(arm=="early")c$R1 else c$R0
   ix <- which(risk&is.finite(mat),arr.ind=TRUE)
   if(!nrow(ix))next
   weights <- mat[ix]; threshold<-as.numeric(quantile(weights,.99))
   rows <- data.table(.uid=c$uid[ix[,1]],.time=ix[,2],weight=weights,tail=weights>=threshold)
   rows <- merge(rows,site[,c(".uid",".time",vars),with=FALSE],by=c(".uid",".time"))
   profiles[[paste(site_name,arm)]] <- rbindlist(lapply(vars,function(v)
    rows[,.(variable=v,mean=mean(get(v)),rows=.N,weight_mean=mean(weight)),by=tail]))[,`:=`(site=site_name,arm=arm)]
  }
 }
 ess <- rbindlist(trajectory,fill=TRUE);write(ess,"ccw_weight_summary_and_ess_by_interval.csv")
 write(rbindlist(clipping),"weight_clipping_thresholds_and_fraction.csv")
 write(rbindlist(risk_counts),"initiation_risk_set_counts.csv")
 write(d[,.(rows=.N,patients=uniqueN(.uid),censor_events=sum(.censor)),by=.site],"natural_censoring_counts.csv")
 # Pooled ESS needs combined weights, not the sum of site ESS values.
 pooled <- rbindlist(lapply(seq_len(analysis$config$horizon),function(m)
  rbindlist(lapply(c("early","later_or_never"),function(arm)
   rbindlist(lapply(c("raw","clipped"),function(version) {
    w <- unlist(lapply(analysis$site_components,function(c) {
     r <- if(arm=="early")c$R1[,m] else c$R0[,m]
     x <- if(version=="raw"){if(arm=="early")c$raw_W1 else c$raw_W0} else {if(arm=="early")c$W1 else c$W0}
     x[r,m]
    }),use.names=FALSE)
    as.data.table(c(list(site="All sites",interval=m,hours=(m-1)*step_hours,arm=arm,version=version),summarize(w)))
   }))))))
 write(pooled,"ccw_pooled_ess_by_interval.csv")
 write(rbindlist(probabilities,fill=TRUE),"initiation_probability_summary.csv")
 write(rbindlist(coefficients,fill=TRUE),"weight_model_coefficients.csv")
 write(rbindlist(fit_checks,fill=TRUE),"weight_model_fit_diagnostics.csv")
 write(rbindlist(profiles,fill=TRUE),"extreme_weight_covariate_profiles.csv")
 landmark <- analysis$fits$landmark_ipw$landmark_diagnostics
 lm_summary <- rbindlist(lapply(names(landmark),function(s) {
  a <- landmark[[s]];fit_audit(a$fit,s,"landmark_propensity")
  rbindlist(lapply(0:1,function(arm)rbindlist(lapply(c("raw","clipped"),function(version) {
   w <- if(version=="raw")a$raw_weights[a$A==arm] else a$weights[a$A==arm]
   as.data.table(c(list(site=s,arm=arm,version=version),summarize(w)))
  }))))
 }))
 write(lm_summary,"landmark_weight_summary_and_ess.csv")
 lm_probs <- rbindlist(lapply(names(landmark),function(s) {
  a <- landmark[[s]]
  data.table(site=s,early=a$A,bin=pmin(49L,floor(a$propensity*50)))[,.(patients=.N),by=.(site,early,bin)]
 }))
 write(lm_probs,"landmark_propensity_bins.csv")
 p <- ggplot(lm_probs,aes((bin+.5)/50,patients,fill=factor(early)))+
  geom_col(position="identity",alpha=.5,width=.02)+facet_wrap(~site,scales="free_y")+
  labs(title="Landmark propensity-score distribution",x="Probability of initiation by grace end",y="Patients",fill="Early initiation")+theme_minimal()
 ggsave(file.path(out,"landmark_propensity_distribution.png"),p,width=12,height=7,dpi=160)
 write(rbindlist(coefficients,fill=TRUE),"weight_model_coefficients.csv")
 write(rbindlist(fit_checks,fill=TRUE),"weight_model_fit_diagnostics.csv")
 saveRDS(list(config=analysis$config,site_components=analysis$site_components,landmark=landmark),
   file.path(sensitive,"weight_diagnostics_patient_level.rds"))
 plot <- ggplot(rbindlist(list(ess,pooled)),aes(hours/24,ess,color=version,linetype=arm))+
   geom_line()+facet_wrap(~site,scales="free_y",ncol=2)+geom_vline(xintercept=analysis$config$tau*step_hours/24,linetype=3)+
   labs(title="CCW effective sample size",subtitle="Exact clone risk sets; raw versus clipped weights",
        x="Days since time zero",y="ESS",color="Weight",linetype="Strategy")+theme_minimal()
 ggsave(file.path(out,"ccw_ess_over_time.png"),plot,width=12,height=10,dpi=160)
 bins <- rbindlist(lapply(names(analysis$site_components),function(s)fread(file.path(out,paste0("initiation_probability_bins_",s,".csv")))))
 plot <- ggplot(bins,aes((bin+.5)/50,records,fill=factor(early)))+
  geom_col(position="identity",alpha=.5,width=.02)+facet_wrap(~site,scales="free_y")+
  labs(title="Conditional initiation probabilities during grace",subtitle="Untreated decision rows; this is a time-specific probability, not one baseline propensity score",
       x="Predicted initiation probability",y="Decision rows",fill="Observed early initiation")+theme_minimal()
 ggsave(file.path(out,"initiation_probability_distribution.png"),plot,width=12,height=7,dpi=160)
 writeLines(c("Raw and clipped summaries use exactly the clone risk sets used for CCW.",
  "ESS across person-intervals is not an independent-patient sample size: use interval-specific ESS.",
  "Weights are clipped at the site/arm 1st and 99th percentiles using the inherited observed-row convention.",
  "Initiation probabilities are conditional six-hour hazards among not-yet-initiated patients during grace.",
  "Extreme-weight profiles compare the top 1% of raw risk-set weights with the remaining rows; ties can include more than 1%.",
  "Covariate means in these profiles are encoded/scaled model values. Inspect ipcw_model_specification.csv for exact terms.",
  "Model nonconvergence, aliased coefficients and very small probabilities require review; no automatic model switching occurs.",
  "Confidence intervals treat fitted nuisance models as fixed, following the inherited five-method implementation."),
  file.path(out,"WEIGHT_DIAGNOSTICS_README.txt"))
}
