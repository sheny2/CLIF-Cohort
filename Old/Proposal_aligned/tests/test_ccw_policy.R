# Verify policy details that differ from the legacy parent's never-treated arm.
self <- sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])
source(file.path(dirname(dirname(normalizePath(self))),'06_generalized_five_methods.R'))
set.seed(9)
N <- 120;M <- 6;tau <- 2
S <- rep(c(1,2,4,Inf),length.out=N)
T <- rep(c(2,5,6,6,6,6,6),length.out=N)
d <- do.call(rbind,lapply(seq_len(N),function(i) data.frame(id=i,m=seq_len(T[i]),
 event=as.integer(seq_len(T[i])==T[i]),censor=0L,init=as.integer(seq_len(T[i])==S[i]),x=i%%3)))
a <- suppressWarnings(run_pneumonia_five_methods(d,id_col='id',time_col='m',event_col='event',
 censor_col='censor',initiation_col='init',grace_period=tau,horizon=M,baseline_covariates='x',
 model_terms=list(initiation_denominator=c('.time','x'),initiation_numerator='.time'),verbose=FALSE))
comp <- a$site_components[[1]]
# At interval 2 (the final grace bin), every observed death remains at risk in
# the treatment clone, even if the patient never initiated before that death.
stopifnot(any(comp$event[,tau]), all(comp$R1[,tau][comp$event[,tau]]))
# The original proposal analysis retains never-initiate-through-followup controls.
late <- comp$S>tau & is.finite(comp$S) & comp$S<=M
stopifnot(any(late),all(!comp$R0[late,4]))
# Landmark propensity terms must be baseline-only in this independent library.
body_text <- paste(deparse(body(g6_landmark_local_stats)),collapse=' ')
stopifnot(grepl('propensity_terms <- cfg\\$baseline',body_text))
stopifnot(length(a$config$lag_names)==0L)
cat('PASS: empty time-varying set, final-grace deaths retained, controls censored at any initiation, baseline-only landmark model.\n')
