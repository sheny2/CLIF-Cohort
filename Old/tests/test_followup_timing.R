suppressPackageStartupMessages({library(data.table);library(ggplot2)})
self <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
root <- dirname(dirname(normalizePath(self)))
source(file.path(root, "followup_timing_helpers.R"))
t0 <- as.POSIXct("2026-01-01", tz = "UTC")
c <- data.table(hospital_block_id=as.character(1:8), respiratory_failure_start=t0,
                event_dttm=t0+c(0,6,24,96,672,673,100,200)*3600,
                discharge_dttm=t0+672*3600)
a <- data.table(hospital_block_id=c('1','2','2','3','4','5','6','7','8','8'),
                admin_dttm=t0+c(0,-1,6,24,96,672,673,101,199,200)*3600)
x <- build_followup_timing(c,a,t0+672*3600,1L)
stopifnot(nrow(x$intervals)==224,
          x$grace_comparison[event_type=='Death' & grace_hours==24,n_reportable]==3,
          x$grace_comparison[event_type=='Death' & grace_hours==96,n_reportable]==4,
          x$grace_comparison[event_type!='Death' & grace_hours==24,n_reportable]==3,
          x$grace_comparison[event_type!='Death' & grace_hours==96,n_reportable]==4,
          x$intervals[event_type=='Death' & interval==1,n_reportable]==2,
          x$intervals[event_type=='Death' & interval==112,n_reportable]==1,
          x$intervals[event_type!='Death',sum(n_reportable)]==6)
y <- build_followup_timing(c,a,t0+25*3600,1L)
stopifnot(y$grace_comparison[event_type=='Death' & grace_hours==96,n_reportable]==3)
empty <- build_followup_timing(c[0],a[0],t0+672*3600)
stopifnot(nrow(empty$intervals)==224,all(empty$intervals$n_reportable==0))
# Small cells must remain hidden; no raw event counts in either export.
z <- build_followup_timing(c,a,t0+672*3600,10L)
stopifnot(!'n' %in% names(z$intervals),!'n' %in% names(z$grace_comparison),
          all(is.na(z$intervals[n_reportable!=0 | is.na(n_reportable),n_reportable])))
# Complementary suppression when a window has exactly one small bin.
big <- c[rep(3L,12L)];big[,hospital_block_id:=as.character(seq_len(.N))]
big[1,event_dttm:=t0+6*3600]
b <- build_followup_timing(big,a[0],t0+672*3600,10L)
stopifnot(sum(b$intervals[event_type=='Death' & interval_end_hours<=24,count_suppressed])==2)
# A single small 24-96h cell cannot be recovered by differencing totals.
single <- c[rep(3L,12L)]; single[,hospital_block_id:=as.character(seq_len(.N))]
single[1,event_dttm:=t0+30*3600]
u <- build_followup_timing(single,a[0],t0+672*3600,10L)
stopifnot(is.na(u$grace_comparison[event_type=='Death' & grace_hours==96,n_reportable]),
          u$grace_comparison[event_type=='Death' & grace_hours==24,n_reportable]==11)
# Smoke-render the actual plotting function extracted from script 05.
program <- parse(file.path(root,'05_cohort_summary_eda.R'))
expr <- Filter(function(e) is.call(e) && identical(e[[1]],as.name('<-')) &&
                 identical(e[[2]],as.name('make_timing_plot')),as.list(program))
stopifnot(length(expr)==1)
eval(expr[[1]])
cutoff_lines <- data.table(hours=c(24,96),label=c('24 hours','96 hours'))
min_cell_n <- 1L
png <- tempfile(fileext='.png', tmpdir='/tmp')
ggsave(png,make_timing_plot(x$intervals),width=12,height=8,dpi=100)
stopifnot(file.info(png)$size>0)
cat('Plot preview:',png,'\n')
cat('PASS: six-hour boundaries, 24/96 totals, first dose, cutoff/death limits, empty data, suppression and plot rendering.\n')
