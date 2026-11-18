#!/usr/bin/env Rscript
# Rscript 11_prior_utilization_eda.R RUN_ROOT [EXPANDED_HOSPITALIZATION_PARQUET] [HISTORY_COVERAGE_START_DATE]
suppressPackageStartupMessages({library(arrow);library(data.table);library(dplyr);library(ggplot2)})
setDTthreads(1L)
a <- commandArgs(TRUE)
if(!length(a))stop("Supply a completed 07 RUN_ROOT")
f <- sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1L])
root <- dirname(normalizePath(gsub("~+~"," ",f,fixed=TRUE)))
for(h in c("hospital_block_helpers.R","site_helpers.R","history_eda_helpers.R"))source(file.path(root,h))
run <- normalizePath(a[1]);expanded <- if(length(a)>=2L)a[2] else "/projects/emulate_icu/rclif_2_1/clif_hospitalization.parquet"
coverage <- if(length(a)>=3L)as.Date(a[3]) else as.Date(NA)
if(length(a)>=3L&&is.na(coverage))stop("Coverage start must be YYYY-MM-DD")
if(!file.exists(file.path(run,"aggregate","analysis","COMPLETE.txt")))stop("Complete 07 first; 09/10 are not required")
inp <- file.path(run,"sensitive","analysis")
c <- as.data.table(read_parquet(file.path(inp,"pneumonia_ccw_cohort.parquet")))
cw <- as.data.table(read_parquet(file.path(inp,"sampled_hospital_block_crosswalk.parquet")))
design <- readRDS(file.path(inp,"analysis_design.rds"))
assert_mapped_sites(c$site)
if(!all(c("patient_id","hospital_block_id","site","admission_dttm","S")%in%names(c)))stop("Prepared cohort lacks required fields")
c[,treatment_group:=ifelse(S<=design$tau_intervals,"Initiation within grace","Later or never initiation")]
ds <- open_dataset(expanded)
required <- c("patient_id","hospitalization_id","admission_dttm","discharge_dttm","admission_type_category")
if(!all(required%in%ds$schema$names))stop("Expanded hospitalization table lacks: ",paste(setdiff(required,ds$schema$names),collapse=", "))
ids <- as.character(c$patient_id)
h <- ds |> filter(.data$patient_id%in%ids) |> select(all_of(required)) |> collect() |> as.data.table()
x <- build_history_eda(c,cw,h,coverage)
out <- file.path(run,"aggregate","prior_utilization_eda");sensitive <- file.path(run,"sensitive","prior_utilization_eda")
dir.create(out,recursive=TRUE,showWarnings=FALSE);dir.create(sensitive,recursive=TRUE,showWarnings=FALSE)
unlink(file.path(out,"COMPLETE.txt"))
write <- function(z,name) {
  if(any(c("patient_id","hospitalization_id","hospital_block_id","expanded_episode_id")%in%names(z)))stop("Identifier in aggregate output")
  fwrite(z,file.path(out,name),na="")
}
write(x$audit,"history_source_audit.csv");write(x$mapping,"admission_type_mapping.csv")
write_parquet(x$patient,file.path(sensitive,"patient_prior_utilization.parquet"))
write_parquet(x$prior,file.path(sensitive,"completed_prior_episodes.parquet"))
write_parquet(x$invalid,file.path(sensitive,"invalid_history_rows.parquet"))
write_parquet(x$overlap,file.path(sensitive,"prior_year_overlap_days.parquet"))
# All patients, including observed zero counts, enter the patient summaries.
scopes <- rbindlist(list(copy(x$patient)[,scope:="Overall"],
  copy(x$patient)[,scope:=paste0("Site: ",site)],
  copy(x$patient)[,scope:=paste0("Treatment: ",treatment_group)]))
availability <- scopes[,.(patients=.N,with_observed_prior_episode=sum(observed_episode_count>0),
  percent_with_observed_prior_episode=100*mean(observed_episode_count>0),
  patients_with_invalid_possible_history=sum(invalid_possible_prior_rows>0),
  patients_with_unknown_ed_origin=sum(unknown_ed_episodes>0),
  verified_year_coverage=sum(year_coverage_verified==TRUE,na.rm=TRUE),
  short_year_coverage=sum(year_coverage_verified==FALSE,na.rm=TRUE),
  unverified_year_coverage=sum(is.na(year_coverage_verified))),by=.(scope,window,encounter_type)]
write(availability,"history_availability_and_prevalence.csv")
metrics <- c("observed_episode_count","total_los_days","mean_los_days","longest_los_days","days_since_last_discharge")
summary <- rbindlist(lapply(metrics,function(metric)scopes[,c(list(metric=metric),history_summary(get(metric))),by=.(scope,window,encounter_type)]))
write(summary,"patient_count_and_duration_summary.csv")
bins <- scopes[,.(patients=.N),by=.(scope,window,encounter_type,count_category=ifelse(observed_episode_count>=5,"5+",as.character(observed_episode_count)))]
bins[,percent:=100*patients/sum(patients),by=.(scope,window,encounter_type)]
write(bins,"prior_episode_count_distribution.csv")
episodes <- merge(x$prior,c[,.(patient_id,site,treatment_group)],by="patient_id")
ep <- rbindlist(lapply(c("All available prior history","Previous 365 days"),function(window) {
  z <- if(window=="Previous 365 days")episodes[within_year==TRUE] else copy(episodes)
  rbind(copy(z)[,`:=`(window=window,encounter_type="All hospitalizations")],
    copy(z[ed_origin==TRUE])[,`:=`(window=window,encounter_type="ED-origin hospitalizations")])
}))
epscopes <- rbindlist(list(copy(ep)[,scope:="Overall"],copy(ep)[,scope:=paste0("Site: ",site)],
  copy(ep)[,scope:=paste0("Treatment: ",treatment_group)]))
write(epscopes[,history_summary(los_days),by=.(scope,window,encounter_type)],"episode_length_of_stay_summary.csv")
overlap <- merge(CJ(patient_id=c$patient_id,encounter_type=c("All hospitalizations","ED-origin hospitalizations")),x$overlap,
  by=c("patient_id","encounter_type"),all.x=TRUE)
overlap[is.na(observed_days_in_previous_year),observed_days_in_previous_year:=0]
write(overlap[,history_summary(observed_days_in_previous_year),by=encounter_type],"hospital_days_overlapping_previous_year.csv")
saveplot <- function(p,name,height=7)ggsave(file.path(out,name),p,width=12,height=height,dpi=180)
bins[,count_category:=factor(count_category,levels=c("0","1","2","3","4","5+"))]
saveplot(ggplot(bins[scope=="Overall"],aes(count_category,percent))+geom_col(fill="#2878A5")+
  facet_grid(encounter_type~window)+theme_minimal()+labs(title="Observed prior hospitalization counts",x="Episodes per patient",y="Patients (%)"),"prior_counts_overall.png")
saveplot(ggplot(bins[grepl("^Treatment:",scope)],aes(count_category,percent,fill=scope))+
  geom_col(position="dodge")+facet_grid(encounter_type~window)+theme_minimal()+
  labs(title="Prior utilization by observed initiation group",subtitle="Unweighted descriptive comparison",x="Episodes per patient",y="Patients (%)",fill=NULL)+
  theme(legend.position="bottom"),"prior_counts_by_treatment.png")
saveplot(ggplot(bins[grepl("^Site:",scope)],aes(count_category,percent,fill=scope))+
  geom_col(position="dodge")+facet_grid(encounter_type~window)+theme_minimal()+
  labs(title="Prior utilization by index hospital site",x="Episodes per patient",y="Patients (%)",fill=NULL)+
  theme(legend.position="bottom"),"prior_counts_by_site.png")
if(nrow(ep)) {
  duration_plot <- ggplot(ep,aes(log1p(los_days)))+geom_histogram(bins=40,fill="#2878A5")+
    facet_grid(encounter_type~window,scales="free_y")+theme_minimal()+
    labs(title="Length of completed prior hospitalization episodes",x="log(1 + episode length in days)",y="Episodes")
} else duration_plot <- ggplot()+annotate("text",x=0,y=0,label="No completed prior episodes observed")+theme_void()
saveplot(duration_plot,"prior_episode_duration.png")
saveplot(ggplot(x$patient,aes(treatment_group,log1p(total_los_days),fill=treatment_group))+
  geom_boxplot(outlier.size=.6)+facet_grid(encounter_type~window)+theme_minimal()+
  labs(title="Cumulative prior episode duration per patient",x=NULL,y="log(1 + total episode days)")+theme(legend.position="none"),"prior_total_duration_by_treatment.png")
writeLines(c("All available history means all observed completed episodes before INDEX BLOCK ADMISSION, never future encounters.",
  "Previous year: episode START within 365 days before index admission, matching 09/10; full episode duration is used.",
  "Episodes overlap/merge when gaps <=6 hours. Duration is episode admission-to-discharge span, including bridging gaps; not exact occupied bed-days.",
  "The entire expanded episode containing any selected index encounter is excluded, including connected transfers.",
  "ED-origin means the FIRST component's admission_type_category is ed. This is a hospital stay, not a standalone ED visit or ED dwell duration.",
  "Other/unknown/missing admission categories have unknown ED origin. ED counts are confirmed observed counts; unknowns are separately reported.",
  "Invalid-date/ID records are excluded from episode construction and audited; observed counts may underestimate utilization.",
  "Zero means no qualifying episode observed, not independently verified absence. No patient is excluded here for short/unverified coverage.",
  paste("Explicit system coverage start:",coverage,". All-history comparisons have unequal observation opportunities."),
  "Patient summaries include observed zeros for count, total duration, and longest duration; mean duration and recency are undefined without an episode.",
  "Episode summaries condition on observed episodes. Counts and durations are shown overall, by index site, and by observed treatment group.",
  "Separate overlap-year days clip all completed prior episodes to the calendar lookback, including episodes that started earlier.",
  "No weights, model fitting, trimming or changes to 01--10. Plots use log1p for duration displays; summary CSVs retain original units and maxima.",
  "Patient-level files are sensitive; aggregate counts are unsuppressed."),file.path(out,"METHOD_NOTES.txt"))
writeLines(capture.output(sessionInfo()),file.path(out,"SESSION_INFO.txt"))
writeLines("Complete",file.path(out,"COMPLETE.txt"))
print(availability[scope=="Overall"]);cat("\nResults:",out,"\n")
