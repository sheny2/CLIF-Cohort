# Shared helpers, local to this standalone analysis.
suppressPackageStartupMessages({library(data.table); library(arrow); library(dplyr)})
cape_norm <- function(x) gsub("[^a-z0-9]", "", tolower(as.character(x)))
cape_read <- function(path, columns, ids = NULL) {
  ds <- open_dataset(path)
  if (!all(columns %in% ds$schema$names)) stop("Missing fields in ", basename(path), ": ", paste(setdiff(columns, ds$schema$names), collapse=", "))
  q <- ds |> select(all_of(columns))
  if (!is.null(ids)) q <- q |> filter(.data$hospitalization_id %in% ids)
  as.data.table(collect(q))
}
cape_steroid_records <- function(data_dir, ids) {
  out <- lapply(c("intermittent", "continuous"), function(source) {
    path <- file.path(data_dir,paste0("clif_medication_admin_",source,".parquet"))
    schema <- open_dataset(path)$schema$names
    cols <- c("hospitalization_id","admin_dttm","med_category","med_dose","med_dose_unit","med_route_category")
    if (source == "continuous" && !"med_route_category" %in% schema) cols <- setdiff(cols,"med_route_category")
    x <- cape_read(path,unique(c(cols, intersect("mar_action_category",schema))),ids)
    if (!"med_route_category" %in% names(x)) x[,med_route_category:="iv"]
    x[, `:=`(source_table=source, drug=cape_norm(med_category),
             route=cape_norm(med_route_category), dose=suppressWarnings(as.numeric(med_dose)),
             unit=cape_norm(gsub("μ","u",as.character(med_dose_unit),fixed=TRUE)),
             admin_dttm=as.POSIXct(admin_dttm,tz="UTC"),hospitalization_id=as.character(hospitalization_id))]
    # Match systemic drugs, never inhaled agents; retain all positive dose levels.
    x <- x[grepl("hydrocortisone|dexamethasone|prednisone|prednisolone|methylprednisolone|betamethasone|cortisoneacetate",drug) &
             route %in% c("iv","intravenous","enteral","oral","po") & is.finite(dose) & dose>0]
    if ("mar_action_category" %in% names(x)) x <- x[!cape_norm(mar_action_category) %in%
      c("held","hold","notgiven","stopped","stop","discontinued","canceled","cancelled","ordered")]
    x[, is_hydrocortisone:=grepl("hydrocortisone",drug)]
    x[, dose_mg:=fcase(unit %in% c("mg","milligram","milligrams"),dose,
                       unit %in% c("mcg","ug"),dose/1000,unit=="g",dose*1000,default=NA_real_)]
    x[, rate_mg_day:=fcase(unit %in% c("mghr","mgh","mgperhour"),dose*24,
                           unit %in% c("mgday","mgd","mg24hr"),dose,
                           unit %in% c("mcgmin","ugmin"),dose*1440/1000,default=NA_real_)]
    x[, qualifies_hydrocortisone:=as.integer(is_hydrocortisone & route %in% c("iv","intravenous") &
      ((source=="continuous" & is.finite(rate_mg_day) & rate_mg_day>=180 & rate_mg_day<=240) |
       (source=="intermittent" & is.finite(dose_mg) & dose_mg>=50 & dose_mg<=100)))]
    if (Sys.getenv("CAPE_HYDRO_MODE","pragmatic_iv") == "continuous_only") x[source_table=="intermittent",qualifies_hydrocortisone:=0L]
    x
  })
  rbindlist(out,fill=TRUE)
}
cape_last <- function(times, values, at, lookback=4) {
  keep <- which(!is.na(times) & times<=at & times>=at-lookback*3600 & is.finite(values))
  if (!length(keep)) return(NA_real_)
  values[keep[which.max(times[keep])]]
}
cape_respiratory_events <- function(events,data_dir,crosswalk,resp,labs,vitals) {
  mode <- Sys.getenv("CAPE_SEVERITY_MODE","pf")
  if (!mode %in% c("pf","sf_proxy")) stop("CAPE_SEVERITY_MODE must be pf or sf_proxy")
  ids <- unique(as.character(crosswalk$hospitalization_id))
  a <- cape_read(file.path(data_dir,"clif_adt.parquet"),c("hospitalization_id","in_dttm","out_dttm","location_category"),ids)
  a[,hospitalization_id:=as.character(hospitalization_id)]
  a <- merge(a,crosswalk[,.(hospitalization_id,hospital_block_id)],by="hospitalization_id")
  a[,`:=`(in_dttm=as.POSIXct(in_dttm,tz="UTC"),out_dttm=as.POSIXct(out_dttm,tz="UTC"))]
  path <- file.path(data_dir,"clif_respiratory_support.parquet")
  pc <- Sys.getenv("CAPE_PEEP_COLUMN","peep_set")
  has_peep <- pc %in% open_dataset(path)$schema$names
  if (!has_peep) stop("PEEP field ",pc," missing. Set CAPE_PEEP_COLUMN to the actual source field; missing PEEP is not assumed normal.")
  p <- cape_read(path,c("hospitalization_id","recorded_dttm",pc),ids)
  p[,hospitalization_id:=as.character(hospitalization_id)]
  p <- merge(p,crosswalk[,.(hospitalization_id,hospital_block_id)],by="hospitalization_id")
  p[,recorded_dttm:=as.POSIXct(recorded_dttm,tz="UTC")]
  po2 <- strsplit(Sys.getenv("CAPE_PAO2_CATEGORIES","po2_arterial;pao2;pao2_arterial;arterial_po2"),";",fixed=TRUE)[[1]]
  events[,`:=`(baseline_spo2=NA_real_,baseline_fio2=NA_real_,baseline_sf_ratio=NA_real_,
                baseline_pf_ratio=NA_real_,baseline_peep=NA_real_,icu_at_time_zero=0L,
                respiratory_initiation_eligible=0L,severity_mode=mode)]
  for (z in list(a,p,resp,labs,vitals)) setkeyv(z,"hospital_block_id")
  # Records at/preceding eligibility only; no future oxygenation or ICU information.
  for (i in seq_len(nrow(events))) {
    id <- events$hospital_block_id[i]; t <- events$eligibility_dttm[i]
    icu <- a[.(id)][cape_norm(location_category)=="icu" & in_dttm<=t & (is.na(out_dttm)|out_dttm>t),.N]>0
    pe <- p[.(id)]
    rr <- resp[.(id)]; ll <- labs[.(id)]; vv <- vitals[.(id)]
    peep <- cape_last(pe$recorded_dttm,as.numeric(pe[[pc]]),t)
    f <- cape_last(rr$recorded_dttm,rr$fio2,t)
    v <- vv[cape_norm(vital_category)=="spo2"]
    sp <- cape_last(v$recorded_dttm,as.numeric(v$vital_value),t)
    l <- ll[tolower(lab_category) %in% tolower(po2)]
    pa <- cape_last(l$lab_collect_dttm,as.numeric(l$lab_value_numeric),t)
    pf <- pa/f; sf <- sp/f
    oxygen <- is.finite(f) && f>=0.5 && if(mode=="pf") is.finite(pf)&&pf<300 else is.finite(sf)&&sf<315&&sp>=80&&sp<=97
    eligible <- icu && if(events$respiratory_support_category[i] %in% c("imv","nippv")) is.finite(peep)&&peep>=5 else oxygen
    set(events,i,c("baseline_spo2","baseline_fio2","baseline_sf_ratio","baseline_pf_ratio","baseline_peep","icu_at_time_zero","respiratory_initiation_eligible"),
        list(sp,f,sf,pf,peep,as.integer(icu),as.integer(eligible)))
  }
  events
}
cape_baseline_exclusions <- function(blocks,data_dir,crosswalk,diagnoses,output_dir) {
  cw <- crosswalk[hospital_block_id %in% blocks$hospital_block_id,.(hospitalization_id,hospital_block_id)]
  med <- cape_steroid_records(data_dir,unique(as.character(cw$hospitalization_id)))
  med <- merge(med,cw,by="hospitalization_id")
  med <- merge(med,blocks[,.(hospital_block_id,admission_dttm,respiratory_failure_start)],by="hospital_block_id")
  prior <- unique(med[admin_dttm>=admission_dttm & admin_dttm<=respiratory_failure_start,hospital_block_id])
  dx <- diagnoses[hospital_block_id %in% blocks$hospital_block_id & poa_present_binary==1L]
  # Observable baseline exclusions; not an assertion that ICD codes recover clinical diagnoses.
  excluded <- unique(dx[grepl("^(U071|J1282|J690|E84|A15|A16|A17|A18|A19|B44|B45|B46|B59)",diagnosis_code_clean),hospital_block_id])
  audit <- blocks[,.(hospital_block_id,patient_id,respiratory_failure_start)]
  audit[,`:=`(exclude_observed_prebaseline_steroid=as.integer(hospital_block_id %in% prior),
                exclude_other_baseline_diagnosis=as.integer(hospital_block_id %in% excluded))]
  write_parquet(audit,file.path(output_dir,"baseline_exclusion_audit.parquet"))
  keep <- audit[exclude_observed_prebaseline_steroid==0L & exclude_other_baseline_diagnosis==0L,hospital_block_id]
  blocks[hospital_block_id %in% keep]
}
