suppressPackageStartupMessages({library(arrow);library(data.table)})
self <- sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])
source(file.path(dirname(dirname(normalizePath(self))),'cape_cod_helpers.R'))
tmp <- tempfile();dir.create(tmp)
x <- data.table(hospitalization_id=as.character(1:6),admin_dttm=as.POSIXct('2025-01-01',tz='UTC'),
 med_category=c('hydrocortisone','hydrocortisone','hydrocortisone','prednisone','hydrocortisone','hydrocortisone'),
 med_dose=c(50,25,50,5,50,50),med_dose_unit='mg',med_route_category=c('iv','iv','iv','enteral','enteral','iv'),
 mar_action_category=c('given','given','held','given','given','ordered'))
y <- data.table(hospitalization_id=c('7','8','9'),admin_dttm=as.POSIXct('2025-01-01',tz='UTC'),
 med_category='hydrocortisone',med_dose=c(200/24,200,5),med_dose_unit=c('mg/hr','mg/day','ml/hr'),mar_action_category='started')
write_parquet(x,file.path(tmp,'clif_medication_admin_intermittent.parquet'))
write_parquet(y,file.path(tmp,'clif_medication_admin_continuous.parquet'))
Sys.setenv(CAPE_HYDRO_MODE='pragmatic_iv')
z <- cape_steroid_records(tmp,as.character(1:9))
stopifnot(setequal(z[qualifies_hydrocortisone==1,hospitalization_id],c('1','7','8')),
          !any(z$hospitalization_id %in% c('3','6')),
          all(c('2','4','5','9') %in% z[qualifies_hydrocortisone==0,hospitalization_id]))
Sys.setenv(CAPE_HYDRO_MODE='continuous_only')
z <- cape_steroid_records(tmp,as.character(1:9))
stopifnot(setequal(z[qualifies_hydrocortisone==1,hospitalization_id],c('7','8')))
unlink(tmp,recursive=TRUE)
cat('PASS: IV/enteral routes, dose and rate units, stopped/held/ordered actions, unsupported-unit deviations and continuous-only mode.\n')
