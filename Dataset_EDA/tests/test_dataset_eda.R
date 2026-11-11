suppressPackageStartupMessages({library(arrow);library(data.table)})
self<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]);folder<-dirname(dirname(normalizePath(self)))
source(file.path(folder,'eda_rclif_dataset.R'))
tmp<-tempfile('dataset_eda_',tmpdir='/tmp');dir.create(tmp);dir.create(file.path(tmp,'data'))
t<-as.POSIXct('2025-01-01',tz='UTC')
put<-function(x,n)write_parquet(x,file.path(tmp,'data',paste0('clif_',n,'.parquet')))
put(data.table(patient_id=c('p1','p2','p3'),sex_category=c('female','male',''),death_dttm=t+c(NA,86400,NA)),'patient')
put(data.table(hospitalization_id=c('h1','h2','h3'),patient_id=c('p1','p1','p2'),admission_dttm=t+c(0,86400,172800),discharge_dttm=t+c(86400,3*86400,4*86400),age_at_admission=c(30,31,60),admission_type_category='ed'),'hospitalization')
put(data.table(hospitalization_id=c('h1','h2','h3'),hospital_id=c('JHH','JHH','SMH'),in_dttm=t),'adt')
put(data.table(hospitalization_id=c('h1','h1','h2','orphan'),vital_category=c('temperature','temperature','spo2','spo2'),vital_value=c(36,38,95,NA),recorded_dttm=t),'vitals')
put(data.table(hospitalization_id='h1',lab_category=c('glucose','glucose','albumin'),lab_value_unit=c('mg/dL','mg/dL','g/dL'),lab_value_numeric=c(100,120,3),lab_collect_dttm=t),'labs')
# High-cardinality category contents must be omitted rather than misleadingly truncated.
put(data.table(hospitalization_id='h1',diagnosis_code=c('J18.9','J15.9','J13','J12.9')),'hospital_diagnosis')
ok<-eda_profile(file.path(tmp,'data'),file.path(tmp,'out'),batch_rows=2,reservoir_n=100,max_categories=3,min_cell=1)
stopifnot(ok,file.exists(file.path(tmp,'out','COMPLETE.txt')))
i<-fread(file.path(tmp,'out','table_inventory.csv'))
stopifnot(i[table=='clif_hospitalization',rows]==3,i[table=='clif_hospitalization',unique_patients_direct]==2,
 i[table=='clif_vitals',unique_hospitalizations]==3,i[table=='clif_vitals',unique_patients_linked_from_hospitalization]==1,
 i[table=='clif_patient',unique_patients_direct]==3,all(i$parquet_fragments==1))
n<-fread(file.path(tmp,'out','numeric_distributions.csv'))
stopifnot(n[table=='clif_vitals'&category=='temperature',mean]==37,
 n[table=='clif_labs'&category=='glucose',median]==110)
m<-fread(file.path(tmp,'out','column_missingness.csv'))
stopifnot(m[table=='clif_patient'&variable=='death_dttm',n_na]==2,
 m[table=='clif_patient'&variable=='sex_category',n_blank_string]==1)
d<-fread(file.path(tmp,'out','timestamp_ranges.csv'),colClasses='character')
stopifnot(d[table=='clif_hospitalization'&variable=='admission_dttm',maximum_utc]=='2025-01-03 00:00:00')
c<-fread(file.path(tmp,'out','category_frequencies.csv'))
stopifnot(c[table=='clif_hospital_diagnosis',status]=='cardinality limit exceeded; frequencies omitted')
s<-fread(file.path(tmp,'out','hospital_sizes_first_adt.csv'))
stopifnot(s[site=='JHH',hospitalizations]==2,s[site=='JHH',linked_unique_patients]==1)
for(p in list.files(file.path(tmp,'out'),pattern='csv$',full.names=TRUE)) {
 z<-fread(p,showProgress=FALSE)
 stopifnot(!any(c('patient_id','hospitalization_id') %in% names(z)))
}
cat('PASS: exact counts/linkage/ranges/missingness, categorized numeric distributions, category overflow and no identifier columns. Outputs:',file.path(tmp,'out'),'\n')

# Streaming numeric moments and bounded samples span later batches too.
state<-eda_num_state();set.seed(19)
for(z in split(1:20000,rep(1:20,each=1000)))eda_add_numeric(state,z,100)
stopifnot(state$n==20000,length(state$sample)==100,abs(state$mean-mean(1:20000))<1e-8,
 abs(sqrt(state$m2/(state$n-1))-sd(1:20000))<1e-8,any(state$sample>1000))
times<-eda_timestamp(c('2025-01-01','2025-01-01T03:00:00+03:00','2025-01-01T00:00:00Z','invalid'))
stopifnot(all(times[1:3]==as.numeric(t)),is.na(times[4]))
# A cardinality-limit failure must be reported as partial, not as complete data.
partial<-suppressMessages(eda_profile(file.path(tmp,'data'),file.path(tmp,'partial'),batch_rows=2,max_ids=1,min_cell=1))
stopifnot(!partial,file.exists(file.path(tmp,'partial','PARTIAL.txt')),
 !file.exists(file.path(tmp,'partial','COMPLETE.txt')))
cat('PASS: full-stream moments, bounded reservoir, timestamp offsets and partial-run reporting.\n')
