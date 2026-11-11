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
ok<-eda_profile(file.path(tmp,'data'),file.path(tmp,'out'),batch_rows=2,min_cell=1)
stopifnot(ok,file.exists(file.path(tmp,'out','COMPLETE.txt')))
i<-fread(file.path(tmp,'out','table_inventory.csv'))
stopifnot(i[table=='clif_hospitalization',rows]==3,i[table=='clif_vitals',rows]==4)
s<-fread(file.path(tmp,'out','sample_size_summary.csv'))
stopifnot(s[quantity=='unique_patients_in_patient_table',count]==3,
 s[quantity=='unique_patients_in_hospitalization_table',count]==2,
 s[quantity=='unique_hospitalizations',count]==3)
a<-fread(file.path(tmp,'out','table_availability.csv'))
stopifnot(a[table=='clif_vitals',source_hospitalizations_with_records]==2,
 a[table=='clif_vitals',orphan_hospitalization_ids]==1)
d<-fread(file.path(tmp,'out','time_span.csv'))
stopifnot(d[timestamp=='admission_dttm',first_month]=='2025-01',d[timestamp=='admission_dttm',last_month]=='2025-01')
sites<-fread(file.path(tmp,'out','sample_size_by_site.csv'))
stopifnot(sites[site=='JHH',hospitalizations]==2,sites[site=='JHH',unique_patients]==1)
mo<-fread(file.path(tmp,'out','admissions_by_month.csv'))
stopifnot(mo[month=='2025-01',hospitalizations]==3,mo[month=='2025-01',unique_patients]==2)
md<-fread(file.path(tmp,'out','measurement_types_available.csv'))
stopifnot(setequal(md[table=='clif_vitals',available_measurement_type],c('temperature','spo2')))
stopifnot(length(list.files(file.path(tmp,'out','plots'),pattern='[.]png$'))==4L,
 !file.exists(file.path(tmp,'out','numeric_distributions.csv')),
 !file.exists(file.path(tmp,'out','category_frequencies.csv')))
for(p in list.files(file.path(tmp,'out'),pattern='csv$',full.names=TRUE)){
 z<-fread(p,showProgress=FALSE);stopifnot(!any(c('patient_id','hospitalization_id')%in%names(z)))
}
# Small site/month counts are omitted rather than displayed as zero.
ok<-eda_profile(file.path(tmp,'data'),file.path(tmp,'suppressed'),batch_rows=2,min_cell=10)
s<-fread(file.path(tmp,'suppressed','sample_size_by_site.csv'))
stopifnot(all(is.na(s$hospitalizations)),all(is.na(s$unique_patients)))
cat('PASS: exact overall/site/month counts, first-ADT attribution, table coverage, variable inventory, four overview plots and small-cell suppression. Outputs:',file.path(tmp,'out'),'\n')
