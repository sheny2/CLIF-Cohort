#!/usr/bin/env Rscript
# Basic overview of the delivered CLIF extract. No analytic cohort filters.
Sys.setenv(ARROW_NUM_THREADS='1',OMP_NUM_THREADS='1')
suppressPackageStartupMessages({library(arrow);library(data.table);library(dplyr);library(ggplot2)})
arrow::set_cpu_count(1L);data.table::setDTthreads(1L)
eda_blank <- function(x) is.na(x)|(is.character(x)&!is.na(x)&trimws(x)=='')
eda_id_state <- function(){e<-new.env(parent=emptyenv());e$values<-new.env(hash=TRUE,parent=emptyenv());e$n<-0L;e}
eda_add_ids <- function(e,x,limit){
 for(v in unique(as.character(x[!eda_blank(x)])))if(!exists(v,e$values,inherits=FALSE)){
  if(e$n>=limit)stop('EDA_MAX_UNIQUE_IDS exceeded; raise only if memory permits.')
  assign(v,TRUE,e$values);e$n<-e$n+1L
 };invisible(NULL)
}
eda_write <- function(x,path){x<-as.data.table(x);if(!ncol(x))x<-data.table(table=character());fwrite(x,path,na='')}
eda_bind <- function(x)if(length(x))rbindlist(x,fill=TRUE)else data.table()
eda_scan <- function(path,columns,callback,batch_rows){
 ds<-open_dataset(path);q<-ds |> select(all_of(columns))
 r<-Scanner$create(q,use_threads=FALSE,batch_size=batch_rows)$ToRecordBatchReader()
 on.exit(r$Close(),add=TRUE)
 repeat{b<-r$read_next_batch();if(is.null(b))break;callback(as.data.table(as.data.frame(b)))}
 invisible(NULL)
}
eda_timestamp <- function(x) {
 if(inherits(x,'POSIXt')) return(as.numeric(x))
 if(inherits(x,'Date')) return(as.numeric(as.POSIXct(x,tz='UTC')))
 if(is.character(x)) {
  x<-gsub('T',' ',x,fixed=TRUE);x<-sub('Z$','+0000',x)
  x<-sub('([+-][0-9]{2}):([0-9]{2})$','\\1\\2',x)
  out<-rep(NA_real_,length(x));present<-which(!is.na(x)&nzchar(trimws(x)))
  zoned<-present[grepl('[+-][0-9]{4}$',x[present])]
  dateonly<-present[grepl('^[0-9]{4}-[0-9]{2}-[0-9]{2}$',x[present])]
  naive<-setdiff(present,c(zoned,dateonly))
  if(length(zoned))out[zoned]<-as.numeric(as.POSIXct(strptime(x[zoned],'%Y-%m-%d %H:%M:%OS%z',tz='UTC')))
  if(length(dateonly))out[dateonly]<-as.numeric(as.POSIXct(x[dateonly],tz='UTC',format='%Y-%m-%d'))
  if(length(naive))out[naive]<-as.numeric(as.POSIXct(x[naive],tz='UTC',format='%Y-%m-%d %H:%M:%OS'))
  return(out)
 }
 rep(NA_real_,length(x))
}
eda_profile <- function(data_dir,output_dir,batch_rows=65536L,max_ids=2000000L,min_cell=10L){
 stopifnot(dir.exists(data_dir),is.finite(batch_rows),batch_rows>0,is.finite(max_ids),max_ids>0,is.finite(min_cell),min_cell>=1)
 if(dir.exists(output_dir))stop('OUTPUT_DIR exists; choose a fresh directory.')
 paths<-list.files(data_dir,pattern='^clif_.*[.]parquet$',full.names=TRUE)
 if(!length(paths))stop('No clif_*.parquet tables found.')
 names(paths)<-sub('[.]parquet$','',basename(paths))
 if(!'clif_hospitalization'%in%names(paths))stop('clif_hospitalization.parquet is required for sample sizes and coverage.')
 dir.create(file.path(output_dir,'plots'),recursive=TRUE)
 dict<-list();inventory<-list();coverage<-list();measurements<-list();errors<-list()
 # Metadata/row counts only: do not load numeric measurements or medication doses.
 for(tab in names(paths))tryCatch({
  message('Inventory: ',tab)
  ds<-open_dataset(paths[[tab]]);cols<-ds$schema$names
  types<-vapply(ds$schema$fields,function(f)f$type$ToString(),character(1))
  dict[[tab]]<-data.table(table=tab,variable=cols,arrow_type=types)
  count<-ds |> summarise(rows=n()) |> collect()
  inventory[[tab]]<-data.table(table=tab,rows=as.numeric(count$rows),columns=length(cols),status='complete')
 },error=function(e){errors[[tab]]<<-data.table(table=tab,error=conditionMessage(e));inventory[[tab]]<<-data.table(table=tab,rows=NA_real_,columns=NA_integer_,status='failed')})
 eda_write(eda_bind(dict),file.path(output_dir,'variable_dictionary.csv'))
 eda_write(eda_bind(inventory),file.path(output_dir,'table_inventory.csv'))
 hp<-paths[['clif_hospitalization']];hc<-open_dataset(hp)$schema$names
 required<-c('hospitalization_id','patient_id','admission_dttm','discharge_dttm')
 if(!all(required%in%hc))stop('Hospitalization table missing: ',paste(setdiff(required,hc),collapse=', '))
 site_col<-intersect(c('admission_hospital_id','hospital_id','site_id','admission_hospital_name','hospital_name'),hc)
 site_col<-if(length(site_col))site_col[1]else character()
 chunks<-list();k<-0L
 eda_scan(hp,c(required,site_col),function(x){k<<-k+1L;chunks[[k]]<<-x},batch_rows)
 h<-if(length(chunks))rbindlist(chunks)else as.data.table(setNames(rep(list(logical()),length(c(required,site_col))),c(required,site_col)));rm(chunks)
 if(nrow(h)>max_ids)stop('Hospitalization row safety limit exceeded; increase EDA_MAX_UNIQUE_IDS if appropriate.')
 h[,`:=`(hospitalization_id=as.character(hospitalization_id),patient_id=as.character(patient_id),
          admission_seconds=eda_timestamp(admission_dttm),discharge_seconds=eda_timestamp(discharge_dttm))]
 missing_hospital_ids<-sum(eda_blank(h$hospitalization_id));missing_patient_ids<-sum(eda_blank(h$patient_id))
 hospital_patients<-unique(h$patient_id[!eda_blank(h$patient_id)])
 raw_hospital_rows<-nrow(h)
 # Date ranges are coarsened to months in the export; exact timestamps stay in memory.
 spans<-rbindlist(lapply(c('admission','discharge'),function(v){
  z<-h[[paste0(v,'_seconds')]];ok<-is.finite(z)
  data.table(timestamp=paste0(v,'_dttm'),n_with_valid_date=sum(ok),n_missing_or_invalid=sum(!ok),
   first_month=if(any(ok))format(as.POSIXct(min(z[ok]),origin='1970-01-01',tz='UTC'),'%Y-%m')else NA_character_,
   last_month=if(any(ok))format(as.POSIXct(max(z[ok]),origin='1970-01-01',tz='UTC'),'%Y-%m')else NA_character_)
 }))
 eda_write(spans,file.path(output_dir,'time_span.csv'))
 conflicts<-h[!eda_blank(hospitalization_id)&!eda_blank(patient_id),.(n=uniqueN(patient_id)),by=hospitalization_id][n>1,hospitalization_id]
 duplicate_hospital_rows<-h[!eda_blank(hospitalization_id),.N]-uniqueN(h$hospitalization_id[!eda_blank(h$hospitalization_id)])
 setorder(h,hospitalization_id,admission_seconds,na.last=TRUE)
 h<-h[!eda_blank(hospitalization_id)&!duplicated(hospitalization_id)]
 h[hospitalization_id%in%conflicts,patient_id:=NA_character_]
 h[,site:=if(length(site_col))as.character(get(site_col))else NA_character_]
 # Earliest observed ADT site is the preferred attribution, with source-field fallback.
 attribution<-'hospitalization field or Unknown'
 if('clif_adt'%in%names(paths))tryCatch({
  ac<-open_dataset(paths[['clif_adt']])$schema$names
  if(all(c('hospitalization_id','hospital_id','in_dttm')%in%ac)){
   first<-new.env(hash=TRUE,parent=emptyenv());nfirst<-0L
   eda_scan(paths[['clif_adt']],c('hospitalization_id','hospital_id','in_dttm'),function(x){
    x[,`:=`(hospitalization_id=as.character(hospitalization_id),seconds=eda_timestamp(in_dttm))]
    x<-x[hospitalization_id%in%h$hospitalization_id&!eda_blank(hospital_id)&is.finite(seconds)]
    setorder(x,hospitalization_id,seconds);x<-x[!duplicated(hospitalization_id)]
    for(i in seq_len(nrow(x))){id<-x$hospitalization_id[i];z<-list(time=x$seconds[i],site=as.character(x$hospital_id[i]))
     if(!exists(id,first,inherits=FALSE)){assign(id,z,first);nfirst<<-nfirst+1L}
     else if(z$time<get(id,first)$time)assign(id,z,first)
    }
   },batch_rows)
   assigned<-vapply(h$hospitalization_id,function(id)if(exists(id,first,inherits=FALSE))get(id,first)$site else NA_character_,character(1))
   h[!is.na(assigned),site:=assigned[!is.na(assigned)]]
   attribution<-'earliest valid recorded ADT hospital, with hospitalization-field fallback'
  }
 },error=function(e){errors[['site_attribution']]<<-data.table(table='site_attribution',error=conditionMessage(e))})
 h[is.na(site)|trimws(site)=='',site:='Unknown']
 patient_table_ids<-eda_id_state()
 # Scan identifier columns only for table-level availability; category LABELS only
 # identify measurement types. No lab/vital values or drug doses are read.
 for(tab in names(paths))tryCatch({
  message('Availability: ',tab)
  ds<-open_dataset(paths[[tab]]);cols<-ds$schema$names
  idcols<-intersect(c('patient_id','hospitalization_id'),cols)
  typecol<-intersect(c('lab_category','vital_category','assessment_category'),cols)
  typecol<-if(length(typecol))typecol[1]else character()
  if(!length(c(idcols,typecol)))return_value<-NULL
  else{
   pset<-eda_id_state();hset<-eda_id_state();typeset<-eda_id_state()
   eda_scan(paths[[tab]],c(idcols,typecol),function(x){
    if('patient_id'%in%idcols)eda_add_ids(pset,x$patient_id,max_ids)
    if('hospitalization_id'%in%idcols)eda_add_ids(hset,x$hospitalization_id,max_ids)
    if(length(typecol))eda_add_ids(typeset,x[[typecol]],max_ids)
   },batch_rows)
   if(tab=='clif_patient')patient_table_ids<-pset
   represented<-ls(hset$values,all.names=TRUE);known<-intersect(represented,h$hospitalization_id)
   linked<-h[hospitalization_id%in%known,uniqueN(patient_id[!is.na(patient_id)])]
   coverage[[tab]]<-data.table(table=tab,
    unique_patients_direct=if('patient_id'%in%idcols)pset$n else NA_real_,
    unique_hospitalizations_direct=if('hospitalization_id'%in%idcols)hset$n else NA_real_,
    source_hospitalizations_with_records=if('hospitalization_id'%in%idcols)length(known)else NA_real_,
    hospitalization_coverage_percent=if('hospitalization_id'%in%idcols&&nrow(h)>0)100*length(known)/nrow(h)else NA_real_,
    linked_patients_with_records=if('hospitalization_id'%in%idcols)linked else NA_real_,
    orphan_hospitalization_ids=if('hospitalization_id'%in%idcols)length(setdiff(represented,h$hospitalization_id))else NA_real_)
   if(length(typecol))measurements[[tab]]<-data.table(table=tab,category_field=typecol,
      available_measurement_type=ls(typeset$values,all.names=TRUE))
  }
 },error=function(e){errors[[paste0(tab,'_coverage')]]<<-data.table(table=tab,error=conditionMessage(e))})
 summary<-data.table(quantity=c('unique_patients_in_patient_table','unique_patients_in_hospitalization_table',
   'unique_hospitalizations','source_hospitalization_rows','duplicate_hospitalization_rows',
   'hospitalization_rows_missing_id','hospitalization_rows_missing_patient_id','hospitalization_ids_with_conflicting_patient_links'),
  count=c(if('clif_patient'%in%names(coverage))patient_table_ids$n else NA_real_,length(hospital_patients),nrow(h),raw_hospital_rows,
          duplicate_hospital_rows,missing_hospital_ids,missing_patient_ids,length(conflicts)))
 eda_write(summary,file.path(output_dir,'sample_size_summary.csv'))
 eda_write(eda_bind(coverage),file.path(output_dir,'table_availability.csv'))
 eda_write(eda_bind(measurements),file.path(output_dir,'measurement_types_available.csv'))
 site<-h[,.(hospitalizations=.N,unique_patients=uniqueN(patient_id[!is.na(patient_id)]),
             admissions_with_unknown_patient=sum(is.na(patient_id))),by=site]
 site[,`:=`(hospitalizations=ifelse(hospitalizations<min_cell,NA_real_,hospitalizations),
             unique_patients=ifelse(unique_patients<min_cell,NA_real_,unique_patients))]
 eda_write(site,file.path(output_dir,'sample_size_by_site.csv'))
 calendar<-h[is.finite(admission_seconds),.(site,patient_id,
    month=format(as.POSIXct(admission_seconds,origin='1970-01-01',tz='UTC'),'%Y-%m'))]
 bytime<-calendar[,.(hospitalizations=.N,unique_patients=uniqueN(patient_id[!is.na(patient_id)])),by=.(site,month)]
 alltime<-calendar[,.(hospitalizations=.N,unique_patients=uniqueN(patient_id[!is.na(patient_id)])),by=month]
 bytime[,`:=`(hospitalizations=ifelse(hospitalizations<min_cell,NA_real_,hospitalizations),
               unique_patients=ifelse(unique_patients<min_cell,NA_real_,unique_patients))]
 alltime[,`:=`(hospitalizations=ifelse(hospitalizations<min_cell,NA_real_,hospitalizations),
                unique_patients=ifelse(unique_patients<min_cell,NA_real_,unique_patients))]
 eda_write(bytime,file.path(output_dir,'admissions_by_site_and_month.csv'))
 eda_write(alltime,file.path(output_dir,'admissions_by_month.csv'))
 theme_set(theme_minimal(base_size=12))
 sites_long<-melt(site,id.vars='site',measure.vars=c('hospitalizations','unique_patients'),variable.name='measure',value.name='n')
 sites_long[,measure:=factor(measure,levels=c('hospitalizations','unique_patients'),labels=c('Hospitalizations','Unique patients'))]
 if(nrow(sites_long)&&any(is.finite(sites_long$n))){
  pp<-ggplot(sites_long,aes(site,n))+geom_col(fill='#2878A8',na.rm=TRUE)+
   geom_text(aes(label=format(n,big.mark=',',trim=TRUE)),vjust=-.3,na.rm=TRUE)+
   scale_y_continuous(expand=expansion(mult=c(0,.12)))+facet_wrap(~measure,scales='free_y')+
   labs(title='Sample size by hospital',subtitle='Patient counts are unique within each site; patients may appear at multiple sites.',x='Hospital',y='Count')
  ggsave(file.path(output_dir,'plots','sample_size_by_site.png'),pp,width=10,height=5,dpi=160)
 }
 if(nrow(alltime)&&any(is.finite(alltime$hospitalizations))){
  monthly_plot<-rbindlist(list(data.table(site='All sites',alltime),bytime),fill=TRUE)
  monthly_plot[,date:=as.Date(paste0(month,'-01'))]
  pp<-ggplot(monthly_plot,aes(date,hospitalizations))+geom_col(fill='#2878A8',na.rm=TRUE)+facet_wrap(~site,scales='free_y',ncol=2)+
   labs(title='Admissions across the observed time span',subtitle='Unique hospitalization IDs per month; small counts omitted.',x='Admission month (UTC)',y='Hospitalizations')
  ggsave(file.path(output_dir,'plots','admissions_over_time.png'),pp,width=12,height=max(5,2.7*ceiling(uniqueN(monthly_plot$site)/2)),dpi=160)
 }
 inv<-eda_bind(inventory)
 if(any(is.finite(inv$rows)&inv$rows>0)){
  pp<-ggplot(inv[is.finite(rows)&rows>0],aes(reorder(table,rows),rows))+geom_col(fill='#4A8C61')+coord_flip()+scale_y_log10()+
   labs(title='Source table sizes',subtitle='Repeated measurement rows are not independent patients.',x=NULL,y='Rows (log scale)')
  ggsave(file.path(output_dir,'plots','table_sizes.png'),pp,width=10,height=max(5,.3*nrow(inv)),dpi=160)
 }
 cv<-eda_bind(coverage)
 if(nrow(cv)&&any(is.finite(cv$hospitalization_coverage_percent))){
  pp<-ggplot(cv[is.finite(hospitalization_coverage_percent)],aes(reorder(table,hospitalization_coverage_percent),hospitalization_coverage_percent))+
   geom_col(fill='#756BB1')+coord_flip()+labs(title='Availability of table records',x=NULL,y='% of source hospitalizations with at least one record')
  ggsave(file.path(output_dir,'plots','table_availability.png'),pp,width=10,height=max(5,.3*nrow(cv)),dpi=160)
 }
 eda_write(if(length(errors))eda_bind(errors)else data.table(table=character(),error=character()),file.path(output_dir,'table_errors.csv'))
 writeLines(c('# Basic rclif_2_1 overview',
  'Start with sample_size_summary.csv, sample_size_by_site.csv, admissions_by_month.csv and time_span.csv.',
  'The four overview plots show site size, admissions over time, table row counts and record availability.',
  'variable_dictionary.csv lists every column/type. measurement_types_available.csv lists lab/vital/assessment labels only, without individual values or distributions.',
  'No numeric vital/lab distributions, medication-dose plots, exact timestamps, free text or identifier values are exported.',
  paste('Source:',normalizePath(data_dir)),paste('Site attribution:',attribution),paste('Small-cell threshold:',min_cell),
  'Patients can have multiple hospitalizations and can contribute to more than one site/month. Do not add site/month patient counts to infer unique overall patients.',
  'Table coverage means at least one row exists for an admission; it does not establish complete measurements or follow-up.',
  'Observed date ranges are exported by month; they do not establish continuous coverage or a data refresh date.',
  'These are working summaries of the delivered filtered extract, not the whole Hopkins EHR. Outputs require normal project disclosure review.'),file.path(output_dir,'START_HERE.md'))
 writeLines(c(paste('created_utc:',format(Sys.time(),tz='UTC')),paste('batch_rows:',batch_rows),
  paste('max_unique_ids:',max_ids),paste('tables_with_errors:',length(errors))),file.path(output_dir,'run_notes.txt'))
 marker<-if(length(errors))'PARTIAL.txt'else'COMPLETE.txt'
 writeLines(if(length(errors))'Inspect table_errors.csv: some summaries are incomplete.'else'Basic overview completed.',file.path(output_dir,marker))
 invisible(!length(errors))
}
if(sys.nframe()==0L){
 a<-commandArgs(trailingOnly=TRUE);if(!length(a))stop('Usage: Rscript Dataset_EDA/eda_rclif_dataset.R DATA_DIR [NEW_OUTPUT_DIR]')
 self<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]);folder<-dirname(normalizePath(self))
 out<-if(length(a)>=2L)a[2]else file.path(folder,'results',paste0('basic_run_',format(Sys.time(),'%Y%m%d_%H%M%S')))
 ok<-eda_profile(a[1],out,batch_rows=as.integer(Sys.getenv('EDA_BATCH_ROWS','65536')),
  max_ids=as.integer(Sys.getenv('EDA_MAX_UNIQUE_IDS','2000000')),min_cell=as.integer(Sys.getenv('SUMMARY_SMALL_CELL_THRESHOLD','10')))
 cat('Basic EDA written to ',out,'\n',sep='');if(!ok)quit(status=1L)
}
