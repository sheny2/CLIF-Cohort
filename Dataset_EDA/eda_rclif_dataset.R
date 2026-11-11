#!/usr/bin/env Rscript
# Profile every CLIF Parquet table in the delivered extract, without cohort filters.
# No individual records or patient/encounter ID values are exported.
Sys.setenv(ARROW_NUM_THREADS='1',OMP_NUM_THREADS='1')
suppressPackageStartupMessages({library(arrow);library(data.table);library(dplyr);library(ggplot2)})
arrow::set_cpu_count(1L);data.table::setDTthreads(1L)

eda_blank <- function(x) is.na(x) | (is.character(x) & !is.na(x) & trimws(x)=='')
eda_id_state <- function() {e<-new.env(parent=emptyenv());e$values<-new.env(hash=TRUE,parent=emptyenv());e$n<-0L;e}
eda_add_ids <- function(e,x,limit) {
 x<-unique(as.character(x[!eda_blank(x)]))
 for(v in x) if(!exists(v,e$values,inherits=FALSE)) {
  if(e$n>=limit) stop('Exact-ID limit exceeded; raise EDA_MAX_UNIQUE_IDS if memory permits.')
  assign(v,TRUE,e$values);e$n<-e$n+1L
 }
 invisible(NULL)
}
eda_num_state <- function() {e<-new.env(parent=emptyenv());e$n<-0;e$mean<-0;e$m2<-0;e$min<-Inf;e$max<--Inf;e$sample<-numeric();e$priority<-numeric();e}
eda_add_numeric <- function(e,x,k) {
 x<-as.numeric(x);x<-x[is.finite(x)];n<-length(x)
 if(!n) return(invisible(NULL))
 mu<-mean(x);delta<-mu-e$mean;total<-e$n+n
 e$m2<-e$m2+sum((x-mu)^2)+delta^2*e$n*n/total
 e$mean<-e$mean+delta*n/total;e$n<-total;e$min<-min(e$min,min(x));e$max<-max(e$max,max(x))
 # Uniform bounded priority reservoir across ALL batches, not just the first rows.
 priority<-runif(n);take<-order(priority)[seq_len(min(n,k))]
 values<-c(e$sample,x[take]);priority<-c(e$priority,priority[take]);keep<-order(priority)[seq_len(min(length(values),k))]
 e$sample<-values[keep];e$priority<-priority[keep]
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
eda_has_time <- function(name,type) grepl('timestamp|date32|date64',type,ignore.case=TRUE)||grepl('(_dttm|_datetime|_date)$',name)
eda_is_id <- function(name) grepl('(^|_)(patient|hospitalization|hospital_block|encounter|visit|record|row)_id$|^id$',name)
eda_is_category <- function(name) !eda_is_id(name) && (grepl('(_category|_group|_code|_format|_unit|_units)$',name)||name %in% c('hospital_id','site_id','admission_hospital_id','poa_present','diagnosis_primary','tracheostomy'))
eda_is_numeric <- function(name,type) !eda_is_id(name)&&!eda_has_time(name,type)&&grepl('int|float|double|decimal|bool',type,ignore.case=TRUE)
eda_write <- function(x,path) {x<-as.data.table(x);if(!ncol(x))x<-data.table(table=character());fwrite(x,path,na='')}
eda_bind <- function(x) if(length(x)) rbindlist(x,fill=TRUE) else data.table()

eda_profile <- function(data_dir,output_dir,batch_rows=65536L,reservoir_n=2000L,
                        max_categories=1000L,max_groups=2000L,max_ids=2000000L,min_cell=10L,seed=20261110L) {
 stopifnot(dir.exists(data_dir),batch_rows>0,reservoir_n>0,max_categories>0,max_groups>0,max_ids>0,min_cell>=1)
 if(dir.exists(output_dir)) stop('OUTPUT_DIR already exists; choose a new directory to avoid mixing runs.')
 dir.create(output_dir,recursive=TRUE);dir.create(file.path(output_dir,'plots'))
 set.seed(seed)
 paths<-list.files(data_dir,pattern='^clif_.*[.]parquet$',full.names=TRUE)
 if(!length(paths)) stop('No clif_*.parquet files or directories were found in ',data_dir)
 # Hospitalization first permits linked patient coverage for child tables.
 paths<-paths[order(basename(paths)!='clif_hospitalization.parquet',basename(paths))]
 roster<-new.env(hash=TRUE,parent=emptyenv());roster_conflicts<-eda_id_state();roster_complete<-FALSE
 totals<-list();schema_all<-list();missing_all<-list();date_all<-list();calendar_all<-list()
 categories_all<-list();numeric_all<-list();hist_all<-list();domains_all<-list();errors<-list();links_all<-list()
 put_outputs<-function() {
  for(nm in c('totals','schema_all','missing_all','date_all','calendar_all','categories_all','numeric_all','hist_all','domains_all','links_all','errors')) {
   filename<-c(totals='table_inventory.csv',schema_all='variable_dictionary.csv',missing_all='column_missingness.csv',
    date_all='timestamp_ranges.csv',calendar_all='timestamp_counts_by_month.csv',categories_all='category_frequencies.csv',
    numeric_all='numeric_distributions.csv',hist_all='sampled_numeric_histograms.csv',domains_all='measurement_groups.csv',
    links_all='table_linkage_coverage.csv',errors='table_errors.csv')[[nm]]
   z<-eda_bind(get(nm,inherits=TRUE))
   if(nm=='errors'&&!ncol(z))z<-data.table(table=character(),error=character())
   eda_write(z,file.path(output_dir,filename))
  }
 }
 for(path in paths) {
  table<-sub('[.]parquet$','',basename(path));message('Profiling ',table,' ...')
  tryCatch({
   ds<-open_dataset(path);cols<-ds$schema$names
   types<-vapply(ds$schema$fields,function(f)f$type$ToString(),character(1))
   names(types)<-cols
   schema_all[[table]]<-data.table(table=table,variable=cols,arrow_type=types,
    role=ifelse(vapply(cols,eda_is_id,logical(1)),'identifier',ifelse(mapply(eda_has_time,cols,types),'timestamp',
     ifelse(vapply(cols,eda_is_category,logical(1)),'coded/category',ifelse(mapply(eda_is_numeric,cols,types),'numeric','other/text')))))
   ids<-lapply(intersect(c('patient_id','hospitalization_id'),cols),function(v)eda_id_state());names(ids)<-intersect(c('patient_id','hospitalization_id'),cols)
   categories<-new.env(hash=TRUE,parent=emptyenv());numbers<-new.env(hash=TRUE,parent=emptyenv())
   numeric_meta<-list();num_count<-0L;groups_omitted<-0L
   times<-list();monthly<-list();missing<-setNames(numeric(length(cols)),cols);blanks<-missing;nonfinite<-missing
   row_n<-0;n_batch<-0;hospital_orphans<-eda_id_state();patient_orphans<-NULL
   time_cols<-cols[mapply(eda_has_time,cols,types)]
   cat_cols<-cols[vapply(cols,eda_is_category,logical(1))]
   num_cols<-cols[mapply(eda_is_numeric,cols,types)]
   for(v in cat_cols) {e<-new.env(parent=emptyenv());e$counts<-integer();e$overflow<-FALSE;assign(v,e,categories)}
   # Domain-specific numeric profiles separate laboratory/vital/medication types and units.
   category_col<-intersect(c('lab_category','vital_category','med_category','assessment_category'),cols)
   category_col<-if(length(category_col)) category_col[1] else NULL
   unit_col<-intersect(c('lab_value_unit','lab_unit','vital_unit','vital_value_unit','med_dose_unit','unit'),cols)
   unit_col<-if(length(unit_col)) unit_col[1] else NULL
   domain_values<-intersect(c('lab_value_numeric','vital_value','med_dose','numerical_value'),num_cols)
   reader<-Scanner$create(ds,use_threads=FALSE,batch_size=batch_rows)$ToRecordBatchReader()
   repeat {
    b<-reader$read_next_batch();if(is.null(b)) break
    x<-as.data.table(as.data.frame(b))
    for(v in names(x)[vapply(x,is.factor,logical(1))]) set(x,j=v,value=as.character(x[[v]]))
    row_n<-row_n+nrow(x);n_batch<-n_batch+1L
    for(v in cols) {
     z<-x[[v]];missing[v]<-missing[v]+sum(is.na(z))
     if(is.character(z)) blanks[v]<-blanks[v]+sum(!is.na(z)&trimws(z)=='')
     if(is.numeric(z)) nonfinite[v]<-nonfinite[v]+sum(!is.na(z)&!is.finite(z))
    }
    for(v in names(ids)) eda_add_ids(ids[[v]],x[[v]],max_ids)
    if(table=='clif_hospitalization' && all(c('patient_id','hospitalization_id') %in% cols)) {
     pairs<-unique(x[!eda_blank(patient_id)&!eda_blank(hospitalization_id),.(hospitalization_id=as.character(hospitalization_id),patient_id=as.character(patient_id))])
     for(j in seq_len(nrow(pairs))) {
      h<-pairs$hospitalization_id[j];p<-pairs$patient_id[j]
      if(exists(h,roster,inherits=FALSE)&&!identical(get(h,roster),p)) eda_add_ids(roster_conflicts,h,max_ids)
      else assign(h,p,roster)
     }
    }
    for(v in time_cols) {
     t<-eda_timestamp(x[[v]]);ok<-is.finite(t)
     if(is.null(times[[v]])) times[[v]]<-c(n=0,invalid=0,min=Inf,max=-Inf)
     times[[v]]['invalid']<-times[[v]]['invalid']+sum(!eda_blank(x[[v]])&!ok)
     if(any(ok)) {
      times[[v]]['n']<-times[[v]]['n']+sum(ok)
      times[[v]]['min']<-min(times[[v]]['min'],min(t[ok]));times[[v]]['max']<-max(times[[v]]['max'],max(t[ok]))
      mo<-format(as.POSIXct(t[ok],origin='1970-01-01',tz='UTC'),'%Y-%m')
      cc<-data.table(month=mo)[,.N,by=month]
      monthly[[v]]<-rbindlist(list(monthly[[v]],cc))[,.(N=sum(N)),by=month]
     }
    }
    for(v in cat_cols) {
     e<-get(v,categories);if(e$overflow) next
     z<-x[[v]];cc<-data.table(value=as.character(z[!eda_blank(z)]))[,.N,by=value]
     values<-union(names(e$counts),cc$value)
     if(length(values)>max_categories) {e$overflow<-TRUE;e$counts<-integer();next}
     old<-setNames(rep(0,length(values)),values);old[names(e$counts)]<-e$counts
     old[cc$value]<-old[cc$value]+cc$N;e$counts<-old
    }
    for(v in num_cols) {
     groups<-if(!is.null(category_col)&&v %in% domain_values) as.character(x[[category_col]]) else rep('ALL',nrow(x))
     units<-if(!is.null(unit_col)&&v %in% domain_values) as.character(x[[unit_col]]) else rep('not supplied/applicable',nrow(x))
     groups[is.na(groups)|groups=='']<-'Unknown category';units[is.na(units)|units=='']<-'Unknown unit'
     chunks<-split(seq_len(nrow(x)),paste(groups,units,sep='\034'))
     for(ii in chunks) {
      group<-groups[ii[1]];unit<-units[ii[1]];key<-paste(v,group,unit,sep='\035')
      if(!exists(key,numbers,inherits=FALSE)) {
       if(num_count>=max_groups) {groups_omitted<-groups_omitted+1L;next}
       num_count<-num_count+1L;assign(key,eda_num_state(),numbers)
       numeric_meta[[key]]<-list(variable=v,category=group,unit=unit)
      }
      eda_add_numeric(get(key,numbers),x[[v]][ii],reservoir_n)
     }
    }
    rm(x,b);if(n_batch%%25L==0L) {message(table,': ',format(row_n,big.mark=','),' rows');invisible(gc(verbose=FALSE))}
   }
   reader$Close()
   # Derived linkage is exact for IDs represented in the supplied hospitalization roster.
   if(table=='clif_hospitalization') roster_complete<-TRUE
   linked<-eda_id_state()
   hosp_ids<-if('hospitalization_id' %in% names(ids)) ls(ids$hospitalization_id$values,all.names=TRUE) else character()
   if(roster_complete) for(h in hosp_ids) {
    if(exists(h,roster_conflicts$values,inherits=FALSE)) next
    if(exists(h,roster,inherits=FALSE)) eda_add_ids(linked,get(h,roster),max_ids) else eda_add_ids(hospital_orphans,h,max_ids)
   }
   totals[[table]]<-data.table(table=table,rows=row_n,columns=length(cols),parquet_fragments=if(dir.exists(path))length(list.files(path,pattern='[.]parquet$',recursive=TRUE))else 1L,
    unique_patients_direct=if('patient_id' %in% names(ids))ids$patient_id$n else NA_real_,
    unique_hospitalizations=if('hospitalization_id' %in% names(ids))ids$hospitalization_id$n else NA_real_,
    unique_patients_linked_from_hospitalization=if(length(hosp_ids)&&roster_complete)linked$n else NA_real_,
    numeric_group_limit_reached=groups_omitted>0,status='complete')
   links_all[[table]]<-data.table(table=table,unique_hospitalization_ids_absent_from_roster=if(roster_complete)hospital_orphans$n else NA_real_,roster_complete=roster_complete,
    roster_conflicting_patient_links=roster_conflicts$n)
   missing_all[[table]]<-data.table(table=table,variable=cols,rows=row_n,n_na=missing,n_blank_string=blanks,
    missing_or_blank_percent=if(row_n)100*(missing+blanks)/row_n else NA_real_,n_nonfinite_nonmissing=nonfinite)
   date_all[[table]]<-eda_bind(lapply(names(times),function(v) {
    t<-times[[v]];data.table(table=table,variable=v,n_valid=t['n'],n_unparseable=t['invalid'],
     minimum_utc=if(t['n']>0)format(as.POSIXct(t['min'],origin='1970-01-01',tz='UTC'),'%Y-%m-%d %H:%M:%S')else NA_character_,
     maximum_utc=if(t['n']>0)format(as.POSIXct(t['max'],origin='1970-01-01',tz='UTC'),'%Y-%m-%d %H:%M:%S')else NA_character_)
   }))
   calendar_all[[table]]<-eda_bind(lapply(names(monthly),function(v) data.table(table=table,variable=v,monthly[[v]])))
   categories_all[[table]]<-eda_bind(lapply(cat_cols,function(v) {
    e<-get(v,categories)
    if(e$overflow) return(data.table(table=table,variable=v,value=NA_character_,n=NA_real_,percent_of_table=NA_real_,status='cardinality limit exceeded; frequencies omitted'))
    cc<-data.table(value=names(e$counts),n=as.numeric(e$counts));setorder(cc,-n,value)
    cc<-head(cc,50L);cc[,status:=ifelse(n<min_cell,'small cell suppressed','exact top frequency')]
    cc[n<min_cell,`:=`(value=NA_character_,n=NA_real_)]
    cc[,`:=`(table=table,variable=v,percent_of_table=100*n/max(1,row_n))];cc
   }))
   numeric_all[[table]]<-eda_bind(lapply(names(numeric_meta),function(key) {
    e<-get(key,numbers);meta<-numeric_meta[[key]];q<-if(e$n>=min_cell)quantile(e$sample,c(.01,.05,.25,.5,.75,.95,.99))else rep(NA_real_,7)
    data.table(table=table,variable=meta$variable,category=meta$category,unit=meta$unit,n_finite=e$n,
     mean=if(e$n>=min_cell)e$mean else NA_real_,sd=if(e$n>=max(2,min_cell))sqrt(e$m2/(e$n-1))else NA_real_,
     minimum=if(e$n>=min_cell)e$min else NA_real_,maximum=if(e$n>=min_cell)e$max else NA_real_,
     p01=q[1],p05=q[2],p25=q[3],median=q[4],p75=q[5],p95=q[6],p99=q[7],sample_n=length(e$sample),
     quantiles=if(e$n<=reservoir_n)'exact (all finite values retained)'else'uniform reservoir estimate')
   }))
   domains_all[[table]]<-copy(numeric_all[[table]])[,intersect(c('table','variable','category','unit','n_finite'),names(numeric_all[[table]])),with=FALSE]
   hist_rows<-list()
   # Only the largest 12 numeric groups per table are plotted/exported as histograms.
   keys<-names(numeric_meta);if(length(keys))keys<-keys[order(vapply(keys,function(k)get(k,numbers)$n,numeric(1)),decreasing=TRUE)]
   for(key in head(keys,12)) {
    e<-get(key,numbers);meta<-numeric_meta[[key]]
    if(length(e$sample)<max(2,min_cell))next
    h<-hist(e$sample,breaks=30,plot=FALSE)
    hist_rows[[key]]<-data.table(table=table,variable=meta$variable,category=meta$category,unit=meta$unit,
     left=head(h$breaks,-1),right=tail(h$breaks,-1),sample_count=ifelse(h$counts>0&h$counts<min_cell,NA_real_,h$counts),
     sample_n=length(e$sample),population_finite=e$n)
   }
   hist_all[[table]]<-eda_bind(hist_rows)
   if(nrow(hist_all[[table]])) {
    hh<-copy(hist_all[[table]]);hh[,panel:=paste(variable,category,unit,sep=' | ')]
    pp<-ggplot(hh,aes(x=(left+right)/2,y=sample_count))+geom_col(aes(width=right-left),fill='#2878A8',na.rm=TRUE)+
     facet_wrap(~panel,scales='free',ncol=3)+theme_minimal(base_size=10)+
     labs(title=paste(table,'numeric distributions'),subtitle='Uniform reservoir samples; panels may have different units. Small histogram cells omitted.',
          x='Value',y='Sampled observations (not patient counts)')
    ggsave(file.path(output_dir,'plots',paste0(table,'_numeric.png')),pp,width=14,height=max(4,3*ceiling(uniqueN(hh$panel)/3)),dpi=150,limitsize=FALSE)
   }
   message('Completed ',table,': ',format(row_n,big.mark=','),' rows.')
   rm(ds,reader,ids,categories,numbers,linked);invisible(gc(verbose=FALSE))
  },error=function(e) {
   message('FAILED ',table,': ',conditionMessage(e))
   errors[[table]]<<-data.table(table=table,error=conditionMessage(e))
   totals[[table]]<<-data.table(table=table,rows=NA_real_,columns=NA_integer_,unique_patients_direct=NA_real_,unique_hospitalizations=NA_real_,unique_patients_linked_from_hospitalization=NA_real_,status='failed')
  })
  put_outputs()
 }
 # Exact encounter-level derived LOS (numeric distribution needs no large full-table read).
 # Available directly from hospitalization timestamps; handled in a separate streamed pass.
 hp<-paths[basename(paths)=='clif_hospitalization.parquet']
 if(length(hp)) tryCatch({
  ds<-open_dataset(hp)
  if(all(c('admission_dttm','discharge_dttm') %in% ds$schema$names)) {
   rr<-Scanner$create(ds |> select(admission_dttm,discharge_dttm),use_threads=FALSE,batch_size=batch_rows)$ToRecordBatchReader()
   state<-eda_num_state();negative<-0
   repeat {b<-rr$read_next_batch();if(is.null(b))break;x<-as.data.frame(b)
    los<-(eda_timestamp(x$discharge_dttm)-eda_timestamp(x$admission_dttm))/86400
    negative<-negative+sum(is.finite(los)&los<0);eda_add_numeric(state,los[is.finite(los)&los>=0],reservoir_n)
   }
   rr$Close();q<-if(length(state$sample))quantile(state$sample,c(.25,.5,.75,.95))else rep(NA_real_,4)
   eda_write(data.table(valid_nonnegative_encounter_rows=state$n,negative_los_rows=negative,
    mean_days=if(state$n)state$mean else NA_real_,q25_days=q[1],median_days=q[2],q75_days=q[3],q95_days=q[4],
    sample_n=length(state$sample),unit='days; encounter rows, not merged hospital blocks'),file.path(output_dir,'hospital_length_of_stay.csv'))
  }
 },error=function(e){errors[['length_of_stay']]<<-data.table(table='length_of_stay',error=conditionMessage(e));put_outputs()})
 inventory<-eda_bind(totals);dates<-eda_bind(date_all);monthly_counts<-eda_bind(calendar_all)
 eda_write(data.table(quantity=c('unique_patients_in_patient_table','patient_rows','unique_patients_in_hospitalization_table','unique_hospitalizations','hospitalization_rows'),
  value=c(inventory[table=='clif_patient',unique_patients_direct][1],inventory[table=='clif_patient',rows][1],
          inventory[table=='clif_hospitalization',unique_patients_direct][1],inventory[table=='clif_hospitalization',unique_hospitalizations][1],
          inventory[table=='clif_hospitalization',rows][1])),file.path(output_dir,'sample_size_summary.csv'))
 if(nrow(monthly_counts)) {
  mm<-monthly_counts[table %in% c('clif_hospitalization','clif_patient') & variable %in% c('admission_dttm','discharge_dttm','death_dttm')]
  if(nrow(mm)) {
   mm[,date:=as.Date(paste0(month,'-01'))]
   pp<-ggplot(mm,aes(date,N))+geom_col(fill='#2878A8')+facet_wrap(~variable,scales='free_y',ncol=1)+theme_minimal()+
    labs(title='Observed records by calendar month',x='Month (UTC)',y='Timestamped rows',subtitle='Record dates do not establish complete coverage or death ascertainment.')
   ggsave(file.path(output_dir,'plots','calendar_coverage.png'),pp,width=12,height=8,dpi=150)
  }
 }
 # Attribute encounters to the earliest recorded ADT hospital, not to every
 # location visited; linked patient denominators use unambiguous roster links.
 ap<-paths[basename(paths)=='clif_adt.parquet']
 if(length(ap)) tryCatch({
  ds<-open_dataset(ap)
  if(all(c('hospitalization_id','hospital_id','in_dttm') %in% ds$schema$names)) {
   rr<-Scanner$create(ds |> select(hospitalization_id,hospital_id,in_dttm),use_threads=FALSE,batch_size=batch_rows)$ToRecordBatchReader()
   first<-new.env(hash=TRUE,parent=emptyenv());first_n<-0L
   repeat {
    b<-rr$read_next_batch();if(is.null(b))break
    x<-as.data.table(as.data.frame(b));x[,seconds:=eda_timestamp(in_dttm)]
    x<-x[!eda_blank(hospitalization_id)&!eda_blank(hospital_id)&is.finite(seconds)]
    x[,hospitalization_id:=as.character(hospitalization_id)];setorder(x,hospitalization_id,seconds)
    x<-x[!duplicated(hospitalization_id)]
    for(i in seq_len(nrow(x))) {
     h<-x$hospitalization_id[i];z<-list(time=x$seconds[i],site=as.character(x$hospital_id[i]))
     if(!exists(h,first,inherits=FALSE)) {if(first_n>=max_ids)stop('ADT exact-ID limit exceeded');assign(h,z,first);first_n<-first_n+1L}
     else if(z$time<get(h,first)$time) assign(h,z,first)
    }
   }
   rr$Close()
   h<-ls(first,all.names=TRUE)
   if(length(h)) {
    ss<-data.table(site=vapply(h,function(id)get(id,first)$site,character(1)),
      patient=vapply(h,function(id)if(roster_complete&&exists(id,roster,inherits=FALSE)&&!exists(id,roster_conflicts$values,inherits=FALSE))get(id,roster)else NA_character_,character(1)))
    sites<-ss[,.(hospitalizations=.N,linked_unique_patients=if(roster_complete)uniqueN(patient[!is.na(patient)])else NA_integer_,
                  hospitalizations_without_unambiguous_patient_link=sum(is.na(patient))),by=site]
    eda_write(sites,file.path(output_dir,'hospital_sizes_first_adt.csv'))
   }
  }
 },error=function(e){errors[['hospital_sizes']]<<-data.table(table='hospital_sizes',error=conditionMessage(e));put_outputs()})
 writeLines(c('# Whole-extract CLIF EDA',
  'Start with sample_size_summary.csv, table_inventory.csv, timestamp_ranges.csv and variable_dictionary.csv.',
  'These summaries describe the delivered extract BEFORE any analytic cohort filters.',
  'category_frequencies.csv counts source rows, not independent patients. hospital_sizes_first_adt.csv attributes each encounter to its earliest recorded ADT hospital.',
  'column_missingness.csv reports NA/blank cells; missing death timestamps do not by themselves mean unknown vital status.',
  'numeric_distributions.csv separates measurement categories and units, with exact finite-value moments and sampled quantiles.',
  'Quantile/histogram samples are uniform reservoirs drawn across the entire scan; they are not first-row samples.',
  paste('Tables successfully profiled:',sum(inventory$status=='complete')),paste('Reported errors:',length(errors)),
  'All outputs should remain on Discovery until disclosure review. No raw identifiers or free-text contents are exported.'),file.path(output_dir,'START_HERE.md'))
 writeLines(c(paste('source:',normalizePath(data_dir)),paste('created_utc:',format(Sys.time(),tz='UTC')),
  paste('batch_rows:',batch_rows),paste('numeric_reservoir_per_group:',reservoir_n),paste('seed:',seed),
  'Row counts, ID cardinalities, missingness, timestamp ranges and retained category counts are exact for completed scans.',
  'Numeric means/SD/min/max use all finite values; quantiles/histograms use bounded uniform reservoir samples.',
  'Outputs describe the delivered filtered extract, NOT the complete Hopkins EHR or an analytic cohort.',
  'Admission/discharge/measurement/death date ranges do not establish complete longitudinal follow-up.',
  'Identifier values, free text and row-level samples are never exported. All outputs remain sensitive pending review.',
  'Exact distinct-ID tracking uses memory proportional to distinct IDs; scan batches and numeric reservoirs are bounded.',
  paste('tables_with_errors:',length(errors))),file.path(output_dir,'run_notes.txt'))
 if(length(errors)) {writeLines('Some tables failed. Inspect table_errors.csv before interpreting coverage.',file.path(output_dir,'PARTIAL.txt'))
  return(invisible(FALSE))}
 writeLines('All discovered CLIF tables profiled.',file.path(output_dir,'COMPLETE.txt'))
 invisible(TRUE)
}

if(sys.nframe()==0L) {
 a<-commandArgs(trailingOnly=TRUE)
 if(!length(a))stop('Usage: Rscript Dataset_EDA/eda_rclif_dataset.R DATA_DIR [NEW_OUTPUT_DIR]')
 self<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]);folder<-dirname(normalizePath(self))
 out<-if(length(a)>=2L)a[2]else file.path(folder,'results',paste0('run_',format(Sys.time(),'%Y%m%d_%H%M%S')))
 env_int<-function(n,default)as.integer(Sys.getenv(n,as.character(default)))
 ok<-eda_profile(a[1],out,batch_rows=env_int('EDA_BATCH_ROWS',65536),reservoir_n=env_int('EDA_RESERVOIR_N',2000),
  max_categories=env_int('EDA_MAX_CATEGORIES',1000),max_groups=env_int('EDA_MAX_NUMERIC_GROUPS',2000),
  max_ids=env_int('EDA_MAX_UNIQUE_IDS',2000000),min_cell=env_int('SUMMARY_SMALL_CELL_THRESHOLD',10))
 cat('Dataset EDA written to ',out,'\n',sep='')
 if(!ok)quit(status=1L)
}
