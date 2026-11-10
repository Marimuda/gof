#!/usr/bin/env Rscript
# Audit and export a completed development screen; never overwrite evidence.
# Rscript tools/export-clmm-inference.R AUDIT_DIRECTORY NEW_OUTPUT
args<-commandArgs(TRUE);stopifnot(length(args)==2L)
root<-normalizePath(args[1],mustWork=TRUE);out<-args[2]
if(file.exists(out))stop('Choose a new output directory')
dir.create(out,recursive=TRUE)
stopifnot(file.copy(file.path(root,'pilot-plan.json'),file.path(out,'pilot-plan.json')),
 file.copy(file.path(root,'boundary-profiles','profiles.csv'),file.path(out,'boundary-profiles.csv')))
source('tools/clmm-study.R')
plan<-jsonlite::fromJSON(file.path(root,'pilot-plan.json'));cells<-gof_clmm_study_cells('full')
files<-list.files(file.path(root,'pilot'),full.names=TRUE,pattern='[.]rds$');stopifnot(length(files)==464L)
x<-lapply(files,readRDS)
keys<-vapply(x,function(z)paste(z$id,z$replicate),'');expected<-expand.grid(id=cells$id,replicate=1:16)
stopifnot(!anyDuplicated(keys),identical(sort(keys),sort(paste(expected$id,expected$replicate))))
field<-function(z,k,default)if(is.null(z[[k]]))default else z[[k]]
a<-do.call(rbind,lapply(x,function(z) {
 stopifnot(z$seed==180000000L+match(z$id,cells$id)*100000L+z$replicate)
 do.call(rbind,lapply(plan$methods,function(method) {
  m<-if(z$status=='accepted')z$methods[[match(method,plan$methods)]]else NULL
  if(!is.null(m))stopifnot(m$method==method,isTRUE(m$informative),all(is.finite(c(m$Sup,m$L2))))
  data.frame(id=z$id,replicate=z$replicate,seed=z$seed,method=method,status=z$status,
    stage=field(z,'stage',''),failure_code=field(z$condition,'code',''),
    condition_stage=field(z$condition,'stage',''),message=field(z,'message',''),Sup=if(is.null(m))NA else m$Sup,
    L2=if(is.null(m))NA else m$L2,warnings=paste(z$warnings,collapse=' | '))
 }))
}))
rows<-list()
for(id in cells$id)for(method in plan$methods)for(statistic in c('Sup','L2')) {
 group<-a[a$id==id & a$method==method,]
 values<-group[[statistic]][group$status=='accepted'];n<-length(values);k<-sum(values<=.05)
 ci<-if(n)binom.test(k,n)$conf.int else c(NA,NA)
 rows[[length(rows)+1L]]<-data.frame(id=id,method=method,statistic=statistic,
   attempted=nrow(group),accepted=n,rejections=k,rate=if(n)k/n else NA,
   lower95=ci[1],upper95=ci[2],refusal=1-n/nrow(group))
}
s<-do.call(rbind,rows)
stopifnot(nrow(a)==1392L,nrow(s)==174L)
con<-gzfile(file.path(out,'attempts.csv.gz'),'wt');write.csv(a,con,row.names=FALSE);close(con)
write.csv(s,file.path(out,'summary.csv'),row.names=FALSE)
write.csv(cells,file.path(out,'designs.csv'),row.names=FALSE)
raw_hash<-vapply(files,function(p)digest::digest(file=p,algo='sha256'),'')
write.csv(data.frame(file=basename(files),sha256=raw_hash),file.path(out,'raw-records-sha256.csv'),row.names=FALSE)
source_files<-names(plan$source_sha256)
stopifnot(all(vapply(source_files,function(p)digest::digest(file=p,algo='sha256'),'')==unlist(plan$source_sha256)))
metadata<-list(scope='Development screen only;16 attempts per design cannot qualify size or refusal rates',
 fitting='Default interior fitter, refit=FALSE; Cauchit finite-start nonadaptive input and nodes=21 explicitly predeclared',
 attempts=464L,accepted=sum(vapply(x,function(z)z$status=='accepted',TRUE)),
 failure_stages=as.list(table(vapply(x,function(z)field(z,'stage','accepted'),''))),
 failure_codes=as.list(table(vapply(x,function(z)field(z$condition,'code','accepted'),''))),
 correction_failure_semantics='A correction-stage error refuses the common comparison; inspect stage rather than counting it as a fitting failure',
 source='pilot-plan.json contains the pre-outcome numerical source hashes',
 tooling='An extra closing parenthesis affected only the original post-run summary. Atomic numerical records were retained; this exporter independently reconstructs every endpoint from those records.')
metadata$original_runner_sha256<-digest::digest(file=file.path(root,'pilot-original.R'),algo='sha256')
metadata$exporter_sha256<-digest::digest(file='tools/export-clmm-inference.R',algo='sha256')
names<-c('pilot-plan.json','boundary-profiles.csv','attempts.csv.gz','summary.csv','designs.csv','raw-records-sha256.csv')
metadata$sha256<-as.list(setNames(vapply(file.path(out,names),function(p)digest::digest(file=p,algo='sha256'),''),names))
jsonlite::write_json(metadata,file.path(out,'PROVENANCE.json'),auto_unbox=TRUE,pretty=TRUE)
print(metadata[c('attempts','accepted','failure_stages')]);print(subset(s,id%in%c('core02','core05','core11','core12','link5')))
