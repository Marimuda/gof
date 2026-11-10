#!/usr/bin/env Rscript
# Replay the sealed 464-case development screen, never release confirmation.
# Select an isolated candidate build and run from the repository root:
# Rscript tools/check-clmm-inference.R NEW_OUTPUT [WORKERS [CELL REPLICATE]]
args<-commandArgs(TRUE)
if(!length(args) || !length(args)%in%c(1L,2L,4L) || file.exists(args[1]))
    stop('Supply NEW_OUTPUT [WORKERS [CELL REPLICATE]]')
library(gof)
if(!'correction'%in%names(formals(gof_clmm_control)))stop('Select a candidate development build')
for(name in c('gof_clmm_problem','gof_polish','gof_select_integrator','gof_clmm_marginal',
             'gof_clmm_ordering','gof_clmm_correction','gof_clmm_informative',
             'gof_ordinal_engine','gof_ordinal_norm'))assign(name,getFromNamespace(name,'gof'))
source('tools/clmm-study.R')
root<-args[1];dir.create(root,recursive=TRUE)
plan<-jsonlite::fromJSON('inst/extdata/clmm-inference/pilot-plan.json')
plan$workers<-if(length(args)>=2L)as.integer(args[2])else 1L
stopifnot(!is.na(plan$workers),plan$workers>=1L)
capture.output(sessionInfo(),file=file.path(root,'session.txt'))
saveRDS(list(plan=plan,args=args,installed_R_md5=tools::md5sum(
    list.files(file.path(find.package('gof'),'R'),full.names=TRUE))),file.path(root,'invocation.rds'))
cells<-gof_clmm_study_cells('full');out<-file.path(root,'pilot');dir.create(out)
jobs<-expand.grid(cell=1:nrow(cells),rep=1:plan$replicates)
if(length(args)==4L) {
  index<-match(args[3],cells$id);replicate<-as.integer(args[4])
  stopifnot(!is.na(index),!is.na(replicate),replicate%in%seq_len(plan$replicates))
  jobs<-jobs[jobs$cell==index & jobs$rep==replicate,,drop=FALSE]
}
one<-function(i) {
  cell<-cells[jobs$cell[i],];rep<-jobs$rep[i];seed<-as.integer(plan$seed_base+jobs$cell[i]*100000L+rep)
  warnings<-character();stage<-'input'
  value<-tryCatch(withCallingHandlers({
    d<-gof_clmm_study_data(cell,seed)
    formula<-switch(cell$structure,intercept=y~x+(1|g),slope=y~x+(1+z|g),nested=y~x+(1|g/h))
    ctl<-gof_clmm_control()
    if(cell$link=='cauchit') {
      pr<-head(cumsum(tabulate(as.integer(d$y),4))/nrow(d),3)
      model<-ordinal::clmm(formula,data=d,link=cell$link,start=c(qcauchy(pr),0,log(.5)),nAGQ=-21L)
      ctl<-gof_clmm_control(nodes=21L,max.refresh=20L)
    } else model<-ordinal::clmm(formula,data=d,link=cell$link)
    stage<-'polish';p<-gof_clmm_problem(model,ctl);fit<-gof_polish(p,gof_select_integrator(p,ctl),ctl)
    marginal<-gof_clmm_marginal(p,fit$theta,ctl$marginal.nodes)
    r<-(outer(p$y,1:3,'<=')-marginal$cdf)*p$weights
    gradient<-marginal$gradient*rep(p$weights,3)
    ordering<-gof_clmm_ordering(p,'x',fit$theta)
    stage<-'correction'
    methods<-lapply(plan$methods,function(method) {
      c<-gof_clmm_correction(p,fit,r,gradient,method)
      informative<-gof_clmm_informative(r,c$gradient,c$multiplier,ordering)
      set.seed(seed+50000L)
      answer<-gof_ordinal_engine(r,c$gradient,c$multiplier,ordering,gof_ordinal_norm,'omnibus',
        informative,plan$R,0L,'x','clmm',list())
      list(method=method,Sup=answer$KS,L2=answer$CvM,informative=informative,metadata=c$metadata)
    })
    list(status='accepted',methods=methods,boundary=fit$boundary)
  },warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),
  error=function(e)list(status='refused',stage=stage,message=conditionMessage(e),condition=e))
  value$id<-cell$id;value$replicate<-rep;value$seed<-seed;value$warnings<-unique(warnings)
  path<-file.path(out,sprintf('%s-%03d.rds',cell$id,rep));saveRDS(value,paste0(path,'.tmp'));file.rename(paste0(path,'.tmp'),path)
}
for(batch in split(1:nrow(jobs),ceiling((1:nrow(jobs))/24))) {
 stopifnot(all(unlist(parallel::mclapply(batch,one,mc.cores=plan$workers,mc.preschedule=FALSE,mc.set.seed=FALSE))))
 cat('Saved',max(batch),'/',nrow(jobs),'\n');flush.console()
}
x<-lapply(list.files(out,full.names=TRUE,pattern='[.]rds$'),readRDS)
s<-do.call(rbind,lapply(cells$id,function(id)do.call(rbind,lapply(plan$methods,function(method) {
 group<-x[vapply(x,function(z)z$id==id,TRUE)];ok<-group[vapply(group,function(z)z$status=='accepted',TRUE)]
 m<-lapply(ok,function(z)z$methods[[match(method,plan$methods)]])
 data.frame(id=id,method=method,attempted=length(group),accepted=length(ok),
 Sup=if(length(m))mean(vapply(m,`[[`,0,'Sup')<=.05)else NA,
 L2=if(length(m))mean(vapply(m,`[[`,0,'L2')<=.05)else NA)
}))))
write.csv(s,file.path(root,'pilot-summary.csv'),row.names=FALSE)
