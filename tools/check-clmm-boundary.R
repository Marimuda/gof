#!/usr/bin/env Rscript
# Independent constrained likelihood profiles for two retained failures.
# Diagnostic evidence only; neither global optimality nor boundary p-values.
# Run from repository root: Rscript tools/check-clmm-boundary.R NEW_OUTPUT
args<-commandArgs(TRUE)
if(length(args)!=1L || file.exists(args[1]))stop('Supply one new output directory')
source('tools/clmm-study.R');source('tools/clmm-reference.R')
output<-args[1];dir.create(output,recursive=TRUE)
results<-list()
for(case in c('core11','core14')) {
  cells<-gof_clmm_study_cells('full');cell<-cells[cells$id==case,]
  replicate<-if(case=='core11')19L else 30L
  seed<-140000000L+match(case,cells$id)*100000L+replicate
  data<-gof_clmm_study_data(cell,seed)
  start<-reference_start(cell)
  reference<-reference_fit(data,cell,list(start),Q=31L,refine=61L)
  phi<-reference$theta;j<-length(phi)
  lower<-c(rep(-Inf,4),if(case=='core11')c(0,-Inf,0) else c(0,0))
  rows<-list()
  for(tau in c(0,.0025,.01,.04)) {
    expand<-function(nuisance){par<-phi;par[-j]<-nuisance;par[j]<-sqrt(tau);par}
    objective<-function(par){v<-reference_loglik(data,cell,expand(par),31L);if(is.finite(v))-v else 1e100}
    fit<-nlminb(phi[-j],objective,lower=lower[-j],control=list(iter.max=300L,eval.max=1000L,rel.tol=1e-11))
    direct61<-reference_loglik(data,cell,expand(fit$par),61L)
    rows[[length(rows)+1L]]<-data.frame(id=case,seed=seed,conditional_variance=tau,
      convergence=fit$convergence,loglik31=-fit$objective,loglik61=direct61,
      quadrature_difference=direct61+fit$objective)
    saveRDS(list(reference=reference,profiles=rows),file.path(output,paste0(case,'.rds')))
  }
  tab<-do.call(rbind,rows);tab$profile_change<-tab$loglik61-tab$loglik61[1]
  results[[case]]<-tab;print(tab);flush.console()
}
write.csv(do.call(rbind,results),file.path(output,'profiles.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(output,'session.txt'))
