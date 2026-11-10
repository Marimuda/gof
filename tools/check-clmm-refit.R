#!/usr/bin/env Rscript
# Deterministic numerical checks, not calibration or release qualification.
# Select an isolated development gof library, then run from the package root:
# Rscript tools/check-clmm-refit.R NEW_OUTPUT_DIRECTORY
args <- commandArgs(TRUE)
if (length(args)!=1L || file.exists(args[1])) stop("Supply one NEW output directory")
for (p in c("gof","ordinal","numDeriv"))
    if (!requireNamespace(p,quietly=TRUE)) stop("Select a development library containing ",p)
if (utils::packageVersion("gof")<"1.4.0") stop("Select an isolated development gof >= 1.4.0 library")
if (!"refit" %in% names(formals(gof::gof_clmm_control))) stop("Select the numerical-repair development build")
source("tools/clmm-study.R")
output <- args[1];dir.create(output,recursive=TRUE)
records <- list(); warnings <- list()
capture <- function(label,expr) withCallingHandlers(expr,warning=function(w) {
    warnings[[label]] <<- c(warnings[[label]],conditionMessage(w));invokeRestart("muffleWarning")
})

## Direct Gaussian-prior integration; neither adaptive nodes nor gof's
## likelihood, score, marginal probabilities or information are used here.
rule <- function(Q) {
    J <- matrix(0,Q,Q);J[cbind(1:(Q-1L),2:Q)]<-sqrt(1:(Q-1L));J<-J+t(J)
    e<-eigen(J,symmetric=TRUE);o<-order(e$values)
    list(z=e$values[o],w=e$vectors[1,o]^2)
}
direct_slope <- function(data,theta,Q) {
    g<-rule(Q);index<-expand.grid(a=seq_len(Q),b=seq_len(Q))
    stopifnot(length(theta)==7L)
    L<-matrix(c(exp(theta[5]),theta[6],0,exp(theta[7])),2,2)
    b0<-L[1,1]*g$z[index$a]
    b1<-L[2,1]*g$z[index$a]+L[2,2]*g$z[index$b]
    w<-g$w[index$a]*g$w[index$b];d<-data
    eta<-outer(theta[4]*d$x,rep(1,length(w)))+outer(rep(1,nrow(d)),b0)+outer(d$z,b1)
    a<-c(-Inf,theta[1:3],Inf);y<-as.integer(d$y)
    probability<-plogis(a[y+1L]-eta)-plogis(a[y]-eta)
    lp<-rowsum(log(probability),d$g,reorder=FALSE)
    sum(apply(lp,1,function(x){mx<-max(x);mx+log(sum(exp(x-mx)*w))}))
}
direct_scalar <- function(data,theta) {
    alpha<-c(-Inf,theta[1:3],Inf);sd<-exp(theta[5])
    sum(vapply(split(seq_len(nrow(data)),data$g),function(rows) {
        y<-as.integer(data$y[rows]);eta<-theta[4]*data$x[rows]
        kernel<-function(b)vapply(b,function(z)
            prod(pcauchy(alpha[y+1L]-eta-z)-pcauchy(alpha[y]-eta-z))*dnorm(z,sd=sd),0)
        log(integrate(kernel,-Inf,Inf,rel.tol=1e-10,abs.tol=1e-300,subdivisions=500L)$value)
    },0))
}

cell<-subset(gof_clmm_study_cells("full"),id=="core08")
d<-gof_clmm_study_data(cell,140800004L)
model<-capture("slope-input",ordinal::clmm(y~x+(1+z|g),data=d))
model$ST[[1L]][2L,2L]<-0 # Platform-independent singular initialization.
set.seed(2741)
answer<-capture("slope-refit",gof::cumres(model,variable="x",R=3,plots=0,
    control=gof::gof_clmm_control(refit=TRUE)))
direct51<-direct_slope(d,answer$fit_check$theta,51L)
direct81<-direct_slope(d,answer$fit_check$theta,81L)
stopifnot(abs(direct51-direct81)<1e-6,abs(direct81-answer$fit_check$loglik)<1e-4)
records[["slope"]]<-list(model=model,answer=answer,independent51=direct51,independent81=direct81)
table<-data.frame(case="sparse slope",seed=140800004L,reported=answer$fit_check$loglik,
                  independent=direct81,error=abs(direct81-answer$fit_check$loglik))

cell<-subset(gof_clmm_study_cells("full"),id=="link5")
for (replicate in c(2L,5L,7L)) {
    seed<-142400000L+replicate;d<-gof_clmm_study_data(cell,seed)
    probability<-head(cumsum(tabulate(as.integer(d$y),4))/nrow(d),3)
    # A finite marginal start and a non-adaptive input rule avoid the
    # failing Laplace mode search on these fixtures. This is a different
    # input workflow from the sealed original calibration cohorts.
    model<-capture(paste0("cauchit-input-",replicate),ordinal::clmm(y~x+(1|g),data=d,
        link="cauchit",start=c(qcauchy(probability),0,log(.5)),nAGQ=-21L))
    set.seed(2741L+replicate)
    answer<-capture(paste0("cauchit-diagnostic-",replicate),gof::cumres(model,variable="x",R=3,plots=0,
        control=gof::gof_clmm_control(nodes=21L,max.refresh=20L)))
    direct<-direct_scalar(d,answer$fit_check$theta)
    stopifnot(abs(direct-answer$fit_check$loglik)<1e-4)
    records[[paste0("cauchit",replicate)]]<-list(model=model,answer=answer,independent=direct)
    table<-rbind(table,data.frame(case=paste("cauchit",replicate),seed=seed,
        reported=answer$fit_check$loglik,independent=direct,error=abs(direct-answer$fit_check$loglik)))
}
saveRDS(list(records=records,warnings=warnings),file.path(output,"numerical-checks.rds"))
write.csv(table,file.path(output,"numerical-checks.csv"),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(output,"session.txt"))
writeLines("Numerical development checks only. Reused failed fixtures, not independent calibration. Original failed cohorts and release hold remain unchanged.",
           file.path(output,"scope.txt"))
print(table)
