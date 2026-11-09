#!/usr/bin/env Rscript
# Reproduce the asymmetric input-likelihood audit without changing any package.
# This is a deterministic software check, not a calibration study.
args <- commandArgs(TRUE)
if (length(args)!=1L) stop("Supply one new output directory")
output <- args[1]
if (file.exists(output)) stop("Refusing to overwrite existing audit output")
for (package in c("gof","ordinal","numDeriv"))
    if (!requireNamespace(package,quietly=TRUE)) stop("Install development dependency: ",package)
if (utils::packageVersion("gof")<"1.4.0") stop("Select an isolated gof >= 1.4.0 library")
dir.create(output,recursive=TRUE)
problem <- getFromNamespace("gof_clmm_problem","gof")
integrator <- getFromNamespace("gof_integrator_nested_agq","gof")
F_link <- function(link) if (link=="cloglog") function(x) -expm1(-exp(x)) else function(x) exp(-exp(-x))
make_data <- function(link) {
    set.seed(3)
    m <- 80L; nc <- 8L; n <- m*nc
    x <- rnorm(n); z <- rnorm(n); group <- rep(seq_len(m),each=nc)
    eta <- .8*x+rnorm(m)[group]
    cdf <- F_link(link)(outer(-eta,c(-1.5,0,1.5),"+"))
    y <- 1L+rowSums(runif(n)>cdf)
    data.frame(y=ordered(y,levels=1:4),x=x,z=z,g=factor(group))
}
theta_input <- function(fit) c(fit$alpha,fit$beta,log(fit$ST$g[1,1]))
direct_loglik <- function(data,theta,link) {
    alpha <- c(-Inf,theta[1:3],Inf)
    beta <- theta[4]; sd <- exp(theta[5]); F <- F_link(link)
    sum(vapply(split(seq_len(nrow(data)),data$g),function(rows) {
        y <- as.integer(data$y[rows]); eta <- beta*data$x[rows]
        kernel <- function(b) vapply(b,function(value) {
            probability <- F(alpha[y+1L]-eta-value)-F(alpha[y]-eta-value)
            prod(probability)*dnorm(value,sd=sd)
        },0)
        log(integrate(kernel,-Inf,Inf,rel.tol=1e-10,abs.tol=1e-300,subdivisions=500L)$value)
    },0))
}
warnings <- list()
capture <- function(label,action) withCallingHandlers(action(),warning=function(w) {
    warnings[[label]] <<- c(warnings[[label]],conditionMessage(w));invokeRestart("muffleWarning")
})
records <- list()
for (link in c("cloglog","loglog")) {
    data <- make_data(link)
    scalar <- capture(paste0(link,"-scalar"),function()
        ordinal::clmm(y~x+(1|g),data=data,link=link,nAGQ=21L))
    p <- problem(scalar)
    ig <- integrator(p,21L)
    gof_at_input <- sum(ig$loglik(p$theta,ig$refresh(p$theta)))
    declared <- direct_loglik(data,theta_input(scalar),link)
    opposite <- direct_loglik(data,theta_input(scalar),if (link=="cloglog") "loglog" else "cloglog")
    set.seed(1940)
    diagnostic <- capture(paste0(link,"-gof"),function()
        gof::cumres(scalar,variable="x",R=3L,plots=0L,
                   control=gof::gof_clmm_control(nodes=21L,max.refresh=10L)))
    polished_direct <- direct_loglik(data,diagnostic$fit_check$theta,link)
    matrix <- capture(paste0(link,"-matrix"),function()
        ordinal::clmm(y~x+(1|g),data=data,link=link,
                     control=ordinal::clmm.control(useMatrix=TRUE,gradTol=1e-10)))
    mp <- problem(matrix); mi <- integrator(mp,1L)
    matrix_laplace <- sum(mi$loglik(mp$theta,mi$refresh(mp$theta)))
    stopifnot(abs(gof_at_input-declared)<1e-7,
              abs(diagnostic$fit_check$loglik-polished_direct)<1e-7,
              abs(matrix_laplace-as.numeric(logLik(matrix)))<1e-7)
    records[[link]] <- data.frame(link=link,
        ordinal_version=as.character(utils::packageVersion("ordinal")),
        reported_scalar_loglik=as.numeric(logLik(scalar)),
        direct_declared_loglik=declared,direct_opposite_loglik=opposite,
        gof_at_input=gof_at_input,gof_at_polished=diagnostic$fit_check$loglik,
        direct_at_polished=polished_direct,
        max_polishing_shift_se=max(abs(diagnostic$fit_check$shift$shift_se)),
        reported_matrix_loglik=as.numeric(logLik(matrix)),gof_matrix_laplace=matrix_laplace,
        scalar_matches_declared=abs(as.numeric(logLik(scalar))-declared)<1e-7,
        scalar_matches_opposite=abs(as.numeric(logLik(scalar))-opposite)<1e-7)
}
table <- do.call(rbind,records)
write.csv(table,file.path(output,"likelihood-checks.csv"),row.names=FALSE)
saveRDS(list(scope="Deterministic input-likelihood check, not calibration evidence",
    results=table,warnings=warnings,seed=3L,m=80L,nc=8L,
    package_versions=vapply(c("gof","ordinal","numDeriv"),function(x)
        as.character(utils::packageVersion(x)),"")),file.path(output,"audit.rds"))
capture.output(sessionInfo(),file=file.path(output,"session.txt"))
print(table)
cat("Declared-link scalar integration and matrix-backend consistency checks passed.\n")
