#!/usr/bin/env Rscript
# Development comparison, NOT a calibrated replacement for cumres.clmm.
# Canonical replay uses an isolated install of ff80632 (see evidence README).
# From the repository root:
# Rscript tools/check-clmm-correction.R NEW_OUTPUT [WORKERS [CELL REPLICATE]]
# With CELL/REPLICATE only that sealed case is replayed; otherwise all 1500.
args <- commandArgs(TRUE)
if (!length(args) || !length(args) %in% c(1L,2L,4L) || file.exists(args[1]))
    stop("Supply NEW_OUTPUT [WORKERS [CELL REPLICATE]]")
for (package in c("gof","ordinal","numDeriv"))
    if (!requireNamespace(package,quietly=TRUE)) stop("Select an isolated library containing ",package)
if (utils::packageVersion("gof")<"1.4.0") stop("Select the isolated CLMM development library")
workers <- if (length(args)>=2L) as.integer(args[2]) else 1L
stopifnot(length(workers)==1L,!is.na(workers),workers>=1L)
source("tools/clmm-study.R")
cells <- gof_clmm_study_cells("full")
cells <- cells[match(c("core05","core11","core12"),cells$id),]
jobs <- expand.grid(cell=seq_len(nrow(cells)),replicate=seq_len(500L))
if (length(args)==4L) {
    id <- match(args[3],cells$id);replicate <- as.integer(args[4])
    stopifnot(!is.na(id),!is.na(replicate),replicate %in% 1:500)
    jobs <- jobs[jobs$cell==id & jobs$replicate==replicate,,drop=FALSE]
}
output <- args[1];dir.create(output,recursive=TRUE)
capture.output(sessionInfo(),file=file.path(output,"session.txt"))
saveRDS(list(scope="Sealed candidate comparison; not full qualification",jobs=jobs,
    baseline_commit="ff80632ae542b4de6719282d21df9b144fe3a15c",
    library=find.package("gof"),installed_R_md5=tools::md5sum(
        list.files(file.path(find.package("gof"),"R"),full.names=TRUE))),file.path(output,"invocation.rds"))

## The three sealed designs have four categories, unit weights and ordering x.
## Construct one full-grid contribution per independent block, independently
## of the native multiplier engine. The observed statistic is unchanged.
process <- function(problem,fit) {
    p <- problem;n <- p$n;m <- p$blocks$m;theta <- fit$theta
    stopifnot(p$K==4L,all(p$weights==1))
    o <- order(p$X[,"x"]);x <- p$X[o,"x"];index <- p$blocks$id[o]
    marginal <- getFromNamespace("gof_clmm_marginal","gof")(p,theta)
    residual <- outer(p$y,1:3,"<=")-marginal$cdf
    raw <- matrix(0,m,3*n);D <- matrix(0,3*n,length(theta));observed <- matrix(0,n,3)
    for (k in 1:3) {
        R <- matrix(0,n,m);R[cbind(seq_len(n),index)] <- residual[o,k]
        raw[,seq_len(n)+(k-1L)*n] <- t(apply(R,2,cumsum))
        D[seq_len(n)+(k-1L)*n,] <- apply(marginal$gradient[o+(k-1L)*n,,drop=FALSE],2,cumsum)
        observed[,k] <- cumsum(residual[o,k])/(n/sqrt(m))
    }
    integrator <- getFromNamespace("gof_integrator_nested_agq","gof")(p)
    score <- integrator$score(theta,integrator$refresh(theta))[[1L]]
    M <- raw+t(D %*% solve(fit$information,t(score)))
    decomposition <- qr(score,tol=1e-10)
    if (decomposition$rank!=ncol(score)) stop("Score projection rank deficient")
    Q <- qr.Q(decomposition);h <- rowSums(Q^2)
    if (any(h>=1)) stop("HC2 leverage is at least one")
    projected <- raw-Q %*% crossprod(Q,raw)
    list(observed=observed,x=x,baseline=M,projection=projected,
         HC2=projected/sqrt(1-h),maximum_leverage=max(h))
}

## Finite-grid Gaussian L2 reference via its quadratic-form eigenvalues.
## Check its algebra against the package's independently implemented L2 norm.
quadratic <- function(process,M,seed) {
    n <- nrow(process$observed);m <- nrow(M);dt <- c(diff(process$x),0)
    observed <- sum(rowSums(process$observed^2)*dt)
    weighted <- sweep(M,2,sqrt(rep(dt,3)),"*")/(n/sqrt(m))
    B <- tcrossprod(weighted)
    eigenvalues <- eigen(B,symmetric=TRUE,only.values=TRUE)$values
    stopifnot(min(eigenvalues)> -1e-10*max(eigenvalues))
    eigenvalues <- pmax(eigenvalues,0)
    set.seed(seed)
    simulation <- drop(crossprod(eigenvalues,matrix(rchisq(length(eigenvalues)*9999L,df=1),length(eigenvalues))))
    set.seed(seed+1L);g <- rnorm(m)
    direct <- matrix(drop(crossprod(g,M))/(n/sqrt(m)),n,3)
    native <- getFromNamespace("L2Test","gof")(sqrt(rowSums(direct^2)),process$x)^2
    identity <- drop(crossprod(g,B %*% g))
    stopifnot(isTRUE(all.equal(native,identity,tolerance=1e-10)))
    list(observed_squared=observed,trace=sum(eigenvalues),squared_trace=sum(eigenvalues^2),
         p=mean(simulation>=observed),quantile95=unname(quantile(simulation,.95)),
         quadratic_identity_error=abs(native-identity))
}

attempt <- function(i) {
    cell <- cells[jobs$cell[i],];replicate <- jobs$replicate[i]
    seed <- as.integer(160000000L+match(cell$id,gof_clmm_study_cells("full")$id)*100000L+replicate)
    warnings <- character();stage <- "input"
    result <- tryCatch(withCallingHandlers({
        data <- gof_clmm_study_data(cell,seed)
        formula <- if (cell$structure=="slope") y~x+(1+z|g) else y~x+(1|g)
        model <- ordinal::clmm(formula,data=data,link=cell$link)
        stage <- "polish";p <- getFromNamespace("gof_clmm_problem","gof")(model)
        fit <- getFromNamespace("gof_polish","gof")(p,getFromNamespace("gof_select_integrator","gof")(p))
        stage <- "gaussian_reference";z <- process(p,fit)
        list(status="accepted",baseline=quadratic(z,z$baseline,seed+50000L),
            projection=quadratic(z,z$projection,seed+60000L),HC2=quadratic(z,z$HC2,seed+70000L),
            theta=fit$theta,boundary=fit$boundary,node_check=fit$node_check,
            maximum_leverage=z$maximum_leverage)
    },warning=function(w) {warnings <<- c(warnings,conditionMessage(w));invokeRestart("muffleWarning")}),
    error=function(e) list(status="refused",stage=stage,message=conditionMessage(e)))
    result$id <- cell$id;result$replicate <- replicate;result$seed <- seed;result$warnings <- unique(warnings)
    path <- file.path(output,sprintf("%s-%04d.rds",cell$id,replicate))
    saveRDS(result,paste0(path,".tmp"));stopifnot(file.rename(paste0(path,".tmp"),path))
    invisible(TRUE)
}
for (batch in split(seq_len(nrow(jobs)),ceiling(seq_len(nrow(jobs))/40L))) {
    stopifnot(all(unlist(parallel::mclapply(batch,attempt,mc.cores=workers,mc.preschedule=FALSE,mc.set.seed=FALSE))))
    cat("Completed",max(batch),"/",nrow(jobs),"sealed cases\n");flush.console()
}
