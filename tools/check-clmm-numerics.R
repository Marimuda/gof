#!/usr/bin/env Rscript
# Reproduce the 96-case numerical cost screen, not a calibration study.
# R_LIBS=ISOLATED_REPAIR_LIBRARY Rscript tools/check-clmm-numerics.R NEW_OUTPUT [WORKERS]
args <- commandArgs(TRUE)
source("tools/clmm-study.R")
seed_for <- function(id,replicate)
    as.integer(170000000L+match(id,gof_clmm_study_cells("full")$id)*100000L+replicate)

if (length(args)==5L && args[1]=="--case") {
    output <- args[2];id <- args[3];replicate <- as.integer(args[4])
    .libPaths(c(args[5],.libPaths()))
    if (!"refit" %in% names(formals(gof::gof_clmm_control))) stop("Select a numerical-repair build")
    cell <- subset(gof_clmm_study_cells("full"),id==args[3]);seed <- seed_for(id,replicate)
    warnings <- character();stage <- "input_fit";started <- proc.time()[3]
    cat("stage:",stage,"seed:",seed,"\n");flush.console()
    result <- tryCatch(withCallingHandlers({
        data <- gof_clmm_study_data(cell,seed)
        if (cell$link=="cauchit") {
            probability <- head(cumsum(tabulate(as.integer(data$y),4))/nrow(data),3)
            model <- ordinal::clmm(y~x+(1|g),data=data,link="cauchit",
                start=c(qcauchy(probability),0,log(.5)),nAGQ=-21L)
            control <- gof::gof_clmm_control(refit=TRUE,nodes=21L,max.refresh=20L)
        } else {
            formula <- if (cell$structure=="slope") y~x+(1+z|g) else y~x+(1|g/h)
            model <- ordinal::clmm(formula,data=data)
            control <- gof::gof_clmm_control(refit=TRUE)
        }
        stage <- "diagnostic";cat("stage:",stage,"\n");flush.console()
        set.seed(seed+50000L)
        answer <- gof::cumres(model,variable="x",R=3,plots=0,control=control)
        list(status="accepted",theta=answer$fit_check$theta,fit_check=answer$fit_check,
             finite_output=all(is.finite(c(answer$KS,answer$CvM))))
    },warning=function(w) {warnings <<- c(warnings,conditionMessage(w));invokeRestart("muffleWarning")}),
    error=function(e) list(status="refused",stage=stage,message=conditionMessage(e),
        code=if (inherits(e,"gof_clmm_error")) e$code else "input_fit_failure",condition=e))
    result$id <- id;result$replicate <- replicate;result$seed <- seed;result$warnings <- unique(warnings)
    result$seconds <- unname(proc.time()[3]-started)
    path <- file.path(output,"records",sprintf("%s-%04d.rds",id,replicate))
    saveRDS(result,paste0(path,".tmp"));stopifnot(file.rename(paste0(path,".tmp"),path))
} else {
    if (!length(args) %in% 1:2 || file.exists(args[1])) stop("Supply NEW_OUTPUT [WORKERS]")
    workers <- if (length(args)==2L) as.integer(args[2]) else 1L
    stopifnot(length(workers)==1L,!is.na(workers),workers>=1L)
    for (package in c("gof","ordinal","numDeriv"))
        if (!requireNamespace(package,quietly=TRUE)) stop("Select an isolated library containing ",package)
    if (!"refit" %in% names(formals(gof::gof_clmm_control))) stop("Select a numerical-repair build")
    output <- args[1];dir.create(output,recursive=TRUE);output <- normalizePath(output)
    dir.create(file.path(output,"records"));dir.create(file.path(output,"logs"))
    capture.output(sessionInfo(),file=file.path(output,"session.txt"))
    jobs <- expand.grid(id=c("core08","core14","link5"),replicate=1:32,stringsAsFactors=FALSE)
    saveRDS(list(scope="Numerical cost/availability only; no calibration claim",jobs=jobs,
        timeout_seconds=60L,R=3L,library=find.package("gof")),file.path(output,"invocation.rds"))
    attempt <- function(i) {
        job <- jobs[i,];name <- sprintf("%s-%04d",job$id,job$replicate)
        path <- file.path(output,"records",paste0(name,".rds"))
        command <- c(normalizePath("tools/check-clmm-numerics.R"),"--case",output,
            job$id,as.character(job$replicate),dirname(find.package("gof")))
        code <- suppressWarnings(system2(file.path(R.home("bin"),"Rscript"),vapply(command,shQuote,""),
            stdout=file.path(output,"logs",paste0(name,".log")),
            stderr=file.path(output,"logs",paste0(name,".log")),timeout=60L))
        if (!file.exists(path)) {
            saveRDS(list(id=job$id,replicate=job$replicate,seed=seed_for(job$id,job$replicate),
                status=if (code==124L) "timed_out" else "process_failure",exit_code=code,
                scope="Execution-budget outcome, not a covariance-boundary conclusion"),paste0(path,".tmp"))
            stopifnot(file.rename(paste0(path,".tmp"),path))
        }
        invisible(TRUE)
    }
    for (batch in split(seq_len(nrow(jobs)),ceiling(seq_len(nrow(jobs))/12L))) {
        stopifnot(all(unlist(parallel::mclapply(batch,attempt,mc.cores=workers,mc.preschedule=FALSE,mc.set.seed=FALSE))))
        cat("Completed",max(batch),"/",nrow(jobs),"numerical cases\n");flush.console()
    }
}
