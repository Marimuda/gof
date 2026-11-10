#!/usr/bin/env Rscript
# Audit compact development evidence. Never overwrite an existing export.
# Usage: Rscript tools/export-clmm-repair-evidence.R REPAIR_DIRECTORY NEW_OUTPUT
args <- commandArgs(TRUE)
if (length(args)!=2L || file.exists(args[2])) stop("Supply REPAIR_DIRECTORY NEW_OUTPUT")
for (package in c("jsonlite","digest"))
    if (!requireNamespace(package,quietly=TRUE)) stop("Install development dependency: ",package)
root <- normalizePath(args[1],mustWork=TRUE);output <- args[2]
source("tools/clmm-study.R")
field <- function(x,name,default) if (is.null(x[[name]])) default else x[[name]]
load <- function(directory,count) {
    paths <- list.files(file.path(root,directory),pattern="[.]rds$",full.names=TRUE)
    stopifnot(length(paths)==count)
    lapply(paths,readRDS)
}
prospective <- load("prospective-screen",1500L)
numerical <- load("numerical-screen",96L)
check_roster <- function(records,ids,n,base) {
    keys <- vapply(records,function(x) paste(x$id,x$replicate,sep="/"),"")
    expected <- expand.grid(id=ids,replicate=seq_len(n),stringsAsFactors=FALSE)
    stopifnot(!anyDuplicated(keys),identical(sort(keys),sort(paste(expected$id,expected$replicate,sep="/"))))
    for (x in records) stopifnot(x$seed==base+match(x$id,gof_clmm_study_cells("full")$id)*100000L+x$replicate)
}
check_roster(prospective,c("core05","core11","core12"),500L,160000000L)
check_roster(numerical,c("core08","core14","link5"),32L,170000000L)
attempts <- do.call(rbind,lapply(prospective,function(x) {
    accepted <- identical(x$status,"accepted")
    stopifnot(x$status %in% c("accepted","refused"))
    row <- data.frame(id=x$id,replicate=x$replicate,seed=x$seed,status=x$status,
        stage=field(x,"stage",""),message=field(x,"message",""),
        warnings=paste(x$warnings,collapse=" | "),maximum_leverage=field(x,"maximum_leverage",NA_real_))
    for (method in c("baseline","projection","HC2")) {
        for (name in c("p","observed_squared","trace","squared_trace","quantile95","quadratic_identity_error"))
            row[[paste(method,name,sep="_")]] <- if (accepted) x[[method]][[name]] else NA_real_
        if (accepted) stopifnot(is.finite(x[[method]]$p),x[[method]]$p>=0,x[[method]]$p<=1,
            x[[method]]$quadratic_identity_error<1e-10)
    }
    row
}))
summary <- do.call(rbind,lapply(c("core05","core11","core12"),function(id)
    do.call(rbind,lapply(c("baseline","projection","HC2"),function(method) {
        group <- attempts[attempts$id==id,];p <- group[[paste0(method,"_p")]][group$status=="accepted"]
        k <- sum(p<=.05);n <- length(p)
        interval <- c(if (!k) 0 else qbeta(.025,k,n-k+1),if (k==n) 1 else qbeta(.975,k+1,n-k))
        data.frame(id=id,method=method,attempted=nrow(group),accepted=n,rejections=k,rate=k/n,
                   lower95=interval[1],upper95=interval[2],refusal=1-n/nrow(group))
    }))))
stopifnot(isTRUE(all.equal(summary,read.csv(file.path(root,"prospective-screen-summary.csv")),
    tolerance=1e-12,check.attributes=FALSE)))
numerical_attempts <- do.call(rbind,lapply(numerical,function(x) {
    stopifnot(x$status %in% c("accepted","refused","timed_out","process_failure"))
    if (x$status=="accepted") stopifnot(isTRUE(x$finite_output))
    data.frame(id=x$id,replicate=x$replicate,seed=x$seed,status=x$status,
        stage=field(x,"stage",""),code=field(x,"code",""),message=field(x,"message",""),
        seconds=field(x,"seconds",NA_real_),exit_code=field(x,"exit_code",NA_integer_),
        finite_output=field(x,"finite_output",NA),warnings=paste(x$warnings,collapse=" | "))
}))
numerical_summary <- do.call(rbind,lapply(c("core08","core14","link5"),function(id) {
    group <- numerical_attempts[numerical_attempts$id==id,]
    counts <- table(factor(group$status,levels=c("accepted","refused","timed_out","process_failure")))
    data.frame(id=id,attempted=nrow(group),as.list(counts),check.names=FALSE)
}))
stopifnot(isTRUE(all.equal(numerical_summary,read.csv(file.path(root,"numerical-screen-summary.csv")),
    tolerance=1e-12,check.attributes=FALSE)))
dir.create(output,recursive=TRUE)
write_gzip <- function(data,name) {
    connection <- gzfile(file.path(output,name),"wt");on.exit(close(connection))
    write.csv(data,connection,row.names=FALSE)
}
write_gzip(attempts,"prospective-attempts.csv.gz")
write_gzip(numerical_attempts,"numerical-attempts.csv.gz")
write.csv(summary,file.path(output,"prospective-summary.csv"),row.names=FALSE)
write.csv(numerical_summary,file.path(output,"numerical-summary.csv"),row.names=FALSE)
for (name in c("prospective-screen-plan.json","numerical-screen-plan.json","numerical-source-sha256.json"))
    stopifnot(file.copy(file.path(root,name),file.path(output,name)))
stopifnot(file.copy(file.path(root,"public-numerical-audit-v2","numerical-checks.csv"),output),
    file.copy(file.path(root,"prospective-session.txt"),file.path(output,"session-info.txt")))
inputs <- c(list.files(file.path(root,"prospective-screen"),pattern="[.]rds$",full.names=TRUE),
    list.files(file.path(root,"numerical-screen"),pattern="[.]rds$",full.names=TRUE),
    file.path(root,c("prospective-screen.R","l2-screen.R","numerical-worker.R","numerical-screen.R")))
hashes <- function(paths) stats::setNames(vapply(paths,function(p)digest::digest(file=p,algo="sha256"),""),
                                        substring(paths,nchar(root)+2L))
input_hashes <- hashes(inputs)
write_gzip(data.frame(path=names(input_hashes),sha256=unname(input_hashes)),"raw-inputs-sha256.csv.gz")
files <- list.files(output,full.names=TRUE)
provenance <- list(schema=1L,status="Development evidence; original release qualification remains failed",
    baseline_commit="ff80632ae542b4de6719282d21df9b144fe3a15c",
    source_scope="Prospective comparison uses unchanged baseline fitting. Numerical screen uses the separately frozen repair source manifest; subsequent explicit stats:: quantile prefixes do not change its mathematics.",
    scope="L2 only, three designs, not Sup/power/full-scope validation. Numerical screen checks finite API output under a 60-second budget, not p-value calibration or global-optimum certification.",
    selection="Candidates compared after completion; none adopted as a public inference default. The HC2 candidate still misses the original refusal criterion in core11.",
    prospective_failure_stages=as.list(table(attempts$stage[attempts$status!="accepted"])),
    source_inputs_sha256="raw-inputs-sha256.csv.gz",
    exported_files_sha256=as.list(stats::setNames(vapply(files,function(p)digest::digest(file=p,algo="sha256"),""),basename(files))))
jsonlite::write_json(provenance,file.path(output,"PROVENANCE.json"),auto_unbox=TRUE,pretty=TRUE)
print(summary);print(numerical_summary)
