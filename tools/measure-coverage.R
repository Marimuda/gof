#!/usr/bin/env Rscript
# Rscript tools/measure-coverage.R [new-output-directory]
# Install covr separately in a development library. No coverage is uploaded.
invocation <- commandArgs(FALSE)
script <- sub("^--file=","",invocation[grepl("^--file=",invocation)])
if (length(script)!=1L) stop("Run this file with Rscript")
package <- normalizePath(file.path(dirname(script),".."),winslash="/",mustWork=TRUE)
args <- commandArgs(TRUE)
if (length(args)>1L) stop("Supply at most one output directory")
output <- if (length(args)) args[1] else file.path(dirname(package),"gof-coverage",
    paste0(format(Sys.time(),"%Y%m%dT%H%M%S"),"-",Sys.getpid()))
if (file.exists(output)) stop("Choose a new output directory")
if (!requireNamespace("covr",quietly=TRUE)) stop("Install covr in a development library first")
dir.create(output,recursive=TRUE)
output <- normalizePath(output,winslash="/",mustWork=TRUE)
if (startsWith(paste0(output,"/"),paste0(package,"/"))) stop("Keep coverage output outside the source package")
snapshot <- file.path(output,"source")
dir.create(snapshot)
entries <- list.files(package,all.files=TRUE,full.names=TRUE,no..=TRUE)
entries <- entries[!basename(entries) %in% c(".git",".Rproj.user",".Rhistory")]
if (!all(file.copy(entries,snapshot,recursive=TRUE,copy.mode=TRUE,copy.date=TRUE)))
    stop("Could not create source snapshot")
cpp <- Sys.getenv("GOF_CXX_STD","")
if (nzchar(cpp)) {
    if (!cpp %in% c("CXX14","CXX17","CXX20","CXX23"))
        stop("GOF_CXX_STD must be CXX14, CXX17, CXX20 or CXX23")
    cat("\nCXX_STD = ",cpp,"\n",file=file.path(snapshot,"src","Makevars"),append=TRUE,sep="")
    descriptor <- file.path(snapshot,"DESCRIPTION")
    fields <- read.dcf(descriptor)
    target <- sub("^CXX","C++",cpp)
    requirements <- if ("SystemRequirements" %in% colnames(fields)) fields[1,"SystemRequirements"] else ""
    requirements <- if (grepl("C\\+\\+ ?[0-9]+",requirements)) {
        sub("C\\+\\+ ?[0-9]+",target,requirements)
    } else paste(c(target,requirements[nzchar(requirements)]),collapse=", ")
    if (!"SystemRequirements" %in% colnames(fields)) fields <- cbind(fields,SystemRequirements=requirements)
    else fields[1,"SystemRequirements"] <- requirements
    write.dcf(fields,descriptor)
}
capture.output(sessionInfo(),file=file.path(output,"session-info.txt"))
Sys.setenv(GOF_CALIBRATION_REPORT=file.path(output,"calibration-rates.rds"),
           GOF_SEM_CALIBRATION_REPORT=file.path(output,"calibration-sem-ordinal.rds"),
           GOF_WEIGHT_CALIBRATION_REPORT=file.path(output,"calibration-weights.rds"),
           GOF_JOINT_CALIBRATION_REPORT=file.path(output,"calibration-joint.rds"),
           GOF_NORMALITY_CALIBRATION_REPORT=file.path(output,"calibration-normality.rds"),
           GOF_BOOTSTRAP_CALIBRATION_REPORT=file.path(output,"calibration-bootstrap.rds"))
client <- paste0("source(",encodeString(file.path(snapshot,"tools","check-client.R"),quote='"'),")")
coverage <- covr::package_coverage(snapshot,type=c("tests","examples"),code=client,
                                  quiet=FALSE,clean=FALSE)
saveRDS(coverage,file.path(output,"coverage.rds"))
summary <- covr::coverage_to_list(coverage)
write.csv(data.frame(file=names(summary$filecoverage),percent=unname(summary$filecoverage)),
          file.path(output,"files.csv"),row.names=FALSE)
write.csv(covr::zero_coverage(coverage),file.path(output,"uncovered.csv"),row.names=FALSE)
frame <- as.data.frame(coverage)
write.csv(frame,file.path(output,"lines.csv"),row.names=FALSE)
r <- frame[grepl("^R/",frame$filename),,drop=FALSE]
functions <- aggregate(value~filename+functions,r,
                       FUN=function(x)c(covered=sum(x>0),total=length(x)))
write.csv(functions,file.path(output,"R-functions.csv"),row.names=FALSE)
cat("Measured coverage:",summary$totalcoverage,"percent\n")
cat("Artifacts:",output,"\n")
cat("Generated wrappers and native functions unused by R tests are included.\n")
cat("Combine with the separate native harness; coverage is not a correctness proof.\n")
