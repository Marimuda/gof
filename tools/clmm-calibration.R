#!/usr/bin/env Rscript
# Usage: Rscript tools/clmm-calibration.R OUTPUT [full|slow|pilot] [N] [workers]
args <- commandArgs(TRUE)
if (length(args)<1L || length(args)>4L) stop("Supply OUTPUT [full|slow|pilot] [N] [workers]")
invocation <- commandArgs(FALSE)
script <- sub("^--file=","",invocation[grepl("^--file=",invocation)])
source(file.path(dirname(script),"clmm-study.R"))
mode <- if (length(args)>=2L) args[2L] else "full"
count <- if (length(args)>=3L) as.integer(args[3L]) else if (mode=="slow") 500L else 1000L
workers <- if (length(args)>=4L) as.integer(args[4L]) else 1L
stopifnot(!is.na(count),count>0L,!is.na(workers),workers>0L)
gof_clmm_study_run(args[1L],mode,count,R=199L,workers=workers)
