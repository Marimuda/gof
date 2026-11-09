#!/usr/bin/env Rscript
# Independently audit a completed cohort before installing its compact evidence.
# Usage: Rscript tools/export-clmm-evidence.R AUDIT_DIRECTORY DESTINATION
args <- commandArgs(TRUE)
if (length(args)!=2L) stop("Supply AUDIT_DIRECTORY DESTINATION")
audit <- normalizePath(args[1],mustWork=TRUE)
output <- args[2]
if (file.exists(output)) stop("Refusing to overwrite existing installed evidence")
if (!requireNamespace("jsonlite",quietly=TRUE)) stop("Install jsonlite in a development library to export JSON provenance")
if (!file.exists(file.path(audit,"confirmation-completed.txt"))) stop("The confirmation cohort is not complete")
study <- file.path(audit,"confirmation")
protocol <- readRDS(file.path(study,"protocol.rds"))
records <- readRDS(file.path(study,"records.rds"))
summary <- readRDS(file.path(study,"summary.rds"))
stopifnot(protocol$mode=="full",protocol$replicates==1000L,protocol$R==199L,
          nrow(protocol$cells)==29L,length(records)==29000L)
md5 <- protocol$installed_R_md5
stopifnot(identical(unname(tools::md5sum(names(md5))),unname(md5)))

field <- function(x,name,default) if (is.null(x[[name]])) default else x[[name]]
attempts <- do.call(rbind,lapply(records,function(x) data.frame(
    id=x$id,replicate=x$replicate,seed=x$seed,status=x$status,
    p_sup=x$p_sup,p_l2=x$p_l2,boundary_z=x$boundary_z,shift_se=x$shift_se,
    node_score=x$node_score,seconds=x$seconds,
    message=field(x,"message",""),warnings=paste(x$warnings,collapse=" | "),stringsAsFactors=FALSE)))
stopifnot(!anyDuplicated(paste(attempts$id,attempts$replicate)),
          all(attempts$id %in% protocol$cells$id),
          all(attempts$status %in% c("accepted","boundary","nodes","shift","convergence","failure","indeterminate")),
          all(is.finite(attempts$seconds) & attempts$seconds>=0))
for (id in protocol$cells$id) stopifnot(identical(sort(attempts$replicate[attempts$id==id]),seq_len(1000L)))
for (record in records) {
    retained <- file.path(study,"records",sprintf("%s-%04d.rds",record$id,record$replicate))
    stopifnot(file.exists(retained),identical(readRDS(retained),record))
}
cell_index <- match(attempts$id,protocol$cells$id)
stopifnot(all(attempts$seed==protocol$seed_base+cell_index*100000L+attempts$replicate))
accepted <- attempts$status=="accepted"
stopifnot(all(is.finite(attempts$p_sup[accepted]) & is.finite(attempts$p_l2[accepted])),
          all(attempts$p_sup[accepted]>=0 & attempts$p_sup[accepted]<=1),
          all(attempts$p_l2[accepted]>=0 & attempts$p_l2[accepted]<=1))

equal <- function(x,y) isTRUE(all.equal(x,y,tolerance=1e-12,check.attributes=FALSE))
# Check the complete table roster, not just the rows that happen to be present.
expected <- expand.grid(id=protocol$cells$id,statistic=c("Sup","L2"),
    alpha=c(.01,.05,.10),stratum=c("all","z<1","z>=1"),stringsAsFactors=FALSE)
summary_key <- function(data) with(data,paste(id,statistic,alpha,stratum,sep="/"))
stopifnot(nrow(summary)==nrow(expected),!anyDuplicated(summary_key(summary)),
    identical(sort(summary_key(summary)),sort(summary_key(expected))),
    equal(summary[,names(protocol$cells)],
          protocol$cells[match(summary$id,protocol$cells$id),]))
for (i in seq_len(nrow(summary))) {
    row <- summary[i,]
    data <- attempts[attempts$id==row$id,]
    use <- data$status=="accepted"
    if (row$stratum=="z<1") use <- use & !is.na(data$boundary_z) & data$boundary_z<1
    if (row$stratum=="z>=1") use <- use & !is.na(data$boundary_z) & data$boundary_z>=1
    p <- data[[if (row$statistic=="Sup") "p_sup" else "p_l2"]][use]
    n <- length(p); k <- sum(p<=row$alpha)
    # Direct beta quantiles independently reproduce binom.test() intervals.
    interval <- if (!n) c(NA_real_,NA_real_) else c(
        if (!k) 0 else qbeta(.025,k,n-k+1),
        if (k==n) 1 else qbeta(.975,k+1,n-k))
    stopifnot(row$attempted==nrow(data),row$accepted==n,row$rejections==k,
        equal(row$rate,if (n) k/n else NA_real_),
        equal(c(row$lower95,row$upper95),interval),
        equal(row$refusal_rate,mean(data$status!="accepted")),
        row$boundary==sum(data$status=="boundary"),
        row$node_failure==sum(data$status=="nodes"),
        row$other_failure==sum(!data$status %in% c("accepted","boundary","nodes")))
    distance <- if (n) {
        ordered <- sort(p)
        max(abs(rep(ordered,each=2L)-as.vector(rbind((seq_len(n)-1)/n,seq_len(n)/n))))
    } else NA_real_
    stopifnot(equal(row$uniformity_D,distance))
}
histograms <- read.csv(file.path(study,"pvalue-histograms.csv"),stringsAsFactors=FALSE)
expected_histogram <- expand.grid(id=protocol$cells$id,statistic=c("p_sup","p_l2"),
    lower=seq(0,.9,.1),stringsAsFactors=FALSE)
histogram_key <- function(data) with(data,paste(id,statistic,sprintf("%.1f",lower),sep="/"))
stopifnot(nrow(histograms)==nrow(expected_histogram),
    !anyDuplicated(histogram_key(histograms)),
    identical(sort(histogram_key(histograms)),sort(histogram_key(expected_histogram))),
    equal(histograms$upper,histograms$lower+.1))
for (i in seq_len(nrow(histograms))) {
    row <- histograms[i,]
    data <- attempts[attempts$id==row$id & attempts$status=="accepted",]
    p <- data[[row$statistic]]
    count <- sum(p<=row$upper & if (row$lower==0) p>=0 else p>row$lower)
    stopifnot(row$count==count)
}
primary <- subset(summary,primary & alpha==.05 & stratum=="all")
criterion <- with(primary,!is.na(rate) & rate>=.03 & rate<=.075 & refusal_rate<.05)
stopifnot(nrow(primary)==36L,identical(primary$criterion_met,criterion))

parent <- dirname(output)
dir.create(parent,recursive=TRUE,showWarnings=FALSE)
staging <- tempfile("clmm-evidence-",tmpdir=parent)
dir.create(staging)
# A failed staging directory is retained for diagnosis by this CLI.
summary$rows_per_block <- ifelse(summary$structure=="nested",3L*summary$nc,summary$nc)
write.csv(summary,file.path(staging,"summary.csv"),row.names=FALSE)
connection <- gzfile(file.path(staging,"attempts.csv.gz"),"wt")
write.csv(attempts,connection,row.names=FALSE)
close(connection)
stopifnot(file.copy(file.path(study,"pvalue-histograms.csv"),staging))
provenance <- list(schema=1L,package="gof",version=protocol$package_version,
    scientific_commit="5bc9e4d8b0f39724a1df931bcfb011d331a3534e",
    original_feature_commit=trimws(readLines(file.path(audit,"source-base-commit.txt"))),
    created_utc=format(Sys.time(),tz="UTC",usetz=TRUE),
    attempts=length(records),cells=protocol$cells,replicates=protocol$replicates,R=protocol$R,
    seed_base=protocol$seed_base,seed_rule="seed_base + cell_index * 100000 + replicate; multiplier seed = data seed + 50000",
    ordering=protocol$ordering,aggregation=protocol$aggregation,nested_design=protocol$nested,
    numerical_controls=protocol$nodes,acceptance=protocol$acceptance,
    intervals="Pointwise exact 95% binomial intervals among accepted fits; no simultaneous or unconditional guarantee",
    project_primary_rule=protocol$primary_rule,primary_cells=nrow(primary),primary_passes=sum(criterion),
    empirical_release_gate_passed=all(criterion),
    review_status="Model-assisted assessment; no independent human statistical review or expert endorsement",
    integrity_audit="Every planned key, seed, status, count, rate, refusal decomposition and binomial interval audited against retained attempt records",
    trajectory_limit="Response datasets and multiplier-statistic matrices were not retained; replay requires the frozen source and seed roster",
    installed_R_md5=as.list(stats::setNames(unname(md5),basename(names(md5)))),
    file_md5=as.list(tools::md5sum(file.path(staging,c("summary.csv","attempts.csv.gz","pvalue-histograms.csv")))))
names(provenance$file_md5) <- basename(names(provenance$file_md5))
jsonlite::write_json(provenance,file.path(staging,"PROVENANCE.json"),pretty=TRUE,auto_unbox=TRUE,na="null")
writeLines(c("# Ordinal mixed-model calibration evidence","",
    "Synthetic confirmation data for the declared CLMM protocol. See the installed ordinal-mixed-models guide for interpretation.","",
    "- summary.csv: conditional rejection rates, pointwise exact intervals, refusal counts and covariance-z strata.",
    "- attempts.csv.gz: all 29,000 planned samples, seeds, outcomes, numerical checks and warning/failure messages.",
    "- pvalue-histograms.csv: descriptive ten-bin distributions among accepted fits.",
    "- PROVENANCE.json: source binding, design roster, seed rule, integrity checks and the unchanged project release criterion.","",
    "All refusals remain in refusal-rate denominators. Rejection rates condition on accepted fits. Failed project criteria remain failed; this panel does not prove universal calibration.",
    "Response datasets and multiplier-statistic matrices were not retained. The frozen source and recorded seeds permit replay; no audit of unretained trajectories is claimed."),
    file.path(staging,"README.md"))
if (!file.rename(staging,output)) stop("Could not atomically install the audited evidence directory")
cat("Audited",length(records),"attempts and",nrow(summary),"summary rows. Primary project gates:",
    sum(criterion),"/",length(criterion),"passed. Evidence:",output,"\n")
