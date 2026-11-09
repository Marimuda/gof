#!/usr/bin/env Rscript
# Export the sealed, completed four-cell CI slow cohort, including failed gates.
# Usage: Rscript tools/export-clmm-slow-evidence.R CALIBRATION_ZIP NEW_DESTINATION
main <- function() {
args <- commandArgs(TRUE)
if(length(args)!=2L) stop("Supply CALIBRATION_ZIP NEW_DESTINATION")
for(package in c("digest","jsonlite"))
    if(!requireNamespace(package,quietly=TRUE)) stop("Install development dependency: ",package)
archive <- normalizePath(args[1],mustWork=TRUE); output <- args[2]
if(file.exists(output)) stop("Refusing to overwrite existing slow-cohort evidence")
expected_sha256 <- "c48d391e4618724bd54330bd0b1bca002da2267b565a1a83c58ab7cdd9b1b9ad"
if(digest::digest(archive,algo="sha256",file=TRUE)!=expected_sha256)
    stop("This is not the sealed calibration artifact from CI run 42621994375")
root <- tempfile("gof-ci-slow-artifact-");dir.create(root)
on.exit(unlink(root,recursive=TRUE),add=TRUE)
entries <- utils::unzip(archive,list=TRUE)$Name
if(anyDuplicated(entries) || any(grepl("(^/|^[A-Za-z]:|(^|/)[.][.](/|$)|\\\\)",entries)))
    stop("Unsafe or duplicate archive paths")
selected <- entries[grepl("^calibration-clmm[.]rds($|[.]d/)|^lib/gof/R/",entries)]
utils::unzip(archive,files=selected,exdir=root)
report <- readRDS(file.path(root,'calibration-clmm.rds'))
protocol <- report$protocol; records <- report$records; summary <- report$summary
stopifnot(protocol$mode=='slow',protocol$replicates==500L,protocol$R==199L,
          protocol$seed_base==120000000L,length(records)==2000L,
          identical(protocol$cells$id,c('core02','core05','core11','link2')),
          identical(records,readRDS(file.path(root,'calibration-clmm.rds.d/records.rds'))))
md5 <- protocol$installed_R_md5
stopifnot(identical(unname(tools::md5sum(file.path(root,'lib/gof/R',basename(names(md5))))),unname(md5)))
keys <- vapply(records,function(x) paste(x$id,x$replicate,sep='/'),'')
stopifnot(!anyDuplicated(keys))
indices <- c(core02=2L,core05=5L,core11=11L,link2=21L)
stopifnot(all(protocol$cells$m==50L),identical(as.integer(protocol$cells$nc),c(3L,10L,10L,10L)),
    identical(protocol$cells$structure,c('intercept','intercept','slope','intercept')),
    identical(protocol$cells$link,c('logit','logit','logit','probit')),
    all(protocol$cells$sd==1),all(protocol$cells$alternative=='null'),all(protocol$cells$primary))
for (id in protocol$cells$id) {
    group <- records[vapply(records,function(x) x$id==id,TRUE)]
    stopifnot(identical(sort(vapply(group,`[[`,0L,'replicate')),seq_len(500L)))
}
for (record in records) {
    original <- file.path(root,'calibration-clmm.rds.d/records',sprintf('%s-%04d.rds',record$id,record$replicate))
    stopifnot(identical(record,readRDS(original)),
        record$seed==protocol$seed_base+indices[[record$id]]*100000L+record$replicate,
        record$status %in% c('accepted','boundary','nodes','shift','convergence','failure','indeterminate'),
        is.finite(record$seconds),record$seconds>=0)
    if(record$status=='accepted') stopifnot(all(is.finite(c(record$p_sup,record$p_l2))),
        all(c(record$p_sup,record$p_l2)>=0),all(c(record$p_sup,record$p_l2)<=1),
        max(abs(199*c(record$p_sup,record$p_l2)-round(199*c(record$p_sup,record$p_l2))))<1e-10)
}
expected <- expand.grid(id=protocol$cells$id,statistic=c('Sup','L2'),alpha=c(.01,.05,.1),
    stratum=c('all','z<1','z>=1'),stringsAsFactors=FALSE)
key <- function(data) with(data,paste(id,statistic,alpha,stratum,sep='/'))
stopifnot(nrow(summary)==72L,!anyDuplicated(key(summary)),identical(sort(key(summary)),sort(key(expected))))
equal <- function(x,y) isTRUE(all.equal(x,y,tolerance=1e-12,check.attributes=FALSE))
stopifnot(equal(summary[,names(protocol$cells)],protocol$cells[match(summary$id,protocol$cells$id),]))
for (i in seq_len(nrow(summary))) {
    row <- summary[i,]; subset <- records[vapply(records,function(x) x$id==row$id,TRUE)]
    accepted <- vapply(subset,function(x) x$status=='accepted',TRUE)
    z <- vapply(subset,`[[`,0,'boundary_z')
    use <- accepted & switch(row$stratum,all=rep(TRUE,length(z)),`z<1`=!is.na(z)&z<1,`z>=1`=!is.na(z)&z>=1)
    p <- vapply(subset,`[[`,0,if(row$statistic=='Sup')'p_sup' else 'p_l2')[use]
    n <- length(p); k <- sum(p<=row$alpha)
    ci <- if(!n)c(NA_real_,NA_real_) else c(if(!k)0 else qbeta(.025,k,n-k+1),if(k==n)1 else qbeta(.975,k+1,n-k))
    sorted <- sort(p)
    distance <- if(n)max(abs(rep(sorted,each=2L)-as.vector(rbind((seq_len(n)-1)/n,seq_len(n)/n)))) else NA_real_
    status <- vapply(subset,`[[`,'','status')
    stopifnot(row$attempted==500L,row$accepted==n,row$rejections==k,equal(row$rate,if(n)k/n else NA_real_),
        equal(c(row$lower95,row$upper95),ci),equal(row$refusal_rate,mean(!accepted)),
        row$boundary==sum(status=='boundary'),row$node_failure==sum(status=='nodes'),
        row$other_failure==sum(!status %in% c('accepted','boundary','nodes')),equal(row$uniformity_D,distance))
}
histograms <- read.csv(file.path(root,'calibration-clmm.rds.d/pvalue-histograms.csv'))
expected_histogram <- expand.grid(id=protocol$cells$id,statistic=c('p_sup','p_l2'),lower=seq(0,.9,.1),stringsAsFactors=FALSE)
hist_key <- function(data) with(data,paste(id,statistic,sprintf('%.1f',lower),sep='/'))
stopifnot(nrow(histograms)==80L,!anyDuplicated(hist_key(histograms)),
    identical(sort(hist_key(histograms)),sort(hist_key(expected_histogram))),
    equal(histograms$lower,expected_histogram$lower[match(hist_key(histograms),hist_key(expected_histogram))]),
    equal(histograms$upper,histograms$lower+.1))
for(i in seq_len(nrow(histograms))) {
    row <- histograms[i,]
    subset <- records[vapply(records,function(x)x$id==row$id&&x$status=='accepted',TRUE)]
    p <- vapply(subset,`[[`,0,row$statistic)
    stopifnot(row$count==sum(p<=row$upper & if(row$lower==0)p>=0 else p>row$lower))
}
primary <- subset(summary,alpha==.05 & stratum=='all')
criterion <- with(primary,!is.na(rate)&rate>=.03&rate<=.075&refusal_rate<.05)
stopifnot(nrow(primary)==8L,identical(primary$criterion_met,criterion))

parent <- dirname(output);dir.create(parent,recursive=TRUE,showWarnings=FALSE)
staging <- tempfile("clmm-slow-evidence-",tmpdir=parent);dir.create(staging)
summary$rows_per_block <- summary$nc
write.csv(summary,file.path(staging,"summary.csv"),row.names=FALSE)
field <- function(x,name,default) if(is.null(x[[name]]))default else x[[name]]
attempts <- do.call(rbind,lapply(records,function(x)data.frame(
    id=x$id,replicate=x$replicate,seed=x$seed,status=x$status,p_sup=x$p_sup,p_l2=x$p_l2,
    boundary_z=x$boundary_z,shift_se=x$shift_se,node_score=x$node_score,seconds=x$seconds,
    message=field(x,"message",""),warnings=paste(x$warnings,collapse=" | "),stringsAsFactors=FALSE)))
connection <- gzfile(file.path(staging,"attempts.csv.gz"),"wt")
write.csv(attempts,connection,row.names=FALSE);close(connection)
write.csv(histograms,file.path(staging,"pvalue-histograms.csv"),row.names=FALSE)
provenance <- list(schema=1L,package="gof",version=protocol$package_version,cohort="slow",
    source_commit="ea75b3b5279274afbfdf5c46a01eb2cfc2db9244",
    ci_run="https://github.com/Marimuda/gof/actions/runs/42621994375",
    ci_artifact=13680316496,artifact_sha256=expected_sha256,
    software_versions=list(R="4.6.2",ordinal="2026.7-26",numDeriv="2016.8-1.1"),
    software_version_source="Recorded CI job session log",
    created_utc=format(Sys.time(),tz="UTC",usetz=TRUE),attempts=length(records),
    cells=protocol$cells,replicates=500L,R=199L,seed_base=protocol$seed_base,
    full_cell_indices=as.list(indices),
    seed_rule="seed_base + full_cell_index * 100000 + replicate; multiplier seed = data seed + 50000",
    ordering=protocol$ordering,aggregation=protocol$aggregation,numerical_controls=protocol$nodes,
    primary_cells=4L,primary_endpoints=8L,primary_passes=sum(criterion),
    project_criteria=protocol$primary_rule,slow_cohort_criteria_met=all(criterion),
    cohort_relationship="Disjoint seeds from the full confirmation cohort. A local run with these same slow seeds is a replay, not additional independent samples; do not pool its attempts with these.",
    intervals="Pointwise exact 95% binomial intervals conditional on accepted fits; no simultaneous or unconditional guarantee",
    integrity_audit="Official artifact SHA256, installed R hashes, every individual attempt/key/seed, all 72 summary rows, independent beta-quantile intervals, refusal decomposition, uniformity distances and all 80 histogram rows audited",
    trajectory_limit="Response datasets and multiplier-statistic matrices were not retained; replay requires the recorded source and seed roster",
    review_status="Model-assisted assessment; no independent human statistical review or expert endorsement",
    installed_R_md5=as.list(stats::setNames(unname(md5),basename(names(md5)))),
    file_md5=as.list(tools::md5sum(file.path(staging,c("summary.csv","attempts.csv.gz","pvalue-histograms.csv")))))
names(provenance$file_md5) <- basename(names(provenance$file_md5))
jsonlite::write_json(provenance,file.path(staging,"PROVENANCE.json"),auto_unbox=TRUE,pretty=TRUE,na="null")
writeLines(c("# Ordinal mixed-model slow-cohort evidence","",
    "Completed synthetic four-cell cohort from the recorded CI source and disjoint slow-study seed roster.","",
    "- summary.csv: all conditional rates, exact intervals, refusal counts and covariance-z strata.",
    "- attempts.csv.gz: all 2,000 attempts, including every refused fit and warning/failure message.",
    "- pvalue-histograms.csv: descriptive ten-bin distributions among accepted fits.",
    "- PROVENANCE.json: immutable artifact digest, source, design, seed rule and unchanged project criteria.","",
    paste0(sum(criterion)," of 8 primary statistic endpoints met the project criteria. Failed endpoints remain failed. This cohort does not qualify broad mixed-model inference."),
    "A local run using these same slow-study seeds is a replay, not another 2,000 independent samples. It must not be pooled with this cohort.",
    "Response datasets and multiplier-statistic matrices were not retained. The public study generator and recorded seeds permit replay; no audit of unretained trajectories is claimed."),
    file.path(staging,"README.md"))
if(!file.rename(staging,output))stop("Could not atomically install audited slow-cohort evidence")
cat("Audited",length(records),"slow-cohort attempts and all72 summary rows; unchanged criteria:",sum(criterion),"/8 passed. Evidence:",output,"\n")
}
main()
