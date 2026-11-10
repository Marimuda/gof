# CLMM simulation protocol and resumable execution helpers. Source after
# loading gof; kept separate from the CLI so opt-in tests use the same design.

gof_clmm_study_cells <- function(mode=c("full","slow","pilot")) {
    mode <- match.arg(mode)
    core <- expand.grid(m=c(20L,50L,150L),nc=c(3L,10L),
        structure=c("intercept","slope","nested"),stringsAsFactors=FALSE)
    core$link <- "logit"; core$sd <- 1; core$alternative <- "null"
    core <- rbind(core,data.frame(m=50L,nc=10L,structure="intercept",link="logit",sd=.5,alternative="null"))
    core$id <- sprintf("core%02d",seq_len(nrow(core)))
    links <- data.frame(m=50L,nc=10L,structure="intercept",
        link=c("logit","probit","cloglog","loglog","cauchit"),sd=1,alternative="null",
        id=paste0("link",1:5),stringsAsFactors=FALSE)
    power <- data.frame(m=50L,nc=10L,structure="intercept",link="logit",sd=1,
        alternative=c("quadratic","wrong-cloglog","wrong-cauchit","omitted-slope","nonproportional"),
        id=paste0("power",1:5),stringsAsFactors=FALSE)
    cells <- rbind(core,links,power)
    cells$primary <- cells$alternative=="null" & cells$m>=50
    if (mode=="slow") cells <- cells[match(c("core02","core05","core11","link2"),cells$id),]
    if (mode=="pilot") cells <- cells[match(c("core02","core11","core17","link3","link5"),cells$id),]
    rownames(cells) <- NULL
    cells
}

gof_clmm_study_data <- function(cell,seed) {
    set.seed(seed)
    child <- if (cell$structure=="nested") 3L else 1L
    n <- cell$m*cell$nc*child
    g <- rep(seq_len(cell$m),each=cell$nc*child)
    h <- rep(rep(seq_len(child),each=cell$nc),cell$m)
    x <- rnorm(n); z <- rnorm(n)
    eta <- .8*x
    if (cell$structure=="slope" || cell$alternative=="omitted-slope") {
        covariance <- matrix(c(cell$sd^2,.3*cell$sd*.5,.3*cell$sd*.5,.5^2),2)
        b <- matrix(rnorm(2*cell$m),cell$m) %*% chol(covariance)
        eta <- eta+b[g,1]+b[g,2]*z
    } else eta <- eta+rnorm(cell$m,sd=cell$sd)[g]
    if (child>1) eta <- eta+rnorm(cell$m*child,sd=.5)[(g-1L)*child+h]
    if (cell$alternative=="quadratic") eta <- eta+.8*(x^2-1)
    link <- if (startsWith(cell$alternative,"wrong-")) sub("wrong-","",cell$alternative) else cell$link
    F <- switch(link,logit=plogis,probit=pnorm,cloglog=function(a) -expm1(-exp(a)),
        loglog=function(a) exp(-exp(-a)),cauchit=pcauchy)
    cutpoints <- c(-1.5,0,1.5)
    if (cell$alternative=="nonproportional") {
        ## Ordered cutpoints for every x: a smooth positive spacing creates
        ## nonparallel conditional logits without repairing crossed CDFs.
        spacing <- 1.5*exp(.45*tanh(x))
        cdf <- cbind(F(-spacing-eta),F(-eta),F(spacing-eta))
    } else cdf <- sapply(cutpoints,function(a) F(a-eta))
    u <- runif(n)
    y <- 1L+rowSums(u>cdf)
    data.frame(y=ordered(y,levels=1:4),x=x,z=z,g=factor(g),h=factor(h))
}

gof_clmm_study_attempt <- function(cell,replicate,R=199L,seed_base=140000000L) {
    ## The seed roster is independent of worker count, run order and retries.
    all_ids <- gof_clmm_study_cells("full")$id
    index <- match(cell$id,all_ids)
    seed <- as.integer(seed_base+index*100000L+replicate)
    started <- proc.time()[3]
    warnings <- character()
    stage <- "data_generation"
    result <- tryCatch(withCallingHandlers({
        data <- gof_clmm_study_data(cell,seed)
        formula <- switch(cell$structure,intercept=y~x+(1|g),slope=y~x+(1+z|g),nested=y~x+(1|g/h))
        stage <- "input_fit"
        model <- ordinal::clmm(formula,data=data,link=cell$link)
        set.seed(seed+50000L)
        stage <- "diagnostic"
        diagnostic <- gof::cumres(model,variable="x",R=R,plots=0L)
        list(status=if (anyNA(diagnostic$KS) || anyNA(diagnostic$CvM)) "indeterminate" else "accepted",
             p_sup=diagnostic$KS[1],p_l2=diagnostic$CvM[1],
             boundary_z=min(diagnostic$fit_check$boundary$z),
             shift_se=max(abs(diagnostic$fit_check$shift$shift_se),na.rm=TRUE),
             shift_comparable=all(is.finite(diagnostic$fit_check$shift$shift_se)),
             node_score=diagnostic$fit_check$node_check$standardized_score)
    },warning=function(w) {warnings <<- c(warnings,conditionMessage(w)); invokeRestart("muffleWarning")}),
    error=function(e) {
        message <- conditionMessage(e)
        code <- if (inherits(e,"gof_clmm_error")) e$code else paste0(stage,"_failure")
        failure_stage <- if (inherits(e,"gof_clmm_error")) e$stage else stage
        status <- if (inherits(e,"gof_clmm_error")) switch(code,
            covariance_boundary="boundary",quadrature_unstable="nodes",parameter_shift="shift",
            input_nonconvergence="convergence",input_gradient="convergence",
            polish_nonconvergence="convergence",refit_failed="convergence","failure") else
            if (grepl("boundary|singular covariance",message)) "boundary" else
            if (grepl("node check",message)) "nodes" else
            if (grepl("shift.tol",message)) "shift" else
            if (grepl("converge|gradient|line search|iteration",message)) "convergence" else "failure"
        list(status=status,message=message,p_sup=NA_real_,p_l2=NA_real_,
             boundary_z=NA_real_,shift_se=NA_real_,node_score=NA_real_,
             failure_code=code,failure_stage=failure_stage,
             failure_diagnostics=if (inherits(e,"gof_clmm_error")) e$diagnostics else list())
    })
    result$id <- cell$id; result$replicate <- replicate; result$seed <- seed
    result$warnings <- unique(warnings); result$seconds <- unname(proc.time()[3]-started)
    result
}

gof_clmm_study_summary <- function(records,cells) {
    rows <- list(); j <- 0L
    for (id in unique(vapply(records,`[[`,"","id"))) {
        subset <- records[vapply(records,function(x) x$id==id,TRUE)]
        cell <- cells[match(id,cells$id),]
        accepted <- vapply(subset,function(x) x$status=="accepted",TRUE)
        for (statistic in c("Sup","L2")) for (alpha in c(.01,.05,.10)) {
            p <- vapply(subset,`[[`,0,if (statistic=="Sup") "p_sup" else "p_l2")
            for (stratum in c("all","z<1","z>=1")) {
                z <- vapply(subset,`[[`,0,"boundary_z")
                use <- accepted & switch(stratum,all=rep(TRUE,length(z)),`z<1`=!is.na(z)&z<1,`z>=1`=!is.na(z)&z>=1)
                n <- sum(use); rejected <- sum(p[use]<=alpha)
                ci <- if (n) binom.test(rejected,n)$conf.int else c(NA_real_,NA_real_)
                j <- j+1L
                rows[[j]] <- data.frame(cell,statistic=statistic,alpha=alpha,stratum=stratum,
                    attempted=length(subset),accepted=n,rejections=rejected,
                    rate=if (n) rejected/n else NA_real_,lower95=ci[1],upper95=ci[2],
                    refusal_rate=mean(!accepted),boundary=sum(vapply(subset,function(x) x$status=="boundary",TRUE)),
                    node_failure=sum(vapply(subset,function(x) x$status=="nodes",TRUE)),
                    other_failure=sum(!accepted & vapply(subset,function(x) !x$status %in% c("boundary","nodes"),TRUE)),
                    seconds=sum(vapply(subset,`[[`,0,"seconds")),
                    uniformity_D=if (n) {
                        ordered <- sort(p[use]); max(seq_len(n)/n-ordered,ordered-(seq_len(n)-1)/n)
                    } else NA_real_,row.names=NULL)
            }
        }
    }
    do.call(rbind,rows)
}

gof_clmm_study_run <- function(output,mode="full",replicates=1000L,R=199L,workers=1L,
                               seed_base=if (mode=="pilot") 110000000L else if (mode=="slow") 120000000L else 140000000L) {
    stopifnot(requireNamespace("ordinal",quietly=TRUE),requireNamespace("numDeriv",quietly=TRUE))
    cells <- gof_clmm_study_cells(mode)
    dir.create(output,recursive=TRUE,showWarnings=FALSE)
    dir.create(file.path(output,"records"),showWarnings=FALSE)
    protocol <- list(schema=2L,mode=mode,cells=cells,replicates=replicates,R=R,seed_base=seed_base,
        failure_accounting="Structured failure code/stage/state; nonconvergence is not classified as a covariance boundary",
        ordering="x",aggregation="omnibus",nodes="package defaults; Q+4 check",nested="three children per top-level cluster, nc rows per child",
        acceptance="conditional on an accepted numerical fit; all refused attempts remain in the denominator of refusal rates",
        primary_rule="m >= 50 null cells: 5% rejection rate in [0.03,0.075] and refusal < 0.05; m=20 descriptive",
        package_version=as.character(utils::packageVersion("gof")),
        installed_R_md5=tools::md5sum(list.files(file.path(find.package("gof"),"R"),full.names=TRUE)))
    seal <- file.path(output,"protocol.rds")
    if (file.exists(seal)) stopifnot(identical(readRDS(seal),protocol)) else {
        saveRDS(protocol,seal); capture.output(dput(protocol),file=file.path(output,"protocol.txt"))
    }
    tasks <- expand.grid(cell=seq_len(nrow(cells)),replicate=seq_len(replicates))
    file <- function(i) file.path(output,"records",sprintf("%s-%04d.rds",cells$id[tasks$cell[i]],tasks$replicate[i]))
    pending <- which(!vapply(seq_len(nrow(tasks)),function(i) file.exists(file(i)),TRUE))
    ## Small resumable batches keep progress durable and worker failures visible.
    for (batch in split(pending,ceiling(seq_along(pending)/max(workers*50L,1L)))) {
        results <- parallel::mclapply(batch,function(i) {
            value <- gof_clmm_study_attempt(cells[tasks$cell[i],],tasks$replicate[i],R,seed_base)
            target <- file(i); temporary <- paste0(target,".tmp-",Sys.getpid())
            saveRDS(value,temporary)
            if (!file.rename(temporary,target)) stop("Could not retain attempt ",target)
            TRUE
        },mc.cores=workers,mc.preschedule=FALSE,mc.set.seed=FALSE)
        if (!all(vapply(results,isTRUE,TRUE))) stop("A worker failed; completed attempts remain resumable")
        cat(format(Sys.time())," completed ",sum(vapply(seq_len(nrow(tasks)),function(i) file.exists(file(i)),TRUE)),"/",nrow(tasks),"\n",sep="")
        flush.console()
    }
    records <- lapply(seq_len(nrow(tasks)),function(i) readRDS(file(i)))
    saveRDS(records,file.path(output,"records.rds"))
    summary <- gof_clmm_study_summary(records,cells)
    summary$criterion_met <- with(summary,ifelse(primary & alpha==.05 & stratum=="all",
        !is.na(rate) & rate>=.03 & rate<=.075 & refusal_rate<.05,NA))
    write.csv(summary,file.path(output,"summary.csv"),row.names=FALSE)
    saveRDS(summary,file.path(output,"summary.rds"))
    capture.output(sessionInfo(),file=file.path(output,"session.txt"))
    histogram <- do.call(rbind,lapply(cells$id,function(id) {
        subset <- records[vapply(records,function(x) x$id==id && x$status=="accepted",TRUE)]
        do.call(rbind,lapply(c("p_sup","p_l2"),function(statistic) {
            p <- vapply(subset,`[[`,0,statistic)
            data.frame(id=id,statistic=statistic,lower=seq(0,.9,.1),upper=seq(.1,1,.1),
                count=hist(p,breaks=seq(0,1,.1),plot=FALSE,include.lowest=TRUE)$counts)
        }))
    }))
    write.csv(histogram,file.path(output,"pvalue-histograms.csv"),row.names=FALSE)
    invisible(list(protocol=protocol,records=records,summary=summary))
}
