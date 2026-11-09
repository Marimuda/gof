## Pure-R reference implementation of the cumulative residual GoF tests.
##
## Written for clarity, not speed, following Lin, Wei & Ying (2002):
##
##   W(t)    = n^{-1/2} sum_i 1{x_i <= t} r_i                    (observed)
##   W*(t)   = n^{-1/2} [ sum_i 1{x_i <= t} r_i g_i
##                        + eta(t)' sum_i ic_i g_i ]             (under H0)
##   eta(t)  = sum_i 1{x_i <= t} dr_i   (derivative of residuals wrt beta)
##
## with g_i iid N(0,1). The p-values are P(stat(W*) >= stat(W)) estimated from
## R realisations, with stat = sup|W| or the L2 norm sqrt(int W(t)^2 dt).
##
## It consumes R's random number stream in exactly the same order as
## gof:::cumres.default (reserved uniform draws, then the R simulated processes
## sharing one g per realisation across variables, then the plotted processes),
## so with the same seed the results must agree with the C++ engine.
##
## For n > subset.max the processes are only evaluated at a subset of the
## (1-based) positions in the sorted order of each variable.

ref_cumres <- function(r, dr, ic, variable, R = 100, plots = 5,
                       breakties = 0, subset.max = 500L, b = 0) {
    r <- unname(r); dr <- unname(as.matrix(dr)); ic <- unname(as.matrix(ic))
    n <- length(r)
    variable <- unname(as.matrix(variable))
    invisible(runif(n))
    p <- ncol(variable)
    ord <- lapply(seq_len(p), function(j) order(variable[, j]))
    ## Independent grid oracle: unique values, brute-force nearest neighbours,
    ## quantile endpoints and the same explicitly reserved random fill-up.
    idx <- lapply(seq_len(p),function(j) {
        s <- sort(variable[,j]); values <- unique(s)
        ends <- which(!duplicated(s,fromLast=TRUE))
        if (b>0) {
            values <- sort(unique(c(s,s+b)))
            ends <- seq_along(values)
        }
        if (is.null(subset.max) || length(values)<=subset.max) return(ends)
        if (subset.max==1) return(tail(ends,1))
        nq <- floor(subset.max/2); ne <- ceiling(subset.max/2)
        eq <- if (ne) vapply(seq(values[1],tail(values,1),length.out=ne),function(x) {
            d <- abs(values-x); tail(which(d==min(d)),1)
        },integer(1)) else integer()
        qt <- if (nq) vapply(quantile(values,seq(0,1,length.out=nq)),function(x)
            max(which(values<=x)),integer(1)) else integer()
        selected <- sort(unique(c(eq,qt,length(values))))
        if (length(selected)<subset.max) {
            remaining <- setdiff(seq_along(values),selected)
            selected <- sort(c(selected,remaining[sample.int(length(remaining),subset.max-length(selected))]))
        }
        ends[selected]
    })

    one_process <- function(o, g) {
        w1 <- cumsum(r[o] * g)
        eta <- apply(dr[o, , drop = FALSE], 2, cumsum)
        B <- colSums(ic[o, , drop = FALSE] * g)
        (w1 + drop(eta %*% B)) / sqrt(n)
    }
    sup_stat <- function(w) max(abs(w))
    l2_stat <- function(w, t) sqrt(sum(c(diff(t), 0) * w^2))

    t0 <- lapply(seq_len(p), function(j) {
        values <- if (b==0) sort(variable[,j]) else sort(unique(c(variable[,j],variable[,j]+b)))
        values[idx[[j]]]
    })
    W <- lapply(seq_len(p), function(j) {
        if (b==0) return((cumsum(r[ord[[j]]])/sqrt(n))[idx[[j]]])
        vapply(t0[[j]],function(t) sum(r*(variable[,j]<=t & variable[,j]+b>t))/sqrt(n),numeric(1))
    })
    sample_at <- function(j,g) {
        if (b==0) return(one_process(ord[[j]],g)[idx[[j]]])
        o <- ord[[j]]
        vapply(t0[[j]],function(t) {
            in_window <- variable[o,j]<=t & variable[o,j]+b>t
            (sum(r[o]*g*in_window)+sum(colSums(dr[o,,drop=FALSE]*in_window)*colSums(ic[o,,drop=FALSE]*g)))/sqrt(n)
        },numeric(1))
    }

    sims <- matrix(0, R, 2 * p)
    for (k in seq_len(R)) {
        g <- rnorm(n)
        for (j in seq_len(p)) {
            w <- sample_at(j,g)
            sims[k, 2 * j - 1] <- sup_stat(w)
            sims[k, 2 * j] <- l2_stat(w, t0[[j]])
        }
    }
    What <- lapply(seq_len(p), function(j) matrix(0, length(idx[[j]]), plots))
    for (k in seq_len(plots)) {
        g <- rnorm(n)
        for (j in seq_len(p))
            What[[j]][, k] <- sample_at(j,g)
    }
    KS <- sapply(seq_len(p), function(j) mean(sims[, 2 * j - 1] >= sup_stat(W[[j]])))
    CvM <- sapply(seq_len(p), function(j) mean(sims[, 2 * j] >= l2_stat(W[[j]], t0[[j]])))
    list(W = W, What = What, x = t0, KS = KS, CvM = CvM)
}

## The inputs that cumres.glm should pass on to cumres.default (ungrouped
## glm). The derivative of the residual r_i = y_i - mu_i wrt beta is
## -mu.eta(eta_i) x_i, where eta_i is the full linear predictor, offsets
## included.
glm_inputs <- function(model) {
    X <- model.matrix(model)
    eta <- predict(model, type = "link")
    list(r = residuals(model, type = "response"),
         dr = -(as.numeric(family(model)$mu.eta(eta)) * X),
         ic = lava::iid(model))
}

fx_poisson_offset <- function(n = 300, seed = 7) {
    set.seed(seed)
    d <- data.frame(x = runif(n, 0, 2), exposure = runif(n, 0.5, 20))
    d$y <- rpois(n, d$exposure * exp(-1 + 0.5 * d$x))
    glm(y ~ x + offset(log(exposure)), data = d, family = poisson())
}

## Verbatim copy of cumres.ordreg from gof 1.3.0 (commit cecbcba; unchanged
## at a3c5731), kept as the bit-identity oracle for the shared ordinal engine.
## Do not edit: the refactored method must reproduce it exactly, including
## the random number stream.
ref_cumres_ordreg_130 <- function(model,variable="predicted",R=100,plot=FALSE,
                          aggregate=c("omnibus","score","cutpoints"),plots=min(R,50L),...) {
    if (any(c("resolution","on_unresolved") %in% names(list(...))))
        stop("resolution and on_unresolved are not implemented for ordinal diagnostics",call.=FALSE)
    if ("simulation" %in% names(list(...)))
        stop("Ordinal diagnostics use shared cutpoint multipliers for one ordering; simulation is a default/linear/GLM option")
    R <- gof_count(R,"R",1L)
    plots <- gof_count(plots,"plots")
    if (!is.logical(plot) || length(plot)!=1L || is.na(plot)) stop("plot must be TRUE or FALSE")
    if (!requireNamespace("numDeriv",quietly=TRUE)) stop("package 'numDeriv' not available")
    if (is.character(aggregate)) {
        aggregation <- match.arg(aggregate,c("omnibus","score","cutpoints"))
        aggregate <- switch(aggregation,omnibus=gof_ordinal_norm,score=sum,cutpoints=identity)
    } else {
        if (!is.function(aggregate)) stop("aggregate must be a function or a named aggregation mode")
        aggregation <- if (identical(aggregate,sum)) "score" else
            if (identical(aggregate,identity)) "cutpoints" else "custom"
    }
    env <- model$up
    if (!is.environment(env) || any(!is.finite(coef(model))) ||
        (!is.null(model$opt$convergence) && !isTRUE(model$opt$convergence==0L)))
        stop("Supply a finite, converged ordinal fit")
    if (length(env$n)!=1L || !is.finite(env$n) || env$n<1 ||
        length(env$K)!=1L || !is.finite(env$K) || env$K<2 ||
        !is.numeric(env$y) || length(env$y)!=env$n || any(!is.finite(env$y)) ||
        !is.matrix(env$X) || nrow(env$X)!=env$n)
        stop("Incompatible ordinal fitted response and predictor dimensions")
    ignored <- intersect(c("offset","subset","weights","na.action"),names(as.list(model$call)))
    if (length(ignored)) stop("lava::ordreg ignores these fitting options: ",paste(ignored,collapse=", "),
                              ". Fit explicitly preprocessed data; ordinal offsets are not supported by lava.")
    ## lava::predict.ordreg mutates this shared environment. Select predictors
    ## at the fitted coefficients and restore the caller's state on all exits.
    saved <- as.list(env,all.names=TRUE)
    on.exit({
        added <- setdiff(ls(env,all.names=TRUE),names(saved))
        if (length(added)) rm(list=added,envir=env)
        list2env(saved,envir=env)
    },add=TRUE)
    cpr <- as.matrix(predict(model,p=coef(model),type="cumulative"))
    if (length(variable)!=1L || !is.character(variable) || is.na(variable)) stop("Select one ordinal covariate")
    if (variable=="predicted") x <- env$lp[,1] else {
        if (!variable %in% colnames(env$X)) stop("Variable not found: ",variable)
        x <- env$X[,variable]
    }
    n <- length(x)
    k <- length(env$alpha)
    if (!n || !k || any(!is.finite(x))) stop("Ordinal covariates must be finite")
    if (nrow(cpr)!=n || ncol(cpr)!=k || length(env$y)!=n ||
        any(!is.finite(env$y)) || any(env$y!=floor(env$y) | env$y<1 | env$y>k+1L))
        stop("Incompatible ordinal fitted response and prediction dimensions")
    informative <- gof_ordinal_information(model,x)
    Y <- outer(env$y,seq_len(k),"<=")*1
    residual <- Y-cpr
    gradient <- -numDeriv::jacobian(function(p) predict(model,p=p,type="cumulative"),coef(model))
    ic <- as.matrix(lava::iid(model))
    if (!identical(dim(gradient),c(n*k,length(coef(model)))) ||
        !identical(dim(ic),c(n,length(coef(model)))))
        stop("Incompatible ordinal gradient or influence dimensions")
    o <- order(x)
    x <- x[o]
    ends <- cumsum(rle(x)$lengths)
    x <- x[ends]
    residual <- residual[o,,drop=FALSE]
    eta <- lapply(seq_len(k),function(j) {
        z <- gradient[seq_len(n)+(j-1L)*n,,drop=FALSE][o,,drop=FALSE]
        apply(z,2,cumsum)
    })
    ic <- ic[o,,drop=FALSE]
    if (any(!is.finite(residual)) || any(!is.finite(gradient)) || any(!is.finite(ic)))
        stop("Ordinal residuals, gradients and influences must be finite")
    combine <- function(w) {
        rows <- lapply(seq_len(nrow(w)),function(i) aggregate(w[i,]))
        sizes <- lengths(rows)
        if (!length(sizes) || !sizes[1] || any(sizes!=sizes[1]) ||
            any(!vapply(rows,function(z) is.numeric(z) && !is.complex(z) && is.null(dim(z)),logical(1))) ||
            any(!is.finite(unlist(rows))))
            stop("aggregate must return a fixed-length finite numeric vector")
        do.call(rbind,rows)
    }
    observed_cdf <- apply(residual,2,cumsum)[ends,,drop=FALSE]/sqrt(n)
    if (identical(informative,FALSE)) observed_cdf[] <- 0
    observed <- combine(observed_cdf)
    q <- ncol(observed)
    What <- lapply(seq_len(q),function(i) matrix(0,length(ends),plots))
    obs.ks <- apply(abs(observed),2,max)
    obs.l2 <- apply(observed,2,function(w) L2Test(w,x))
    ks <- l2 <- matrix(0,R,q)
    for (j in seq_len(max(R,plots))) {
        g <- rnorm(n)
        B <- colSums(ic*g)
        w <- vapply(seq_len(k),function(i)
            (cumsum(residual[,i]*g)+drop(eta[[i]] %*% B))[ends]/sqrt(n),
            numeric(length(ends)))
        w <- matrix(w,ncol=k)
        if (identical(informative,FALSE)) w[] <- 0
        w <- combine(w)
        if (ncol(w)!=q) stop("aggregate returned a different number of processes")
        if (j<=R) {
            ks[j,] <- apply(abs(w),2,max)
            l2[j,] <- apply(w,2,function(z) L2Test(z,x))
        }
        if (j<=plots) for (i in seq_len(q)) What[[i]][,j] <- w[,i]
    }
    labels <- if (q==1L) variable else paste(variable,seq_len(q),sep=".")
    result <- structure(list(W=lapply(seq_len(q),function(i) observed[,i]),
                             What=What,x=rep(list(x),q),
                             KS=colMeans(sweep(ks,2,obs.ks,">=")),
                             CvM=colMeans(sweep(l2,2,obs.l2,">=")),
                             n=n,R=R,type=rep(if (aggregation=="omnibus") "ordinal-norm" else "normal",q),
                             aggregation=aggregation,informative=rep(informative,q),sd=NULL,cvalues=NULL,
                             variable=labels,model="ordreg",
                             design_info=gof_design_info(unit="subject",n_rows=n,
                                 fit_parameters=length(coef(model)),estimated_parameters=length(coef(model)),
                                 correction_dimension=ncol(ic),n_components=k)),class="cumres")
    if (is.na(informative)) result$KS[] <- result$CvM[] <- NA_real_
    if (isTRUE(plot)) graphics::plot(result,...)
    result
}
