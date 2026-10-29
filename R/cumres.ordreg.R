##' Cumulative residual diagnostics for an ordinal regression
##'
##' @param model A fitted lava ordinal regression.
##' @param variable Covariate name, or "predicted" for the linear predictor.
##' @param R Number of null processes to simulate.
##' @param plot Display the resulting process.
##' @param aggregate Function combining cumulative residuals across categories.
##' It must return a fixed-length finite real numeric vector for every row.
##' @param plots Number of simulated processes to retain.
##' @param ... Arguments passed to plot when plot=TRUE.
##' @return A cumres object with Sup and L2 tests.
##' @note Predictions use the fitted coefficients and the fit's internal state
##' is restored after diagnostics. The upstream \code{lava::ordreg} fitter
##' ignores \code{offset}, \code{subset}, \code{weights} and \code{na.action}
##' arguments; fits recording these options are rejected to avoid diagnosing
##' an unintended model. Subset and complete the data before fitting. Ordinal
##' offsets require upstream fitting support and cannot be added by diagnostics.
##' @export
cumres.ordreg <- function(model,variable="predicted",R=100,plot=FALSE,
                          aggregate=sum,plots=min(R,50L),...) {
    R <- gof_count(R,"R",1L)
    plots <- gof_count(plots,"plots")
    if (!is.logical(plot) || length(plot)!=1L || is.na(plot)) stop("plot must be TRUE or FALSE")
    if (!requireNamespace("numDeriv",quietly=TRUE)) stop("package 'numDeriv' not available")
    if (!is.function(aggregate)) stop("aggregate must be a function")
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
    observed <- combine(apply(residual,2,cumsum)[ends,,drop=FALSE]/sqrt(n))
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
        w <- combine(matrix(w,ncol=k))
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
                             n=n,R=R,type=rep("normal",q),sd=NULL,cvalues=NULL,
                             variable=labels,model="ordreg"),class="cumres")
    if (isTRUE(plot)) graphics::plot(result,...)
    result
}
