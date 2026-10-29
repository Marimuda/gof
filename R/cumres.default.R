## Validate counts before they reach unsigned C++ arguments or allocation sizes.
gof_count <- function(x, name, minimum=0L, infinite=FALSE) {
    if (infinite && identical(x, Inf)) return(x)
    if (!is.numeric(x) || length(x)!=1L || !is.finite(x) ||
        x < minimum || x != floor(x) || x > .Machine$integer.max)
        stop(name, " must be a whole number >= ", minimum,
             if (infinite) " (or Inf)" else "", call.=FALSE)
    as.integer(x)
}

subset_t <- function(t, n.qt=250L, n.eq=250L) {
    n.qt <- gof_count(n.qt, "n.qt")
    n.eq <- gof_count(n.eq, "n.eq")
    count <- as.double(n.qt)+n.eq
    if (!length(t) || any(!is.finite(t)) || count < 1L)
        stop("Supply finite evaluation points and a positive grid size")
    if (length(t) <= count) return(seq_along(t))
    t <- sort(t)
    if (count==1L) return(length(t))
    eq <- if (n.eq) .fastapprox(t, seq(t[1], tail(t,1), length.out=n.eq))[,1] else integer()
    qt <- if (n.qt) findInterval(quantile(t, seq(0,1,length.out=n.qt)),t) else integer()
    idx <- sort(unique(c(eq, qt, length(t))))
    if (length(idx) < count) {
        remaining <- setdiff(seq_along(t), idx)
        idx <- sort(c(idx, remaining[sample.int(length(remaining),count-length(idx))]))
    }
    idx
}

## A process at a tied value contains ALL residuals at that value,
## irrespective of the order of rows or the units of the covariate.
gof_grid <- function(x, subset.max, subset) {
    values <- sort(unique(x))
    ends <- findInterval(values, sort(x))
    if (!is.null(subset.max) && length(values) > subset.max) {
        idx <- subset(values, n.qt=floor(subset.max/2), n.eq=ceiling(subset.max/2))
        if (!is.numeric(idx) || !length(idx) || any(!is.finite(idx)) ||
            any(idx != floor(idx) | idx < 1 | idx > length(values)) ||
            anyDuplicated(idx) || length(idx) > subset.max)
            stop("subset must return unique valid 1-based positions within subset.max")
        ends <- ends[sort(idx)]
    }
    as.integer(ends)
}

gof_samplestat <- function(cumres_obj, R=1e3L, x=NULL, idx=NULL) {
    R <- gof_count(R,"R",1L)
    if (!is.null(x)) cumres_obj$order(as.matrix(x))
    if (ncol(cumres_obj$inp)!=1L) stop("Supply a single ordering variable")
    t <- cumres_obj$inp[,1]
    if (is.null(idx)) idx <- gof_grid(t,Inf,subset_t)
    if (!length(idx) || any(!is.finite(idx)) || any(idx!=floor(idx) | idx<1 | idx>length(t)))
        stop("idx must contain valid 1-based positions")
    idx <- sort(unique(findInterval(t[idx],t)))
    teststat <- cumres_obj$samplestat(R, cbind(idx-1L), FALSE)
    t0 <- cumres_obj$inp[idx,1]
    W0 <- cumres_obj$obs()[idx,1]
    ks0 <- SupTest(W0)
    cvm0 <- L2Test(W0,t0)
    list(pval.ks=mean(teststat[,1]>=ks0), ks=ks0,
         pval.cvm=mean(teststat[,2]>=cvm0), cvm=cvm0, t=t0, w=W0)
}

##' @export
cumres <- function(model,...) UseMethod("cumres")

##' @export
cumres.default <- function(model, variable, r=stats::residuals(model), dr,
                           R=1e3L, b=0, plots=min(R,50L),
                           coef=stats::coef(model), ic=lava::iid(model), breakties=0,
                           subset=subset_t, subset.max=500L, ...) {
    R <- gof_count(R,"R",1L)
    plots <- gof_count(plots,"plots")
    if (!is.null(subset.max)) subset.max <- gof_count(subset.max,"subset.max",1L,TRUE)
    if (!is.function(subset)) stop("subset must be a function")
    if (!is.numeric(b) || length(b)!=1L || !is.finite(b) || b<0)
        stop("b must be a finite nonnegative bandwidth")
    if (!is.numeric(breakties) || length(breakties)!=1L || !is.finite(breakties) || breakties<0)
        stop("breakties must be a finite nonnegative number")
    if (breakties!=0) warning("breakties is deprecated; tied observations are evaluated together")
    if (anyNA(coef)) stop("Over-parametrized model")
    if (!is.null(coef) && (!is.numeric(coef) || is.complex(coef) || any(!is.finite(coef))))
        stop("Coefficients must be finite numeric values")
    residual_function <- if (is.function(r)) r else NULL
    if (!is.null(residual_function)) r <- residual_function(coef,...)
    if (missing(dr)) {
        dr <- attributes(r)$grad
        if (is.null(dr)) {
            if (is.null(residual_function)) stop("Supply a residual gradient dr (or a residual function)")
            if (!length(coef)) dr <- matrix(0,length(r),0L) else {
                if (!requireNamespace("numDeriv",quietly=TRUE)) stop("Supply gradient")
                dr <- numDeriv::jacobian(function(p) residual_function(p,...),coef)
            }
        }
    } else if (is.function(dr)) dr <- dr(coef,...)
    if (missing(ic) && !length(coef) && ncol(as.matrix(dr))==0L) ic <- matrix(0,length(r),0L)
    if (!is.numeric(r) || is.complex(r) || !is.numeric(dr) || is.complex(dr) ||
        !is.numeric(ic) || is.complex(ic))
        stop("Residuals, gradients and influences must be real numeric values")
    r <- as.numeric(r); dr <- as.matrix(dr); ic <- as.matrix(ic)
    vnames <- names(as.data.frame(variable))
    variable <- as.matrix(variable)
    n <- length(r); p <- ncol(variable)
    if (!n || !p || nrow(variable)!=n || nrow(dr)!=n || nrow(ic)!=n ||
        ncol(dr)!=ncol(ic)) stop("Incompatible residual, variable, gradient or influence dimensions")
    if (!is.numeric(variable) || is.complex(variable) || any(!is.finite(variable)) ||
        any(!is.finite(r)) || any(!is.finite(dr)) || any(!is.finite(ic)))
        stop("Residuals, variables, gradients and influences must be finite numeric values")

    ## Reserve historical uniform draws to preserve seeded Monte Carlo streams
    ## for untied models. Numerical jitter is no longer added to observations.
    invisible(runif(n))
    idx <- lapply(seq_len(p),function(i) gof_grid(variable[,i],subset.max,subset))
    lens <- lengths(idx)
    ## Repeating final endpoints pads unequal grids for C++ without changing
    ## Sup or L2 statistics. Padding is trimmed from the returned processes.
    Idx <- vapply(idx,function(i) c(i,rep(tail(i,1),max(lens)-length(i)))-1L,integer(max(lens)))
    Idx <- matrix(Idx,ncol=p)
    cr <- new(CumRes,r,dr,ic)
    cr$order(variable)
    x <- lapply(seq_len(p),function(i) cr$inp[idx[[i]],i])
    W0 <- cr$obs()
    lower <- if (b>0) lapply(seq_len(p),function(i) findInterval(x[[i]]-b,cr$inp[,i])) else NULL
    project <- function(w) lapply(seq_len(p),function(i) {
        z <- w[idx[[i]],i]
        if (b>0) z <- z-c(0,w[,i])[lower[[i]]+1L]
        z
    })
    W <- project(W0)
    stats <- function(w) unlist(lapply(seq_len(p),function(i) c(SupTest(w[[i]]),L2Test(w[[i]],x[[i]]))))
    if (b==0) {
        teststat <- cr$samplestat(R,Idx,FALSE)
    } else {
        ## A window (t-b,t] is a difference of cumulative processes. Full
        ## prefixes keep lower endpoints exact even with a reduced grid.
        teststat <- t(vapply(seq_len(R),function(k) stats(project(cr$sample1(matrix(integer(),0,0)))),numeric(2*p)))
    }
    observed <- stats(W)
    pv <- colMeans(sweep(teststat,2,observed,">="))
    What <- lapply(lens,function(k) matrix(0,k,plots))
    for (j in seq_len(plots)) {
        sim <- if (b==0) {
            w <- cr$sample1(Idx)
            lapply(seq_len(p),function(i) w[seq_len(lens[i]),i])
        } else project(cr$sample1(matrix(integer(),0,0)))
        for (i in seq_len(p)) What[[i]][,j] <- sim[[i]]
    }
    structure(list(W=W,What=What,x=x,KS=pv[seq(1,2*p,2)],CvM=pv[seq(2,2*p,2)],
                   R=R,n=n,sd=NULL,cvalues=NULL,variable=vnames,
                   type=rep("residual",p),untie=rep(0,n),model=class(model)[1]),class="cumres")
}
