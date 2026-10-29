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

gof_window_grid <- function(x, b, subset.max, subset) {
    exits <- x+b
    if (any(!is.finite(exits)))
        stop("Moving-window exit coordinates must be finite; rescale the covariate and bandwidth")
    if (any(exits<=x))
        stop("Moving-window bandwidth cannot be resolved at the covariate's numeric scale; rescale the covariate and bandwidth")
    # A window changes when a residual enters at x or leaves at x+b. Include
    # the zero-valued final exit so the complete left-step integral is defined.
    values <- sort(unique(c(x,exits)))
    values[gof_grid(values,subset.max,subset)]
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

##' Cumulative residual diagnostics for fitted models
##'
##' Check a fitted mean structure by accumulating residuals over an ordering
##' variable and comparing the observed process with simulated null processes.
##' The fitted model is the first argument, so model objects can be passed with
##' the native pipe \code{|>}.
##'
##' @param model A fitted model. Methods are supplied for \code{lm}, supported
##' \code{glm}, \code{lava::lvmfit} and \code{lava::ordreg} objects. Use
##' \code{NULL} with the default method when supplying residuals and influences.
##' @param ... Arguments for the selected method, including ordering variables,
##' the number of simulations and the number of curves to retain.
##' @details
##' Use \code{cumres.glm} for normal, binomial and Poisson mean models,
##' \code{cumres.lvmfit} for supported Gaussian structural equation models,
##' and \code{cumres.ordreg} for ordinal models fitted by \code{lava}.
##' The default method accepts residuals, their parameter derivatives and
##' coefficient influence contributions for other justified model structures.
##'
##' By default, each ordering defines a separate diagnostic. These marginal tests are not
##' a global omnibus check of every aspect of a fitted model. Saved simulation
##' columns from different orderings must not be combined as a joint null draw:
##' multipliers are paired to sorted ranks separately for each ordering. The
##' ordinal omnibus mode instead combines cutpoints for one ordering using
##' shared subject multipliers. Use \code{\link{cumres_joint}} for simultaneous
##' inference across prespecified orderings of a linear/generalized linear
##' model or a justified default residual process.
##'
##' For cumulative processes, the Sup statistic is the largest absolute
##' observed process value. The L2 statistic is the norm of the step process
##' on its evaluation grid. The historically named \code{KS} and \code{CvM}
##' components contain the corresponding Monte Carlo \emph{p-values}, not
##' those observed statistics. A simulated statistic equal to the observed
##' statistic is included in the tail probability. The legacy marginal p-value
##' resolution is \eqn{1/R}; a simulated p-value of zero is not proof of an
##' impossible null. Joint simulations use the plus-one rule described in
##' \code{\link{cumres_joint}}.
##' Examples use few simulations to keep them fast; increase \code{R} when
##' assessing small tail probabilities or the stability of a reported result.
##'
##' For default/linear/GLM/SEM processes, the default evaluation limit,
##' \code{subset.max=500}, limits distinct grid
##' values, never the fitted observations contributing to a process. Tied
##' observations are accumulated together. Set a seed before calling
##' \code{cumres} to reproduce simulations. Stored curves are intended for
##' visual comparison; they are not simultaneous confidence bands.
##' Positive-bandwidth grids include both entries \eqn{x_i} and exits
##' \eqn{x_i+b}, through the final zero tail. A finite grid cap approximates
##' the full-domain statistic. Prespecify bandwidth, ordering and grid-selection
##' rules for formal testing; response-driven selection requires additional
##' statistical adjustment.
##'
##' This is an asymptotic multiplier approximation. Rare counts, severe class
##' imbalance, high leverage or little effective information can invalidate a
##' Gaussian process approximation even when the nominal sample is large and
##' the fitter reports convergence. No universal finite-sample guarantee or
##' exact small-count calibration is supplied.
##' @returns
##' For one diagnostic, an S3 object of class \code{cumres}, represented by a
##' list with the following common components:
##' \describe{
##'   \item{\code{W}}{A list of observed numeric process vectors, one per
##'   ordering variable or ordinal aggregation component.}
##'   \item{\code{What}}{A list of matrices of saved null processes. Rows
##'   match the corresponding \code{W} and \code{x}; columns are saved draws.
##'   With \code{plots=0}, the matrices have zero columns.}
##'   \item{\code{x}}{A list of sorted numeric evaluation grids. Different
##'   processes can have different numbers of grid values.}
##'   \item{\code{KS}, \code{CvM}}{Numeric vectors of Sup and L2 Monte Carlo
##'   p-values, respectively, in the same order as \code{W}.}
##'   \item{\code{variable}}{Character labels for the processes.}
##'   \item{\code{R}}{Number of null simulations used for the p-values.}
##'   \item{\code{n}}{Number of independent residual units. Grouped binomial
##'   diagnostics count expanded trials, rather than original fitted rows.}
##'   \item{\code{type}, \code{model}}{Labels describing the process and
##'   fitted model adapter.}
##' }
##' Normal/GLM methods additionally provide \code{nobs}, the number of
##' positive-weight fitted rows, and \code{response}, the response expression.
##' Normal/GLM and ordinal methods provide \code{informative}, a logical
##' information flag. \code{FALSE} denotes a structurally recognized
##' zero-information residual/CDF process; p-values are one. Built-in
##' ordinal aggregations then return zero curves, while custom aggregators
##' return their value at the zero vector. \code{NA} denotes numerically
##' unresolved score-space membership: computed curves are retained but
##' p-values are missing. Other methods may omit this metadata. Ordinal
##' objects additionally record their \code{aggregation} mode.
##'
##' Structural equation calls requesting multiple diagnostics can return
##' lists, including nested lists when supplied lists of formulas. Their
##' individual elements follow the same \code{cumres} contract.
##' @seealso \code{\link{cumres.glm}}, \code{\link{cumres.default}},
##' \code{\link{cumres.lvmfit}}, \code{\link{cumres.ordreg}},
##' \code{\link{plot.cumres}}.
##' @references D. Y. Lin, L. J. Wei and Z. Ying (2002).
##' \emph{Model-Checking Techniques Based on Cumulative Residuals}.
##' Biometrics 58, 1--12. \doi{10.1111/j.0006-341X.2002.00001.x}.
##' @examples
##' fit <- lm(mpg ~ wt, data = mtcars)
##' set.seed(11)
##' diagnostic <- fit |> cumres(variable = "wt", R = 99, plots = 5)
##' coef(diagnostic)
##' plot(diagnostic)
##' @export
cumres <- function(model,...) UseMethod("cumres")

##' Cumulative residual processes from residuals and influence contributions
##'
##' Construct a diagnostic from supplied residuals, residual derivatives and
##' coefficient influence contributions. This interface supports custom models
##' when the caller can justify the required conditional-centering and
##' independent-unit influence representation.
##'
##' @param model Model providing the default coefficients, residuals and
##' influences, or \code{NULL} when these quantities are supplied explicitly.
##' @param variable A finite real numeric ordering vector or matrix with one
##' row per independent residual unit. Matrix columns define separate processes.
##' @param r A real residual vector, or a function of \code{coef} and
##' \code{...} returning that vector. A \code{grad} attribute can supply the
##' residual gradient.
##' @param dr Residual gradient matrix, with one row per residual unit and one
##' column per fitted parameter, or a function of \code{coef} and \code{...}.
##' If omitted, use the residual's gradient attribute or differentiate the
##' residual function with \pkg{numDeriv}.
##' @param R Positive whole number of multiplier simulations.
##' @param b Nonnegative moving-window bandwidth. Zero gives cumulative
##' residuals; positive values give sums over the window \eqn{(t-b,t]}.
##' @param plots Nonnegative whole number of null processes to retain for
##' plotting. This is independent of \code{R}; zero retains no simulated curves.
##' @param simulation Null simulation coupling. The default \code{"marginal"}
##' preserves the original separate-ordering simulations. \code{"joint"}
##' shares multipliers by independent subject and adds family statistics and
##' simultaneous null reference envelopes. See \code{\link{cumres_joint}};
##' this option applies to the default and linear/GLM methods.
##' @param coef Finite real coefficient vector used to evaluate residual and
##' gradient callbacks. It may be \code{NULL} when residuals, derivatives and
##' influences are already supplied. A known mean with no fitted parameters
##' uses an empty vector and zero-column \code{dr} and \code{ic} matrices.
##' @param ic Finite real matrix of coefficient influence contributions, with
##' one row per residual unit and the same parameter columns as \code{dr}.
##' These are contributions to coefficient perturbations, not unscaled score
##' vectors. The scaling is
##' \eqn{\hat\beta-\beta_0 = \sum_i IC_i + o_p(n^{-1/2})}; typically
##' \eqn{IC_i=I_n^{-1}s_i}, with total sample information \eqn{I_n}.
##' An unscaled order-one influence function must be divided by \eqn{n}.
##' Matched zero-column matrices are supported; when no coefficients
##' or gradient columns exist, the default influence matrix has zero columns.
##' @param breakties Deprecated jitter width. Nonzero values warn; tied
##' observations are retained and evaluated together.
##' @param subset Function receiving sorted distinct values and arguments
##' \code{n.qt} and \code{n.eq}. It must return unique valid 1-based evaluation
##' positions within \code{subset.max}.
##' @param subset.max Positive whole-number cap on distinct evaluation values
##' per process. \code{Inf} or \code{NULL} evaluates all distinct values.
##' Every residual unit contributes irrespective of this cap.
##' @param ... Additional arguments passed to residual and gradient callbacks.
##' @returns A \code{cumres} object. See \code{\link{cumres}} for the common
##' result components and interpretation of the Monte Carlo p-values.
##' @seealso \code{\link{cumres}}, \code{\link{cumres.glm}}.
##' @examples
##' # Known zero mean: no fitted parameters or influence adjustment.
##' set.seed(12)
##' x <- seq(-2, 2, length.out = 80)
##' r <- rnorm(length(x))
##' diagnostic <- cumres(NULL, variable = cbind(x = x), r = r,
##'   dr = matrix(numeric(), length(x), 0), coef = numeric(),
##'   R = 99, plots = 5)
##' coef(diagnostic)
##' @export
cumres.default <- function(model, variable, r=stats::residuals(model), dr,
                           R=1e3L, b=0, plots=min(R,50L),
                           coef=stats::coef(model), ic=lava::iid(model), breakties=0,
                           subset=subset_t, subset.max=500L,
                           simulation=c("marginal","joint"), ...) {
    simulation <- match.arg(simulation)
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

    if (simulation=="joint")
        return(gof_joint_engine(model,variable,r,dr,ic,R,b,plots,subset,subset.max,vnames))

    ## Reserve historical uniform draws to preserve seeded Monte Carlo streams
    ## for untied models. Numerical jitter is no longer added to observations.
    invisible(runif(n))
    if (b==0) {
        idx <- lapply(seq_len(p),function(i) gof_grid(variable[,i],subset.max,subset))
        lens <- lengths(idx)
        ## Repeating final endpoints pads unequal grids for C++ without changing
        ## Sup or L2 statistics. Padding is trimmed from the returned processes.
        Idx <- vapply(idx,function(i) c(i,rep(tail(i,1),max(lens)-length(i)))-1L,integer(max(lens)))
        Idx <- matrix(Idx,ncol=p)
    } else {
        x <- lapply(seq_len(p),function(i) gof_window_grid(variable[,i],b,subset.max,subset))
        lens <- lengths(x)
    }
    cr <- new(CumRes,r,dr,ic)
    cr$order(variable)
    if (b==0) x <- lapply(seq_len(p),function(i) cr$inp[idx[[i]],i])
    W0 <- cr$obs()
    if (b==0) {
        project <- function(w) lapply(seq_len(p),function(i) w[idx[[i]],i])
    } else {
        upper <- lapply(seq_len(p),function(i) findInterval(x[[i]],cr$inp[,i]))
        lower <- lapply(seq_len(p),function(i) findInterval(x[[i]],cr$inp[,i]+b))
        project <- function(w) lapply(seq_len(p),function(i) {
            prefix <- c(0,w[,i])
            prefix[upper[[i]]+1L]-prefix[lower[[i]]+1L]
        })
    }
    W <- project(W0)
    stats <- function(w) unlist(lapply(seq_len(p),function(i) c(SupTest(w[[i]]),L2Test(w[[i]],x[[i]]))))
    if (b==0) {
        teststat <- cr$samplestat(R,Idx,FALSE)
    } else {
        ## A window is active from entry x through (but not including) exit
        ## x+b. Use full prefixes at both events even on a reduced grid.
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
