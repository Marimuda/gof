##' Cumulative residual diagnostics for an ordinal regression
##'
##' Accumulate ordinal residuals over a fitted covariate or the linear
##' predictor and compare the resulting process with multiplier simulations.
##'
##' @param model A finite, converged \code{lava::ordreg} fit.
##' @param variable Covariate name, or "predicted" for the linear predictor.
##' @param R Positive whole number of null processes to simulate.
##' @param plot Logical; display the resulting process after calculating it.
##' The diagnostic object is returned whether or not a plot is drawn.
##' @param aggregate One of \code{"omnibus"} (default), \code{"score"} or
##' \code{"cutpoints"}, or a custom aggregation function. Omnibus uses the
##' Euclidean norm of the joint cumulative-cutpoint process at each grid value.
##' Score sums its components and checks the integer-coded category-score mean
##' only. Cutpoints returns the signed component processes. A custom function
##' must return a fixed-length finite real vector and be continuous at the
##' limiting null process for the simulated approximation to be justified.
##' @param plots Nonnegative whole number of simulated processes to retain.
##' With zero, the returned simulated-curve matrices have zero columns.
##' @param ... Graphical arguments passed to \code{plot.cumres} when
##' \code{plot=TRUE}. This method uses all distinct ordering values and does
##' not take \code{subset.max} or moving-window arguments through \code{...}.
##' The linear/GLM \code{resolution} and \code{on_unresolved} controls are
##' rejected rather than silently ignored.
##' @returns A \code{cumres} object with one process for every component
##' returned by \code{aggregate}. A scalar aggregator produces one process;
##' vector outputs receive labels such as \code{x.1}, \code{x.2}. See
##' \code{\link{cumres}} for the common observed/simulated curve and p-value
##' components. The \code{aggregation} field records the selected mode.
##' \code{informative} records structurally recognized zero-information or
##' numerically unresolved binary score processes. Recognized raw CDF processes
##' are set to zero before aggregation and have p-values one; a custom function
##' is evaluated at that zero vector and may itself return a nonzero constant.
##' Numerically unresolved cases retain curves and have missing p-values.
##' Fitted-row \code{nobs} metadata are not added by this method.
##' @details Omnibus Sup is \eqn{\sup_t\|W(t)\|_2}; omnibus L2 is
##' \eqn{(\int\sum_k W_k(t)^2\,dt)^{1/2}}, across cumulative cutpoints for the
##' one selected ordering. The same subject multipliers and complete parameter
##' influence vector are used across cutpoints, retaining their covariance.
##' Continuity of the Euclidean norm transfers the joint process approximation
##' to these statistics. This is not a joint test over every model covariate.
##' Omnibus curves are nonnegative norms, not zero-centered signed Gaussian
##' processes. Use \code{aggregate="cutpoints"} for directional component plots.
##' Before version 1.1.0 the default was \code{sum}. Request \code{"score"} or
##' \code{sum} explicitly to reproduce that mean-only test. Mean-preserving
##' changes in category probabilities can be invisible to this score sum.
##' @note Predictions use the fitted coefficients and the fit's internal state
##' is restored after diagnostics. The upstream \code{lava::ordreg} fitter
##' ignores \code{offset}, \code{subset}, \code{weights} and \code{na.action}
##' arguments; fits recording these options are rejected to avoid diagnosing
##' an unintended model. Subset and complete the data before fitting. Ordinal
##' offsets require upstream fitting support and cannot be added by diagnostics.
##' @seealso \code{\link{cumres}}, \code{\link{plot.cumres}}.
##' @examples
##' if (requireNamespace("numDeriv", quietly = TRUE)) {
##'   set.seed(15)
##'   d <- data.frame(x = rnorm(120))
##'   d$y <- ordered(cut(d$x + rnorm(nrow(d)),
##'     breaks = c(-Inf, -0.5, 0.5, Inf)))
##'   fit <- lava::ordreg(y ~ x, data = d)
##'   diagnostic <- cumres(fit, variable = "x", R = 49, plots = 3)
##'   coef(diagnostic)
##'   plot(diagnostic)
##'   # Inspect signed cumulative-cutpoint components as well.
##'   components <- cumres(fit, variable = "x", aggregate = "cutpoints", R = 49)
##'   as.data.frame(components)
##' }
##' @export
cumres.ordreg <- function(model,variable="predicted",R=100,plot=FALSE,
                          aggregate=c("omnibus","score","cutpoints"),plots=min(R,50L),...) {
    gof_ordinal_check_dots(...)
    R <- gof_count(R,"R",1L)
    plots <- gof_count(plots,"plots")
    if (!is.logical(plot) || length(plot)!=1L || is.na(plot)) stop("plot must be TRUE or FALSE")
    if (!requireNamespace("numDeriv",quietly=TRUE)) stop("package 'numDeriv' not available")
    aggregation <- gof_ordinal_aggregation(aggregate)
    aggregate <- aggregation$fn
    aggregation <- aggregation$aggregation
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
    if (any(!is.finite(residual)) || any(!is.finite(gradient)) || any(!is.finite(ic)))
        stop("Ordinal residuals, gradients and influences must be finite")
    ## One multiplier per fitted row; the engine sorts rows by the ordering.
    multiplier <- gof_multiplier_design(list(list(name="rows",index=seq_len(n),influence=ic)),
                                        n_units=n,divisor=sqrt(n))
    result <- gof_ordinal_engine(residual,gradient,multiplier,ordering=x,aggregate=aggregate,
                                 aggregation=aggregation,informative=informative,R=R,plots=plots,
                                 variable=variable,model="ordreg",
                                 design_info=gof_design_info(unit="subject",n_rows=n,
                                     fit_parameters=length(coef(model)),estimated_parameters=length(coef(model)),
                                     correction_dimension=ncol(ic),n_components=k))
    if (isTRUE(plot)) graphics::plot(result,...)
    result
}
