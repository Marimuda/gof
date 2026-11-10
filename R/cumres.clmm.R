#' Cumulative residual diagnostics for ordinal mixed models
#'
#' Check marginal cumulative-category residuals from an
#' [ordinal::clmm()] fit using one Gaussian multiplier per independent
#' top-level cluster.
#'
#' @param model A converged `ordinal::clmm` fit retaining its model frame.
#' @param variable One numeric ordering variable: `"predicted"` for the
#'   marginal fixed predictor, a fixed-design column, or a numeric covariate
#'   retained in the model frame. Conditional predictions and BLUPs are not
#'   supported orderings.
#' @param R Number of multiplier simulations for the Sup and L2 tests.
#' @param plot Whether to plot the result immediately.
#' @param aggregate Ordinal aggregation: `"omnibus"` (Euclidean norm across
#'   cumulative cutpoints), `"score"` (sum), `"cutpoints"` (signed components),
#'   or a function returning a fixed-length finite numeric vector.
#' @param plots Number of simulated curves to retain, possibly larger than `R`.
#' @param control Numerical settings from [gof_clmm_control()].
#' @param ... Arguments for [plot.cumres()] when `plot = TRUE`.
#' @returns A [cumres-object] with `multiplier`, `random_effects` and
#'   `fit_check` metadata. `n` counts independent blocks, not data rows.
#' @details
#' Single grouping factors and nesting chains, random intercepts and slopes,
#' all five cumulative links, and all `clmm` threshold structures are supported.
#' Crossed effects are rejected. Weights multiply each likelihood and residual
#' contribution; they are frequency/likelihood weights, not a survey-design
#' correction. Formula offsets must occur before random-effect terms because
#' some `ordinal` versions read later offsets from the wrong frame column.
#'
#' The supplied estimates initialize adaptive quadrature. The diagnostic uses
#' an estimating-equation root with quadrature nodes frozen at their final
#' adaptive locations. A second calculation with four more nodes checks
#' numerical sensitivity. This finite-quadrature root is not a claim of an
#' exact integrated-likelihood solution. The input fit is not modified.
#' Parameter shifts, node checks and covariance-boundary diagnostics are
#' retained in `fit_check`; actual boundary fits and failed checks are refused.
#' Weakly estimated interior covariances and large shifts warn and are reported.
#' Numerical refusals inherit `gof_clmm_error`. The condition's `code` and
#' `stage` distinguish input convergence, unsuccessful polishing, covariance
#' boundaries, unidentified information and quadrature failures. Its
#' `diagnostics` retain the available last parameter, covariance, information
#' and score state. Unsuccessful polishing does not certify a boundary
#' optimum. The boundary table reports `at_boundary` separately from
#' `identified`; no ordinary p-values are produced for these refusals.
#' This development method is experimental: the completed slow cohort failed
#' its prespecified size/refusal criteria. Broad mixed-model release
#' qualification is on hold; see the mixed-model guide for the retained results.
#' The diagnosed model uses the declared link. In ordinal 2026.7-26, its
#' scalar backend interchanges the two asymmetric Gumbel links; see
#' `vignette("ordinal-mixed-models", package = "gof")` for the verified input
#' refit and independent likelihood checks. The upstream reported objective
#' in `fit_check$loglik_clmm` need not be the same likelihood as gof's
#' `fit_check$loglik`; their difference is not a likelihood-ratio statistic.
#'
#' For cutpoint k, the residual is `w_i * (I(Y_i <= k) - P_ik)`, with P
#' integrated over the random effects. All rows in a top-level cluster share
#' a multiplier; a score/information correction accounts for estimated
#' parameters. The process scale is sqrt(m)/sum(w), where m is the number of
#' independent clusters. This is an asymptotic procedure: small cluster counts,
#' weakly identified covariances and heavy-tailed links can affect calibration.
#' The score aggregation can miss departures that preserve the ordinal mean.
#'
#' Ties enter together and the full distinct ordering grid is used. Fitted
#' predictor ordering is response-dependent and retains numerical and
#' statistical ordering limitations. Prefer a prespecified observed covariate
#' for formal checks. A numerically vanishing multiplier covariance withholds
#' p-values; it does not certify that a model is adequate.
#' @seealso [cumres.ordreg()], [gof_clmm_control()], [cumres_joint()].
#' @references Pan, Z. and Lin, D. Y. (2005). Goodness-of-Fit Methods for
#'   Generalized Linear Mixed Models. Biometrics 61, 1000--1009.
#'   \doi{10.1111/j.1541-0420.2005.00365.x}.
#' @examples
#' if (requireNamespace("ordinal", quietly = TRUE) &&
#'     requireNamespace("numDeriv", quietly = TRUE)) {
#'   data("wine", package = "ordinal")
#'   fit <- ordinal::clmm(rating ~ temp + contact + (1 | judge), data = wine)
#'   diagnostic <- suppressWarnings(cumres(fit, variable = "tempwarm",
#'                                         R = 19, plots = 2))
#'   as.data.frame(diagnostic)
#'   diagnostic$fit_check$shift
#' }
#' @export
cumres.clmm <- function(model,variable="predicted",R=100,plot=FALSE,
                        aggregate=c("omnibus","score","cutpoints"),plots=min(R,50L),
                        control=gof_clmm_control(),...) {
    gof_ordinal_check_dots(...)
    if (!requireNamespace("ordinal",quietly=TRUE))
        stop("Install the optional 'ordinal' package to diagnose clmm fits",call.=FALSE)
    if (!requireNamespace("numDeriv",quietly=TRUE))
        stop("Install the optional 'numDeriv' package to diagnose clmm fits",call.=FALSE)
    R <- gof_count(R,"R",1L)
    plots <- gof_count(plots,"plots")
    if (!is.logical(plot) || length(plot)!=1L || is.na(plot)) stop("plot must be TRUE or FALSE",call.=FALSE)
    aggregation <- gof_ordinal_aggregation(aggregate)
    problem <- gof_clmm_problem(model,control)
    if (!identical(problem$structure,"nesting-chain"))
        stop("Crossed or non-chain random effects are not supported by cumres.clmm; ",
             "use a single grouping factor or a nesting chain",call.=FALSE)
    gof_clmm_separation(problem)
    ## Resolve named orderings before the expensive numerical fit.
    gof_clmm_ordering(problem,variable,problem$theta)
    fit <- gof_polish(problem,gof_select_integrator(problem,control),control)
    ordering <- gof_clmm_ordering(problem,variable,fit$theta)
    marginal <- gof_clmm_marginal(problem,fit$theta,control$marginal.nodes)
    residual <- (outer(problem$y,seq_len(problem$K-1L),"<=")-marginal$cdf)*problem$weights
    gradient <- marginal$gradient*rep(problem$weights,problem$K-1L)
    multiplier <- gof_clmm_multiplier(problem,fit)
    informative <- gof_clmm_informative(residual,gradient,multiplier,ordering)
    re <- gof_clmm_unpack(problem,fit$theta)$re
    random_effects <- lapply(seq_along(problem$terms),function(i) {
        tm <- problem$terms[[i]]
        list(term=tm$label,group=tm$group,q=tm$q,levels=tm$n_levels,
             relation=problem$structure,covariance=re[[i]]$Sigma)
    })
    result <- gof_ordinal_engine(residual,gradient,multiplier,ordering,
        aggregate=aggregation$fn,aggregation=aggregation$aggregation,informative=informative,
        R=R,plots=plots,variable=variable,model="clmm",
        design_info=gof_design_info(unit="block",n_rows=problem$n,
            fit_parameters=length(fit$theta),estimated_parameters=length(fit$theta),
            correction_dimension=length(fit$theta),n_components=problem$K-1L),
        metadata=list(random_effects=random_effects,
            fit_check=fit[c("theta","loglik","standardized_score","integrator","nodes",
                "iterations","refreshes","early_refreshes","trace","shift","node_check",
                "boundary","diagnostics","initialization")]))
    result$fit_check$loglik_clmm <- as.numeric(stats::logLik(model))
    result$fit_check$ordering <- if (variable=="predicted") "fitted fixed predictor" else "observed covariate"
    if (isTRUE(plot)) graphics::plot(result,...)
    result
}

gof_clmm_ordering <- function(problem,variable,theta) {
    if (!is.character(variable) || length(variable)!=1L || is.na(variable) || !nzchar(variable))
        stop("Select one numeric clmm ordering variable",call.=FALSE)
    if (variable %in% c("conditional","blup","BLUP","ranef","fitted"))
        stop("Conditional predictions and BLUP orderings are response-dependent and unsupported",call.=FALSE)
    if (variable=="predicted") {
        up <- gof_clmm_unpack(problem,theta)
        x <- drop(problem$X %*% up$beta)+problem$offset
    } else if (variable %in% colnames(problem$X)) {
        x <- problem$X[,variable]
    } else {
        frame <- problem$source$data
        response <- names(frame)[1L]
        if (variable==response || !variable %in% names(frame))
            stop("Variable not found among the retained clmm covariates: ",variable,call.=FALSE)
        x <- frame[[variable]][problem$source$rows]
    }
    if (!is.numeric(x) || is.complex(x) || !is.null(dim(x)) || length(x)!=problem$n || any(!is.finite(x)))
        stop("clmm ordering covariates must be finite numeric vectors; use a named design column for factors",call.=FALSE)
    as.numeric(x)
}

## A fixed-effect direction that improves every observed category also
## improves each conditional likelihood, hence the integrated likelihood.
gof_clmm_separation <- function(problem) {
    y <- problem$y
    upper <- which(y<problem$K)
    lower <- which(y>1L)
    r <- ncol(problem$tJac)
    p <- ncol(problem$X)
    D <- matrix(0,length(upper)+length(lower),r+p)
    iu <- seq_along(upper)
    il <- length(upper)+seq_along(lower)
    D[iu,seq_len(r)] <- problem$tJac[y[upper],,drop=FALSE]
    D[il,seq_len(r)] <- -problem$tJac[y[lower]-1L,,drop=FALSE]
    if (p) {
        D[iu,r+seq_len(p)] <- -problem$X[upper,,drop=FALSE]
        D[il,r+seq_len(p)] <- problem$X[lower,,drop=FALSE]
    }
    if (gof_mle_direction(D,seq_len(nrow(D))))
        stop("Ordinal separation prevents a regular interior mixed-model fit",call.=FALSE)
    invisible(NULL)
}

## Evaluate the exact covariance of the linear multiplier process, up to
## roundoff. Zero covariance is indeterminate, never a structural adequacy
## certificate. Work one cutpoint at a time to bound memory.
gof_clmm_informative <- function(residual,gradient,multiplier,ordering) {
    n <- nrow(residual)
    o <- order(ordering)
    ends <- cumsum(rle(ordering[o])$lengths)
    f <- multiplier$factors[[1L]]
    m <- nrow(f$influence)
    index <- f$index[o]
    for (k in seq_len(ncol(residual))) {
        partial <- matrix(0,n,m)
        partial[cbind(seq_len(n),index)] <- residual[o,k]
        partial <- apply(partial,2L,cumsum)[ends,,drop=FALSE]
        G <- gradient[seq_len(n)+(k-1L)*n,,drop=FALSE][o,,drop=FALSE]
        G <- apply(G,2L,cumsum)[ends,,drop=FALSE]
        correction <- G %*% t(f$influence)
        scale <- pmax(rowSums(abs(partial)+abs(correction)),.Machine$double.xmin)
        error <- 64*.Machine$double.eps*(n+ncol(G)+m)*scale
        if (any(sqrt(rowSums((partial+correction)^2))>error)) return(TRUE)
    }
    NA
}
