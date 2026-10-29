# Information in the unaggregated ordinal CDF residual process. This helper
# uses fitted design/score structure, never the magnitude of observed curves
# or p-values. Root aggregation must happen after a certified zero vector is
# enforced, because a custom aggregator need not map zero to zero.
gof_ordinal_information <- function(model,ordering,grid=sort(unique(ordering))) {
    env <- model$up
    X <- env$X
    n <- nrow(X)
    p <- ncol(X)
    K <- env$K
    if (length(ordering)!=n || !is.numeric(ordering) || any(!is.finite(ordering)) ||
        !is.numeric(grid) || !length(grid) || any(!is.finite(grid)))
        stop("Ordinal information checks require finite fitted ordering values and a nonempty grid",call.=FALSE)
    if (any(tabulate(env$y,nbins=K)==0L))
        stop("A regular interior ordinal fit requires every modelled response category to be observed",call.=FALSE)
    beta <- stats::coef(model)
    k <- K-1L
    if (length(beta)!=k+p)
        stop("Incompatible ordinal coefficient and design dimensions",call.=FALSE)
    # Fresh raw cutpoints and fitted linear predictors: env$lp and env$alpha
    # can retain a perturbed state from prediction or numerical derivatives.
    alpha <- env$threshold(beta,K)
    lp <- matrix(rep(alpha,each=n),nrow=n)-drop(X%*%beta[k+seq_len(p)])
    q <- env$h(lp)
    derivative <- env$dh(lp)
    if (length(q)!=n*k || length(derivative)!=n*k)
        stop("Incompatible ordinal probability and link derivative dimensions",call.=FALSE)
    q <- matrix(q,nrow=n)
    derivative <- matrix(derivative,nrow=n)
    upper_rows <- which(env$y<K)
    lower_rows <- which(env$y>1L)
    upper_cdf <- q[cbind(upper_rows,env$y[upper_rows])]
    lower_cdf <- q[cbind(lower_rows,env$y[lower_rows]-1L)]
    upper_density <- derivative[cbind(upper_rows,env$y[upper_rows])]
    lower_density <- derivative[cbind(lower_rows,env$y[lower_rows]-1L)]
    active_probability <- rep(1,n)
    active_probability[upper_rows] <- upper_cdf
    active_probability[lower_rows] <- active_probability[lower_rows]-lower_cdf
    if (any(!is.finite(c(upper_cdf,lower_cdf,upper_density,lower_density))) ||
        any(c(upper_cdf,lower_cdf)<=0 | c(upper_cdf,lower_cdf)>=1) ||
        any(c(upper_density,lower_density)<=0) ||
        any(!is.finite(active_probability) | active_probability<=0))
        stop("Ordinal diagnostics require interior observed-category probabilities and positive link derivatives",call.=FALSE)

    # A sufficient numerical separation check for a strictly increasing CDF:
    # move an observation's upper cutpoint predictor up and its lower one down.
    # Every probability then has a nonnegative directional derivative; one
    # strict bound makes the likelihood derivative positive. Such a direction
    # rules out an interior score root. LP infeasibility is not a proof of
    # global MLE existence for every link or ill-conditioned design.
    D <- matrix(0,length(upper_rows)+length(lower_rows),k+p)
    ur <- seq_along(upper_rows)
    lr <- length(upper_rows)+seq_along(lower_rows)
    D[cbind(ur,env$y[upper_rows])] <- 1
    D[cbind(lr,env$y[lower_rows]-1L)] <- -1
    if (p) {
        D[ur,k+seq_len(p)] <- -X[upper_rows,,drop=FALSE]
        D[lr,k+seq_len(p)] <- X[lower_rows,,drop=FALSE]
    }
    if (gof_mle_direction(D,seq_len(nrow(D))))
        stop("Ordinal separation prevents a regular interior maximum-likelihood fit",call.=FALSE)

    if (!p) {
        # K-1 free cutpoints saturate the K marginal category probabilities.
        # Each fitted marginal CDF equals its empirical CDF, so its complete
        # residual sum (and an empty sum) is identically zero. This certificate
        # applies only to all-row/empty evaluation indicators, not arbitrary
        # within-sample orderings supplied to this lower-level helper.
        trivial <- vapply(grid,function(t) {
            selected <- ordering<=t
            all(selected) || !any(selected)
        },logical(1))
        return(!all(trivial))
    }
    if (K!=2L) return(TRUE)

    # With two categories the ordinal model is a binomial mean model for
    # 1{Y<=1}, with linear predictor alpha_1-X*beta and design [1,-X].
    # Reconstruct from the stored coefficients: predict()/gradient evaluation
    # can mutate env$lp and env$alpha, which are not reliable fitted state.
    Z <- cbind(1,-X)
    precision <- as.numeric(derivative/(q*(1-q)))
    gof_informative_processes(Z,precision,cbind(ordering=ordering),list(grid))[1L]
}
