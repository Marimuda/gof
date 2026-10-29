## Pure-R reference implementation of the cumulative residual GoF tests.
##
## Written for clarity, not speed, following Lin, Wei & Ying (2002):
##
##   W(t)    = n^{-1/2} sum_i 1{x_i <= t} r_i                    (observed)
##   W*(t)   = n^{-1/2} [ sum_i 1{x_i <= t} r_i g_i
##                        + eta(t)' sum_i ic_i g_i ]             (under H0)
##   eta(t)  = sum_i 1{x_i <= t} dr_i   (derivative of residuals wrt beta)
##
## with g_i iid N(0,1). The p-values are P(stat(W*) > stat(W)) estimated from
## R realisations, with stat = sup|W| or the L2 norm sqrt(int W(t)^2 dt).
##
## It consumes R's random number stream in exactly the same order as
## gof:::cumres.default (tie-breaking jitter, then the R simulated processes
## sharing one g per realisation across variables, then the plotted processes),
## so with the same seed the results must agree with the C++ engine.
##
## For n > subset.max the processes are only evaluated at a subset of the
## (1-based) positions in the sorted order of each variable.

ref_cumres <- function(r, dr, ic, variable, R = 100, plots = 5,
                       breakties = 1e-12, subset.max = 500L) {
    r <- unname(r); dr <- unname(as.matrix(dr)); ic <- unname(as.matrix(ic))
    n <- length(r)
    variable <- unname(as.matrix(variable))
    untie <- runif(n, 0, breakties)
    variable <- variable + untie
    p <- ncol(variable)
    ord <- lapply(seq_len(p), function(j) order(variable[, j]))
    idx <- rep(list(seq_len(n)), p)
    if (n > subset.max) {
        ## Evaluation points chosen by the package's own subset_t (it may
        ## draw random fill-up points, so it must be called at this point of
        ## the random stream). subset_t itself is tested in test-subset.R.
        ns <- ceiling(subset.max / 2)
        idx <- lapply(seq_len(p), function(j)
            gof:::subset_t(variable[, j], n.qt = ns, n.eq = ns))
    }

    one_process <- function(o, g) {
        w1 <- cumsum(r[o] * g)
        eta <- apply(dr[o, , drop = FALSE], 2, cumsum)
        B <- colSums(ic[o, , drop = FALSE] * g)
        (w1 + drop(eta %*% B)) / sqrt(n)
    }
    sup_stat <- function(w) max(abs(w))
    l2_stat <- function(w, t) sqrt(sum(c(diff(t), 0) * w^2))

    t0 <- lapply(seq_len(p), function(j) sort(variable[, j])[idx[[j]]])
    W <- lapply(seq_len(p), function(j) (cumsum(r[ord[[j]]]) / sqrt(n))[idx[[j]]])

    sims <- matrix(0, R, 2 * p)
    for (k in seq_len(R)) {
        g <- rnorm(n)
        for (j in seq_len(p)) {
            w <- one_process(ord[[j]], g)[idx[[j]]]
            sims[k, 2 * j - 1] <- sup_stat(w)
            sims[k, 2 * j] <- l2_stat(w, t0[[j]])
        }
    }
    What <- lapply(seq_len(p), function(j) matrix(0, length(idx[[j]]), plots))
    for (k in seq_len(plots)) {
        g <- rnorm(n)
        for (j in seq_len(p))
            What[[j]][, k] <- one_process(ord[[j]], g)[idx[[j]]]
    }
    KS <- sapply(seq_len(p), function(j) mean(sims[, 2 * j - 1] > sup_stat(W[[j]])))
    CvM <- sapply(seq_len(p), function(j) mean(sims[, 2 * j] > l2_stat(W[[j]], t0[[j]])))
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
