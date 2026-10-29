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

    t0 <- lapply(seq_len(p), function(j) sort(variable[, j])[idx[[j]]])
    W <- lapply(seq_len(p), function(j) {
        if (b==0) return((cumsum(r[ord[[j]]])/sqrt(n))[idx[[j]]])
        vapply(t0[[j]],function(t) sum(r*(variable[,j]>t-b & variable[,j]<=t))/sqrt(n),numeric(1))
    })
    sample_at <- function(j,g) {
        if (b==0) return(one_process(ord[[j]],g)[idx[[j]]])
        o <- ord[[j]]
        vapply(t0[[j]],function(t) {
            in_window <- variable[o,j]>t-b & variable[o,j]<=t
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
