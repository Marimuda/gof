## Evaluation on a subset of points when n > subset.max (default 500).
##
## The subset is only a computational device: the process is evaluated at
## fewer points, but it must be the *same* process. So on the subset path the
## observed process must equal the full process at the selected points.

test_that("cumres works for n > 500 (GitHub issue #1)", {
    set.seed(31)
    expect_no_error(g <- cumres(fx_large_ties(), R = 50))
    expect_s3_class(g, "cumres")
    expect_equal(g$variable, c("predicted", "group", "z"))
    for (i in seq_along(g$W)) {
        expect_length(g$x[[i]], 500)
        expect_length(g$W[[i]], 500)
        expect_equal(nrow(g$What[[i]]), 500)
    }
    expect_true(all(g$KS >= 0 & g$KS <= 1))
    expect_true(all(g$CvM >= 0 & g$CvM <= 1))
})

test_that("subset process equals the full process at the selected points", {
    m <- fx_large_ties()
    set.seed(32); sub <- cumres(m, R = 10, plots = 1)
    set.seed(32); full <- cumres(m, R = 10, plots = 1, subset.max = Inf)
    for (i in seq_along(sub$W)) {
        pos <- match(sub$x[[i]], full$x[[i]])
        expect_false(anyNA(pos))
        ## The 1e-12 tie-breaking jitter does not separate every tied value
        ## (e.g. 'group' keeps a few exact duplicates), and matching by value
        ## is ambiguous there. Compare at the values that are unique.
        dup <- full$x[[i]] %in% full$x[[i]][duplicated(full$x[[i]])]
        ok <- !dup[pos]
        expect_gt(mean(ok), 0.8)
        expect_equal(sub$W[[i]][ok], full$W[[i]][pos][ok])
    }
})

test_that("the observed process is the cumulated residuals at the x-values", {
    m <- fx_large_ties()
    set.seed(33)
    g <- cumres(m, R = 10, plots = 1)
    r <- residuals(m, type = "response")
    z <- model.matrix(m)[, "z"]
    ## z is continuous: no ties, the jitter does not change the order
    Wz <- unname(cumsum(r[order(z)]) / sqrt(length(r)))
    pos <- match(g$x[[3]], sort(z + g$untie))
    expect_false(anyNA(pos))
    expect_equal(g$W[[3]], Wz[pos])
})

test_that("plot() works on the subset path", {
    set.seed(34)
    g <- cumres(fx_large_ties(), R = 20)
    pdf(NULL); on.exit(dev.off())
    expect_no_error(plot(g))
})

## --- subset_t: selection of evaluation points ---------------------------

test_that("subset_t returns all points when n is small", {
    expect_equal(gof:::subset_t(runif(100), n.qt = 60, n.eq = 60), 1:100)
})

test_that("subset_t returns sorted, unique, valid 1-based indices", {
    set.seed(35)
    t <- rexp(2000)
    idx <- gof:::subset_t(t, n.qt = 50, n.eq = 50)
    expect_length(idx, 100)
    expect_false(anyNA(idx))
    expect_false(is.unsorted(idx, strictly = TRUE))
    expect_true(all(idx >= 1 & idx <= length(t)))
    expect_true(length(t) %in% idx)   # the process is always evaluated at the end
})

test_that("subset_t includes the point nearest to each equidistant grid point", {
    set.seed(36)
    t <- sort(rexp(2000))
    idx <- gof:::subset_t(t, n.qt = 50, n.eq = 50)
    grid <- seq(t[1], t[length(t)], length.out = 50)
    nearest <- vapply(grid, function(g) which.min(abs(t - g)), integer(1))
    expect_true(all(nearest %in% idx))
})

test_that("subset_t includes the empirical quantiles", {
    set.seed(37)
    t <- sort(rnorm(3000))
    idx <- gof:::subset_t(t, n.qt = 40, n.eq = 40)
    q <- quantile(t, seq(0, 1, length.out = 40), type = 7)
    ## a quantile between two observations is represented by the last
    ## observation <= the quantile (the process is right-continuous)
    last_le <- vapply(q, function(v) max(which(t <= v)), integer(1))
    expect_true(all(last_le %in% idx))
})
