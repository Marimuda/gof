## Grouped binomial responses, cbind(successes, failures).
##
## cumres expands grouped data to one 0/1 row per trial. The result must be
## identical to fitting the same model to the expanded data directly.

## Separate grouped/expanded fits must converge to the same numerical MLE.
## Diagnostics now preserve the original fit rather than silently refitting.

test_that("grouped data give the same result as the expanded 0/1 data", {
    d <- fx_grouped_data()
    m_grp <- glm(control=glm.control(epsilon=1e-14,maxit=100), cbind(s, f) ~ x, data = d, family = binomial())
    m_ind <- glm(control=glm.control(epsilon=1e-14,maxit=100), y ~ x, data = expand_grouped(d), family = binomial())
    set.seed(41); a <- cumres(m_grp, R = 100)
    set.seed(41); b <- cumres(m_ind, R = 100)
    expect_same_cumres(a, b)
})

test_that("grouped data with equal group sizes work", {
    ## apply() simplifies to a matrix when all groups have the same size
    d <- fx_grouped_data(sizes = rep(50, 8))
    m_grp <- glm(control=glm.control(epsilon=1e-14,maxit=100), cbind(s, f) ~ x, data = d, family = binomial())
    m_ind <- glm(control=glm.control(epsilon=1e-14,maxit=100), y ~ x, data = expand_grouped(d), family = binomial())
    set.seed(42); expect_no_error(a <- cumres(m_grp, R = 100))
    set.seed(42); b <- cumres(m_ind, R = 100)
    expect_same_cumres(a, b)
})

test_that("variables are reported and selected by their original names", {
    m <- fx_grouped()
    set.seed(43)
    g <- cumres(m, R = 50)
    expect_equal(g$variable, c("predicted", "x"))
    set.seed(43)
    g <- cumres(m, variable = "x", R = 50)
    expect_equal(g$variable, "x")
})

test_that("grouped model fitted with subset= works", {
    d <- fx_grouped_data()
    d$keep <- d$x != 3
    m <- glm(control=glm.control(epsilon=1e-14,maxit=100), cbind(s, f) ~ x, data = d, family = binomial(), subset = keep)
    m_ind <- glm(control=glm.control(epsilon=1e-14,maxit=100), y ~ x, data = expand_grouped(d[d$keep, ]), family = binomial())
    set.seed(44); expect_no_error(a <- cumres(m, R = 50))
    set.seed(44); b <- cumres(m_ind, R = 50)
    expect_same_cumres(a, b)
})

test_that("grouped model with an offset works", {
    d <- fx_grouped_data()
    d$o <- log(d$size) / 10
    m <- glm(control=glm.control(epsilon=1e-14,maxit=100), cbind(s, f) ~ x + offset(o), data = d, family = binomial())
    di <- expand_grouped(d)
    di$o <- rep(d$o, d$s + d$f)
    m_ind <- glm(control=glm.control(epsilon=1e-14,maxit=100), y ~ x + offset(o), data = di, family = binomial())
    set.seed(45); expect_no_error(a <- cumres(m, R = 50))
    set.seed(45); b <- cumres(m_ind, R = 50)
    expect_same_cumres(a, b)
})

test_that("cumres works when the model is fitted inside a function", {
    f <- function(dd) {
        mm <- glm(control=glm.control(epsilon=1e-14,maxit=100), cbind(s, f) ~ x, data = dd, family = binomial())
        cumres(mm, R = 50)
    }
    set.seed(46)
    expect_no_error(g <- f(fx_grouped_data(sizes = rep(30, 6))))
    expect_s3_class(g, "cumres")
})

test_that("non-canonical link (probit) with grouped data works", {
    d <- fx_grouped_data()
    m <- glm(control=glm.control(epsilon=1e-14,maxit=100), cbind(s, f) ~ x, data = d, family = binomial(link = "probit"))
    m_ind <- glm(control=glm.control(epsilon=1e-14,maxit=100), y ~ x, data = expand_grouped(d), family = binomial(link = "probit"))
    set.seed(47); a <- cumres(m, R = 50)
    set.seed(47); b <- cumres(m_ind, R = 50)
    expect_same_cumres(a, b)
})

test_that("prior weights are rejected with an informative error", {
    d <- fx_grouped_data()
    m <- glm(control=glm.control(epsilon=1e-14,maxit=100), s / size ~ x, data = d, weights = size, family = binomial())
    expect_error(cumres(m, R = 10), "Weight argument not supported")
})
