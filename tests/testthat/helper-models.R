## Shared model fixtures for the test suite.
##
## Each fixture returns a fitted model built from seeded simulated data, so the
## same models can be used for golden-value tests, reference-implementation
## tests and invariance tests.

fx_lm <- function(n = 100, seed = 1) {
    set.seed(seed)
    d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
    d$y <- 1 + d$x1 + 0.5 * d$x2 + rnorm(n)
    lm(y ~ x1 + x2, data = d)
}

fx_logistic <- function(n = 300, seed = 2) {
    set.seed(seed)
    d <- data.frame(x = rnorm(n))
    d$y <- rbinom(n, 1, plogis(-0.5 + d$x))
    glm(y ~ x, data = d, family = binomial())
}

fx_probit <- function(n = 300, seed = 3) {
    set.seed(seed)
    d <- data.frame(x = rnorm(n))
    d$y <- rbinom(n, 1, pnorm(0.2 + 0.8 * d$x))
    glm(y ~ x, data = d, family = binomial(link = "probit"))
}

fx_poisson <- function(n = 200, seed = 4) {
    set.seed(seed)
    d <- data.frame(x = runif(n, 0, 2))
    d$y <- rpois(n, exp(0.5 + 0.7 * d$x))
    glm(y ~ x, data = d, family = poisson())
}

## Grouped binomial data: cbind(successes, failures) with *unequal* group
## sizes (this layout works in gof 1.0.1, so it can be pinned by golden values).
fx_grouped_data <- function(sizes = c(40, 45, 50, 55, 60, 65, 70, 75), seed = 5) {
    set.seed(seed)
    x <- seq_along(sizes)
    s <- rbinom(length(sizes), sizes, plogis(-2 + 0.5 * x))
    data.frame(x = x, s = s, f = sizes - s, size = sizes)
}

fx_grouped <- function(...) {
    d <- fx_grouped_data(...)
    glm(cbind(s, f) ~ x, data = d, family = binomial())
}

## Logistic model with n > 500 and heavy ties in the covariate (like the
## Company_default example from the course, which triggered GitHub issue #1).
fx_large_ties <- function(n = 982, seed = 6) {
    set.seed(seed)
    d <- data.frame(group = sample(1:5, n, replace = TRUE), z = rnorm(n))
    d$y <- rbinom(n, 1, plogis(-2 + 0.3 * d$group + 0.5 * d$z))
    glm(y ~ group + z, data = d, family = binomial())
}

expand_grouped <- function(d) {
    data.frame(x=rep(d$x,d$s+d$f),
               y=unlist(lapply(seq_len(nrow(d)),function(i) rep(c(1,0),c(d$s[i],d$f[i])))))
}

expect_same_cumres <- function(a,b) {
    expect_equal(a$KS,b$KS)
    expect_equal(a$CvM,b$CvM)
    expect_equal(a$W,b$W)
    expect_equal(a$x,b$x)
    expect_equal(a$variable,b$variable)
}
