## Statistical calibration: under a correct model the p-values should be
## approximately uniform (5% rejection at the 5% level), and a misspecified
## model should be detected. These are slow Monte Carlo checks, so they only
## run when GOF_SLOW_TESTS=true, e.g.
##   GOF_SLOW_TESTS=true Rscript -e 'testthat::test_local()'

skip_if_not(identical(Sys.getenv("GOF_SLOW_TESTS"), "true"),
            "slow tests: set GOF_SLOW_TESTS=true")

reject_rate <- function(sim_fit, nsim = 200, R = 200) {
    p <- replicate(nsim, {
        g <- cumres(sim_fit(), R = R, plots = 1)
        c(g$KS, g$CvM)
    })
    rowMeans(p < 0.05)
}

## With nsim = 200 the rejection rate under H0 lies in [0.015, 0.09] with
## probability > 99%.
expect_calibrated <- function(rate) {
    expect_true(all(rate > 0.015 & rate < 0.09), info = paste(round(rate, 3), collapse = " "))
}

test_that("type I error is controlled, n <= 500 (all points)", {
    set.seed(61)
    rate <- reject_rate(function() {
        x <- rnorm(300); y <- rbinom(300, 1, plogis(-0.5 + x))
        glm(y ~ x, family = binomial())
    })
    expect_calibrated(rate)
})

test_that("type I error is controlled, n > 500 (subset path)", {
    set.seed(62)
    rate <- reject_rate(function() {
        x <- rnorm(1500); y <- rbinom(1500, 1, plogis(-0.5 + x))
        glm(y ~ x, family = binomial())
    }, nsim = 200)
    expect_calibrated(rate)
})

test_that("type I error is controlled for grouped binomial data", {
    set.seed(63)
    rate <- reject_rate(function() {
        x <- 1:10; s <- rbinom(10, 40, plogis(-2 + 0.4 * x))
        glm(cbind(s, 40 - s) ~ x, family = binomial())
    })
    expect_calibrated(rate)
})

test_that("a missing quadratic term is detected, n > 500 (subset path)", {
    set.seed(64)
    rate <- reject_rate(function() {
        x <- rnorm(1500); y <- rbinom(1500, 1, plogis(-1 + x + 0.5 * x^2))
        glm(y ~ x, family = binomial())
    }, nsim = 50)
    expect_true(all(rate > 0.8), info = paste(round(rate, 3), collapse = " "))
})
