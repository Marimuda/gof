## Statistical calibration: under a correct model the p-values should be
## approximately uniform (5% rejection at the 5% level), and a misspecified
## model should be detected. These are slow Monte Carlo checks, so they only
## run when GOF_SLOW_TESTS=true, e.g.
##   GOF_SLOW_TESTS=true Rscript -e 'testthat::test_local()'

skip_if_not(identical(Sys.getenv("GOF_SLOW_TESTS"), "true"),
            "slow tests: set GOF_SLOW_TESTS=true")

reject_rate <- function(sim_fit, nsim = 400, R = 200, ...) {
    extra <- list(...)
    p <- replicate(nsim, {
        ## replicate() itself introduces ..., containing its loop index.
        ## Capture the caller's options before entering that scope.
        g <- do.call(cumres,c(list(model=sim_fit(),R=R,plots=0),extra))
        c(g$KS, g$CvM)
    })
    result <- rowMeans(p < 0.05)
    message("Rejection rates (",nsim," datasets): ",paste(round(result,3),collapse=" "))
    result
}

## Each test checks 4 rejection rates (Sup and L2 test for 2 variables).
## With nsim = 400 a single rate under H0 has SE 0.011 and lies in
## [0.01, 0.10] with probability > 99.9%, so all rates of all tests pass
## jointly with high probability, while gross miscalibration is caught.
expect_calibrated <- function(rate) {
    calibration_records[[length(calibration_records)+1L]] <<- rate
    expect_true(all(rate >= 0.01 & rate <= 0.10), info = paste(round(rate, 3), collapse = " "))
}
calibration_records <- list()

test_that("type I error is controlled, n <= 500 (all points)", {
    set.seed(61)
    rate <- reject_rate(function() {
        x <- rnorm(300); y <- rbinom(300, 1, plogis(-0.5 + x))
        glm(y ~ x, family = binomial())
    })
    expect_calibrated(rate)
})

test_that("grouped calibration is invariant to predictor units", {
    p <- matrix(0,200,4)
    for (i in 1:200) {
        set.seed(7000+i)
        d <- data.frame(x=1:10)
        d$s <- rbinom(10,40,plogis(-2+.4*d$x)); d$f <- 40-d$s; d$large <- d$x*1e6
        a <- glm(cbind(s,f)~x,d,family=binomial())
        b <- glm(cbind(s,f)~large,d,family=binomial())
        set.seed(9000+i); ga <- cumres(a,variable="x",R=200,plots=0)
        set.seed(9000+i); gb <- cumres(b,variable="large",R=200,plots=0)
        p[i,] <- c(ga$KS,ga$CvM,gb$KS,gb$CvM)
    }
    expect_equal(p[,1:2],p[,3:4])
    rates <- colMeans(p<.05)
    message("Grouped unit-invariance rejection rates: ",paste(rates,collapse=" "))
    expect_calibrated(rates)
})

test_that("type I error is controlled with all points above 500", {
    set.seed(66)
    rate <- reject_rate(function() {
        x <- rnorm(650); y <- rbinom(650,1,plogis(-.5+x))
        glm(y~x,family=binomial())
    },nsim=200,subset.max=Inf)
    expect_calibrated(rate)
})

test_that("type I error is controlled for probit and normal models", {
    set.seed(67)
    probit <- reject_rate(function() {
        x <- rnorm(300); y <- rbinom(300,1,pnorm(.2+.8*x))
        glm(y~x,family=binomial(link="probit"))
    },nsim=200)
    expect_calibrated(probit)
    normal <- reject_rate(function() {
        x <- rnorm(300); y <- 1+.5*x+rnorm(300)
        lm(y~x)
    },nsim=200)
    expect_calibrated(normal)
})

test_that("type I error is controlled for moving windows", {
    set.seed(68)
    rate <- reject_rate(function() {
        x <- rnorm(300); y <- rbinom(300,1,plogis(-.5+x))
        glm(y~x,family=binomial())
    },nsim=200,b=.7,subset.max=99)
    expect_calibrated(rate)
})

test_that("type I error is controlled for ordinal processes", {
    set.seed(69)
    p <- replicate(200, {
        x <- rnorm(300)
        y <- ordered(cut(x+rnorm(300),breaks=c(-Inf,-.5,.5,Inf)))
        m <- lava::ordreg(y~x)
        g <- cumres(m,variable="x",R=200,plots=0)
        c(g$KS,g$CvM)
    })
    rate <- rowMeans(p<.05)
    message("Ordinal rejection rates: ",paste(rate,collapse=" "))
    expect_calibrated(rate)
})

test_that("type I error is controlled for structural-equation residuals", {
    set.seed(70)
    m <- lava::lvm(y~x)
    p <- replicate(200, {
        d <- lava::sim(m,300)
        e <- lava::estimate(m,d)
        g <- cumres(e,"y",x="x",R=200,plots=0)
        c(g$KS,g$CvM)
    })
    rate <- rowMeans(p<.05)
    message("SEM rejection rates: ",paste(rate,collapse=" "))
    expect_calibrated(rate)
})

test_that("type I error is controlled, n > 500 (subset path)", {
    set.seed(62)
    rate <- reject_rate(function() {
        x <- rnorm(1500); y <- rbinom(1500, 1, plogis(-0.5 + x))
        glm(y ~ x, family = binomial())
    })
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

test_that("type I error is controlled for poisson with an offset", {
    set.seed(65)
    rate <- reject_rate(function() {
        x <- runif(300, 0, 2); e <- runif(300, 0.5, 20)
        y <- rpois(300, e * exp(-1 + 0.5 * x))
        glm(y ~ x + offset(log(e)), family = poisson())
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

if (nzchar(Sys.getenv("GOF_CALIBRATION_REPORT")))
    saveRDS(calibration_records,Sys.getenv("GOF_CALIBRATION_REPORT"))
