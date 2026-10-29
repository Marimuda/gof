skip_if_not(identical(Sys.getenv("GOF_SLOW_TESTS"),"true"),
            "slow tests: set GOF_SLOW_TESTS=true")

bootstrap_calibration_records <- list()
bootstrap_record_rate <- function(label,p,upper=.12,lower=.005) {
    rates <- colMeans(p<=.05)
    bootstrap_calibration_records[[label]] <<- list(rates=rates,p=p)
    message(label," rejection rates: ",paste(round(rates,4),collapse=" "))
    expect_true(all(is.finite(p)))
    expect_true(all(rates>=lower & rates<=upper),info=paste(rates,collapse=" "))
    invisible(rates)
}

test_that("specified-null Monte Carlo is calibrated for sparse Poisson and grouped binomial", {
    set.seed(4201)
    n <- 80
    x <- seq(-1,1,length.out=n)
    mu <- exp(-3+.25*x)
    p <- replicate(200,{
        y <- rpois(n,mu)
        # Empty mean model is just a stored design/container; the specified
        # means, not an estimate from y, define both residuals and simulation.
        fit <- glm(y~0+offset(log(mu)),family=poisson())
        g <- cumres_bootstrap(fit,"x",data=data.frame(x=x),R=99,plots=0,
                              null="specified",mean=mu)
        c(g$KS,g$CvM)
    })
    bootstrap_record_rate("specified-sparse-poisson",t(p),lower=0)
    trials <- rep(c(2,5,10,20),10)
    x <- seq(-1,1,length.out=length(trials))
    probability <- plogis(-3+.5*x)
    p <- replicate(200,{
        s <- rbinom(length(x),trials,probability)
        fit <- glm(cbind(s,trials-s)~0+offset(qlogis(probability)),family=binomial())
        g <- cumres_bootstrap(fit,"x",data=data.frame(x=x),R=99,plots=0,
                              null="specified",mean=probability)
        c(g$KS,g$CvM)
    })
    bootstrap_record_rate("specified-sparse-grouped-binomial",t(p),lower=0)
})

test_that("fitted bootstrap is calibrated in bounded low-count designs", {
    set.seed(4202)
    n <- 100
    x <- seq(-1,1,length.out=n)
    exposure <- exp(.3*x)
    p <- replicate(120,{
        y <- rpois(n,exposure*.08)
        fit <- glm(y~1+offset(log(exposure)),family=poisson(),
                   control=glm.control(epsilon=1e-10,maxit=100))
        g <- cumres_bootstrap(fit,"x",data=data.frame(x=x),R=99,plots=0)
        c(g$KS,g$CvM)
    })
    bootstrap_record_rate("fitted-low-count-poisson",t(p),lower=0)
    trials <- rep(5,50)
    x <- seq(-1,1,length.out=length(trials))
    p <- replicate(120,{
        s <- rbinom(length(x),trials,.04)
        fit <- glm(cbind(s,trials-s)~1,family=binomial(),
                   control=glm.control(epsilon=1e-10,maxit=100))
        g <- cumres_bootstrap(fit,"x",data=data.frame(x=x),R=99,plots=0)
        c(g$KS,g$CvM)
    })
    bootstrap_record_rate("fitted-low-count-binomial",t(p),lower=0)
})

test_that("response bootstrap detects a missing mean term", {
    set.seed(4203)
    x <- seq(-1,1,length.out=150)
    p <- replicate(60,{
        y <- rpois(length(x),exp(-1+3*x))
        fit <- glm(y~1,family=poisson())
        g <- cumres_bootstrap(fit,"x",data=data.frame(x=x),R=99,plots=0)
        c(g$KS,g$CvM)
    })
    rates <- rowMeans(p<=.05)
    bootstrap_calibration_records[["missing-mean-term"]] <- list(rates=rates,p=t(p))
    message("Bootstrap missing-term power: ",paste(rates,collapse=" "))
    expect_true(all(rates>.8))
})

if (nzchar(Sys.getenv("GOF_BOOTSTRAP_CALIBRATION_REPORT")))
    saveRDS(bootstrap_calibration_records,Sys.getenv("GOF_BOOTSTRAP_CALIBRATION_REPORT"))
