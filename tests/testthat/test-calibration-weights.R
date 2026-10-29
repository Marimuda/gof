skip_if_not(identical(Sys.getenv("GOF_SLOW_TESTS"),"true"),"set GOF_SLOW_TESTS=true")
weight_calibration <- list()
weight_rate <- function(generate,seed) {
    p <- matrix(NA_real_,200,4)
    for (i in seq_len(nrow(p))) {
        set.seed(seed+i)
        model <- generate()
        g <- cumres(model,R=199,plots=0)
        p[i,] <- c(g$KS,g$CvM)
    }
    rates <- colMeans(p<.05)
    weight_calibration[[as.character(seed)]] <<- list(pvalues=p,rejection=rates)
    expect_true(all(rates>=.005 & rates<=.12),info=paste(round(rates,4),collapse=" "))
}

test_that("inverse-variance weighted normal diagnostics are calibrated", {
    weight_rate(function() {
        d <- data.frame(x=rnorm(300),w=runif(300,.4,4))
        d$y <- 1+.5*d$x+rnorm(300)/sqrt(d$w)
        lm(y~x,d,weights=w)
    },9100)
})

test_that("Poisson response-mean precision weights are calibrated", {
    weight_rate(function() {
        d <- data.frame(x=runif(300,-1,1),exposure=runif(300,.8,8))
        d$count <- rpois(300,d$exposure*exp(.2+.5*d$x)); d$rate <- d$count/d$exposure
        suppressWarnings(glm(rate~x,d,family=poisson(),weights=exposure))
    },9400)
})

test_that("binomial proportions with trial counts and a probit link are calibrated", {
    weight_rate(function() {
        d <- data.frame(x=rnorm(60),trials=sample(8:20,60,TRUE))
        d$s <- rbinom(60,d$trials,pnorm(.2+.5*d$x))
        glm(s/trials~x,d,weights=trials,family=binomial("probit"))
    },9700)
})

if (nzchar(Sys.getenv("GOF_WEIGHT_CALIBRATION_REPORT")))
    saveRDS(weight_calibration,Sys.getenv("GOF_WEIGHT_CALIBRATION_REPORT"))
