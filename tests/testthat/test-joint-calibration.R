skip_if_not(identical(Sys.getenv("GOF_SLOW_TESTS"),"true"),
            "slow tests: set GOF_SLOW_TESTS=true")

joint_calibration_records <- list()

test_that("joint family tests and Sup envelopes calibrate for normal models", {
    set.seed(88001)
    nsim <- 300L
    results <- replicate(nsim,{
        n <- 180
        d <- data.frame(x=rnorm(n),z=rnorm(n))
        d$y <- 1+d$x-.4*d$z+rnorm(n)
        g <- cumres_joint(lm(y~x+z,d),variable=c("x","z"),R=199,plots=0)
        c(g$joint$p.value,g$joint$observed[1L]>g$joint$critical)
    })
    rate <- rowMeans(results[1L:2L,,drop=FALSE]<=.05)
    envelope_rate <- mean(results[3L,])
    joint_calibration_records$normal <<- list(pvalues=t(results[1L:2L,,drop=FALSE]),
        rejection=rate,envelope.crossing=as.logical(results[3L,]),envelope.rate=envelope_rate)
    message("Joint normal rejection rates Sup/L2.range/envelope: ",
            paste(c(rate,envelope_rate),collapse=" "))
    expect_true(all(rate>=.01 & rate<=.10))
    expect_true(envelope_rate>=.01 && envelope_rate<=.10)
    expect_equal(results[1L,]<=.05,as.logical(results[3L,]))
})

test_that("joint family tests calibrate for binomial models and detect curvature", {
    set.seed(88002)
    nsim <- 200L
    null <- alternative <- matrix(0,2L,nsim)
    for (k in seq_len(nsim)) {
        n <- 250
        d <- data.frame(x=rnorm(n),z=rnorm(n))
        d$y <- rbinom(n,1,plogis(-.3+.5*d$x+.4*d$z))
        a <- cumres_joint(glm(y~x+z,d,family=binomial()),variable=c("x","z"),R=199,plots=0)
        null[,k] <- a$joint$p.value
        d$y <- rbinom(n,1,plogis(-.3+.5*d$x+.4*d$z+1.5*(d$x^2-1)))
        b <- cumres_joint(glm(y~x+z,d,family=binomial()),variable=c("x","z"),R=199,plots=0)
        alternative[,k] <- b$joint$p.value
    }
    rates <- rowMeans(null<=.05)
    power <- rowMeans(alternative<=.05)
    joint_calibration_records$binomial <<- list(pvalues=t(null),rejection=rates,
        alternative.pvalues=t(alternative),power=power)
    message("Joint binomial rejection rates: ",paste(rates,collapse=" "),
            "; curvature power: ",paste(power,collapse=" "))
    expect_true(all(rates>=.005 & rates<=.105))
    expect_true(all(power>.8))
})

if (nzchar(Sys.getenv("GOF_JOINT_CALIBRATION_REPORT")))
    saveRDS(joint_calibration_records,Sys.getenv("GOF_JOINT_CALIBRATION_REPORT"))
