skip_if_not(identical(Sys.getenv("GOF_SLOW_TESTS"),"true"),"set GOF_SLOW_TESTS=true")

sem_ordinal_calibration <- list()
check_sem_ordinal_calibration <- function(rate) {
    ## 160 independent datasets: this catches gross distortion while allowing
    ## Monte Carlo uncertainty around a 5% nominal rejection probability.
    expect_true(all(rate>=.005 & rate<=.12),info=paste(round(rate,4),collapse=" "))
}

test_that("full SEM projections of latent and observed covariates are calibrated", {
    m <- lava::lvm(list(c(y1,y2,y3)~eta,eta~x)); lava::latent(m) <- ~eta
    pvalues <- matrix(NA_real_,160,4)
    for (i in seq_len(nrow(pvalues))) {
        set.seed(8100+i)
        d <- lava::sim(m,400)
        e <- lava::estimate(m,d)
        set.seed(8200+i); a <- cumres(e,"y1",x="eta",full=TRUE,R=199,plots=0)
        set.seed(8300+i); b <- cumres(e,"y1",x="y2",full=TRUE,R=199,plots=0)
        pvalues[i,] <- c(a$KS,a$CvM,b$KS,b$CvM)
    }
    rates <- colMeans(pvalues<.05)
    sem_ordinal_calibration$full_sem <<- list(pvalues=pvalues,rejection=rates)
    check_sem_ordinal_calibration(rates)
})

test_that("logit ordinal diagnostics with tied and predicted covariates are calibrated", {
    pvalues <- matrix(NA_real_,160,4)
    for (i in seq_len(nrow(pvalues))) {
        set.seed(8500+i)
        d <- data.frame(x=rnorm(300),g=factor(rep(c("A","B"),150)))
        d$y <- ordered(cut(.6*d$x+.3*(d$g=="B")+rlogis(300),c(-Inf,-.6,.6,Inf)))
        m <- lava::ordreg(y~x+g,d,family=binomial("logit"))
        set.seed(8600+i); a <- cumres(m,variable="gB",R=199,plots=0)
        set.seed(8700+i); b <- cumres(m,variable="predicted",R=199,plots=0)
        pvalues[i,] <- c(a$KS,a$CvM,b$KS,b$CvM)
    }
    rates <- colMeans(pvalues<.05)
    sem_ordinal_calibration$ordinal_logit <<- list(pvalues=pvalues,rejection=rates)
    check_sem_ordinal_calibration(rates)
})

if (nzchar(Sys.getenv("GOF_SEM_CALIBRATION_REPORT")))
    saveRDS(sem_ordinal_calibration,Sys.getenv("GOF_SEM_CALIBRATION_REPORT"))
