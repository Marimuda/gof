skip_if_not(identical(Sys.getenv("GOF_SLOW_TESTS"),"true"),"set GOF_SLOW_TESTS=true")

normality_calibration_records <- list()

test_that("one-step SEM normality has controlled size for marginal, conditional and latent processes", {
    model <- lava::lvm(list(c(y1,y2,y3)~eta,eta~x))
    lava::latent(model) <- ~eta
    pvalues <- matrix(NA_real_,200L,6L)
    grid <- seq(.001,.999,length.out=201L)
    for (i in seq_len(nrow(pvalues))) {
        set.seed(3400+i)
        data <- lava::sim(model,300L)
        fit <- lava::estimate(model,data)
        for (j in seq_len(3L)) {
            set.seed(4400+3L*i+j)
            result <- cumres_normality(fit,variable=list(NULL,"y1","eta")[[j]],
                grid=grid,R=199L,plots=0L)
            pvalues[i,2L*j-c(1L,0L)] <- c(result$KS,result$CvM)
        }
    }
    rates <- colMeans(pvalues<.05)
    limits <- stats::qbinom(c(.0005,.9995),nrow(pvalues),.05)/nrow(pvalues)
    normality_calibration_records$size <<- list(pvalues=pvalues,rates=rates,
        nominal=.05,acceptance=limits)
    message("SEM normality null rejection rates: ",paste(rates,collapse=" "))
    expect_true(all(rates>=limits[1L] & rates<=limits[2L]),
        info=paste(rates,collapse=" "))
})

test_that("projected SEM normality detects pronounced nonnormal structural errors", {
    model <- lava::lvm(list(c(y1,y2,y3)~eta,eta~x))
    lava::latent(model) <- ~eta
    pvalues <- matrix(NA_real_,50L,4L)
    for (i in seq_len(nrow(pvalues))) {
        for (alternative in 1:2) {
            set.seed(5400+2L*i+alternative)
            x <- stats::rnorm(400)
            latent.error <- if (alternative==1L) (stats::rchisq(400,1)-1)/sqrt(2) else stats::rnorm(400)
            eta <- .7*x+latent.error
            data <- data.frame(x=x,y1=eta+.4*stats::rnorm(400),
                y2=1.2*eta+.4*stats::rnorm(400),y3=.8*eta+.4*stats::rnorm(400))
            if (alternative==2L) data$y1 <- eta+.8*(stats::rchisq(400,1)-1)/sqrt(2)
            fit <- lava::estimate(model,data)
            result <- cumres_normality(fit,variable=c("eta","y1")[alternative],
                R=199L,plots=0L,grid=seq(.001,.999,length.out=201L))
            pvalues[i,2L*alternative-c(1L,0L)] <- c(result$KS,result$CvM)
        }
    }
    rates <- colMeans(pvalues<.05)
    normality_calibration_records$power <<- list(pvalues=pvalues,rates=rates)
    message("SEM normality alternative rejection rates: ",paste(rates,collapse=" "))
    expect_true(all(rates>.8),info=paste(rates,collapse=" "))
})

if (nzchar(Sys.getenv("GOF_NORMALITY_CALIBRATION_REPORT")))
    saveRDS(normality_calibration_records,Sys.getenv("GOF_NORMALITY_CALIBRATION_REPORT"))
