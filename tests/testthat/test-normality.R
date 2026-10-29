normality_regression_fixture <- function(seed=3301,n=140L) {
    set.seed(seed)
    data <- data.frame(x=seq(-1,2,length.out=n))
    data$y <- .3+.7*data$x+stats::rnorm(n,sd=1.2)
    lava::estimate(lava::lvm(y~x),data)
}

test_that("normality location-scale correction matches an independent Gaussian regression derivation", {
    fit <- normality_regression_fixture()
    grid <- seq(.01,.99,length.out=51)
    setup <- gof:::gof_normality_setup(fit,NULL,grid,1e-3)
    data <- model.frame(fit)
    variance <- unname(lava::pars(fit)["y~~y"])
    z <- stats::qnorm(grid)
    expected <- cbind(stats::dnorm(z)/sqrt(variance),
        stats::dnorm(z)*mean(data$x)/sqrt(variance),
        stats::dnorm(z)*z/(2*variance))
    expect_equal(setup$delta,expected,tolerance=1e-7,ignore_attr=TRUE)
    X <- cbind(1,data$x)
    information <- matrix(0,3,3)
    information[1:2,1:2] <- crossprod(X)/variance
    information[3,3] <- nrow(data)/(2*variance^2)
    expect_equal(setup$information,information,tolerance=1e-7)
    expect_lt(setup$score.norm,1e-5)
})

test_that("normality simulations and tests match an independent score-corrected ECDF oracle", {
    fit <- normality_regression_fixture()
    data <- model.frame(fit)
    n <- nrow(data); grid <- seq(.01,.99,length.out=51)
    variance <- unname(lava::pars(fit)["y~~y"])
    X <- cbind(1,data$x); quantiles <- stats::qnorm(grid)
    delta <- cbind(stats::dnorm(quantiles)/sqrt(variance),
        stats::dnorm(quantiles)*mean(data$x)/sqrt(variance),
        stats::dnorm(quantiles)*quantiles/(2*variance))
    information <- matrix(0,3,3)
    information[1:2,1:2] <- crossprod(X)/variance
    information[3,3] <- n/(2*variance^2)
    ecdf <- function(z) vapply(quantiles,function(t) mean(z<=t),numeric(1L))
    stat <- function(w) c(max(abs(w)),sqrt(sum(diff(grid)*head(w,-1L)^2)))
    set.seed(3302)
    oracle <- replicate(64L,{
        e <- stats::rnorm(n)*sqrt(variance)
        score <- c(drop(crossprod(X,e))/variance,sum(e^2-variance)/(2*variance^2))
        sqrt(n)*(ecdf(e/sqrt(variance))-grid+drop(delta %*% solve(information,score)))
    })
    fitted.residual <- data$y-drop(X %*% lava::pars(fit)[1:2])
    expected.observed <- sqrt(n)*(ecdf(fitted.residual/sqrt(variance))-grid)
    null.stats <- t(apply(oracle[,1:59,drop=FALSE],2,stat))
    expected.p <- (1+colSums(sweep(null.stats,2,stat(expected.observed),">=")))/60
    set.seed(3302)
    actual <- cumres_normality(fit,R=59,plots=5,grid=grid)
    expect_equal(actual$W[[1L]],expected.observed,tolerance=1e-8)
    expect_equal(actual$What[[1L]],oracle[,60:64,drop=FALSE],tolerance=1e-7)
    expect_equal(c(actual$KS,actual$CvM),expected.p,tolerance=1e-8,ignore_attr=TRUE)
})

test_that("latent and conditional projections match Gaussian covariance algebra", {
    set.seed(3303)
    model <- lava::lvm(list(c(y1,y2,y3)~eta,eta~x))
    lava::latent(model) <- ~eta
    fit <- lava::estimate(model,lava::sim(model,220))
    data <- model.frame(fit); ii <- lava::index(fit)
    moments <- lava::moments(fit,data=data,conditional=TRUE)
    Y <- ii$endogenous
    for (variable in c("eta","y1")) {
        result <- cumres_normality(fit,variable=variable,R=19,plots=0)
        raw <- moments$P[match(variable,ii$vars),,drop=FALSE] %*%
            t(moments$IAi[match(Y,ii$vars),,drop=FALSE]) %*% solve(moments$C)
        expected <- raw/sqrt(drop(raw %*% moments$C %*% t(raw)))
        expect_equal(result$normality$projection,expected,ignore_attr=TRUE)
        expect_equal(result$normality$residuals,
            (as.matrix(data[Y])-moments$xi) %*% t(expected),ignore_attr=TRUE)
        expect_equal(drop(expected %*% moments$C %*% t(expected)),1)
    }
    setup <- gof:::gof_normality_setup(fit,"eta",seq(.1,.9,length.out=9),1e-3)
    expect_equal(setup$information,lava::information(fit),tolerance=1e-6,ignore_attr=TRUE)
    expected.score <- as.numeric(lava::score(fit,p=lava::pars(fit)))
    actual.score <- gof:::gof_normality_score(as.matrix(data[Y])-moments$xi,setup)
    expect_lt(max(abs(actual.score-expected.score)),1e-6)
})

test_that("normality preserves its calibration and meaning through result interfaces", {
    set.seed(3304)
    result <- normality_regression_fixture() |> cumres_normality(R=29,plots=0)
    expect_s3_class(result,c("cumres_normality","cumres"))
    before <- .Random.seed
    table <- as.data.frame(result)
    expect_equal(table$calibration,rep("gaussian-one-step",2L))
    expect_equal(table$diagnostic,rep("normality",2L))
    expect_equal(table$residual,rep("marginal",2L))
    expect_equal(as.data.frame(summary(result)),table)
    expect_equal(.Random.seed,before)
    expect_output(print(result),"normality diagnostic")
    expect_output(print(summary(result)),"One-step Gaussian calibration")
    path <- tempfile(fileext=".pdf")
    grDevices::pdf(path)
    tryCatch(expect_identical(invisible(plot(result,legend="none")),result),
        finally={grDevices::dev.off(); unlink(path)})
})

test_that("normality simulation count and retained curves are independent", {
    fit <- normality_regression_fixture()
    set.seed(3305); a <- cumres_normality(fit,R=39,plots=0)
    set.seed(3305); b <- cumres_normality(fit,R=39,plots=4)
    expect_equal(c(a$KS,a$CvM),c(b$KS,b$CvM))
    expect_equal(dim(a$What[[1L]]),c(401L,0L))
    expect_equal(dim(b$What[[1L]]),c(401L,4L))
})

test_that("normality rejects unsupported data and invalid numerical contracts", {
    fit <- normality_regression_fixture()
    expect_error(cumres_normality(stats::lm(y~x,model.frame(fit))),"lava Gaussian")
    expect_error(cumres_normality(fit,variable="x"),"endogenous or latent")
    expect_error(cumres_normality(fit,variable=c("x","y")),"endogenous or latent")
    for (grid in list(c(0,.5),c(.5,1),c(.5,.2),c(.5,.5),NA_real_,matrix(c(.1,.9))))
        expect_error(cumres_normality(fit,grid=grid),"strictly increasing probabilities")
    expect_error(cumres_normality(fit,R=0),"whole number")
    expect_error(cumres_normality(fit,plots=-1),"whole number")
    expect_error(cumres_normality(fit,score.tol=0),"positive finite")
    bad <- fit; bad$weights <- rep(1,nrow(model.frame(fit)))
    expect_error(cumres_normality(bad),"Weighted")
    bad <- fit; bad$cluster <- seq_len(nrow(model.frame(fit)))
    expect_error(cumres_normality(bad),"independent subjects")
    bad <- fit; bad$opt$convergence <- 1L
    expect_error(cumres_normality(bad),"not converged")
    bad <- fit; bad$data$model.frame$y[1L] <- NA_real_
    expect_error(cumres_normality(bad),"complete finite")
    bad <- fit; bad$estimator <- "glm"
    expect_error(cumres_normality(bad),"Gaussian maximum-likelihood")
})

test_that("normality refuses a stale optimum or boundary structural variance", {
    fit <- normality_regression_fixture()
    bad <- fit
    bad$opt$par[1L] <- bad$opt$estimate[1L] <- bad$opt$par[1L]+.2
    expect_error(cumres_normality(bad),"interior root")
    bad <- fit
    bad$opt$par[3L] <- bad$opt$estimate[3L] <- 0
    expect_error(cumres_normality(bad),"positive-definite observed covariance")
})

test_that("normality detects pronounced skewed residual distributions", {
    set.seed(3306)
    data <- data.frame(x=stats::rnorm(500))
    data$y <- .3+.7*data$x+(stats::rchisq(500,df=1)-1)/sqrt(2)
    fit <- lava::estimate(lava::lvm(y~x),data)
    result <- cumres_normality(fit,R=199,plots=0)
    expect_true(all(c(result$KS,result$CvM)<.01))
})

test_that("known Gaussian structural parameters need no nuisance correction", {
    model <- lava::lvm(y~x)
    lava::regression(model,y~x) <- .5
    lava::intercept(model,~y) <- 0
    lava::covariance(model,y~y) <- 1
    set.seed(3307)
    data <- lava::sim(model,80L)
    fit <- structure(list(model=model,model0=model,
        data=list(model.frame=data,n=nrow(data)),estimator="gaussian",
        opt=list(estimate=numeric(),par=numeric(),convergence=0L)),class="lvmfit")
    result <- cumres_normality(fit,R=19,plots=2)
    expect_equal(result$normality$information,matrix(0,0,0))
    expect_equal(result$normality$delta,matrix(0,401,0))
    expect_equal(result$normality$residuals,matrix(data$y-.5*data$x,ncol=1L),ignore_attr=TRUE)
    expect_equal(result$normality$score.norm,0)
})
