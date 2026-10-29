ordinal_saturated_binary <- function(link="logit") {
    d <- data.frame(x=rep(0:1,each=100),
                    y=ordered(c(rep(1,80),rep(2,20),rep(1,30),rep(2,70))))
    list(data=d,model=lava::ordreg(y~x,d,family=binomial(link)))
}

test_that("free ordinal cutpoints leave no intercept-only CDF information", {
    set.seed(6073)
    d <- data.frame(y=ordered(sample(1:3,150,replace=TRUE)))
    model <- lava::ordreg(y~1,d)
    expect_identical(gof:::gof_ordinal_information(model,rep(0,150)),FALSE)
    expect_identical(gof:::gof_ordinal_information(model,seq_len(150)),TRUE)
    for (aggregation in c("omnibus","score","cutpoints")) {
        set.seed(6074)
        g <- cumres(model,R=29,plots=2,aggregate=aggregation)
        expect_true(all(g$informative %in% FALSE))
        expect_true(all(unlist(g$W)==0))
        expect_true(all(unlist(g$What)==0))
        expect_true(all(c(g$KS,g$CvM)==1))
    }
})

test_that("custom ordinal aggregation acts on the certified zero CDF vector", {
    set.seed(6073)
    d <- data.frame(y=ordered(sample(1:3,150,replace=TRUE)))
    model <- lava::ordreg(y~1,d)
    aggregate <- function(z) c(2,7+sum(z^2))
    g <- cumres(model,R=11,plots=3,aggregate=aggregate)
    expect_identical(g$informative,c(FALSE,FALSE))
    expect_equal(unlist(g$W),c(2,7))
    expect_true(all(g$What[[1]]==2))
    expect_true(all(g$What[[2]]==7))
    expect_equal(c(g$KS,g$CvM),rep(1,4))
})

test_that("binary ordinal saturated design cells are exactly uninformative", {
    for (link in c("logit","probit")) {
        f <- ordinal_saturated_binary(link)
        model <- f$model
        expect_identical(gof:::gof_ordinal_information(model,f$data$x),FALSE)
        for (aggregation in c("omnibus","score","cutpoints")) {
            set.seed(6074)
            g <- cumres(model,variable="x",R=29,plots=2,aggregate=aggregation)
            expect_identical(g$informative,FALSE)
            expect_true(all(g$W[[1]]==0))
            expect_true(all(g$What[[1]]==0))
            expect_equal(c(g$KS,g$CvM),c(1,1))
        }
    }
})

test_that("ordinal information checks use fitted coefficients without mutating shared state", {
    f <- ordinal_saturated_binary("logit")
    model <- f$model
    altered <- coef(model); altered[1] <- altered[1]+.7
    invisible(predict(model,p=altered,type="cumulative"))
    saved <- as.list(model$up,all.names=TRUE)
    expect_identical(gof:::gof_ordinal_information(model,f$data$x),FALSE)
    expect_identical(as.list(model$up,all.names=TRUE),saved)
})

test_that("exact ordinal no-information handling preserves multiplier draw counts", {
    set.seed(6073)
    d <- data.frame(y=ordered(sample(1:3,150,replace=TRUE)))
    model <- lava::ordreg(y~1,d)
    set.seed(6074)
    g <- cumres(model,R=7,plots=9)
    after <- .Random.seed
    set.seed(6074)
    invisible(rnorm(150*9))
    expect_identical(.Random.seed,after)
    expect_equal(ncol(g$What[[1]]),9L)
})

test_that("binary ordinal existence checks reject a separating direction", {
    f <- ordinal_saturated_binary("logit")
    model <- f$model
    model$up$y <- ifelse(f$data$x==0,1L,2L)
    expect_error(gof:::gof_ordinal_information(model,f$data$x),"separation")
})

test_that("general ordinal category-bound separation is a sufficient rejection certificate", {
    set.seed(6075)
    d <- data.frame(x=rep(0:2,each=60))
    d$y <- ordered(rep(rep(1:3,each=20),3))
    model <- lava::ordreg(y~x,d,family=binomial("logit"))
    expect_identical(gof:::gof_ordinal_information(model,d$x),TRUE)
    ## Test the likelihood geometry on the same full-rank design with an
    ## exactly ordered response. This bypasses the upstream optimizer so the
    ## validity helper cannot rely on its convergence flag.
    model$up$y <- d$x+1L
    expect_error(gof:::gof_ordinal_information(model,d$x),"Ordinal separation")
})
