ordinal_options_fixture <- function(n=150L,link="probit") {
    set.seed(510)
    d <- data.frame(x=rnorm(n),g=factor(rep(c("A","B"),length.out=n)))
    latent <- .6*d$x+.3*(d$g=="B")+if (link=="probit") rnorm(n) else rlogis(n)
    d$y <- ordered(cut(latent,c(-Inf,-.6,.6,Inf)))
    list(data=d,model=lava::ordreg(y~x+g,d,family=binomial(link)))
}

test_that("ordinal named and fitted predictors retain complete tie groups", {
    f <- ordinal_options_fixture(); m <- f$model
    cp <- predict(m,p=coef(m),type="cumulative")
    residual <- outer(m$up$y,seq_len(ncol(cp)),"<=")-cp
    x <- m$up$X[,"gB"]
    g <- cumres(m,variable="gB",R=20,plots=3)
    expected <- vapply(sort(unique(x)),function(t) sum(residual[x<=t,,drop=FALSE])/sqrt(length(x)),numeric(1))
    expect_equal(g$W[[1]],expected)
    expect_equal(unname(g$x[[1]]),c(0,1))
    expect_equal(g$n,length(x))
    pdf(NULL); on.exit(dev.off())
    expect_no_error(cumres(m,variable="x",R=5,plots=2,plot=TRUE,col=NULL))
})

test_that("ordinal prediction/gradient evaluation preserves the supplied fit state", {
    f <- ordinal_options_fixture(); m <- f$model
    p <- coef(m); p[length(p)] <- p[length(p)]*.6
    invisible(predict(m,p=p,type="cumulative"))
    saved <- as.list(m$up,all.names=TRUE)
    set.seed(511); a <- cumres(m,R=20,plots=3)
    expect_identical(as.list(m$up,all.names=TRUE),saved)
    expected <- predict(m,p=coef(m),type="cumulative")
    x <- m$up$lp[,1]
    set.seed(511); b <- cumres(m,R=20,plots=3)
    expect_equal(a$x[[1]],sort(unique(x)))
    expect_equal(a$W,b$W)
    expect_equal(a$What,b$What)
    expect_equal(a$KS,b$KS)
    expect_equal(a$CvM,b$CvM)
    saved <- as.list(m$up,all.names=TRUE)
    expect_error(cumres(m,variable=NA_character_,R=2),"one ordinal")
    expect_identical(as.list(m$up,all.names=TRUE),saved)
})

test_that("ordinal binary, logit, intercept-only and pre-subset fits run", {
    f <- ordinal_options_fixture(link="logit")
    expect_no_error(cumres(f$model,R=5,plots=0))
    d <- f$data[f$data$x>-.4,,drop=FALSE]
    d$y <- droplevels(d$y)
    m <- lava::ordreg(y~x,d,family=binomial("logit"))
    g <- cumres(m,variable="x",R=5,plots=8)
    expect_equal(g$n,nrow(d))
    expect_equal(ncol(g$What[[1]]),8L)
    m0 <- lava::ordreg(y~1,d)
    g0 <- cumres(m0,R=5,plots=0)
    expect_length(g0$W[[1]],1L)
    expect_equal(g0$CvM,1)
})

test_that("invalid ordinal controls and aggregation outputs are explicit errors", {
    f <- ordinal_options_fixture(); m <- f$model
    expect_error(cumres(m,plot=NA,R=2),"plot must")
    expect_error(cumres(m,variable="unknown",R=2),"Variable not found")
    expect_error(cumres(m,R=2,aggregate=NULL),"aggregate must be")
    for (aggregate in list(function(x) numeric(),function(x) NA_real_,function(x) matrix(x,ncol=1),function(x) x+1i))
        expect_error(cumres(m,R=2,aggregate=aggregate),"fixed-length finite numeric vector")
    calls <- 0L
    changing <- function(x) {calls <<- calls+1L; if (calls==1L) 1 else 1:2}
    expect_error(cumres(m,R=2,aggregate=changing),"fixed-length")
    bad <- m; bad$opt$convergence <- 1L
    expect_error(cumres(bad,R=2),"converged ordinal")
    bad <- m; bad$up <- new.env(parent=parent.env(m$up))
    for (nm in ls(m$up,all.names=TRUE)) assign(nm,get(nm,m$up),bad$up)
    bad$up$y[1] <- NA_real_
    expect_error(cumres(bad,R=2),"fitted response")
})

test_that("ignored upstream ordinal fit arguments cannot silently validate a different model", {
    f <- ordinal_options_fixture(); d <- f$data
    offset_fit <- lava::ordreg(y~x,d,offset=rep(.2,nrow(d)))
    subset_fit <- lava::ordreg(y~x,d,subset=x>0)
    expect_error(cumres(offset_fit,R=2),"ignores these fitting options: offset")
    expect_error(cumres(subset_fit,R=2),"ignores these fitting options: subset")
})
