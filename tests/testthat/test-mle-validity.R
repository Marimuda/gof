test_that("finite IRLS floors do not qualify no-event or one-class estimated models", {
    d <- data.frame(x=seq(-1,1,length.out=100),y=0)
    for (f in list(binomial(),binomial("probit"),binomial("cloglog"),poisson())) {
        m <- suppressWarnings(glm(y~1,d,family=f,control=glm.control(maxit=100)))
        expect_true(m$converged)
        expect_true(all(is.finite(coef(m))))
        expect_error(cumres(m,variable="x",data=d,R=5,plots=0),
                     "separation|mean-domain boundary")
    }
    d$y <- 1
    m <- suppressWarnings(glm(y~1,d,family=binomial(),control=glm.control(maxit=100)))
    expect_true(m$converged)
    expect_error(cumres(m,variable="x",data=d,R=5,plots=0),"separation")
})

test_that("complete and quasi separation are determined from design geometry", {
    d <- data.frame(x=seq(-1,1,length.out=100))
    d$y <- as.numeric(d$x>0)
    m <- suppressWarnings(glm(y~x,d,family=binomial(),control=glm.control(maxit=100)))
    expect_true(m$converged)
    expect_error(cumres(m,R=5,plots=0),"separation")
    q <- data.frame(x=rep(c(-1,0,0,1),10),y=rep(c(0,0,1,1),10),o=rep(.2,40))
    m <- suppressWarnings(glm(y~x+offset(o),q,family=binomial(),control=glm.control(maxit=100)))
    expect_error(cumres(m,R=5,plots=0),"separation")
    q$w <- rep(c(1,2,3,4),10)
    m <- suppressWarnings(glm(y~x+offset(o),q,weights=w,family=binomial(),control=glm.control(maxit=100)))
    expect_error(cumres(m,R=5,plots=0),"separation")
    X <- cbind(1,c(-1,0,0,1)); y <- c(0,0,1,1)
    expect_true(gof:::gof_mle_direction(X,c(which(y>0),-which(y<1))))
    expect_true(gof:::gof_mle_direction(sweep(X,2,c(1,1e6),"*"),c(which(y>0),-which(y<1))))
})

test_that("mixed grouped responses and zero-weight rows have the right separation geometry", {
    d <- data.frame(x=c(-1,0,1),s=c(0,5,10),f=c(10,5,0),o=c(.1,-.2,.3))
    m <- suppressWarnings(glm(cbind(s,f)~x+offset(o),d,family=binomial(),
                              control=glm.control(maxit=100)))
    expect_error(cumres(m,variable="x",R=5,plots=0),"separation")
    d$s <- c(1,5,9); d$f <- 10-d$s
    m <- glm(cbind(s,f)~x+offset(o),d,family=binomial())
    expect_no_error(cumres(m,variable="x",R=5,plots=0))
    m <- glm(s/10~x+offset(o),d,weights=rep(10,3),family=binomial())
    expect_no_error(cumres(m,variable="x",R=5,plots=0))
    d <- rbind(d,data.frame(x=1e6,s=0,f=0,o=0))
    m <- suppressWarnings(glm(cbind(s,f)~x+offset(o),d,family=binomial()))
    expect_no_error(cumres(m,variable="x",R=5,plots=0))
    e <- data.frame(x=rep(-1:1,4),y=rep(c(0,0,1,1),3),w=1)
    e <- rbind(e,data.frame(x=100,y=1,w=0))
    m <- glm(y~x,e,weights=w,family=binomial())
    expect_no_error(cumres(m,variable="x",R=5,plots=0))
})

test_that("constrained-link finite boundary optima are not mistaken for regular fits", {
    d <- data.frame(x=0:2,y=c(0,1,1),w=c(20,1,1))
    m <- suppressWarnings(glm(y~x,d,weights=w,family=poisson("identity"),start=c(.01,1),
                              control=glm.control(maxit=200,epsilon=1e-12)))
    expect_true(m$converged)
    expect_error(cumres(m,variable="x",R=5,plots=0),"interior maximum-likelihood mean fit")
    d <- data.frame(x=0:2,y=c(1,0,1),w=c(1,1,20))
    m <- suppressWarnings(glm(y~x,d,weights=w,family=binomial("log"),start=c(-2,.5),
                              control=glm.control(maxit=200,epsilon=1e-12)))
    expect_true(m$converged)
    expect_error(cumres(m,variable="x",R=5,plots=0),"interior maximum-likelihood mean fit")
})

test_that("a boundary mean root can have zero score and positive observed curvature", {
    X <- cbind(1,0:3); y <- c(0,2,2,1)
    beta <- c(0,5/6)
    mu <- drop(X%*%beta)
    row_score <- c(-1,y[-1]/mu[-1]-1)
    H <- crossprod(X[-1,,drop=FALSE],X[-1,,drop=FALSE]*as.numeric(y[-1]/mu[-1]^2))
    expect_equal(drop(crossprod(X,row_score)),c(0,0),tolerance=1e-14)
    expect_true(all(eigen(H,symmetric=TRUE,only.values=TRUE)$values>0))
    # A numerical floor just inside the domain must still fail qualification.
    expect_error(gof:::gof_mle_interior(X,c(1e-14,5/6),rep(0,4),
                    "poisson.identity",y,NULL,rep(1,4)),"interior maximum-likelihood mean fit")
})

test_that("noncanonical Gaussian means require a finite locally regular root", {
    d <- data.frame(x=seq(-1,1,length.out=100),y=seq(-2,-1,length.out=100))
    m <- suppressWarnings(glm(y~1,d,family=gaussian("log"),start=0,
                              control=glm.control(maxit=100)))
    expect_true(m$converged)
    expect_error(cumres(m,variable="x",data=d,R=5,plots=0),"interior maximum-likelihood mean fit")
    d$y <- rep(c(-1,1),50)
    m <- suppressWarnings(glm(y~1,d,family=gaussian("inverse"),start=1,
                              control=glm.control(maxit=100)))
    expect_true(m$converged)
    expect_error(cumres(m,variable="x",data=d,R=5,plots=0),"interior maximum-likelihood mean fit")
    # Dispersion cancels from the linear mean influence: preserve exact fits.
    d$y <- 2+d$x
    expect_no_error(cumres(lm(y~x,d),variable="x",R=5,plots=0))
})

test_that("ordinary interior fits and known means retain their supported paths", {
    families <- c(lapply(c("identity","log","inverse"),gaussian),
                  lapply(c("logit","probit","cauchit","cloglog","log"),binomial),
                  lapply(c("log","identity","sqrt"),poisson))
    for (f in families) {
        set.seed(13001)
        d <- data.frame(x=runif(180,.2,.8))
        eta <- if (f$family=="binomial" && f$link=="log") -1+.2*d$x else 1+.2*d$x
        mu <- f$linkinv(eta)
        d$y <- switch(f$family,gaussian=rnorm(nrow(d),mu,.05),
                      binomial=rbinom(nrow(d),1,mu),poisson=rpois(nrow(d),mu))
        m <- glm(y~x,d,family=f,control=glm.control(epsilon=1e-12,maxit=100))
        expect_no_error(cumres(m,variable="x",R=5,plots=0))
    }
    d <- data.frame(y=rep(0,20),x=seq_len(20),o=rep(0,20))
    known <- glm(y~0+offset(o),d,family=poisson())
    expect_true(known$boundary)
    expect_no_error(cumres(known,variable="x",data=d,R=5,plots=0))
    d$o <- rep(qlogis(.3),20)
    known <- glm(y~0+offset(o),d,family=binomial())
    expect_no_error(cumres(known,variable="x",data=d,R=5,plots=0))
})
