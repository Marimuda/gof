bootstrap_oracle <- function(y,mean,ordering,trials=rep(1,length(y)),grid=NULL,b=0) {
    ordering <- as.matrix(ordering)
    if (is.null(grid)) grid <- lapply(seq_len(ncol(ordering)),function(j)
        sort(unique(c(ordering[,j],if (b>0) ordering[,j]+b))))
    n <- if (all(trials==1)) length(y) else sum(trials)
    W <- lapply(seq_len(ncol(ordering)),function(j)
        vapply(grid[[j]],function(t) {
            inside <- ordering[,j]<=t
            if (b>0) inside <- inside & ordering[,j]+b>t
            sum((y-trials*mean)[inside])/sqrt(n)
        },numeric(1)))
    statistic <- unlist(lapply(seq_along(W),function(j) {
        w <- W[[j]]; x <- grid[[j]]
        c(max(abs(w)),sqrt(sum(c(diff(x),0)*w^2)))
    }),use.names=FALSE)
    list(W=W,statistic=statistic,grid=grid)
}

test_that("specified Poisson simulation matches a direct indicator oracle", {
    d <- data.frame(x=c(0,0,1,2,4),z=c(2,0,1,0,3),y=c(1,0,2,0,1))
    m <- glm(y~x,d,family=poisson())
    mu <- c(.2,.4,.8,1.1,.7)
    for (bandwidth in c(0,1.5)) {
        set.seed(4101)
        g <- cumres_bootstrap(m,c("x","z"),data=d,null="specified",mean=mu,
                              R=29,plots=4,b=bandwidth,subset.max=Inf)
        observed <- bootstrap_oracle(d$y,mu,d[c("x","z")],b=bandwidth)
        expect_equal(g$x,observed$grid)
        expect_equal(g$W,observed$W)
        set.seed(4101)
        simulations <- replicate(29,{
            y <- rpois(nrow(d),mu)
            bootstrap_oracle(y,mu,d[c("x","z")],b=bandwidth)
        },simplify=FALSE)
        statistics <- t(vapply(simulations,function(s) s$statistic,numeric(4)))
        expect_equal(unname(g$bootstrap$statistics),statistics,tolerance=1e-13)
        for (j in 1:2) expect_equal(unname(g$What[[j]]),
            vapply(simulations[1:4],function(s) s$W[[j]],numeric(length(g$x[[j]]))),tolerance=1e-13)
        # Independent summation and L2 reduction can round mathematical ties
        # differently. Bracket only those ties; require every resolved tail
        # comparison to agree rather than accepting a fixed p-value tolerance.
        difference <- sweep(statistics,2,observed$statistic,"-")
        error <- 128*.Machine$double.eps*pmax(1,abs(statistics),
                    matrix(abs(observed$statistic),29,4,byrow=TRUE))
        lower <- (1+colSums(difference>error))/30
        upper <- (1+colSums(difference>= -error))/30
        actual <- as.vector(rbind(g$KS,g$CvM))
        expect_true(all(actual>=lower & actual<=upper))
        expect_equal(g$bootstrap$status,rep("specified",29))
        expect_identical(g$calibration,"specified-null-monte-carlo")
    }
})

test_that("fitted Poisson simulation refits offsets and fitted rows without replaying calls", {
    make_fit <- function() {
        set.seed(4102)
        d <- data.frame(x=seq(-1.5,1.5,length.out=85),exposure=seq(.3,3,length.out=85))
        d$extra <- sin(d$x)
        d$y <- rpois(nrow(d),d$exposure*exp(.1+.4*d$x))
        d$y[7] <- NA
        d$weight <- 1
        d$weight[12] <- 0
        list(model=glm(y~x+offset(log(exposure)),d,subset=x> -1.4,
                       weights=weight,family=poisson()),data=d)
    }
    saved <- make_fit()
    m <- saved$model
    mf <- model.frame(m)
    active <- model.weights(mf)>0
    original <- saved$data[match(rownames(mf),rownames(saved$data)),][active,]
    mu <- as.numeric(fitted(m))[active]
    set.seed(4103)
    g <- cumres_bootstrap(m,c("x","extra"),data=saved$data,R=19,plots=3,subset.max=Inf)
    observed <- bootstrap_oracle(original$y,mu,original[c("x","extra")])
    expect_equal(g$W,observed$W,tolerance=1e-10)
    expect_equal(g$nobs,nrow(original))
    set.seed(4103)
    simulations <- replicate(19,{
        original$yb <- rpois(nrow(original),mu)
        fitted <- glm(yb~x+offset(log(exposure)),original,family=poisson(),
                      control=glm.control(epsilon=1e-10,maxit=100))
        bootstrap_oracle(original$yb,as.numeric(fitted(fitted)),original[c("x","extra")])$statistic
    })
    expect_equal(unname(g$bootstrap$statistics),t(simulations),tolerance=1e-9)
    expect_equal(g$bootstrap$failed,0L)
})

test_that("grouped binomial and weighted proportions preserve trials", {
    set.seed(4104)
    d <- data.frame(x=seq(-1,1,length.out=30),trials=rep(c(10,20,30),10))
    d$s <- rbinom(nrow(d),d$trials,plogis(-.5+.6*d$x))
    d$f <- d$trials-d$s
    d$w <- rep(c(1,2),15)
    grouped <- glm(cbind(s,f)~x,d,weights=w,family=binomial())
    proportion <- glm(I(s/trials)~x,d,weights=trials*w,family=binomial())
    set.seed(4105)
    a <- cumres_bootstrap(grouped,"x",R=29,plots=3)
    set.seed(4105)
    b <- cumres_bootstrap(proportion,"x",R=29,plots=3)
    expect_equal(a$bootstrap$trials,d$trials*d$w)
    expect_equal(a$n,sum(d$trials*d$w))
    expect_equal(a$W,b$W,tolerance=1e-9)
    expect_equal(a$bootstrap$statistics,b$bootstrap$statistics,tolerance=1e-8)
    expect_equal(a$KS,b$KS)
    expect_equal(a$CvM,b$CvM)
    mu <- as.numeric(fitted(grouped))
    set.seed(4105)
    oracle <- replicate(29,{
        d$yb <- rbinom(nrow(d),d$trials*d$w,mu)
        refit <- glm(cbind(yb,trials*w-yb)~x,d,family=binomial(),
                     control=glm.control(epsilon=1e-10,maxit=100))
        bootstrap_oracle(d$yb,as.numeric(fitted(refit)),d["x"],trials=d$trials*d$w)$statistic
    })
    expect_equal(unname(a$bootstrap$statistics),t(oracle),tolerance=1e-9)
})

test_that("group residual curves equal individually expanded Bernoulli curves", {
    d <- data.frame(x=c(0,1,3),s=c(2,0,3),f=c(1,2,2))
    m <- glm(cbind(s,f)~1,d,family=binomial())
    mu <- c(.2,.3,.8)
    g <- cumres_bootstrap(m,"x",data=d,null="specified",mean=mu,R=9,plots=0)
    trials <- d$s+d$f
    y <- unlist(lapply(seq_along(trials),function(i) c(rep(1,d$s[i]),rep(0,d$f[i]))))
    x <- rep(d$x,trials)
    oracle <- bootstrap_oracle(y,rep(mu,trials),matrix(x))
    expect_equal(g$W,oracle$W)
    expect_equal(g$x,oracle$grid)
    expect_equal(g$n,length(y))
})

test_that("analytical boundary fits retain sparse zero and all-success draws", {
    d <- data.frame(x=seq_len(40),y=c(1,rep(0,39)))
    m <- glm(y~1,d,family=poisson())
    set.seed(4106)
    g <- cumres_bootstrap(m,"x",data=d,R=99,plots=4)
    expect_true(any(g$bootstrap$status=="extended-zero"))
    expect_equal(g$bootstrap$failed,0L)
    expect_true(all(is.finite(g$bootstrap$statistics)))
    expect_true(all(g$bootstrap$statistics[g$bootstrap$status=="extended-zero",]==0))
    zero <- glm(I(y*0)~1,d,family=poisson())
    z <- cumres_bootstrap(zero,"x",data=d,R=9,plots=0)
    expect_equal(z$bootstrap$observed.status,"extended-zero")
    expect_equal(z$KS,1,ignore_attr=TRUE)
    expect_equal(z$CvM,1,ignore_attr=TRUE)
    expect_identical(z$informative,FALSE)
    successes <- suppressWarnings(glm(I(rep(1,40))~1,d,family=binomial()))
    one <- cumres_bootstrap(successes,"x",data=d,R=9,plots=0)
    expect_equal(one$bootstrap$observed.status,"extended-one")
    expect_true(all(one$bootstrap$status=="extended-one"))
    expect_equal(one$W,list(rep(0,40)))
})

test_that("unresolved refits are retained with withheld p-values and explicit bounds", {
    set.seed(4107)
    d <- data.frame(x=seq(-1,1,length.out=80))
    d$y <- rpois(nrow(d),exp(.5+.4*d$x))
    m <- glm(y~x,d,family=poisson())
    expect_warning(g <- cumres_bootstrap(m,"x",R=9,plots=2,control=list(maxit=1)),
                   "9 of 9 bootstrap refits unresolved")
    expect_true(all(g$bootstrap$status=="unresolved"))
    expect_true(all(is.na(c(g$KS,g$CvM))))
    expect_equal(length(g$bootstrap$reasons),9L)
    expect_true(all(grepl("converge",unlist(g$bootstrap$reasons))))
    expect_equal(g$bootstrap$tail.bounds$p.lower,rep(.1,2))
    expect_equal(g$bootstrap$tail.bounds$p.upper,rep(1,2))
    expect_equal(ncol(g$What[[1]]),0L)
    table <- as.data.frame(g)
    expect_equal(table$failed,rep(9L,2))
    expect_equal(table$p.lower,rep(.1,2))
    expect_output(print(g),"9 of 9 refits unresolved")
    expect_equal(summary(g)$bootstrap$failed,9L)
})

test_that("storage and extraction do not change response draws", {
    d <- data.frame(x=seq_len(20),y=rep(c(0,1,0,2),5))
    m <- glm(y~x,d,family=poisson())
    set.seed(4108)
    a <- cumres_bootstrap(m,"x",R=19,plots=0)
    seed <- .Random.seed
    set.seed(4108)
    b <- cumres_bootstrap(m,"x",R=19,plots=4)
    expect_identical(.Random.seed,seed)
    expect_identical(a$bootstrap$statistics,b$bootstrap$statistics)
    expect_identical(a$KS,b$KS)
    table <- as.data.frame(b)
    expect_equal(as.data.frame(summary(b)),table)
    expect_identical(.Random.seed,seed)
    expect_equal(table$calibration,rep("parametric-bootstrap",2))
    expect_equal(table$p.lower,table$p.upper)
    expect_equal(table$p.value,table$p.lower)
    expect_true(all(table$p.value>=1/20))
    expect_equal(b$bootstrap$curve.indices,1:4)
})

test_that("fitted saturated cells have structurally zero curves", {
    d <- data.frame(x=rep(1:2,each=20),y=rep(c(1,2,3,2),10))
    m <- glm(y~factor(x),d,family=poisson())
    g <- cumres_bootstrap(m,"x",data=d,R=19,plots=3)
    expect_identical(g$informative,FALSE)
    expect_equal(g$W,list(c(0,0)))
    expect_true(all(g$What[[1]]==0))
    expect_equal(unname(g$KS),1)
    expect_equal(unname(g$CvM),1)
})

test_that("design-certified zero processes stay zero when other fits are unresolved", {
    d <- data.frame(x=rep(1:2,each=10),y=rep(c(1,rep(0,9)),2))
    for (f in list(poisson(),binomial())) {
        m <- glm(y~factor(x),d,family=f)
        set.seed(4110)
        expect_warning(g <- cumres_bootstrap(m,"x",data=d,R=39,plots=2),
                       "informative p-values withheld")
        expect_gt(g$bootstrap$failed,0)
        expect_identical(g$informative,FALSE)
        expect_equal(g$KS,1)
        expect_equal(g$CvM,1)
        expect_equal(g$bootstrap$tail.bounds$p.lower,c(1,1))
        expect_equal(g$bootstrap$tail.bounds$p.upper,c(1,1))
        expect_true(all(g$bootstrap$statistics==0))
        expect_equal(length(g$bootstrap$status),39)
    }
})

test_that("specified-null sparse distributions use response simulation rather than Gaussian noise", {
    d <- data.frame(x=seq_len(1000),y=0)
    m <- suppressWarnings(glm(y~1,d,family=poisson()))
    tiny <- cumres_bootstrap(m,"x",data=d,R=29,plots=0,null="specified",mean=rep(1e-12,1000))
    expect_equal(unname(tiny$KS),1)
    expect_equal(unname(tiny$CvM),1)
    wrong <- cumres_bootstrap(m,"x",data=d,R=99,plots=0,null="specified",mean=rep(1,1000))
    expect_equal(unname(wrong$KS),.01)
    expect_equal(unname(wrong$CvM),.01)
})

test_that("unsupported sampling and response-dependent orderings fail explicitly", {
    d <- data.frame(x=seq_len(20),y=rep(c(0,1,0,2),5))
    m <- glm(y~x,d,family=poisson())
    expect_error(cumres_bootstrap(m,"predicted",R=9),"fixed orderings")
    expect_error(cumres_bootstrap(m,"y",R=9),"tested response")
    expect_error(cumres_bootstrap(m,"absent",R=9),"Variable not found")
    expect_error(cumres_bootstrap(m,"x",R=9,null="specified"),"prespecified mean")
    expect_error(cumres_bootstrap(m,"x",R=9,mean=rep(1,20)),"only used")
    expect_error(cumres_bootstrap(m,"x",R=9,plots=10),"must not exceed")
    expect_error(cumres_bootstrap(m,"x",R=0),"whole number")
    expect_error(cumres_bootstrap(m,"x",R=9,b=-1),"bandwidth")
    expect_error(cumres_bootstrap(m,"x",R=9,subset.max=0),"whole number")
    weighted <- glm(y~x,d,weights=rep(2,20),family=poisson())
    expect_error(cumres_bootstrap(weighted,"x",R=9),"unit prior weights")
    expect_error(cumres_bootstrap(lm(y~x,d),"x",R=9),"glm fitted")
    probit <- glm(I(y>0)~x,d,family=binomial("probit"))
    expect_error(cumres_bootstrap(probit,"x",R=9),"binomial/logit")
    unretained <- glm(y~x,d,family=poisson(),model=FALSE)
    expect_error(cumres_bootstrap(unretained,"x",R=9),"stored model frame")
})

test_that("specified binomial simulation preserves fixed trials and includes exact ties", {
    d <- data.frame(x=c(0,0,1,2),s=c(1,0,1,0),trials=c(2,1,2,1))
    m <- glm(cbind(s,trials-s)~1,d,family=binomial())
    mu <- rep(.5,4)
    set.seed(4109)
    g <- cumres_bootstrap(m,"x",data=d,null="specified",mean=mu,R=39,plots=3)
    observed <- bootstrap_oracle(d$s,mu,d["x"],trials=d$trials)
    set.seed(4109)
    simulations <- replicate(39,{
        y <- rbinom(nrow(d),d$trials,mu)
        bootstrap_oracle(y,mu,d["x"],trials=d$trials)$statistic
    })
    expect_equal(unname(g$bootstrap$statistics),t(simulations))
    expected <- (1+rowSums(simulations>=observed$statistic))/40
    expect_equal(c(g$KS,g$CvM),expected)
    expect_equal(g$bootstrap$tail.bounds$p.lower,expected)
})

test_that("one fixed grid and its response-independent selector serve every bootstrap draw", {
    d <- data.frame(x=seq(-1,1,length.out=60),y=rep(c(0,1,0),20))
    m <- glm(y~x,d,family=poisson())
    calls <- 0L
    selector <- function(t,n.qt,n.eq) {
        calls <<- calls+1L
        c(1L,length(t))
    }
    g <- cumres_bootstrap(m,"x",R=9,plots=1,subset.max=2,subset=selector)
    expect_equal(calls,1L)
    expect_equal(g$x,list(range(d$x)))
    expect_equal(dim(g$What[[1]]),c(2L,1L))
    expect_equal(g$nobs,60L)
    expect_error(cumres_bootstrap(m,"x",R=9,subset.max=2,
                                  subset=function(...) c(0,2)),"valid 1-based")
})
