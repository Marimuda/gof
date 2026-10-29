test_that("tied processes contain all residuals and are invariant to units", {
    d <- fx_grouped_data(sizes=rep(40,10),seed=63)
    d$large <- d$x*1e6
    m <- glm(cbind(s,f)~x,d,family=binomial())
    ml <- glm(cbind(s,f)~large,d,family=binomial())
    set.seed(100); a <- cumres(m,variable="x",R=100,plots=3)
    set.seed(100); b <- cumres(ml,variable="large",R=100,plots=3)
    expected <- cumsum(d$s-(d$s+d$f)*fitted(m))/sqrt(sum(d$s+d$f))
    expect_equal(a$W[[1]],unname(expected),tolerance=1e-10)
    expect_equal(a$W,b$W,tolerance=1e-10)
    expect_equal(a$KS,b$KS)
    expect_equal(a$CvM,b$CvM)
    expect_equal(a$What,b$What,tolerance=1e-10)
    expect_equal(b$x[[1]]/1e6,a$x[[1]])
})

test_that("na.exclude and na.omit use the identical fitted rows", {
    set.seed(1); d <- data.frame(x=rnorm(100),y=rbinom(100,1,.4)); d$x[1] <- NA
    for (type in c("lm","glm")) {
        fit <- if (type=="glm") function(action) glm(y~x,d,family=binomial(),na.action=action) else
            function(action) lm(y~x,d,na.action=action)
        set.seed(10); a <- cumres(fit(na.exclude),R=20,plots=2)
        set.seed(10); b <- cumres(fit(na.omit),R=20,plots=2)
        expect_same_cumres(a,b)
        expect_equal(a$What,b$What)
        expect_equal(a$n,99L)
    }
})

test_that("grouped expansion retains the actual fit and precise controls", {
    d <- fx_grouped_data(sizes=rep(40,10),seed=12)
    count <- 0L
    method <- function(...) {count <<- count+1L; glm.fit(...)}
    control <- glm.control(epsilon=1e-14,maxit=100)
    m <- glm(cbind(s,f)~x,d,family=binomial(link="cloglog"),control=control,method=method)
    mi <- glm(y~x,expand_grouped(d),family=family(m),control=control)
    set.seed(1); a <- cumres(m,R=100,plots=2)
    set.seed(1); b <- cumres(mi,R=100,plots=2)
    expect_equal(a$KS,b$KS)
    expect_equal(a$CvM,b$CvM)
    expect_equal(a$variable,b$variable)
    ## glm's deviance stopping criterion leaves the independent optimizers
    ## about 5e-9 apart in coefficients. Bound absolute process errors rather
    ## than relative errors in small, nearly cancelling residual sums.
    expect_lt(max(abs(unlist(a$W)-unlist(b$W))),1e-8)
    expect_lt(max(abs(unlist(a$What)-unlist(b$What))),1e-7)
    expect_identical(count,1L) # diagnostics do not re-fit the user's model
})

test_that("moving windows match a direct-indicator oracle on full and reduced grids", {
    m <- fx_logistic(n=80)
    inp <- glm_inputs(m)
    v <- cbind(predicted=fitted(m),x=model.matrix(m)[,"x"])
    for (cap in c(Inf,19L)) {
        set.seed(88); a <- cumres(m,R=30,plots=3,b=.7,subset.max=cap)
        set.seed(88); b <- ref_cumres(inp$r,inp$dr,inp$ic,v,R=30,plots=3,b=.7,subset.max=cap)
        expect_equal(a$W,b$W)
        expect_equal(a$What,b$What)
        expect_equal(a$KS,b$KS)
        expect_equal(a$CvM,b$CvM)
    }
})

test_that("zero plots, odd and singleton caps, and custom selectors work", {
    m <- fx_logistic(n=501)
    g <- cumres(m,R=5,plots=0,subset.max=499)
    expect_true(all(lengths(g$W)<=499))
    expect_true(all(vapply(g$What,ncol,integer(1))==0L))
    pdf(NULL); on.exit(dev.off())
    expect_no_error(plot(g,col=NULL))
    expect_no_error(plot(g,idx="x"))
    expect_error(plot(g,idx="unknown"),"Unknown")
    expect_equal(lengths(cumres(m,R=5,plots=0,subset.max=1)$W),c(1L,1L))
    calls <- 0L
    select <- function(t,n.qt,n.eq) {calls <<- calls+1L; c(1L,length(t))}
    g <- cumres(m,R=5,plots=0,subset.max=5,subset=select)
    expect_identical(calls,2L)
    expect_equal(lengths(g$W),c(2L,2L))
    expect_error(cumres(m,R=5,subset.max=5,subset=function(...) c(0,1)),"subset must return")
})

test_that("invalid arguments and mismatched dimensions stop before native allocation", {
    m <- fx_logistic()
    for (cap in list(0,-1,NA,NaN,500.5)) expect_error(cumres(m,subset.max=cap),"subset.max")
    for (r in list(0,-1,NA,1.5)) expect_error(cumres(m,R=r),"R must be")
    expect_error(cumres(m,R=5,plots=-1),"plots")
    expect_error(cumres(m,R=5,b=-.5),"bandwidth")
    expect_error(cumres(m,variable="unknown",R=5),"Variable not found")
    expect_error(gof:::cumres.default(NULL,variable=1:5,r=1:4,dr=matrix(1,4,1),ic=matrix(1,4,1),coef=NULL),"dimensions")
    expect_warning(cumres(m,R=5,plots=0,breakties=1e-12),"deprecated")
})

test_that("zero-length L2 integration has p-value one", {
    m <- glm(y~1,data=data.frame(y=rep(c(0,1),50)),family=binomial())
    g <- cumres(m,R=10,plots=0)
    expect_equal(g$CvM,1)
})

test_that("normal glm matches lm and unsupported fits have informative errors", {
    m <- fx_lm()
    mg <- glm(formula(m),model.frame(m),family=gaussian())
    set.seed(90); a <- cumres(m,R=20,plots=2)
    set.seed(90); b <- cumres(mg,R=20,plots=2)
    expect_same_cumres(a,b)
    expect_equal(a$What,b$What)
    expect_error(cumres(glm(y~x1+x2,transform(model.frame(m),y=exp(y)),family=Gamma()),R=5),"Unsupported family")
    mw <- lm(formula(m),model.frame(m),weights=rep(1,100))
    expect_error(cumres(mw,R=5),"Weight argument")
})
