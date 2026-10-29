sem_options_fixture <- function(n=100L) {
    set.seed(410)
    m <- lava::lvm(list(c(y1,y2,y3)~eta,eta~x)); lava::latent(m) <- ~eta
    d <- lava::sim(m,n)
    list(model=lava::estimate(m,d),data=d)
}

sem_compare <- function(actual,reference) {
    expect_equal(actual$x,reference$x,tolerance=1e-8)
    expect_equal(actual$W,reference$W,tolerance=1e-8)
    expect_equal(actual$What,reference$What,tolerance=1e-8)
    expect_equal(actual$KS,reference$KS)
    expect_equal(actual$CvM,reference$CvM)
}

test_that("full SEM predictor selection uses structural covariance paths", {
    f <- sem_options_fixture()
    expect_setequal(gof:::cumresMaxSet(f$model,"y1"),c("x","y2","y3"))
    m <- lava::lvm(list(y~x,z~y))
    expect_setequal(gof:::cumresMaxSet(m,"y"),"x")
    m <- lava::lvm(list(y~x,z~x))
    expect_setequal(gof:::cumresMaxSet(m,"y"),c("x","z"))
    lava::covariance(m,y~z) <- -.2
    expect_setequal(gof:::cumresMaxSet(m,"y"),"x")
    cyclic <- lava::lvm(list(y~z,z~y,y~x))
    lava::regression(cyclic,y~z) <- .2
    lava::regression(cyclic,z~y) <- -.3
    expect_setequal(gof:::cumresMaxSet(cyclic,"y"),"x")
    ## Generate the stable feedback system analytically; lava::sim itself
    ## only handles acyclic graphs. y=.2*z+x+e_y, z=-.3*y+e_z.
    set.seed(409); d <- data.frame(x=rnorm(120),ey=rnorm(120),ez=rnorm(120))
    d$y <- (d$x+d$ey+.2*d$ez)/1.06
    d$z <- -.3*d$y+d$ez
    fit.cyclic <- lava::estimate(cyclic,d)
    expect_no_error(cumres(fit.cyclic,"y",x="z",full=TRUE,R=5,plots=0))
    expect_error(gof:::cumresMaxSet(m,"unknown"),"named structural residual")
})

test_that("full SEM latent predictions match the Gaussian conditional formula", {
    f <- sem_options_fixture(); e <- f$model; d <- model.frame(e); p <- lava::pars(e)
    mo <- lava::moments(e)
    pp <- mo$P; xx <- lava::index(e)$exo.idx; pp[xx,xx] <- 0
    covariance <- mo$IAi %*% pp %*% t(mo$IAi)
    baseline <- predict(e,x=~1,p=p,data=d)
    used <- c("y2","y3")
    expected <- baseline[,"eta"]+drop((as.matrix(d[used])-baseline[,used]) %*%
        solve(covariance[used,used],covariance[used,"eta",drop=FALSE]))
    rf <- function(q) predict(e,x=lava::vars(e),residual=TRUE,p=q,data=d)[,"y1"]
    dr <- numDeriv::jacobian(rf,p,method=lava::lava.options()$Dmethod)
    set.seed(411); actual <- cumres(e,y1~eta,full=TRUE,R=20,plots=3,subset.max=17)
    set.seed(411); reference <- ref_cumres(rf(p),dr,lava::iid(e),cbind(expected),R=20,plots=3,subset.max=17)
    sem_compare(actual,reference)
    ## y2 is itself allowed in the conditioning set: request it explicitly
    ## rather than taking the default prediction columns (which omit it).
    set.seed(412); actual <- cumres(e,"y1",x="y2",full=TRUE,R=10,plots=2)
    expect_equal(actual$x[[1]],sort(d$y2),tolerance=1e-10)
})

test_that("full SEM handles one endogenous variable and no eligible predictor", {
    set.seed(420); m <- lava::lvm(y~x); e <- lava::estimate(m,lava::sim(m,100))
    expect_setequal(gof:::cumresMaxSet(e,"y"),"x")
    expect_no_error(cumres(e,"y",x="x",full=TRUE,R=5,plots=0))
    expect_no_error(cumres(e,y~1,R=5,plots=0))
    m0 <- lava::lvm(y~1); e0 <- lava::estimate(m0,lava::sim(m0,100))
    expect_length(gof:::cumresMaxSet(e0,"y"),0L)
    expect_no_error(cumres(e0,"y",x="y",full=TRUE,R=5,plots=0))
    ## Ordinary covariance specification reclassifies the correlated parent
    ## as endogenous. Forcing it back into exogeneity contradicts the SEM's
    ## conditional centering: even ~1 then conditions on that unsafe variable.
    set.seed(911); x <- rnorm(300); ey <- .2*x+sqrt(.96)*rnorm(300)
    d <- data.frame(y=x+ey,x=x)
    correlated <- lava::lvm(y~x); lava::covariance(correlated,y~x) <- .2
    regular <- lava::estimate(correlated,d)
    expect_length(lava::exogenous(regular),0L)
    g <- cumres(regular,"y",x="x",full=TRUE,R=5,plots=0)
    expect_length(g$x[[1]],1L)
    lava::exogenous(correlated) <- "x"
    forced <- lava::estimate(correlated,d)
    expect_error(cumres(forced,"y",x="x",full=TRUE,R=5,plots=0),"forced exogeneity")
})

test_that("SEM collections use evaluated objects in local function scopes", {
    f <- sem_options_fixture(); e <- f$model; d <- f$data
    fit_local <- function(ee) {
        xfun <- function(p) d$x
        yfun <- function(p) predict(ee,x=lava::vars(ee),residual=TRUE,p=p)[,"y1"]
        list(
            formulas=cumres(ee,list(y1~eta,y2~x),R=5,plots=0),
            crossed=cumres(ee,c(y1,y2)~c(x,y3),R=5,plots=0),
            functions=cumres(ee,list(yfun,"y2"),x=list(xfun,d$x^2),R=5,plots=0),
            singleton=cumres(ee,list("y1"),x=list("x"),R=5,plots=0))
    }
    result <- fit_local(e)
    expect_length(result$formulas,2L)
    expect_length(result$crossed,4L)
    expect_named(result$crossed,c("y1<-x","y1<-y3","y2<-x","y2<-y3"))
    expect_length(result$functions,4L)
    expect_s3_class(result$singleton,"cumres")
    expect_true(all(vapply(result$functions,inherits,logical(1),"cumres")))
})

test_that("explicit SEM parameter and data options affect every calculation", {
    f <- sem_options_fixture(); e <- f$model; d <- model.frame(e)
    p <- lava::pars(e); p[grep("eta~x",names(p))] <- p[grep("eta~x",names(p))]*.9
    d <- d[rev(seq_len(nrow(d))),,drop=FALSE]
    d$x <- d$x+.1
    rf <- function(q) predict(e,x=lava::vars(e),residual=TRUE,p=q,data=d)[,"y1"]
    dr <- numDeriv::jacobian(rf,p,method=lava::lava.options()$Dmethod)
    score <- lava::score(e,p=p,data=d,indiv=TRUE)
    information <- -numDeriv::jacobian(function(q) colSums(lava::score(e,p=q,data=d,indiv=TRUE)),p,
                                      method=lava::lava.options()$Dmethod)
    influence <- score %*% solve(information)
    variable <- predict(e,x=~1,p=p,data=d)[,"eta",drop=FALSE]
    set.seed(431); actual <- cumres(e,"y1",x="eta",p=rev(p),data=d,R=20,plots=3)
    set.seed(431); reference <- ref_cumres(rf(p),dr,influence,variable,R=20,plots=3)
    sem_compare(actual,reference)
    expect_equal(lava::pars(e),lava::pars(f$model))
    ## A raw numeric predictor is a single process, not a collection.
    expect_s3_class(cumres(e,"y1",x=d$x,p=unname(p),data=d,R=5,plots=0),"cumres")
})

test_that("SEM residual functions can supply an analytic gradient", {
    set.seed(440); m <- lava::lvm(y~x); d <- lava::sim(m,100); e <- lava::estimate(m,d)
    p <- lava::pars(e)
    rf <- function(q) predict(e,x=lava::vars(e),residual=TRUE,p=q)[,"y"]
    grad <- numDeriv::jacobian(rf,p,method=lava::lava.options()$Dmethod)
    with_grad <- function(q) structure(rf(q),grad=grad)
    set.seed(441); a <- cumres(e,with_grad,x=function(q) d$x,R=10,plots=2,b=.5)
    set.seed(441); b <- cumres(e,rf,x=d$x,R=10,plots=2,b=.5)
    sem_compare(a,b)
    with_grad_function <- function(q) structure(rf(q),grad=function(p) grad)
    set.seed(441); c <- cumres(e,with_grad_function,x=d$x,R=10,plots=2,b=.5)
    sem_compare(a,c)
})

test_that("invalid SEM option combinations fail before native processing", {
    f <- sem_options_fixture(); e <- f$model; p <- lava::pars(e); d <- model.frame(e)
    expect_error(cumres(e,"y1",x="unknown",R=2),"Predictor not found")
    expect_error(cumres(e,"y1",x="y1",R=2),"use full=TRUE")
    expect_error(cumres(e,"unknown",x="x",R=2),"endogenous or latent")
    expect_error(cumres(e,"y1",R=2),"nonempty")
    expect_error(cumres(e,~x,R=2),"Formula must name")
    expect_error(cumres(e,character(),x="x",R=2),"one endogenous")
    expect_error(cumres(e,list(),x="x",R=2),"nonempty")
    expect_error(cumres(e,"y1",x=list(),R=2),"nonempty")
    expect_error(cumres(e,"y1",x="x",full=NA,R=2),"full")
    expect_error(cumres(e,"y1",x=rep(1,3),R=2),"per data row")
    expect_error(cumres(e,"y1",x=function(p) matrix(1,nrow(d),2),R=2),"per data row")
    expect_error(cumres(e,function(p) rep(0,nrow(d)),x="eta",full=TRUE,R=2),"named structural residual")
    expect_error(cumres(e,"y1",x="x",p=p[-1],R=2),"one finite parameter")
    expect_error(cumres(e,"y1",x="x",p=matrix(p,ncol=1),R=2),"one finite parameter")
    negative <- p; negative[grep("^y1~~y1$",names(p))] <- -1
    expect_error(cumres(e,"y1",x="x",p=negative,R=2),"error covariance")
    impossible <- e
    specification <- lava::Model(e)
    lava::covariance(specification,y1~y2) <- 100
    impossible$model <- impossible$model0 <- specification
    expect_error(cumres(impossible,"y1",x="x",R=2),"positive semidefinite")
    pn <- p; names(pn)[1] <- "unknown"
    expect_error(cumres(e,"y1",x="x",p=pn,R=2),"p names")
    expect_error(cumres(e,"y1",x="x",data=d[-1],R=2),"all manifest")
    d$y1[1] <- NA
    expect_error(cumres(e,"y1",x="x",data=d,R=2),"complete finite")
    expect_error(cumres(e,function(p) rep(1i,nrow(d)),x="x",R=2),"Each residual")
    bad <- e; bad$opt$convergence <- 1L
    expect_error(cumres(bad,"y1",x="x",R=2),"not converged")
    bad <- e; bad$cluster <- rep(1:2,length.out=nrow(d))
    expect_error(cumres(bad,"y1",x="x",R=2),"independent subjects")
    bad <- e; bad$weights <- rep(1,nrow(d))
    expect_error(cumres(bad,"y1",x="x",R=2),"Weighted structural")
    bad <- e; bad$estimator <- "glm"
    expect_error(cumres(bad,"y1",x="x",R=2),"Gaussian fitted model")
})

test_that("a fixed-parameter SEM needs no artificial optimization or projection", {
    ## lava 1.10 estimate() itself errors for zero free parameters. A fixed
    ## Model already defines the likelihood; this minimal lvmfit is a valid
    ## input-contract fixture with known coefficients and no optimization.
    m <- lava::lvm(y~x)
    lava::regression(m,y~x) <- .5
    lava::intercept(m,~y) <- 0
    lava::covariance(m,y~y) <- 1
    set.seed(460); d <- lava::sim(m,80)
    fixed <- structure(list(model=m,model0=m,data=list(model.frame=d,n=nrow(d)),
                            opt=list(estimate=numeric(),par=numeric(),convergence=0L)),class="lvmfit")
    set.seed(461); g <- cumres(fixed,"y",x="x",R=10,plots=2)
    r <- d$y-.5*d$x; o <- order(d$x)
    expect_equal(g$W[[1]],unname(cumsum(r[o])/sqrt(nrow(d))))
    set.seed(461); invisible(runif(nrow(d)))
    for (i in 1:10) invisible(rnorm(nrow(d)))
    for (i in 1:2) {
        w <- cumsum(r[o]*rnorm(nrow(d)))/sqrt(nrow(d))
        expect_equal(g$What[[1]][,i],unname(w))
    }
    singular <- lava::lvm(list(y~z,z~y))
    lava::regression(singular,y~z) <- 1
    lava::regression(singular,z~y) <- 1
    lava::intercept(singular,~y+z) <- 0
    lava::covariance(singular,y~y) <- 1
    lava::covariance(singular,z~z) <- 1
    invalid <- fixed; invalid$model <- invalid$model0 <- singular
    invalid$data$model.frame <- data.frame(y=d$y,z=d$x)
    expect_error(cumres(invalid,"y",x="z",full=TRUE,R=5),"invertible")
})
