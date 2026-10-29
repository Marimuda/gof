test_that("precision-weighted normal and Poisson processes agree with independent scores", {
    set.seed(501)
    d <- data.frame(x=runif(90,-1,1),w=runif(90,.5,4))
    d$y <- 1+.4*d$x+rnorm(90)/sqrt(d$w)
    normal <- lm(y~x,d,weights=w)
    d$count <- rpois(90,d$w*exp(.2+.5*d$x)); d$rate <- d$count/d$w
    poisson <- suppressWarnings(glm(rate~x,d,weights=w,family=poisson(),control=glm.control(epsilon=1e-13)))
    for (m in list(normal,poisson)) {
        X <- model.matrix(m); y <- model.response(model.frame(m)); mu <- fitted(m)
        if (inherits(m,"glm")) {
            dm <- family(m)$mu.eta(predict(m,type="link")); variance <- family(m)$variance(mu)
        } else {dm <- rep(1,length(y)); variance <- rep(1,length(y))}
        r <- y-mu; dr <- -X*as.numeric(dm)
        information <- crossprod(X,X*as.numeric(d$w*dm^2/variance))
        ic <- (X*as.numeric(r*d$w*dm/variance)) %*% solve(information)
        variable <- cbind(predicted=mu,x=d$x)
        set.seed(502); g <- cumres(m,R=40,plots=3)
        set.seed(502); ref <- ref_cumres(r,dr,ic,variable,R=40,plots=3)
        expect_equal(g$W,ref$W)
        expect_equal(g$What,ref$What,tolerance=1e-7)
        expect_equal(g$KS,ref$KS)
        expect_equal(g$CvM,ref$CvM)
    }
})

test_that("binomial counts, weighted proportions and weighted groups retain the same fit", {
    d <- fx_grouped_data()
    d$weight <- rep(c(0,1,2,3),2)
    m <- glm(cbind(s,f)~x,d,family=binomial(),control=glm.control(epsilon=1e-14))
    mp <- glm(s/size~x,d,weights=size,family=binomial(),control=glm.control(epsilon=1e-14))
    set.seed(503); a <- cumres(m,R=30,plots=3)
    set.seed(503); b <- cumres(mp,R=30,plots=3)
    expect_same_cumres(a,b); expect_equal(a$What,b$What,tolerance=1e-7)
    size_process <- cumres(mp,variable="size",data=d,R=5,plots=0)
    expect_equal(size_process$variable,"size")
    expect_equal(size_process$x[[1]],sort(unique(d$size)))
    mw <- glm(cbind(s,f)~x,d,weights=weight,family=binomial(),control=glm.control(epsilon=1e-14))
    dd <- transform(d,s=s*weight,f=f*weight)
    me <- glm(cbind(s,f)~x,dd,family=binomial(),control=glm.control(epsilon=1e-14))
    set.seed(504); a <- cumres(mw,R=30,plots=2)
    set.seed(504); b <- cumres(me,R=30,plots=2)
    expect_same_cumres(a,b); expect_equal(a$What,b$What,tolerance=1e-7)
    expect_equal(a$n,sum(dd$s+dd$f)); expect_equal(a$nobs,sum(d$weight>0))
})

test_that("zero precision weights and missing rows are excluded consistently", {
    set.seed(505); d <- data.frame(x=rnorm(100),w=runif(100,.5,4))
    d$y <- 1+d$x+rnorm(100)/sqrt(d$w)
    d$w[c(1,4)] <- 0; d$x[8] <- NA
    a <- lm(y~x,d,weights=w,na.action=na.exclude)
    b <- lm(y~x,subset(d,w>0),weights=w,na.action=na.omit)
    set.seed(506); ga <- cumres(a,R=30,plots=2)
    set.seed(506); gb <- cumres(b,R=30,plots=2)
    expect_same_cumres(ga,gb); expect_equal(ga$What,gb$What)
    expect_equal(ga$n,97L); expect_equal(ga$nobs,97L)
})

test_that("external ordering data follows fitted subsets, NA handling and expansion", {
    set.seed(507)
    d <- data.frame(x=seq(-2,2,length.out=90),z=rnorm(90),w=rep(10,90))
    d$y <- rbinom(90,10,plogis(d$x)); d$x[8] <- NA
    m <- glm(y/w~x,d,weights=w,subset=seq_len(90)%%3!=0,family=binomial())
    set.seed(508); g <- cumres(m,variable="z",data=d,R=20,plots=2)
    expected <- sort(unique(d$z[as.integer(rownames(model.frame(m)))]))
    expect_equal(g$x[[1]],expected)
    rownames(d) <- paste0("id",seq_len(90))
    mn <- lm(z~x,d,subset=seq_len(90)%%3!=0)
    expect_error(cumres(mn,variable="new",data=data.frame(new=1:3),R=5),"fitted rows")
    expect_error(cumres(mn,variable="unknown",data=d,R=5),"Variable not found")
    d$ordered <- ordered(rep(letters[1:3],30))
    expect_equal(length(cumres(mn,variable="ordered",data=d,R=5)$W[[1]]),2L)
    d$unordered <- factor(d$ordered,ordered=FALSE)
    expect_error(cumres(mn,variable="unordered",data=d,R=5),"numeric vector or ordered")
})

test_that("explicit binary/constant selections and zero-parameter models have defined processes", {
    d <- data.frame(x=seq(-1,1,length.out=40),group=rep(0:1,20),y=sin(seq_len(40)))
    m <- lm(y~x+group,d)
    expect_equal(length(cumres(m,variable="group",R=5)$W[[1]]),2L)
    expect_equal(length(cumres(m,variable="(Intercept)",R=5)$W[[1]]),1L)
    fixed <- lm(y~0,d)
    set.seed(509); g <- cumres(fixed,variable="x",data=d,R=10,plots=2)
    expect_equal(g$W[[1]],unname(cumsum(d$y[order(d$x)])/sqrt(nrow(d))))
    expect_equal(g$model,"lm"); expect_equal(g$response,"y")
    fixedp <- glm(y~0+offset(o),data=data.frame(y=1:40,o=rep(0,40)),family=poisson())
    expect_no_error(cumres(fixedp,R=5,plots=0))
})

test_that("bad fitted-model state and noninteger trials are rejected explicitly", {
    m <- fx_logistic()
    m$coefficients[1] <- Inf
    expect_error(cumres(m,R=5),"finite numeric")
    m <- fx_logistic(); m$converged <- FALSE
    expect_error(cumres(m,R=5),"not converged")
    m <- fx_logistic(); m$fitted.values <- m$fitted.values[-1]
    expect_error(cumres(m,R=5),"dimensions")
    m <- fx_logistic()
    expect_error(cumres(m,variable=NA_character_,R=5),"names")
    d <- data.frame(x=seq(-1,1,length.out=40),y=rep(.3,40))
    m <- suppressWarnings(glm(y~x,d,family=binomial()))
    expect_error(cumres(m,R=5),"trial counts")
    m <- suppressWarnings(glm(y~x,d,weights=rep(1.5,40),family=binomial()))
    expect_error(cumres(m,R=5),"integer success/failure")
    m <- fx_lm(); m$model$`(weights)` <- rep(-1,nrow(m$model))
    expect_error(cumres(m,R=5),"nonnegative")
})

test_that("predictions retain true design ties and saturated score processes carry no test information", {
    set.seed(607); d <- data.frame(group=rep(0:1,20),y=rnorm(40),z=rnorm(40))
    m <- lm(y~group+z,d)
    g <- cumres(m,variable=c("group","(Intercept)"),R=30,plots=2)
    expect_false(any(g$informative))
    expect_equal(g$KS,c(1,1)); expect_equal(g$CvM,c(1,1))
    expect_true(all(unlist(g$W)==0)); expect_true(all(unlist(g$What)==0))
    cell_model <- lm(y~group,d)
    g <- cumres(cell_model,R=30,plots=2)
    expect_equal(length(g$x[[1]]),2L)
    expect_false(g$informative)
    expect_equal(g$KS,1)
    m <- glm(I(rep(1,40))~z,d,family=poisson())
    expect_no_error(cumres(m,variable="z",R=5))
    expect_error(cumres(lm(y~factor(seq_len(40)),d),R=5),"residual degrees")
})

test_that("all standard supported-family links agree with direct score calculations", {
    families <- c(lapply(c("logit","probit","cauchit","cloglog","log"),binomial),
                  lapply(c("identity","log","inverse"),gaussian),
                  lapply(c("log","identity","sqrt"),poisson))
    for (i in seq_along(families)) {
        f <- families[[i]]
        set.seed(10000+i)
        d <- data.frame(x=runif(180,.2,.8))
        eta <- if (f$family=="binomial" && f$link=="log") -1+.2*d$x else 1+.2*d$x
        mu <- f$linkinv(eta)
        d$y <- switch(f$family,binomial=rbinom(nrow(d),1,mu),
                      gaussian=rnorm(nrow(d),mu,.05),poisson=rpois(nrow(d),mu))
        m <- glm(y~x,d,family=f,control=glm.control(epsilon=1e-12,maxit=100))
        X <- model.matrix(m); fitted <- f$linkinv(drop(X%*%coef(m)))
        derivative <- f$mu.eta(drop(X%*%coef(m)))
        r <- d$y-fitted; dr <- -X*as.numeric(derivative)
        phi <- if (f$family=="gaussian") summary(m)$dispersion else 1
        ic <- (X*as.numeric(r*derivative/f$variance(fitted))) %*% (vcov(m)/phi)
        set.seed(10200+i); g <- cumres(m,R=15,plots=2)
        set.seed(10200+i); ref <- ref_cumres(r,dr,ic,cbind(fitted,d$x),R=15,plots=2)
        expect_equal(g$W,ref$W,tolerance=1e-7,info=paste(f$family,f$link))
        expect_equal(g$What,ref$What,tolerance=1e-6,info=paste(f$family,f$link))
        expect_equal(g$KS,ref$KS)
        expect_equal(g$CvM,ref$CvM)
    }
})

test_that("score-space ambiguity is basis independent without erasing small positive information", {
    x <- rep(1:40,2); group <- rep(0:1,each=40)
    ordering <- cbind(group=group); weights <- rep(1,length(x))
    direct <- cbind(1,x,group)
    hidden <- cbind(1,x,x+group)
    expect_identical(gof:::gof_informative_processes(direct,weights,ordering),FALSE)
    expect_identical(gof:::gof_informative_processes(hidden,weights,ordering),NA)
    ## This is genuine positive test variance, even though the score space is
    ## close to saturation. Its indicator departure exceeds roundoff bounds.
    near <- cbind(1,x,x+group+1e-9*sin(seq_along(x)))
    expect_identical(gof:::gof_informative_processes(near,weights,ordering),TRUE)
    expect_identical(gof:::gof_informative_processes(hidden,weights,ordering,list(c(0,1)),.5),NA)
    expect_identical(gof:::gof_informative_processes(hidden,weights,cbind(x=x),list(max(x))),NA)
})

test_that("indeterminate score processes preserve curves and withhold p-values", {
    set.seed(12101)
    d <- data.frame(x=rep(1:40,2),group=rep(0:1,each=40),y=rnorm(80))
    d$hidden <- d$x+d$group
    d$u <- (d$x-20.5)/40+d$group; d$v <- (d$x-20.5)/40-d$group
    for (m in list(lm(y~x+hidden,d),lm(y~u+v,d))) {
        for (bandwidth in c(0,.5)) {
            set.seed(12102)
            g <- cumres(m,variable="group",data=d,R=15,plots=2,b=bandwidth)
            expect_identical(g$informative,NA)
            expect_true(all(is.na(c(g$KS,g$CvM))))
            X <- model.matrix(m)
            set.seed(12102)
            raw <- gof:::cumres.default(NULL,variable=cbind(group=d$group),
                        r=as.numeric(m$residuals),dr=-X,ic=lava::iid(m),coef=NULL,
                        R=15,plots=2,b=bandwidth)
            expect_identical(g$W,raw$W)
            expect_identical(g$What,raw$What)
        }
    }
    cap <- cumres(lm(y~x,d),variable="x",R=10,plots=0,subset.max=1)
    expect_identical(cap$informative,NA)
    expect_true(all(is.na(c(cap$KS,cap$CvM))))
    structural <- cumres(lm(y~group,d),variable="group",R=10,plots=0)
    expect_identical(structural$informative,FALSE)
    expect_equal(c(structural$KS,structural$CvM),c(1,1))
    expect_true(all(unlist(structural$W)==0))
})

test_that("indeterminate inference is explicit in plot print and summary", {
    set.seed(12103)
    d <- data.frame(x=rep(1:40,2),group=rep(0:1,each=40),y=rnorm(80))
    d$hidden <- d$x+d$group
    g <- cumres(lm(y~x+hidden,d),variable="group",data=d,R=10,plots=2)
    expect_true(is.na(g$informative))
    expect_match(paste(capture.output(print(g)),collapse=" "),"Numerically indeterminate")
    expect_match(paste(capture.output(summary(g)),collapse=" "),"Numerically indeterminate")
    labels <- NULL
    testthat::local_mocked_bindings(legend=function(x,legend,...) {labels <<- legend; invisible(NULL)},
                                   .package="graphics")
    grDevices::pdf(NULL); on.exit(grDevices::dev.off(),add=TRUE)
    expect_no_error(plot(g,legend="type1"))
    expect_true(any(grepl("Numerically indeterminate",labels,fixed=TRUE)))
})
