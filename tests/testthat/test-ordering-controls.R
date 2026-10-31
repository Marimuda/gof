ordering_flat_fit <- function() {
    set.seed(7281)
    y <- unlist(lapply(c(200,187,200),function(s)
        drop(rmultinom(1,s,rep(.01,100)))))
    glm(y~x,data=data.frame(y=y,x=rep(c(-1,0,1),each=100)),
        family=poisson(),control=glm.control(epsilon=1e-12))
}

test_that("uncertain fitted ordering withholds inference without modifying coordinates", {
    m <- ordering_flat_fit()
    expect_warning(g <- cumres(m,variable=c("predicted","x"),R=19,plots=2),"ordering.*predicted")
    expect_identical(g$ordering$status,c("unresolved","fixed"))
    expect_true(all(is.na(c(g$KS[1],g$CvM[1]))))
    expect_true(all(is.finite(c(g$KS[2],g$CvM[2]))))
    expect_equal(g$x[[1]],sort(unique(exp(drop(model.matrix(m)%*%coef(m))))))
    mu <- exp(drop(model.matrix(m)%*%coef(m)))
    residual <- model.response(model.frame(m))-mu
    expected <- max(abs(vapply(sort(unique(mu)),function(t) sum(residual[mu<=t]),numeric(1))))/sqrt(length(mu))
    expect_equal(max(abs(g$W[[1]])),expected,tolerance=1e-10)
    expect_identical(g$ordering$ties,rep("group",2))
    expect_equal(g$n,300)
    # A platform may return an exactly flat numerical fit. A single fitted
    # coordinate then spans only the intercept score; split coordinates do
    # not. Ordering precision and score-space information remain separate.
    expect_identical(g$informative,c(length(g$x[[1]])>1L,TRUE))
})

test_that("strict ordering policy fails before simulation", {
    m <- ordering_flat_fit();set.seed(912);state <- .Random.seed
    expect_error(cumres(m,variable="predicted",R=19,on_unresolved="error"),"ordering.*predicted")
    expect_identical(.Random.seed,state)
    expect_error(cumres(m,variable="x",on_unresolved="ignore"),"on_unresolved")
})

test_that("ordering availability survives all presentation interfaces", {
    m <- ordering_flat_fit()
    g <- suppressWarnings(cumres(m,variable=c("predicted","x"),R=19,plots=0))
    tab <- as.data.frame(g)
    expect_identical(tab$ordering_status,c("unresolved","unresolved","fixed","fixed"))
    expect_true(all(!is.na(tab$ordering_reason[1:2])))
    expect_identical(as.data.frame(summary(g)),tab)
    expect_identical(summary(g)$ordering,g$ordering)
    expect_true(all(is.na(coef(g)[1,])))
    expect_output(print(g),"Fitted ordering unresolved")
    expect_output(print(summary(g)),"Fitted ordering unresolved")
    if (requireNamespace("generics",quietly=TRUE) && requireNamespace("tibble",quietly=TRUE))
        expect_identical(as.data.frame(generics::tidy(g)),tab)
    path <- tempfile(fileext=".pdf");pdf(path);on.exit({dev.off();unlink(path)},add=TRUE)
    expect_no_warning(expect_identical(withVisible(plot(g,idx="predicted"))$value,g))
})

test_that("joint families retain an unresolved member and withhold family inference", {
    m <- ordering_flat_fit()
    expect_warning(g <- cumres_joint(m,variable=c("predicted","x"),R=19,plots=0),"ordering.*predicted")
    expect_identical(g$variable,c("predicted","x"))
    expect_true(all(is.na(g$joint$p.value)))
    expect_true(all(is.na(g$joint$adjusted.p.value)))
    expect_true(is.na(g$joint$critical))
    expect_true(all(vapply(g$joint$envelope,function(z)all(is.na(z$lower)),logical(1))))
    tab <- as.data.frame(g)
    expect_true(all(tab$ordering_status[tab$scope=="family"]=="unresolved"))
    expect_true(all(is.na(tab$p.value[tab$scope=="family"])))
    expect_true(all(is.finite(tab$p.value[tab$variable=="x"])))
    expect_identical(as.data.frame(summary(g)),tab)
})

test_that("declared measurement resolution matches an explicitly transformed fixed ordering", {
    set.seed(450);d <- data.frame(x=seq(-1,1,length.out=70));d$y <- rpois(70,exp(.4+.5*d$x))
    d$rounded <- round(d$x/.25)*.25
    m <- glm(y~x,d,family=poisson())
    for(b in c(0,.5)) {
        set.seed(451);a <- cumres(m,variable="x",resolution=c(x=.25),b=b,R=19,plots=2,subset.max=Inf)
        set.seed(451);z <- cumres(m,variable="rounded",data=d,b=b,R=19,plots=2,subset.max=Inf)
        for(field in c("x","W","What","KS","CvM","n"))expect_equal(a[[field]],z[[field]])
        expect_equal(a$ordering$resolution,.25)
        expect_equal(a$ordering$distinct_values,length(unique(d$rounded)))
        expect_equal(a$ordering$changed_units,sum(d$x!=d$rounded))
        expect_equal(a$ordering$source,"fixed")
    }
    # No default tolerance coarsens observed coordinates.
    literal <- cumres(m,variable="x",R=5,plots=0,subset.max=Inf)
    expect_equal(literal$x[[1]],sort(unique(d$x)))
})

test_that("bootstrap uses the declared fixed coordinate in every draw", {
    set.seed(452);d <- data.frame(x=seq(-1,1,length.out=45));d$y <- rpois(45,exp(.5+.3*d$x))
    d$rounded <- round(d$x/.5)*.5;m <- glm(y~x,d,family=poisson())
    for(null in c("fitted","specified")) {
        args <- list(model=m,R=19,plots=2,null=null)
        if(null=="specified")args$mean <- exp(.5+.3*d$x)
        set.seed(453);a <- do.call(cumres_bootstrap,c(args,list(variable="x",resolution=c(x=.5))))
        set.seed(453);z <- do.call(cumres_bootstrap,c(args,list(variable="rounded",data=d)))
        for(field in c("x","W","What","KS","CvM","n"))expect_equal(a[[field]],z[[field]])
        expect_equal(unname(a$bootstrap$statistics),unname(z$bootstrap$statistics))
        expect_equal(a$ordering$resolution,.5)
        expect_equal(as.data.frame(a)$ordering_status,rep("fixed",2))
    }
})

test_that("resolution controls reject ambiguous or unsupported declarations", {
    m <- fx_logistic()
    for(r in list(.1,c(x=0),c(x=-1),c(x=NA_real_),c(x=Inf),list(x=.1),c(x=.1,x=.2),
                  c(unknown=.1),c(predicted=.1)))
        expect_error(cumres(m,variable="x",resolution=r,R=3,plots=0),"resolution")
    expect_error(cumres(m,variable="x",resolution=c(x=1e-300),R=3,plots=0),"resolution")
    expect_error(cumres(NULL,variable=1:3,r=1:3,dr=matrix(0,3,0),ic=matrix(0,3,0),
                       resolution=c(x=1),R=3,plots=0),"resolution")
})

test_that("structural prediction ties and saturated cells need no arbitrary tie breaking", {
    d <- data.frame(x=rep(0:1,each=15),y=rep(c(0,1,2),10))
    m <- glm(y~x,d,family=poisson())
    expect_no_warning(g <- cumres(m,variable="predicted",R=5,plots=0))
    expect_true(g$ordering$status %in% c("structural","saturated"))
    expect_identical(g$KS,1)
    expect_true(g$ordering$tie_blocks>=1)
    fixed <- glm(y~0+offset(o),data=transform(d,o=0),family=poisson())
    expect_no_warning(z <- cumres(fixed,variable="predicted",R=5,plots=0))
    expect_identical(z$ordering$status,"specified")
})

test_that("well separated fits keep seeded calculations and consume no extra random draws", {
    m <- fx_logistic()
    set.seed(455);a <- cumres(m,variable=c("predicted","x"),R=19,plots=2);state <- .Random.seed
    set.seed(455);z <- do.call(ref_cumres,c(glm_inputs(m),list(variable=cbind(predicted=plogis(drop(model.matrix(m)%*%coef(m))),x=model.matrix(m)[,"x"]),R=19,plots=2)))
    expect_equal(a$KS,z$KS);expect_equal(a$CvM,z$CvM)
    expect_equal(a$W,z$W);expect_equal(a$What,z$What)
    expect_identical(.Random.seed,state)
    expect_identical(a$ordering$status,c("checked","fixed"))
})

test_that("numerically collapsed or reversed predictions are not repaired by merging", {
    for(s in c(-1,1)) {
        d <- data.frame(x=c(-1,0,1),y=c(200,187,200)*1e12+c(0,0,s))
        m <- glm(y~x,d,family=poisson(),control=glm.control(epsilon=1e-12))
        expect_warning(g <- cumres(m,variable="predicted",R=5,plots=0),"ordering.*predicted")
        expect_identical(g$ordering$status,"unresolved")
        expect_true(all(is.na(c(g$KS,g$CvM))))
        expect_identical(g$ordering$resolution,NA_real_)
    }
})

test_that("coarsening retains independent binomial trials and grouped equivalence", {
    d <- fx_grouped_data();control <- glm.control(epsilon=1e-14,maxit=100)
    grouped <- glm(cbind(s,f)~x,d,family=binomial(),control=control)
    expanded <- glm(y~x,expand_grouped(d),family=binomial(),control=control)
    set.seed(456);a <- cumres(grouped,variable="x",resolution=c(x=2),R=19,plots=2)
    set.seed(456);b <- cumres(expanded,variable="x",resolution=c(x=2),R=19,plots=2)
    expect_same_cumres(a,b)
    expect_equal(a$ordering,b$ordering)
    expect_equal(a$n,sum(d$size))
    expect_equal(a$nobs,nrow(d))
})

test_that("measurement metadata follows positive-weight fitted rows", {
    set.seed(457);d <- data.frame(x=seq(-1,1,length.out=50),w=rep(c(0,1,2,1,1),10))
    d$y <- 1+.5*d$x+rnorm(50);d$x[7] <- NA
    m <- lm(y~x,d,weights=w,na.action=na.exclude)
    mf <- model.frame(m);use <- model.weights(mf)>0;x <- model.matrix(m)[use,"x"]
    g <- cumres(m,variable="x",resolution=c(x=.5),R=9,plots=0)
    expect_equal(g$nobs,sum(use))
    expect_equal(g$x[[1]],sort(unique(round(x/.5)*.5)))
    expect_equal(g$ordering$changed_units,sum(x!=round(x/.5)*.5))
    expect_equal(g$ordering$tied_units,sum(use))
    expect_equal(g$ordering$tie_blocks,5)
})

test_that("unsupported adapters reject ordering controls instead of ignoring them", {
    ordinal <- ordinal_fixture()
    expect_error(cumres(ordinal,resolution=c(x=.1),R=3,plots=0),"not implemented")
    expect_error(cumres(ordinal,on_unresolved="error",R=3,plots=0),"not implemented")
    set.seed(458);d <- data.frame(x=rnorm(50));d$y <- 1+d$x+rnorm(50)
    sem <- lava::estimate(lava::lvm(y~x),d)
    expect_error(cumres(sem,y~x,resolution=c(x=.1),R=3,plots=0),"not implemented")
    expect_error(cumres(sem,y~x,on_unresolved="error",R=3,plots=0),"not implemented")
})

test_that("moving-window checks include entry/exit comparisons", {
    d <- data.frame(x=rep(c(-1,0,1),each=20))
    d$y <- 2+.5*d$x+rep(c(.05,-.1,.05),each=20)+rep(c(-1,1),30)
    m <- lm(y~x,d)
    ordinary <- cumres(m,variable="predicted",R=9,plots=0)
    expect_identical(ordinary$ordering$status,"checked")
    mu <- drop(model.matrix(m)%*%coef(m));bandwidth <- mu[21]-mu[1]
    expect_warning(window <- cumres(m,variable="predicted",b=bandwidth,R=9,plots=0),"ordering.*predicted")
    expect_identical(window$ordering$status,"unresolved")
    expect_match(window$ordering$reason,"window")
    expect_gt(window$ordering_check$unresolved_window_comparisons,0)
    expect_true(all(is.na(c(window$KS,window$CvM))))
    expect_no_warning(fixed <- cumres(m,variable="x",b=bandwidth,R=9,plots=0))
    expect_true(all(is.finite(c(fixed$KS,fixed$CvM))))
})

test_that("misspelled controls cannot silently change the intended diagnostic", {
    m <- fx_logistic()
    expect_error(cumres(m,variable="x",resolutions=c(x=.1),R=3,plots=0),"Unused.*resolutions")
    expect_error(cumres(m,variable="x",on_unresolve="error",R=3,plots=0),"Unused.*on_unresolve")
    expect_error(cumres(m,variable="x",ties.method="random",R=3,plots=0),"Unused.*ties.method")
    expect_error(cumres(fx_lm(),variable="x1",resolutions=c(x1=.1),R=3,plots=0),"Unused.*resolutions")
    set.seed(459);a <- cumres(m,variable="x",subset.max=5,R=9,plots=0)
    set.seed(459);b <- cumres(m,variable="x",subset.ma=5,R=9,plots=0)
    expect_same_cumres(a,b)
})

test_that("bootstrap tie metadata counts trials rather than stored binomial groups", {
    d <- fx_grouped_data();m <- glm(cbind(s,f)~x,d,family=binomial())
    for(resolution in list(NULL,c(x=2))) {
        a <- cumres(m,variable="x",resolution=resolution,R=5,plots=0)
        b <- cumres_bootstrap(m,variable="x",resolution=resolution,R=5,plots=0,
                              null="specified",mean=rep(.5,nrow(d)))
        expect_equal(a$ordering,b$ordering)
        expect_equal(b$ordering$tied_units,sum(d$size))
        expect_equal(b$n,sum(d$size))
        expect_equal(b$nobs,nrow(d))
    }
})
