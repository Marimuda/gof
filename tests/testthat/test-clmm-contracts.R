test_that("multiplier designs validate their factors and rules", {
    ic <- matrix(1:6/10,3L)
    d <- gof_multiplier_design(list(list(name="g",index=c(1,1,2,3),influence=ic)),n_units=3,divisor=2)
    expect_s3_class(d,"gof_multiplier_design")
    expect_identical(d$n,4L)
    expect_identical(d$d,2L)
    expect_equal(gof_multiplier_summary(d)$levels,3)
    expect_error(gof_multiplier_design(list(),3,1),"at least one factor")
    expect_error(gof_multiplier_design(list(list(name="g",index=c(1,1,3),influence=ic)),3,1),"used level")
    expect_error(gof_multiplier_design(list(list(name="g",index=c(1,2,3),influence=ic)),3,0),"divisor")
    expect_error(gof_multiplier_design(list(list(name="g",index=c(1,2,3),influence=ic)),3,1,draw="pigeonhole"),
                 "Unknown multiplier draw")
    two <- list(list(name="a",index=1:3,influence=ic),list(name="a",index=c(1,1,1),influence=ic[1,,drop=FALSE]))
    expect_error(gof_multiplier_design(two,3,1),"distinct")
})

test_that("the sum rule reproduces a single legacy row multiplier exactly", {
    g <- list(c(.3,-1,2))
    expect_identical(gof_multiplier_rules$combine$sum(g,list(c(3L,1L,2L)),list(NULL)),c(2,.3,-1))
    two <- gof_multiplier_rules$combine$sum(list(c(1,2),c(10,20,30)),list(c(1L,2L,2L),c(3L,1L,2L)),
                                            list(NULL,c(1,0,1)))
    expect_identical(two,c(31,2,22))
    set.seed(4); a <- gof_multiplier_rules$draw$`independent-normal`(c(3L,2L))
    set.seed(4); expect_identical(a,list(stats::rnorm(3),stats::rnorm(2)))
})

test_that("link table matches ordinal conventions and derivatives", {
    x <- c(-3,-.4,0,1.2,4)
    for (nm in gof_link_names) {
        lk <- gof_ordinal_link(nm)
        expect_equal(lk$f(x),vapply(x,function(t) (lk$F(t+1e-6)-lk$F(t-1e-6))/2e-6,0),tolerance=1e-6)
        expect_equal(lk$df(x),vapply(x,function(t) (lk$f(t+1e-6)-lk$f(t-1e-6))/2e-6,0),tolerance=1e-5)
        expect_equal(lk$logF(x),log(lk$F(x)),tolerance=1e-12)
        naive <- log1p(-lk$F(x))
        ok <- is.finite(naive) & naive > -30
        expect_equal(lk$logS(x)[ok],naive[ok],tolerance=1e-10)
        expect_true(all(is.finite(lk$logS(x)) & is.finite(lk$logF(x))))
    }
    expect_equal(gof_ordinal_link("cloglog")$F(0),1-exp(-1))
    expect_equal(gof_ordinal_link("loglog")$F(0),exp(-1))
    expect_error(gof_ordinal_link("identity"),"Unsupported")
    extreme <- c(-Inf,-1e4,-800,-710,710,800,1e4,Inf)
    for (nm in gof_link_names) {
        lk <- gof_ordinal_link(nm)
        expect_false(anyNA(c(lk$F(extreme),lk$f(extreme),lk$df(extreme))),label=nm)
        expect_equal(lk$f(c(-Inf,Inf)),c(0,0),label=nm)
        expect_equal(lk$df(c(-Inf,Inf)),c(0,0),label=nm)
    }
    expect_gt(gof_ordinal_link("cauchit")$f(800),0)
})

test_that("log-Cholesky maps round-trip and differentiate correctly", {
    S <- matrix(c(1.2,.3,.1,.3,.8,-.2,.1,-.2,.5),3L)
    om <- gof_logchol_omega(S)
    expect_named(om,gof_logchol_names(3L))
    expect_equal(gof_logchol_sigma(om,3L)$Sigma,S)
    ds <- gof_logchol_sigma(om,3L)$dSigma
    for (k in seq_along(om)) {
        h <- replace(numeric(length(om)),k,1e-6)
        num <- (gof_logchol_sigma(om+h,3L)$Sigma-gof_logchol_sigma(om-h,3L)$Sigma)/2e-6
        expect_equal(ds[[k]],num,tolerance=1e-7)
    }
})

test_that("problem fixtures satisfy the contract", {
    for (s in c("intercept","slope","nested")) {
        p <- clmm_problem_fixture(s,m=10L,nc=4L)
        expect_s3_class(p,"gof_clmm_problem")
        u <- gof_clmm_unpack(p)
        expect_length(u$alpha,p$K-1L)
        expect_true(all(diff(u$alpha)>0))
        expect_identical(sum(vapply(p$terms,function(t) length(t$omega),0L))+
                         length(p$index$xi)+length(p$index$beta),length(p$theta))
    }
    p <- clmm_problem_fixture("intercept",m=10L,nc=4L)
    expect_error(new_gof_clmm_problem(y=p$y,K=p$K,X=p$X,offset=p$offset,weights=replace(p$weights,1,0),
        tJac=p$tJac,link=p$link,terms=p$terms,blocks=p$blocks,structure=p$structure,theta=p$theta),
        "positive")
    expect_error(new_gof_clmm_problem(y=p$y,K=p$K,X=p$X,offset=p$offset,weights=p$weights,
        tJac=p$tJac,link=p$link,terms=p$terms,blocks=p$blocks,structure="tree",theta=p$theta),
        "structure")
})

test_that("integrator and control constructors validate inputs", {
    f <- function(...) NULL
    expect_s3_class(gof_integrator("x",f,f,f,f),"gof_integrator")
    expect_error(gof_integrator("x",f,1,f,f),"loglik")
    ctl <- gof_clmm_control()
    expect_s3_class(ctl,"gof_clmm_control")
    expect_identical(gof_clmm_default_nodes(2L,ctl),7L)
    expect_identical(gof_clmm_default_nodes(2L,gof_clmm_control(nodes=15)),15L)
    expect_error(gof_clmm_control(shift.tol=0),"shift.tol")
    expect_silent(gof_clmm_control(shift.tol=Inf))
    expect_identical(ctl$shift.tol,Inf)
    expect_identical(ctl$shift.warn,0.5)
    expect_identical(ctl$boundary.tol,1e-6)
    expect_error(gof_clmm_control(boundary.tol=0),"boundary.tol")
    expect_error(gof_clmm_control(shift.warn=-1),"shift.warn")
})
