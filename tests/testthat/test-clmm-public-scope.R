skip_if_not_installed("ordinal")
skip_if_not_installed("numDeriv")

test_that("the public adapter accepts regular random slopes and nesting chains", {
    for (structure in c("slope","nested")) {
        d <- clmm_simulate(structure,m=100L,nc=8L,seed=1L)
        # A fixed formula environment excludes changing test-local bindings
        # from the serialized model used for the immutability assertion.
        formula <- stats::as.formula(if (structure=="slope") "y~x+(1+z|g)" else "y~x+(1|g/h)",
                                     env=.GlobalEnv)
        fit <- ordinal::clmm(formula,data=d)
        original <- serialize(fit,NULL)
        set.seed(1940)
        diagnostic <- cumres(fit,variable="x",R=3,plots=1)
        expect_s3_class(diagnostic,"cumres")
        expect_identical(diagnostic$n,100L)
        expect_equal(diagnostic$design_info$n_rows,800)
        expect_true(all(is.finite(as.data.frame(diagnostic)$p.value)))
        expect_identical(diagnostic$multiplier$factor,"g")
        expect_length(diagnostic$random_effects,if (structure=="slope") 1L else 2L)
        expect_identical(original,serialize(fit,NULL))
    }
})

test_that("the public adapter integrates all four threshold parametrizations", {
    d <- clmm_simulate("intercept",m=40L,nc=6L,seed=3L)
    for (threshold in c("flexible","symmetric","symmetric2","equidistant")) {
        fit <- ordinal::clmm(y~x+(1|g),data=d,threshold=threshold)
        set.seed(1940)
        diagnostic <- cumres(fit,variable="x",R=3,plots=1)
        expect_equal(diagnostic$design_info$fit_parameters,length(fit$alpha)+length(fit$beta)+1L)
        expect_equal(diagnostic$design_info$n_components,3)
        expect_true(all(is.finite(as.data.frame(diagnostic)$p.value)),info=threshold)
        expect_lt(diagnostic$fit_check$node_check$standardized_score,.05)
    }
})

test_that("the public adapter accepts a regular fitted model for each supported link", {
    for (link in c("logit","probit","cloglog","loglog","cauchit")) {
        # These are deterministic software fixtures, not calibration evidence.
        # The cauchit fixture uses a known converged Laplace input; ordinal's
        # upstream optimizer can fail on other datasets before gof is called.
        if (link=="cauchit") {
            d <- clmm_simulate("intercept",m=50L,nc=10L,seed=113000001L,link=link)
            fit <- ordinal::clmm(y~x+(1|g),data=d,link=link)
        } else {
            d <- clmm_simulate("intercept",m=80L,nc=8L,seed=3L,link=link)
            fit <- ordinal::clmm(y~x+(1|g),data=d,link=link,nAGQ=21L)
        }
        set.seed(1940)
        warnings <- character()
        diagnostic <- withCallingHandlers(cumres(fit,variable="x",R=3,plots=1,
            control=gof_clmm_control(nodes=21L,max.refresh=10L)),warning=function(w) {
                warnings <<- c(warnings,conditionMessage(w))
                invokeRestart("muffleWarning")
            })
        # Large shifts are intentionally reported, rather than refused.
        # They are particularly important for ordinal's asymmetric-link
        # scalar backend; the independent likelihood oracle checks gof's
        # declared-link target without relying on that backend's objective.
        expect_true(all(grepl("^Polished estimates moved",warnings)))
        if (length(warnings)) expect_gt(max(abs(diagnostic$fit_check$shift$shift_se)),.5)
        expect_s3_class(diagnostic,"cumres")
        expect_true(all(is.finite(as.data.frame(diagnostic)$p.value)),info=link)
        expect_lt(diagnostic$fit_check$node_check$standardized_score,.05)
    }
})

test_that("the public wine slope example refuses the covariance boundary", {
    wine <- NULL
    data("wine",package="ordinal",envir=environment())
    fit <- ordinal::clmm(rating~temp+(1+temp|judge),data=wine)
    expect_error(suppressWarnings(cumres(fit,R=3,plots=0)),"at or near its boundary")
})

test_that("the public adapter supports its default three-effect dimension limit", {
    set.seed(1942)
    m <- 100L; nc <- 12L; n <- m*nc
    group <- rep(seq_len(m),each=nc)
    x <- rnorm(n); z1 <- rnorm(n); z2 <- rnorm(n)
    effects <- sweep(matrix(rnorm(m*3L),m,3L),2,c(1.2,.8,.6),"*")
    eta <- .6*x+rowSums(cbind(1,z1,z2)*effects[group,])
    y <- 1L+rowSums(runif(n)>plogis(outer(-eta,c(-1.5,0,1.5),"+")))
    d <- data.frame(y=ordered(y,levels=1:4),x=x,z1=z1,z2=z2,g=factor(group))
    fit <- ordinal::clmm(y~x+(1+z1+z2|g),data=d)
    expect_error(cumres(fit,R=3,plots=0,control=gof_clmm_control(max.dim=2L)),
                 "dimension 3 per level, more than max.dim = 2")
    set.seed(1943)
    diagnostic <- cumres(fit,variable="x",R=3,plots=1)
    expect_identical(diagnostic$random_effects[[1]]$q,3L)
    expect_identical(diagnostic$n,m)
    expect_true(all(is.finite(as.data.frame(diagnostic)$p.value)))
    expect_lt(diagnostic$fit_check$node_check$standardized_score,.05)
})
