test_that("unidentified information does not certify a covariance boundary", {
    p <- clmm_problem_fixture("intercept",m=30L,nc=5L,seed=6)
    H <- diag(length(p$theta))
    H[length(p$theta),length(p$theta)] <- -1
    b <- gof_clmm_boundary(p,p$theta,H)
    expect_gt(b$lambda_min,0.1)
    expect_false(b$identified)
    expect_false(b$at_boundary)
    factory <- function(nodes) {
        ig <- gof_integrator("test-information",function(theta) list(),
            function(theta,state) rep(0,p$blocks$m),
            function(theta,state) list(g=matrix(0,p$blocks$m,length(theta))),
            function(theta,state,check=TRUE) H)
        ig$nodes <- 3L
        ig
    }
    error <- tryCatch(gof_polish(p,factory),error=identity)
    expect_s3_class(error,"gof_clmm_information_unidentified")
    expect_identical(error$stage,"polish")
    expect_identical(error$diagnostics$theta,p$theta)
    expect_equal(as.vector(error$diagnostics$info),as.vector(H))
    expect_identical(attr(error$diagnostics$info,"asymmetry"),0)
    expect_false(error$diagnostics$boundary$at_boundary)
    expect_length(error$diagnostics$covariance,1L)
})

test_that("unsuccessful polishing retains state and has a convergence code", {
    skip_if_not_installed("numDeriv")
    p <- clmm_problem_fixture("intercept",m=50L,nc=6L,seed=2)
    error <- tryCatch(gof_polish(p,gof_select_integrator(p),gof_clmm_control(maxit=1L)),
                      error=identity)
    expect_s3_class(error,"gof_clmm_polish_nonconvergence")
    expect_identical(error$stage,"polish")
    expect_length(error$diagnostics$theta,length(p$theta))
    expect_identical(error$diagnostics$theta0,p$theta)
    expect_true(is.data.frame(error$diagnostics$trace))
    expect_true(is.matrix(error$diagnostics$H))
    expect_length(error$diagnostics$last_frozen_score,length(p$theta))
    expect_false(grepl("approximation is not valid at the boundary",conditionMessage(error)))
    copy <- unserialize(serialize(error,NULL))
    expect_identical(copy$diagnostics,error$diagnostics)
})

test_that("singular input fits carry original covariances and an input stage", {
    skip_if_not_installed("ordinal")
    d <- clmm_simulate("intercept",m=40L,nc=6L,seed=4,sd=0)
    f <- ordinal::clmm(y~x+(1|g),data=d)
    f$ST[[1L]][] <- 0
    p <- gof_clmm_problem(f)
    error <- tryCatch(gof_polish(p,gof_select_integrator(p)),error=identity)
    expect_s3_class(error,"gof_clmm_covariance_boundary")
    expect_identical(error$stage,"input_validation")
    expect_equal(error$diagnostics$input_covariance[[1L]],tcrossprod(unname(f$ST[[1L]])))
    expect_identical(error$diagnostics$theta0,p$theta)
})

test_that("ucminf termination codes follow its own convention and retain gradient checks", {
    skip_if_not_installed("ordinal")
    d <- clmm_simulate("intercept",m=40L,nc=6L,seed=3)
    f <- ordinal::clmm(y~x+(1|g),data=d,
                      control=ordinal::clmm.control(method="ucminf"))
    expect_s3_class(f$optRes,"ucminf")
    expect_true(f$optRes$convergence %in% c(1L,2L))
    expect_s3_class(gof_clmm_problem(f),"gof_clmm_problem")
    for(code in 1:2) {
        other <- f; other$optRes$convergence <- code
        expect_s3_class(gof_clmm_problem(other),"gof_clmm_problem")
    }
    for(code in c(0L,3L,4L,-2L)) {
        other <- f;other$optRes$convergence <- code
        expect_error(gof_clmm_problem(other),class="gof_clmm_input_nonconvergence")
    }
    other <- f;other$gradient[] <- Inf
    expect_error(gof_clmm_problem(other),class="gof_clmm_input_gradient")
    # Requested method is not authoritative when the backend used nlminb.
    ordinary <- ordinal::clmm(y~x+(1|g),data=d)
    ordinary$control$method <- "ucminf"
    expect_s3_class(gof_clmm_problem(ordinary),"gof_clmm_problem")
    ordinary$optRes$convergence <- 1L
    expect_error(gof_clmm_problem(ordinary),class="gof_clmm_input_nonconvergence")
})

test_that("an optional interior refit recovers a singular Laplace start against an independent likelihood oracle", {
    skip_if_not_installed("ordinal")
    skip_if_not_installed("numDeriv")
    source(test_path("..","..","tools","clmm-study.R"),local=TRUE)
    cell <- subset(gof_clmm_study_cells("full"),id=="core08")
    d <- gof_clmm_study_data(cell,140800004L)
    input_environment <- list2env(list(d=d),parent=globalenv())
    f <- evalq(ordinal::clmm(y~x+(1+z|g),data=d),input_environment)
    # Emulate its singular Laplace solution without relying on a platform's
    # tiny returned Cholesky diagonal being exactly zero.
    f$ST[[1L]][2L,2L] <- 0
    before <- serialize(f,NULL)
    expect_error(cumres(f,variable="x",R=3),class="gof_clmm_covariance_boundary")
    messages <- character()
    answer <- withCallingHandlers(cumres(f,variable="x",R=3,
        control=gof_clmm_control(refit=TRUE)),warning=function(w) {
            messages <<- c(messages,conditionMessage(w));invokeRestart("muffleWarning")
        })
    expect_true(any(grepl("interior score root.*refitting",messages)))
    expect_identical(serialize(f,NULL),before)
    expect_true(all(is.finite(c(answer$KS,answer$CvM))))
    expect_identical(answer$fit_check$initialization$strategy,"interior ridge")
    expect_equal(answer$fit_check$nodes,11L)
    expect_gt(answer$fit_check$boundary$lambda_min,.05)
    # Separate direct Gaussian-prior integration and constrained optimization:
    # 31,61,101-point refinements agree on this development fixture.
    expect_lt(abs(answer$fit_check$loglik-(-185.49557338848)),1e-5)
    expect_true(anyNA(answer$fit_check$shift$estimate_clmm))
    expect_true(all(is.finite(answer$fit_check$shift$estimate_start)))
    expect_error(cumres(f,variable="x",R=3,control=gof_clmm_control(refit=TRUE,shift.tol=.5)),
                 class="gof_clmm_parameter_shift_unavailable")
    error <- tryCatch(cumres(f,variable="x",R=3,
        control=gof_clmm_control(refit=TRUE,max.refit.pairs=1L)),error=identity)
    expect_s3_class(error,"gof_clmm_refit_failed")
    expect_length(error$diagnostics$attempts,3L)
})
