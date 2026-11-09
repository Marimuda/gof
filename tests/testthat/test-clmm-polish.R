## Polishing clmm estimates to the root of the AGQ score, acceptance checks
## and the cluster multiplier design.

polish_problem <- function(d,structure,xi,beta,omega,link="logit") {
    n <- nrow(d)
    top <- as.integer(d$g)
    m <- max(top)
    one <- matrix(1,n,1L,dimnames=list(NULL,"(Intercept)"))
    terms <- switch(structure,
        intercept=list(list(label="1 | g",group="g",index=top,n_levels=m,Z=one)),
        slope=list(list(label="1 + z | g",group="g",index=top,n_levels=m,
                        Z=cbind("(Intercept)"=1,z=d$z))))
    theta <- c(stats::setNames(xi,paste0("xi",seq_along(xi))),x=unname(beta),
               stats::setNames(omega,paste0("omega",seq_along(omega))))
    new_gof_clmm_problem(y=as.integer(d$y),K=nlevels(d$y),X=cbind(x=d$x),offset=numeric(n),
                         weights=rep(1,n),tJac=diag(nlevels(d$y)-1L),link=link,terms=terms,
                         blocks=list(id=top,m=m),structure="nesting-chain",chain="g",theta=theta)
}

polish_from_clmm <- function(fit,d) polish_problem(d,"intercept",fit$alpha,fit$beta,log(fit$ST$g[1L]))

test_that("polishing reaches the AGQ root and agrees with clmm(nAGQ = 25)", {
    skip_if_not_installed("ordinal")
    skip_if_not_installed("numDeriv")
    d <- clmm_simulate("intercept",m=40L,nc=6L,seed=3)
    lap <- ordinal::clmm(y~x+(1|g),data=d)
    p <- polish_from_clmm(lap,d)
    fac <- gof_select_integrator(p,gof_clmm_control())
    fit <- gof_polish(p,fac,gof_clmm_control())
    expect_s3_class(fit,"gof_clmm_fit")
    expect_identical(fit$integrator,"nested-agq")
    expect_identical(fit$nodes,11L)
    expect_lt(fit$standardized_score,1e-6)
    integ <- fac()
    start <- sum(integ$loglik(p$theta,integ$refresh(p$theta)))
    expect_gte(fit$loglik,start)
    expect_equal(fit$loglik,sum(integ$loglik(fit$theta,fit$state)))
    ## Exact-likelihood estimates: clmm with 25 adaptive nodes.
    agq <- ordinal::clmm(y~x+(1|g),data=d,nAGQ=25)
    ref <- c(agq$alpha,agq$beta,log(agq$ST$g[1L]))
    expect_lt(max(abs((fit$theta-ref)/fit$shift$se)),0.01)
    expect_named(fit$shift,c("parameter","estimate_clmm","estimate","se","shift_se"))
    expect_equal(fit$shift$estimate_clmm,unname(p$theta))
    expect_lt(fit$node_check$standardized_score,gof_clmm_control()$node.tol)
    expect_identical(fit$node_check$check_nodes,15L)
    expect_named(fit$boundary,c("term","group","lambda_min","se","z","at_boundary"))
    expect_false(any(fit$boundary$at_boundary))
    expect_gt(fit$boundary$z,2)
    expect_gte(fit$refreshes,2L)
    expect_true(all(fit$iterations>=0L))
    ## Influence: per-block score times the inverse information.
    mult <- gof_clmm_multiplier(p,fit)
    expect_s3_class(mult,"gof_multiplier_design")
    expect_identical(mult$n_units,40L)
    expect_equal(mult$divisor,p$N/sqrt(40))
    ic <- mult$factors[[1L]]$influence
    expect_identical(dim(ic),c(40L,length(p$theta)))
    expect_lt(max(abs(colSums(ic))),1e-6*max(abs(ic)))
    expect_identical(mult$factors[[1L]]$name,"g")
    expect_identical(mult$factors[[1L]]$index,p$blocks$id)
})

test_that("polishing works for slopes and nesting chains without ordinal", {
    skip_if_not_installed("numDeriv")
    for (s in c("slope","nested")) {
        p <- clmm_problem_fixture(s,m=100L,nc=8L,seed=1)
        fac <- gof_select_integrator(p)
        fit <- suppressWarnings(gof_polish(p,fac,gof_clmm_control(shift.tol=Inf)))
        expect_lt(fit$standardized_score,1e-6)
        integ <- fac()
        expect_gte(fit$loglik,sum(integ$loglik(p$theta,integ$refresh(p$theta))))
        ## Root of the frozen-node score at the final state.
        S <- colSums(integ$score(fit$theta,fit$state)[[1L]])
        expect_lt(sqrt(drop(crossprod(S,solve(fit$information,S)))),1e-6)
        expect_equal(fit$information,t(fit$information))
        expect_identical(nrow(fit$boundary),length(p$terms))
        mult <- gof_clmm_multiplier(p,fit)
        expect_lt(max(abs(colSums(mult$factors[[1L]]$influence))),1e-6*max(abs(mult$factors[[1L]]$influence)))
    }
})

test_that("covariances at the boundary are rejected", {
    skip_if_not_installed("ordinal")
    skip_if_not_installed("numDeriv")
    ## No between-cluster variation: the variance runs to zero.
    d <- clmm_simulate("intercept",m=40L,nc=6L,seed=4,sd=0)
    p <- polish_from_clmm(ordinal::clmm(y~x+(1|g),data=d),d)
    expect_error(gof_polish(p,gof_select_integrator(p)),"at or near its boundary")
    ## Wine random slope: the covariance is barely identified by 9 judges.
    wine <- NULL
    utils::data("wine",package="ordinal",envir=environment())
    fw <- ordinal::clmm(rating~temp+contact+(1+temp|judge),data=wine)
    n <- nrow(wine)
    X <- stats::model.matrix(~temp+contact,wine)[,-1L,drop=FALSE]
    j <- as.integer(wine$judge)
    om <- gof_logchol_omega(tcrossprod(fw$ST$judge))
    names(om) <- paste0("judge.",names(om))
    pw <- new_gof_clmm_problem(y=as.integer(wine$rating),K=5L,X=X,offset=numeric(n),weights=rep(1,n),
        tJac=diag(4L),link="logit",
        terms=list(list(label="1 + temp | judge",group="judge",index=j,n_levels=9L,
                        Z=cbind("(Intercept)"=1,tempwarm=X[,1L]))),
        blocks=list(id=j,m=9L),structure="nesting-chain",chain="judge",
        theta=c(fw$alpha,fw$beta,om))
    expect_error(gof_polish(pw,gof_select_integrator(pw)),"at or near its boundary")
    ## Structure code flags singular supplied covariances up front.
    pb <- p
    pb$source$boundary <- "1 | g"
    expect_error(gof_polish(pb,function(nodes) stop("not reached")),"at or near its boundary.*1 \\| g")
})

test_that("shift, node and convergence checks are triggerable", {
    skip_if_not_installed("numDeriv")
    p <- clmm_problem_fixture("intercept",m=50L,nc=6L,seed=2)
    fac <- gof_select_integrator(p)
    expect_error(gof_polish(p,fac,gof_clmm_control(shift.tol=1e-6)),"shift.tol=Inf")
    expect_error(suppressWarnings(gof_polish(p,fac,gof_clmm_control(shift.tol=Inf,node.tol=1e-12))),
                 "increase nodes via gof_clmm_control\\(nodes=\\)")
    expect_error(gof_polish(p,fac,gof_clmm_control(shift.tol=Inf,maxit=1L)),"did not converge")
    ## A start far from the root warns about the shift but is accepted by
    ## default (shift.tol = Inf); small shifts are reported silently.
    far <- p
    far$theta[["x"]] <- far$theta[["x"]]+.3
    expect_warning(fit <- gof_polish(far,gof_select_integrator(far)),"fit_check\\$shift")
    expect_gt(max(abs(fit$shift$shift_se)),0.5)
    expect_error(suppressWarnings(gof_polish(far,gof_select_integrator(far),gof_clmm_control(shift.tol=.5))),
                 "shift.tol")
    ## Three nodes leave a quadrature bias of about 0.05 to 0.2 SE that the
    ## Q + 4 node check detects.
    q3 <- function(nodes) gof_integrator_nested_agq(p,if (is.null(nodes)) 3L else nodes)
    expect_error(suppressWarnings(gof_polish(p,q3,gof_clmm_control(shift.tol=Inf,node.tol=.01,max.refresh=20L))),
                 "Quadrature node check failed")
    ## Frozen Laplace nodes do not define an estimator.
    q1 <- function(nodes) gof_integrator_nested_agq(p,if (is.null(nodes)) 1L else nodes)
    expect_error(gof_polish(p,q1,gof_clmm_control(shift.tol=Inf)),"at least 3 quadrature nodes")
})

test_that("multiplier designs are refused for unsupported structures", {
    skip_if_not_installed("numDeriv")
    p <- clmm_problem_fixture("intercept",m=30L,nc=5L,seed=6)
    fit <- suppressWarnings(gof_polish(p,gof_select_integrator(p),gof_clmm_control(shift.tol=Inf)))
    cr <- p
    cr$structure <- "crossed"
    expect_error(gof_clmm_multiplier(cr,fit),"not supported")
    expect_error(gof_clmm_multiplier(p,unclass(fit)),"gof_polish")
    two <- fit
    two$score <- c(fit$score,list(item=fit$score[[1L]]))
    expect_error(gof_clmm_multiplier(p,two),"one score matrix")
})

test_that("imprecise but interior covariances warn and are accepted", {
    skip_if_not_installed("ordinal")
    skip_if_not_installed("numDeriv")
    d <- clmm_simulate("intercept",m=40L,nc=6L,seed=3)
    p <- polish_from_clmm(ordinal::clmm(y~x+(1|g),data=d),d)
    ## Any finite z below boundary.z only warns; the fit is not at its boundary.
    expect_warning(fit <- gof_polish(p,gof_select_integrator(p),gof_clmm_control(boundary.z=1e3)),
                   "imprecisely estimated")
    expect_s3_class(fit,"gof_clmm_fit")
    expect_false(any(fit$boundary$at_boundary))
    ## A covariance numerically at zero is a real boundary under any boundary.z.
    expect_true(gof_clmm_boundary(p,replace(p$theta,length(p$theta),-10),fit$information,1e-6)$at_boundary)
})
