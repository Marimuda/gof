test_that("all five declared-link likelihoods match independent scalar integrals", {
    for (link in c("logit","probit","cloglog","loglog","cauchit")) {
        problem <- clmm_problem_fixture("intercept",m=4L,nc=4L,seed=21L,link=link,sd=.9)
        F <- switch(link,logit=stats::plogis,probit=stats::pnorm,
            cloglog=function(x) -expm1(-exp(x)),loglog=function(x) exp(-exp(-x)),cauchit=stats::pcauchy)
        alpha <- c(-Inf,-1.5,0,1.5,Inf)
        eta <- .8*problem$X[,"x"]
        expected <- vapply(seq_len(4L),function(block) {
            rows <- which(problem$blocks$id==block)
            y <- problem$y[rows]
            integrand <- function(b) vapply(b,function(value) {
                probability <- F(alpha[y+1L]-eta[rows]-value)-F(alpha[y]-eta[rows]-value)
                prod(probability)*stats::dnorm(value,sd=.9)
            },0)
            log(stats::integrate(integrand,-Inf,Inf,rel.tol=1e-11,abs.tol=1e-300)$value)
        },0)
        integrator <- gof_integrator_nested_agq(problem,81L)
        actual <- integrator$loglik(problem$theta,integrator$refresh(problem$theta))
        expect_equal(unname(actual),expected,tolerance=1e-8,info=link)
    }
})
