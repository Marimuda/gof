## Marginal probabilities P(Y <= k) = E F(alpha_k - eta - tau Z) of a clmm
## and their analytic theta-gradients (R/clmm-marginal.R).

## A fixture problem with optional offsets and equidistant thresholds
## (tJac = cbind(1, 0:(K-2)), two threshold parameters).
marginal_problem <- function(structure,link,offset=TRUE,equidistant=TRUE,m=8L,...) {
    p <- clmm_problem_fixture(structure,link=link,m=m,...)
    off <- if (offset) {set.seed(11); stats::rnorm(p$n,0,0.4)} else p$offset
    tJac <- p$tJac; theta <- p$theta
    if (equidistant) {
        tJac <- cbind(1,0:(p$K-2L))
        theta <- c(xi1=-1.1,xi2=0.9,theta[-p$index$xi])
    }
    new_gof_clmm_problem(y=p$y,K=p$K,X=p$X,offset=off,weights=p$weights,tJac=tJac,link=link,
                         terms=p$terms,blocks=p$blocks,structure=p$structure,chain=p$chain,
                         theta=theta,levels=p$levels)
}

## Independent tau_i^2 as a function of theta.
marginal_tau2 <- function(problem,theta=problem$theta) {
    out <- numeric(problem$n)
    for (tm in problem$terms) {
        S <- gof_logchol_sigma(theta[tm$omega],tm$q)$Sigma
        out <- out+rowSums((tm$Z %*% S)*tm$Z)
    }
    out
}

marginal_c <- function(problem,theta=problem$theta) {
    u <- gof_clmm_unpack(problem,theta)
    outer(-(problem$offset+drop(problem$X %*% u$beta)),u$alpha,`+`)
}

## E G(c - tau Z) by adaptive integrate(), the oracle for every link.
marginal_oracle <- function(G,cc,tau)
    mapply(function(a,s) if (s==0) G(a) else
        stats::integrate(function(z) G(a-s*z)*stats::dnorm(z),-Inf,Inf,
                         rel.tol=1e-12,abs.tol=1e-15,subdivisions=2000L)$value,cc,tau)

test_that("probit marginals match the closed form on every fixture structure", {
    for (s in c("intercept","slope","nested")) {
        p <- marginal_problem(s,"probit",m=30L)
        out <- gof_clmm_marginal(p,nodes=30L)
        tau2 <- marginal_tau2(p)
        expect_equal(out$tau2,tau2,tolerance=1e-13)
        C <- marginal_c(p)
        r <- sqrt(1+tau2)
        expect_equal(out$cdf,stats::pnorm(C/r),tolerance=1e-10)
        ## dP/d(tau^2) = -c/2 (1+tau^2)^(-3/2) phi(c/sqrt(1+tau^2)), chained to omega
        ## through an independent numerical d tau^2 / d omega.
        dP <- -C/2*r^-3*stats::dnorm(C/r)
        omega <- unlist(lapply(p$terms,`[[`,"omega"))
        Dtau <- numDeriv::jacobian(function(t) marginal_tau2(p,t),p$theta)[,omega,drop=FALSE]
        expected <- -do.call(rbind,lapply(seq_len(p$K-1L),function(k) dP[,k]*Dtau))
        expect_equal(unname(out$gradient[,omega,drop=FALSE]),expected,tolerance=1e-9)
    }
})

test_that("the generic quadrature forms reproduce the probit closed form", {
    C <- cbind(c(-30,-4,-1,0,0.3,2.5,8,30))
    for (tau in c(0.05,0.5,1,2,3.5,8)) {
        tt <- rep(tau,nrow(C))
        exact <- gof_marginal_probit(C,tt)
        z <- gof_marginal_zform("probit",C,tt,gof_marginal_trapezoid(40*2^ceiling(log2(max(1,5*tau/2)))+1))
        x <- seq(-40,40,by=0.05)
        e <- gof_marginal_eform(C,tt,gof_marginal_lanczos(x,stats::dnorm(x),30L))
        for (v in c("P","dc","dtau2")) {
            expect_equal(z[[v]],exact[[v]],tolerance=1e-12)
            if (tau>=1) expect_equal(e[[v]],exact[[v]],tolerance=1e-10)
        }
    }
})

test_that("link-density Gauss rules have the moments of the link distributions", {
    gam <- -digamma(1)
    r <- gof_marginal_density_rule("logit",30L)
    expect_equal(sum(r$w),1)
    expect_equal(sum(r$w*r$x),0,tolerance=1e-12)
    expect_equal(sum(r$w*r$x^2),pi^2/3,tolerance=1e-12)
    expect_equal(sum(r$w*r$x^4),7*pi^4/15,tolerance=1e-12)
    r <- gof_marginal_density_rule("loglog",30L)
    expect_equal(sum(r$w*r$x),gam,tolerance=1e-12)
    expect_equal(sum(r$w*r$x^2)-gam^2,pi^2/6,tolerance=1e-12)
    r <- gof_marginal_density_rule("cloglog",30L)
    expect_equal(sum(r$w*r$x),-gam,tolerance=1e-12)
})

test_that("every link matches integrate() across tau regimes and far tails", {
    cc <- c(-30,-8,-2,-0.4,0,0.9,3,8,30)
    taus <- c(0,1e-7,0.3,1,1.7,2.5,4.5,9)
    g <- expand.grid(c=cc,tau=taus)
    for (nm in gof_link_names) {
        lk <- gof_ordinal_link(nm)
        I <- gof_marginal_integrals(nm,cbind(g$c),g$tau,nodes=30L)
        expect_true(all(is.finite(unlist(I))),info=nm)
        expect_true(all(I$P>=0 & I$P<=1),info=nm)
        expect_equal(I$P[,1],marginal_oracle(lk$F,g$c,g$tau),tolerance=1e-8,info=nm)
        expect_equal(I$dc[,1],marginal_oracle(lk$f,g$c,g$tau),tolerance=1e-8,info=nm)
        mod <- abs(g$c)<=8
        ## The contract's Gumbel df gives 0 * Inf = NaN far in a tail; its limit is 0.
        df <- function(x) {v <- lk$df(x); v[!is.finite(v)] <- 0; v}
        expect_equal(2*I$dtau2[mod,1],marginal_oracle(df,g$c[mod],g$tau[mod]),
                     tolerance=1e-8,info=nm)
    }
    ## Extreme arguments stay finite and ordered.
    g <- expand.grid(c=c(-1e3,-60,0,60,1e3),tau=c(0,1e-12,1,3,30,1e4))
    for (nm in gof_link_names) {
        I <- gof_marginal_integrals(nm,cbind(g$c,g$c+1),g$tau)
        expect_true(all(is.finite(unlist(I))),info=nm)
        expect_true(all(I$P[,2]>=I$P[,1]),info=nm)
    }
})

test_that("the analytic gradient matches numDeriv for all links and structures", {
    for (s in c("intercept","slope","nested")) for (nm in gof_link_names) {
        p <- marginal_problem(s,nm)
        out <- gof_clmm_marginal(p)
        expect_identical(dim(out$gradient),c(p$n*(p$K-1L),length(p$theta)))
        expect_identical(colnames(out$gradient),names(p$theta))
        J <- -numDeriv::jacobian(function(t) c(gof_clmm_marginal(p,theta=t)$cdf),p$theta)
        expect_equal(unname(out$gradient),J,tolerance=1e-8,info=paste(s,nm))
    }
    ## Large random-effect variances reach the e-form (logit tau >= 4,
    ## Gumbel tau >= 2) and the switch between forms.
    for (nm in c("logit","cloglog","loglog")) {
        p <- marginal_problem("slope",nm,sd=if (nm=="logit") 3 else 1.6,sd2=1.2)
        tau <- sqrt(gof_clmm_marginal(p)$tau2)
        expect_true(any(tau>=gof_marginal_switch[[nm]]) && any(tau<gof_marginal_switch[[nm]]),info=nm)
        J <- -numDeriv::jacobian(function(t) c(gof_clmm_marginal(p,theta=t)$cdf),p$theta)
        expect_equal(unname(gof_clmm_marginal(p)$gradient),J,tolerance=1e-8,info=nm)
    }
})

test_that("vanishing random effects give the conditional cumulative probabilities", {
    for (s in c("intercept","slope","nested")) for (nm in gof_link_names) {
        p <- marginal_problem(s,nm)
        theta <- p$theta
        for (tm in p$terms) {
            om <- theta[tm$omega]
            om[gof_logchol_names(tm$q) %in% paste0("L",seq_len(tm$q),seq_len(tm$q))] <- -30
            om[!gof_logchol_names(tm$q) %in% paste0("L",seq_len(tm$q),seq_len(tm$q))] <- 0
            theta[tm$omega] <- om
        }
        out <- gof_clmm_marginal(p,theta=theta)
        F <- gof_ordinal_link(nm)$F
        cond <- F(marginal_c(p,theta))
        expect_equal(out$cdf,cond,tolerance=1e-12,info=paste(s,nm))
        omega <- unlist(lapply(p$terms,`[[`,"omega"))
        expect_lt(max(abs(out$gradient[,omega])),1e-12)
        ## Fixed-effect columns reduce to the conditional derivatives.
        f <- gof_ordinal_link(nm)$f(marginal_c(p,theta))
        expect_equal(unname(out$gradient[,p$index$beta]),c(f)*rep(p$X[,1],p$K-1L),tolerance=1e-10)
    }
    ## An exactly zero design gives tau = 0 and the exact conditional values.
    p <- marginal_problem("intercept","logit")
    p$terms[[1L]]$Z[] <- 0
    out <- gof_clmm_marginal(p)
    expect_identical(out$tau2,numeric(p$n))
    expect_identical(out$cdf,stats::plogis(marginal_c(p)))
    expect_identical(unname(out$gradient[,p$terms[[1L]]$omega]),numeric(p$n*(p$K-1L)))
})

test_that("marginal cumulative probabilities increase in k", {
    for (s in c("intercept","slope","nested")) for (nm in gof_link_names) {
        out <- gof_clmm_marginal(marginal_problem(s,nm,K=6L,sd=2))
        expect_true(all(out$cdf[,-1L]>=out$cdf[,-ncol(out$cdf)]),info=paste(s,nm))
        expect_true(all(out$cdf>0 & out$cdf<1),info=paste(s,nm))
    }
})

test_that("random-effect variance flattens the marginal toward 1/2", {
    p <- marginal_problem("intercept","logit",offset=FALSE,equidistant=FALSE)
    lo <- gof_clmm_marginal(p,theta=replace(p$theta,p$terms[[1L]]$omega,log(0.2)))$cdf
    hi <- gof_clmm_marginal(p,theta=replace(p$theta,p$terms[[1L]]$omega,log(3)))$cdf
    expect_true(all(abs(hi-0.5)<=abs(lo-0.5)+1e-15))
})

test_that("link conventions agree with ordinal's distribution functions", {
    skip_if_not_installed("ordinal")
    x <- c(-30,-6,-1.3,0,0.4,2,6,30)
    ## ordinal evaluates 1 - exp(-exp(x)) without expm1, hence 1e-12.
    expect_equal(gof_ordinal_link("loglog")$F(x),ordinal::pgumbel(x),tolerance=1e-12)
    expect_equal(gof_ordinal_link("cloglog")$F(x),ordinal::pgumbel(x,max=FALSE),tolerance=1e-12)
    expect_equal(gof_ordinal_link("loglog")$f(x),ordinal::dgumbel(x),tolerance=1e-12)
    expect_equal(gof_ordinal_link("cloglog")$f(x),ordinal::dgumbel(x,max=FALSE),tolerance=1e-12)
    ## Conditional marginals (tau = 0) use exactly these CDFs.
    for (nm in c("cloglog","loglog")) {
        p <- marginal_problem("intercept",nm)
        p$terms[[1L]]$Z[] <- 0
        ref <- ordinal::pgumbel(marginal_c(p),max=nm=="loglog")
        expect_equal(gof_clmm_marginal(p)$cdf,matrix(ref,p$n),tolerance=1e-12,info=nm)
    }
    ## ordinal's clm uses the same parametrisation: P(Y <= k) = F(alpha_k - x'beta).
    d <- clmm_simulate("intercept",m=40L,link="cloglog",seed=5L)
    fit <- ordinal::clm(y~x,data=d,link="cloglog")
    prob <- stats::predict(fit,newdata=d["x"],type="cum.prob")$cprob1
    p <- marginal_problem("intercept","cloglog",m=40L,seed=5L,offset=FALSE,equidistant=FALSE)
    p$terms[[1L]]$Z[] <- 0
    theta <- p$theta
    theta[p$index$xi] <- fit$alpha
    theta[p$index$beta] <- fit$beta
    out <- gof_clmm_marginal(p,theta=theta)
    expect_equal(out$cdf,prob[,seq_len(p$K-1L)],tolerance=1e-10,ignore_attr=TRUE)
})

test_that("gof_clmm_marginal validates its inputs", {
    p <- marginal_problem("intercept","logit")
    expect_error(gof_clmm_marginal(list()),"gof_clmm_problem")
    expect_error(gof_clmm_marginal(p,theta=replace(p$theta,1L,NA)),"finite")
    expect_error(gof_clmm_marginal(p,theta=p$theta[-1L]),"wrong length")
    expect_error(gof_clmm_marginal(p,nodes=1),"nodes")
})
