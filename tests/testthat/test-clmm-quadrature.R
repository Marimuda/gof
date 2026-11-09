## Oracles for the nested adaptive Gauss--Hermite integrator.

## Problem from clmm_simulate() data at given parameters (no ordinal needed).
agq_problem <- function(d,structure,xi,beta,omega,weights=NULL,offset=NULL,link="logit",
                        tJac=NULL,K=nlevels(d$y)) {
    n <- nrow(d)
    top <- as.integer(d$g)
    m <- max(top)
    one <- matrix(1,n,1L,dimnames=list(NULL,"(Intercept)"))
    terms <- switch(structure,
        intercept=list(list(label="1 | g",group="g",index=top,n_levels=m,Z=one)),
        slope=list(list(label="1 + z | g",group="g",index=top,n_levels=m,
                        Z=cbind("(Intercept)"=1,z=d$z))),
        nested={
            hg <- as.integer(interaction(d$h,d$g,drop=TRUE,lex.order=TRUE))
            list(list(label="1 | g",group="g",index=top,n_levels=m,Z=one),
                 list(label="1 | h:g",group="h:g",index=hg,n_levels=max(hg),Z=one))
        })
    if (is.null(tJac)) tJac <- diag(K-1L)
    theta <- c(stats::setNames(xi,paste0("xi",seq_along(xi))),x=unname(beta),
               stats::setNames(omega,paste0("omega",seq_along(omega))))
    new_gof_clmm_problem(y=as.integer(d$y),K=K,X=cbind(x=d$x),
                         offset=if (is.null(offset)) numeric(n) else offset,
                         weights=if (is.null(weights)) rep(1,n) else weights,
                         tJac=tJac,link=link,terms=terms,blocks=list(id=top,m=m),
                         structure="nesting-chain",
                         chain=if (structure=="nested") c("g","h:g") else "g",
                         theta=theta,levels=levels(d$y))
}

agq_problem_from_fit <- function(fit,d,structure) {
    ST <- fit$ST
    om <- switch(structure,intercept=log(ST$g[1L]),
                 slope=gof_logchol_omega(tcrossprod(ST$g)),
                 nested=c(log(ST$g[1L]),log(ST$`h:g`[1L])))
    agq_problem(d,structure,fit$alpha,fit$beta,om,link=fit$link)
}

## Row log-likelihood contributions from the plain CDF (independent oracle).
agq_rowlog <- function(p,theta,b) {
    up <- gof_clmm_unpack(p,theta)
    lk <- gof_ordinal_link(p$link)
    a <- c(-Inf,up$alpha,Inf)
    eta <- drop(p$X %*% up$beta)+p$offset+b
    p$weights*log(lk$F(a[p$y+1L]-eta)-lk$F(a[p$y]-eta))
}

## Log-likelihood of `rows` for a matrix of random linear predictors
## (rows x evaluation points), from the plain CDF.
agq_rowlog_mat <- function(p,theta,rows,B) {
    up <- gof_clmm_unpack(p,theta)
    lk <- gof_ordinal_link(p$link)
    a <- c(-Inf,up$alpha,Inf)
    eta <- drop(p$X[rows,,drop=FALSE] %*% up$beta)+p$offset[rows]+B
    yy <- p$y[rows]
    colSums(p$weights[rows]*log(lk$F(a[yy+1L]-eta)-lk$F(a[yy]-eta)))
}

test_that("Gauss-Hermite rules are exact to degree 2Q-1 and integrate exp(tZ)", {
    dfact <- function(k) if (k<=0) 1 else prod(seq(k-1,1,by=-2))
    for (Q in c(1:12,20L,30L)) {
        gh <- gof_gauss_hermite(Q)
        expect_length(gh$nodes,Q)
        expect_equal(sum(gh$weights),1,tolerance=1e-14)
        expect_true(all(gh$weights>0))
        for (k in 0:(2*Q-1)) {
            truth <- if (k%%2) 0 else dfact(k)
            scale <- sum(gh$weights*abs(gh$nodes)^k)
            expect_lt(abs(sum(gh$weights*gh$nodes^k)-truth),1e-11*max(1,scale),
                      label=paste("Q",Q,"degree",k))
        }
        if (Q>=2) {
            ## Degree 2Q is not exact.
            expect_gt(abs(sum(gh$weights*gh$nodes^(2*Q))-dfact(2*Q)),1e-6)
        }
    }
    gh <- gof_gauss_hermite(20)
    for (t in c(-1,.5,1,2)) expect_equal(sum(gh$weights*exp(t*gh$nodes)),exp(t^2/2),tolerance=1e-12)
    grid <- gof_gh_grid(5,2)
    expect_equal(nrow(grid$z),25L)
    expect_equal(sum(exp(grid$logw)*grid$z[,1L]^2*grid$z[,2L]^4),3,tolerance=1e-12)
})

test_that("the row kernel is stable and differentiates correctly for all links", {
    for (nm in gof_link_names) {
        lk <- gof_clmm_link(nm)
        alpha <- c(-1,.2,1.5)
        y <- rep(1:4,each=5)
        eta <- rep(c(-3,-.5,0,.7,4),4)
        k <- gof_clmm_kernel(eta,y,alpha,4L,lk,2L)
        a <- c(-Inf,alpha,Inf)
        ## Oracle from the tail that avoids cancellation.
        hi <- a[y+1L]-eta; lo <- a[y]-eta
        up <- (hi+lo)>0 & y>1L & y<4L
        naive <- ifelse(up,log(exp(lk$logS(lo))-exp(lk$logS(hi))),log(lk$F(hi)-lk$F(lo)))
        ok <- is.finite(naive)
        expect_true(all(is.finite(k$lp)))
        expect_equal(k$lp[ok],naive[ok],tolerance=1e-12)
        h <- 1e-5
        kp <- gof_clmm_kernel(eta+h,y,alpha,4L,lk,1L)
        km <- gof_clmm_kernel(eta-h,y,alpha,4L,lk,1L)
        expect_equal(k$d1,(kp$lp-km$lp)/(2*h),tolerance=1e-7)
        expect_equal(k$d2,(kp$d1-km$d1)/(2*h),tolerance=1e-6)
        ## Far tails stay finite on the log scale.
        far <- gof_clmm_kernel(c(-40,40),c(4L,1L),alpha,4L,lk,1L)
        expect_true(all(is.finite(far$lp)))
    }
})

test_that("1-D AGQ matches integrate() on tiny clusters", {
    d <- clmm_simulate("intercept",m=6L,nc=3L,seed=11)
    d <- d[-c(1,2,7),]                         # clusters of 1 to 3 rows
    d$g <- droplevels(d$g)
    w <- seq(.5,2,length.out=nrow(d))
    o <- seq(-.4,.4,length.out=nrow(d))
    ## Default 11 nodes: worst errors on 1-3 row clusters are about 2e-7
    ## (logit), 6e-6 (probit) and 4e-4 (cauchit); high Q converges to the oracle.
    for (link in c("logit","probit","cauchit")) {
        p <- agq_problem(d,"intercept",c(-1.2,.1,1.4),.7,log(1.3),weights=w,offset=o,link=link)
        sig <- exp(p$theta[["omega1"]])
        brute <- vapply(seq_len(p$blocks$m),function(c) {
            r <- which(p$blocks$id==c)
            f <- function(b) exp(agq_rowlog_mat(p,p$theta,r,matrix(b,length(r),length(b),byrow=TRUE))+
                                 stats::dnorm(b,0,sig,log=TRUE))
            log(stats::integrate(f,-Inf,Inf,rel.tol=1e-12)$value)
        },0)
        ll <- function(Q) {
            ig <- gof_integrator_nested_agq(p,Q)
            ig$loglik(p$theta,ig$refresh(p$theta))
        }
        tol11 <- c(logit=1e-6,probit=2e-5,cauchit=2e-3)[[link]]
        expect_lt(max(abs(ll(11L)-brute)),tol11,label=link)
        expect_equal(ll(if (link=="cauchit") 81L else 41L),brute,tolerance=1e-9,info=link)
    }
})

test_that("nested AGQ matches nested integrate()", {
    d <- clmm_simulate("nested",m=3L,nc=4L,seed=5)
    p <- agq_problem(d,"nested",c(-1,0,1.2),.6,log(c(1.1,.8)),weights=rep(c(1,2),6),
                     offset=rep(c(.2,-.1),6))
    s1 <- exp(p$theta[["omega1"]]); s2 <- exp(p$theta[["omega2"]])
    hg <- p$terms[[2L]]$index
    brute <- vapply(seq_len(p$blocks$m),function(c) {
        rows <- which(p$blocks$id==c)
        inner <- function(b1) prod(vapply(unique(hg[rows]),function(s) {
            r <- rows[hg[rows]==s]
            f <- function(b2) exp(agq_rowlog_mat(p,p$theta,r,matrix(b1+b2,length(r),length(b2),byrow=TRUE))+
                                  stats::dnorm(b2,0,s2,log=TRUE))
            stats::integrate(f,-Inf,Inf,rel.tol=1e-11)$value
        },0))
        outer <- function(b1) vapply(b1,function(bb) inner(bb)*stats::dnorm(bb,0,s1),0)
        log(stats::integrate(outer,-Inf,Inf,rel.tol=1e-10)$value)
    },0)
    ig <- gof_integrator_nested_agq(p)
    expect_equal(ig$loglik(p$theta,ig$refresh(p$theta)),brute,tolerance=1e-6)
    ig <- gof_integrator_nested_agq(p,c(25L,25L))
    expect_equal(ig$loglik(p$theta,ig$refresh(p$theta)),brute,tolerance=1e-9)
})

test_that("2-D slope AGQ matches nested integrate() cubature", {
    d <- clmm_simulate("slope",m=3L,nc=4L,seed=8)
    S <- matrix(c(1,.3,.3,.5),2L)
    p <- agq_problem(d,"slope",c(-1,.2,1.3),.5,gof_logchol_omega(S))
    Si <- solve(S)
    brute <- vapply(seq_len(p$blocks$m),function(c) {
        r <- which(p$blocks$id==c)
        inner <- function(b1) {
            f <- function(b2) {
                B <- outer(rep(1,length(r)),rep(b1,length(b2)))+outer(d$z[r],b2)
                q <- Si[1,1]*b1^2+2*Si[1,2]*b1*b2+Si[2,2]*b2^2
                exp(agq_rowlog_mat(p,p$theta,r,B)-log(2*pi)-log(det(S))/2-q/2)
            }
            stats::integrate(f,-Inf,Inf,rel.tol=1e-11)$value
        }
        log(stats::integrate(function(b1) vapply(b1,inner,0),-Inf,Inf,rel.tol=1e-10)$value)
    },0)
    ig <- gof_integrator_nested_agq(p,15L)
    expect_equal(ig$loglik(p$theta,ig$refresh(p$theta)),brute,tolerance=1e-8)
    ## Default 7 nodes per dimension is already close.
    ig7 <- gof_integrator_nested_agq(p)
    expect_equal(ig7$loglik(p$theta,ig7$refresh(p$theta)),brute,tolerance=1e-5)
})

test_that("Q = 1 reproduces clmm's Laplace log-likelihood", {
    skip_if_not_installed("ordinal")
    ctl <- ordinal::clmm.control(gradTol=1e-12)
    for (s in c("intercept","slope","nested")) {
        d <- clmm_simulate(s,m=40L,nc=6L,seed=3)
        fm <- switch(s,intercept=y~x+(1|g),slope=y~x+(1+z|g),nested=y~x+(1|g/h))
        ## Tight inner tolerance: clmm's default gradTol = 1e-4 leaves its
        ## Laplace value about 6e-4 off for the scalar-intercept path.
        fit <- suppressWarnings(ordinal::clmm(fm,data=d,control=ctl))
        p <- agq_problem_from_fit(fit,d,s)
        ig <- gof_integrator_nested_agq(p,1L)
        expect_equal(sum(ig$loglik(p$theta,ig$refresh(p$theta))),c(stats::logLik(fit)),
                     tolerance=1e-8,info=s)
    }
})

test_that("AGQ with 10 nodes reproduces clmm(nAGQ = 10)", {
    skip_if_not_installed("ordinal")
    d <- clmm_simulate("intercept",m=40L,nc=6L,seed=3)
    fit <- ordinal::clmm(y~x+(1|g),data=d,nAGQ=10)
    p <- agq_problem_from_fit(fit,d,"intercept")
    ig <- gof_integrator_nested_agq(p,10L)
    expect_equal(sum(ig$loglik(p$theta,ig$refresh(p$theta))),c(stats::logLik(fit)),tolerance=1e-8)
})

test_that("analytic frozen-node scores equal numerical gradients for all links", {
    skip_if_not_installed("numDeriv")
    for (link in gof_link_names) for (s in c("intercept","slope","nested")) {
        d <- clmm_simulate(s,m=5L,nc=4L,seed=21,link=link)
        n <- nrow(d)
        omega <- switch(s,intercept=log(.9),slope=gof_logchol_omega(matrix(c(1,.2,.2,.4),2L)),
                        nested=log(c(.9,.6)))
        equi <- s=="slope"
        p <- agq_problem(d,s,if (equi) c(-1,1.1) else c(-1,.1,1.2),.6,omega,
                         weights=seq(.5,1.5,length.out=n),offset=sin(seq_len(n)),link=link,
                         tJac=if (equi) cbind(1,0:2) else NULL,K=4L)
        ig <- gof_integrator_nested_agq(p)
        st <- ig$refresh(p$theta)
        ## Evaluate away from the refresh point so the frozen state matters.
        th <- p$theta+seq(-.05,.05,length.out=length(p$theta))
        sc <- ig$score(th,st)
        expect_named(sc,"g")
        num <- numDeriv::jacobian(function(t) ig$loglik(t,st),th)
        expect_equal(unname(sc$g),num,tolerance=1e-6,info=paste(link,s))
    }
})

test_that("three-level nesting chains reduce and differentiate correctly", {
    ## Validated offline against a triple nested integrate() (error 9e-11 at
    ## 9 nodes per level); here the faster reductions and the score.
    set.seed(2)
    m <- 3L; n <- m*8L
    g <- rep(seq_len(m),each=8L); h <- rep(seq_len(2L*m),each=4L); k <- rep(seq_len(4L*m),each=2L)
    x <- stats::rnorm(n); y <- sample(1:3,n,TRUE); y[1:3] <- 1:3
    one <- matrix(1,n,1L,dimnames=list(NULL,"(Intercept)"))
    tk <- list(label="1 | k",group="k",index=k,n_levels=4L*m,Z=one)
    tg <- list(label="1 | g",group="g",index=g,n_levels=m,Z=one)
    th <- list(label="1 | h",group="h",index=h,n_levels=2L*m,Z=one)
    mk <- function(terms,chain,theta)
        new_gof_clmm_problem(y=y,K=3L,X=cbind(x=x),offset=numeric(n),weights=rep(1,n),tJac=diag(2L),
                             link="probit",terms=terms,blocks=list(id=g,m=m),structure="nesting-chain",
                             chain=chain,theta=theta)
    ll <- function(p,nodes=NULL) {
        ig <- gof_integrator_nested_agq(p,nodes)
        ig$loglik(p$theta,ig$refresh(p$theta))
    }
    base <- c(a1=-.4,a2=.6,x=.5)
    p3 <- mk(list(tk,tg,th),c("g","h","k"),c(base,k=log(.6),g=log(.9),h=log(.7)))
    ## Tiny finest or middle variance: the two-level chain without it.
    expect_equal(ll(mk(list(tk,tg,th),c("g","h","k"),c(base,k=-20,g=log(.9),h=log(.7)))),
                 ll(mk(list(tg,th),c("g","h"),c(base,g=log(.9),h=log(.7)))),tolerance=1e-12)
    expect_equal(ll(mk(list(tk,tg,th),c("g","h","k"),c(base,k=log(.6),g=log(.9),h=-20))),
                 ll(mk(list(tg,tk),c("g","k"),c(base,g=log(.9),k=log(.6)))),tolerance=1e-12)
    expect_equal(ll(p3,9L),ll(p3,15L),tolerance=1e-8)
    skip_if_not_installed("numDeriv")
    ig <- gof_integrator_nested_agq(p3,c(7L,5L,5L))
    st <- ig$refresh(p3$theta)
    th3 <- p3$theta+.03
    expect_equal(unname(ig$score(th3,st)$g),numDeriv::jacobian(function(t) ig$loglik(t,st),th3),
                 tolerance=1e-6)
    expect_identical(ig$diagnostics(st)$grid,c(7L,5L,5L))
})

test_that("information is the symmetric negative Jacobian of the score", {
    skip_if_not_installed("numDeriv")
    p <- clmm_problem_fixture("nested",m=10L,nc=4L)
    ig <- gof_integrator_nested_agq(p)
    st <- ig$refresh(p$theta)
    H <- ig$information(p$theta,st)
    expect_equal(H,t(H))
    expect_lt(attr(H,"asymmetry"),1e-6)
    hes <- numDeriv::hessian(function(t) sum(ig$loglik(t,st)),p$theta)
    expect_equal(unname(H),-hes,tolerance=1e-5,ignore_attr=TRUE)
    dg <- ig$diagnostics(st)
    expect_identical(dg$grid,c(11L,11L))
    expect_true(all(dg$unconverged==0L))
})

test_that("a vanishing covariance reproduces clm's log-likelihood", {
    skip_if_not_installed("ordinal")
    for (s in c("intercept","nested")) {
        d <- clmm_simulate(s,m=20L,nc=4L,seed=2)
        fit <- ordinal::clm(y~x,data=d)
        p <- agq_problem(d,s,fit$alpha,fit$beta,rep(-20,if (s=="nested") 2L else 1L))
        ig <- gof_integrator_nested_agq(p)
        expect_equal(sum(ig$loglik(p$theta,ig$refresh(p$theta))),c(stats::logLik(fit)),tolerance=1e-10)
    }
})

test_that("duplicating a row equals giving it weight 2", {
    for (s in c("intercept","nested")) {
        d <- clmm_simulate(s,m=6L,nc=4L,seed=9)
        th <- c(-1,0,1.1)
        om <- if (s=="nested") log(c(1,.7)) else log(1.2)
        dup <- d[c(seq_len(nrow(d)),5L),]
        dup <- dup[order(as.integer(dup$g)),]
        p1 <- agq_problem(dup,s,th,.5,om)
        p2 <- agq_problem(d,s,th,.5,om,weights=replace(rep(1,nrow(d)),5L,2))
        i1 <- gof_integrator_nested_agq(p1); s1 <- i1$refresh(p1$theta)
        i2 <- gof_integrator_nested_agq(p2); s2 <- i2$refresh(p2$theta)
        expect_equal(i1$loglik(p1$theta,s1),i2$loglik(p2$theta,s2),tolerance=1e-12)
        expect_equal(i1$score(p1$theta,s1),i2$score(p2$theta,s2),tolerance=1e-10)
    }
})

test_that("term order and shared grouping factors follow problem$chain", {
    d <- clmm_simulate("nested",m=6L,nc=4L,seed=13)
    n <- nrow(d)
    top <- as.integer(d$g)
    hg <- as.integer(interaction(d$h,d$g,drop=TRUE,lex.order=TRUE))
    one <- matrix(1,n,1L,dimnames=list(NULL,"(Intercept)"))
    t_g <- list(label="1 | g",group="g",index=top,n_levels=6L,Z=one)
    t_z <- list(label="0 + z | g",group="g",index=top,n_levels=6L,Z=cbind(z=d$z))
    t_h <- list(label="1 | h:g",group="h:g",index=hg,n_levels=max(hg),Z=one)
    build <- function(terms,omega) {
        theta <- c(xi1=-1,xi2=.1,xi3=1.2,x=.5,omega)
        new_gof_clmm_problem(y=as.integer(d$y),K=4L,X=cbind(x=d$x),offset=numeric(n),weights=rep(1,n),
                             tJac=diag(3L),link="logit",terms=terms,blocks=list(id=top,m=6L),
                             structure="nesting-chain",chain=c("g","h:g"),theta=theta)
    }
    ## clmm order: finest first, then the shared-group pair.
    pa <- build(list(t_h,t_z,t_g),c(h=log(.7),z=log(.5),g=log(1.1)))
    pb <- build(list(t_g,t_z,t_h),c(g=log(1.1),z=log(.5),h=log(.7)))
    ia <- gof_integrator_nested_agq(pa); sa <- ia$refresh(pa$theta)
    ib <- gof_integrator_nested_agq(pb); sb <- ib$refresh(pb$theta)
    expect_equal(ia$loglik(pa$theta,sa),ib$loglik(pb$theta,sb),tolerance=1e-12)
    sca <- ia$score(pa$theta,sa)$g
    scb <- ib$score(pb$theta,sb)$g
    expect_equal(sca[,names(pb$theta)],scb,tolerance=1e-10)
    ## The shared pair stacks into one 2-vector: (1|g)+(0+z|g) equals
    ## (1+z|g) with zero correlation, up to the extra L21 parameter.
    t_s <- list(label="1 + z | g",group="g",index=top,n_levels=6L,Z=cbind("(Intercept)"=1,z=d$z))
    pc <- build(list(t_s,t_h),c(L11=log(1.1),L21=0,L22=log(.5),h=log(.7)))
    ic <- gof_integrator_nested_agq(pc); sc <- ic$refresh(pc$theta)
    expect_equal(ic$loglik(pc$theta,sc),ib$loglik(pb$theta,sb),tolerance=1e-10)
    expect_equal(unname(ic$score(pc$theta,sc)$g[,c("L11","L22","h")]),unname(scb[,c("g","z","h")]),
                 tolerance=1e-8)
})

test_that("cauchit non-PD Hessians fall back to the expected information", {
    d <- clmm_simulate("intercept",m=4L,nc=4L,seed=1,link="cauchit")
    d$x[1:4] <- 0
    brute1 <- function(p) {
        r <- which(p$blocks$id==1L)
        f <- function(b) exp(agq_rowlog_mat(p,p$theta,r,matrix(b,length(r),length(b),byrow=TRUE))+
                             stats::dnorm(b,0,exp(p$theta[["omega1"]]),log=TRUE))
        log(stats::integrate(f,-Inf,Inf,rel.tol=1e-12,subdivisions=1000L)$value)
    }
    ## Indefinite observed curvature on the way, a proper mode at the end.
    d$y[1:4] <- factor(c(1,4,4,4),levels=1:4)
    p <- agq_problem(d,"intercept",c(-.5,0,.5),0,log(3),link="cauchit")
    ig <- gof_integrator_nested_agq(p,41L)
    st <- ig$refresh(p$theta)
    dg <- ig$diagnostics(st)
    expect_gte(dg$expected_steps,1L)
    expect_identical(dg$expected_information,0L)
    expect_equal(ig$loglik(p$theta,st)[1L],brute1(p),tolerance=1e-5)
    ## Symmetric conflicting outcomes: the stationary point is not a maximum
    ## (bimodal posterior); the scale uses the expected information and the
    ## diagnostics flag it.
    d$y[1:4] <- factor(c(1,1,4,4),levels=1:4)
    p <- agq_problem(d,"intercept",c(-.5,0,.5),0,log(30),link="cauchit")
    ig <- gof_integrator_nested_agq(p)
    st <- ig$refresh(p$theta)
    expect_gte(ig$diagnostics(st)$expected_information,1L)
    expect_true(all(is.finite(ig$loglik(p$theta,st))))
})

test_that("the registry selects nested AGQ and rejects crossed structures", {
    p <- clmm_problem_fixture("nested",m=6L,nc=4L)
    fac <- gof_select_integrator(p,gof_clmm_control())
    expect_identical(attr(fac,"integrator"),"nested-agq")
    expect_identical(fac()$nodes,c(11L,11L))
    expect_identical(fac(c(3L,5L))$nodes,c(3L,5L))
    expect_identical(gof_select_integrator(p,gof_clmm_control(nodes=4L))()$nodes,c(4L,4L))
    s <- clmm_problem_fixture("slope",m=6L,nc=4L)
    expect_identical(gof_select_integrator(s)()$nodes,7L)
    cr <- p
    cr$structure <- "crossed"
    expect_error(gof_select_integrator(cr),"Crossed random effects are not supported")
    expect_error(gof_integrator_nested_agq(cr),"nesting-chain")
    bad <- p
    bad$terms[[2L]]$index <- c(bad$terms[[2L]]$index[-1L],bad$terms[[2L]]$index[1L])
    expect_error(gof_integrator_nested_agq(bad),"not nested")
})
