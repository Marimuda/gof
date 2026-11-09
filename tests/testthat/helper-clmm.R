## Fixtures for cumulative link mixed-model components. These build
## gof_clmm_problem objects directly from simulated data at known parameters,
## so integrator/marginal/engine tests do not depend on ordinal::clmm.

clmm_simulate <- function(structure=c("intercept","slope","nested"),m=30L,nc=6L,K=4L,
                          link="logit",seed=1L,sd=1,sd2=0.5,rho=0.3,beta=0.8,
                          weights=NULL,offset=NULL) {
    structure <- match.arg(structure)
    set.seed(seed)
    lk <- gof_ordinal_link(link)
    n <- m*nc
    top <- rep(seq_len(m),each=nc)
    x <- stats::rnorm(n)
    z <- stats::rnorm(n)
    eta <- beta*x
    if (!is.null(offset)) eta <- eta+offset
    if (structure=="intercept") {
        b <- stats::rnorm(m,0,sd)
        eta <- eta+b[top]
    } else if (structure=="slope") {
        S <- matrix(c(sd^2,rho*sd*sd2,rho*sd*sd2,sd2^2),2L)
        B <- matrix(stats::rnorm(2L*m),m) %*% chol(S)
        eta <- eta+B[top,1L]+B[top,2L]*z
    } else {
        stopifnot(nc%%2L==0L)
        sub <- rep(seq_len(2L*m),each=nc/2L)
        eta <- eta+stats::rnorm(m,0,sd)[top]+stats::rnorm(2L*m,0,sd2)[sub]
    }
    cut <- seq(-1.5,1.5,length.out=K-1L)
    u <- stats::runif(n)
    ## Inverse-CDF draw: Y = smallest k with u <= F(alpha_k - eta).
    cdf <- sapply(cut,function(a) lk$F(a-eta))
    y <- 1L+rowSums(u>matrix(cdf,n))
    d <- data.frame(y=factor(y,levels=seq_len(K),ordered=TRUE),x=x,z=z,
                    g=factor(top),stringsAsFactors=FALSE)
    if (structure=="nested") d$h <- factor(rep(rep(1:2,each=nc/2L),m))
    if (!is.null(weights)) d$w <- weights
    if (!is.null(offset)) d$o <- offset
    d
}

## Build a problem directly (no clmm) at the generating parameters.
clmm_problem_fixture <- function(structure=c("intercept","slope","nested"),m=30L,nc=6L,K=4L,
                                 link="logit",seed=1L,sd=1,sd2=0.5,rho=0.3,beta=0.8) {
    structure <- match.arg(structure)
    d <- clmm_simulate(structure,m,nc,K,link,seed,sd,sd2,rho,beta)
    n <- nrow(d)
    top <- as.integer(d$g)
    one <- matrix(1,n,1L,dimnames=list(NULL,"(Intercept)"))
    terms <- switch(structure,
        intercept=list(list(label="1 | g",group="g",index=top,n_levels=m,Z=one)),
        slope=list(list(label="1 + z | g",group="g",index=top,n_levels=m,
                        Z=cbind("(Intercept)"=1,z=d$z))),
        nested=list(list(label="1 | g",group="g",index=top,n_levels=m,Z=one),
                    list(label="1 | h:g",group="h:g",index=as.integer(interaction(d$h,d$g,drop=TRUE,lex.order=TRUE)),
                         n_levels=2L*m,Z=one)))
    cut <- seq(-1.5,1.5,length.out=K-1L)
    omega <- switch(structure,
        intercept=c("1 | g.L11"=log(sd)),
        slope=stats::setNames(gof_logchol_omega(matrix(c(sd^2,rho*sd*sd2,rho*sd*sd2,sd2^2),2L)),
                              paste0("1 + z | g.",gof_logchol_names(2L))),
        nested=c("1 | g.L11"=log(sd),"1 | h:g.L11"=log(sd2)))
    theta <- c(stats::setNames(cut,paste0(seq_len(K-1L),"|",seq_len(K-1L)+1L)),x=beta,omega)
    new_gof_clmm_problem(y=as.integer(d$y),K=K,X=cbind(x=d$x),offset=numeric(n),weights=rep(1,n),
                         tJac=diag(K-1L),link=link,terms=terms,
                         blocks=list(id=top,m=m),structure="nesting-chain",
                         chain=if (structure=="nested") c("g","h:g") else "g",
                         theta=theta,levels=as.character(seq_len(K)),source=list(data=d))
}
