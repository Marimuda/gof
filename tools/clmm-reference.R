# Development oracle: direct integration under a Gaussian prior, independent
# of gof's likelihood, score, adaptive modes, marginal integrals and optimizer.
reference_rule <- local({
    cache <- list()
    function(Q) {
        key <- as.character(Q)
        if (!is.null(cache[[key]])) return(cache[[key]])
        J <- matrix(0,Q,Q)
        J[cbind(1:(Q-1L),2:Q)] <- sqrt(1:(Q-1L))
        J <- J+t(J)
        e <- eigen(J,symmetric=TRUE)
        o <- order(e$values)
        rule <- list(z=e$values[o],w=e$vectors[1,o]^2)
        cache[[key]] <<- rule
        rule
    }
})

reference_logcdf <- function(x,link,lower=TRUE) switch(link,
    logit=plogis(x,lower.tail=lower,log.p=TRUE),
    probit=pnorm(x,lower.tail=lower,log.p=TRUE),
    cauchit=pcauchy(x,lower.tail=lower,log.p=TRUE),
    cloglog=if(lower) ifelse(x< -35,x,log(-expm1(-exp(x)))) else -exp(x),
    loglog=if(lower) -exp(-x) else ifelse(x>35,-x,log(-expm1(-exp(-x)))))

reference_logprob <- function(eta,y,alpha,link) {
    a <- c(-Inf,alpha,Inf)
    hi <- a[y+1L]-eta
    lo <- a[y]-eta
    Fu <- reference_logcdf(hi,link)
    Fl <- reference_logcdf(lo,link)
    use_survival <- Fl>log(.5)
    value <- Fu+log(-expm1(Fl-Fu))
    if(any(use_survival)) {
        Sl <- reference_logcdf(lo[use_survival],link,FALSE)
        Su <- reference_logcdf(hi[use_survival],link,FALSE)
        value[use_survival] <- Sl+log(-expm1(Su-Sl))
    }
    value
}

reference_start <- function(cell,alpha=c(-1.5,0,1.5),beta=.8,ST=NULL) {
    if(is.null(ST)) ST <- switch(cell$structure,
        intercept=cell$sd,slope=c(cell$sd,.15,.5*sqrt(.91)),nested=c(cell$sd,.5))
    c(alpha[1],log(diff(alpha)),beta,ST)
}

reference_unpack <- function(phi) list(alpha=c(phi[1],phi[1]+cumsum(exp(phi[2:3]))),
                                      beta=phi[4],ST=phi[-(1:4)])

reference_loglik <- function(data,cell,phi,Q=51L) {
    up <- reference_unpack(phi)
    if(any(!is.finite(phi)) || any(!is.finite(up$alpha)))return(-Inf)
    rule <- reference_rule(Q)
    lse <- function(x,w) {mx<-max(x); if(!is.finite(mx))return(mx); mx+log(sum(exp(x-mx)*w))}
    y <- as.integer(data$y)
    if(cell$structure!='nested') {
        if(cell$structure=='slope') {
            index <- expand.grid(i=seq_len(Q),j=seq_len(Q))
            b0 <- up$ST[1]*rule$z[index$i]
            b1 <- up$ST[2]*rule$z[index$i]+up$ST[3]*rule$z[index$j]
            w <- rule$w[index$i]*rule$w[index$j]
        } else {b0<-up$ST[1]*rule$z;b1<-0*b0;w<-rule$w}
        eta <- outer(up$beta*data$x,rep(1,length(b0)))+
            outer(rep(1,nrow(data)),b0)+outer(data$z,b1)
        lp <- matrix(reference_logprob(as.vector(eta),rep(y,length(b0)),up$alpha,cell$link),nrow(data))
        blocks <- rowsum(lp,data$g,reorder=FALSE)
        return(sum(apply(blocks,1,lse,w=w)))
    }
    index <- expand.grid(i=seq_len(Q),j=seq_len(Q))
    re <- up$ST[1]*rule$z[index$j]+up$ST[2]*rule$z[index$i]
    total <- 0
    for(rows in split(seq_len(nrow(data)),data$g)) {
        parent <- rep(0,Q)
        for(child in split(rows,data$h[rows])) {
            eta <- outer(up$beta*data$x[child],rep(1,length(re)))+
                matrix(re,length(child),length(re),byrow=TRUE)
            lp <- matrix(reference_logprob(as.vector(eta),rep(y[child],length(re)),up$alpha,cell$link),length(child))
            child_ll <- matrix(colSums(lp),Q,Q)
            parent <- parent+apply(child_ll,2,lse,w=rule$w)
        }
        total <- total+lse(parent,rule$w)
    }
    total
}

reference_fit <- function(data,cell,starts,Q=31L,refine=61L) {
    lower <- c(rep(-Inf,4),switch(cell$structure,intercept=0,slope=c(0,-Inf,0),nested=c(0,0)))
    objective <- function(phi) {
        value <- reference_loglik(data,cell,phi,Q)
        if(is.finite(value)) -value else .Machine$double.xmax/1000
    }
    fits <- lapply(starts,function(phi) tryCatch(nlminb(phi,objective,lower=lower,
        control=list(iter.max=300L,eval.max=1000L,rel.tol=1e-10)),error=function(e)e))
    ok <- which(vapply(fits,function(x)is.list(x)&&!inherits(x,'error')&&is.finite(x$objective),TRUE))
    if(!length(ok))return(list(fits=fits,success=FALSE))
    winner <- ok[which.min(vapply(fits[ok],`[[`,0,'objective'))]
    best <- fits[[winner]]
    list(success=best$convergence==0L,theta=best$par,loglik=-best$objective,
         refined_loglik=reference_loglik(data,cell,best$par,refine),nodes=Q,refinement_nodes=refine,
         fits=fits,gradient=numDeriv::grad(function(phi)reference_loglik(data,cell,phi,refine),best$par))
}
