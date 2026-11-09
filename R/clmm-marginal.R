## Marginal cumulative probabilities of a cumulative link mixed model and
## their analytic derivatives with respect to theta = c(xi, beta, omega).
##
## For one observation the random part sum_t z_it' b_t is N(0, tau_i^2) with
## tau_i^2 = sum_t z_it' Sigma_t z_it, whatever the relation between the
## grouping factors (nested or crossed). With c_ik = alpha_k - eta_i and a
## link variable e ~ F independent of Z ~ N(0, 1),
##   P_ik = P(e + tau_i Z <= c_ik) = E F(c - tau Z)            (Z-form)
##                                = E_e Phi((c - e) / tau)     (e-form),
##   dP/dc = E f(c - tau Z),  dP/d(tau^2) = E f'(c - tau Z) / 2
## (Stein's identity; regular at tau = 0). These are one-dimensional
## integrals; which variable to integrate over depends on tau:
##
## * probit: closed form, P = Phi(c / s), s = sqrt(1 + tau^2).
## * Otherwise the Z-form is integrated by the trapezoidal rule on
##   [-10, 10], which converges exponentially for integrands analytic in a
##   strip: F(c - tau z) is analytic for |Im z| < D/tau (D = pi for logit,
##   pi/2 for cloglog/loglog, 1 for cauchit), so the step is
##   h = min(1/2, 0.18 D/tau), giving errors near 1e-13. Gauss--Hermite in Z
##   is not used: with 30 nodes it is accurate only for tau < ~1 (errors
##   ~1e-3 at tau = 5 for logit).
## * For large tau (logit tau >= 4, cloglog/loglog tau >= 2) the e-form is
##   cheaper: Phi((c - e)/tau) is then smooth on the scale of the link
##   density, and a Gauss rule for the density f itself (Lanczos on a fine
##   discretisation, cached) is accurate to ~1e-12 with 30 nodes. Both rules
##   agree to that accuracy at the switch, so the result is smooth in theta
##   for all practical purposes. The Cauchy density has no Gauss rule; the
##   trapezoidal Z-form is used for every tau.
##
## Rows with tau = 0 get the exact conditional values F(c), f(c), f'(c).

## Gauss rule (normalised weights) from recurrence coefficients a (diagonal)
## and b (off-diagonal) of the Jacobi matrix.
gof_marginal_jacobi <- function(a,b) {
    Q <- length(a)
    J <- diag(a,Q)
    if (Q>1L) {
        J[cbind(seq_len(Q-1L),2:Q)] <- b[seq_len(Q-1L)]
        J[cbind(2:Q,seq_len(Q-1L))] <- b[seq_len(Q-1L)]
    }
    e <- eigen(J,symmetric=TRUE)
    o <- order(e$values)
    w <- e$vectors[1L,o]^2
    list(x=e$values[o],w=w/sum(w))
}

## Q-point Gauss rule for the discrete measure sum_i wt_i delta(x_i) by the
## Lanczos (discretised Stieltjes) recurrence, with a local
## reorthogonalisation against the current vector. For the smooth, densely
## sampled measures used here this agrees with a fully reorthogonalised
## Lanczos to ~1e-13 at a fraction of the cost.
gof_marginal_lanczos <- function(x,wt,Q) {
    q <- sqrt(wt); q <- q/sqrt(sum(q^2))
    qprev <- 0*q; bprev <- 0
    a <- b <- numeric(Q)
    for (j in seq_len(Q)) {
        v <- x*q
        a[j] <- sum(q*v)
        v <- v-a[j]*q-bprev*qprev
        v <- v-sum(v*q)*q
        b[j] <- sqrt(sum(v^2))
        qprev <- q; bprev <- b[j]
        q <- v/b[j]
    }
    gof_marginal_jacobi(a,b)
}

gof_marginal_cache <- new.env(parent=emptyenv())

## @noRd
## Q-point Gauss rule for the density f of the link variable e (logistic,
## Gumbel-min for cloglog, Gumbel-max for loglog). The measure is
## discretised by the trapezoidal rule (spectrally accurate for these smooth
## exponentially decaying densities) over a range that carries all moments
## up to degree 2Q-1; cached per link and Q.
gof_marginal_density_rule <- function(link,Q) {
    key <- paste(link,Q)
    if (!is.null(rule <- gof_marginal_cache[[key]])) return(rule)
    L <- 4*Q+60; h <- 0.1
    rule <- switch(link,
        logit={x <- seq(-L,L,by=h); gof_marginal_lanczos(x,stats::dlogis(x),Q)},
        loglog={x <- seq(-6,L,by=h); gof_marginal_lanczos(x,exp(-x-exp(-x)),Q)},
        cloglog={r <- gof_marginal_density_rule("loglog",Q); list(x=-rev(r$x),w=rev(r$w))},
        stop("No density rule for link ",link,call.=FALSE))
    assign(key,rule,envir=gof_marginal_cache)
    rule
}

## Largest tau handled by the trapezoidal Z-form, and the analyticity
## half-width D of F(x) around the real axis (|Im x| < D).
gof_marginal_switch <- c(logit=4,cloglog=2,loglog=2,cauchit=Inf)
gof_marginal_strip <- c(logit=pi,cloglog=pi/2,loglog=pi/2,cauchit=1)

## The integrators below take C (m x (K-1) values c) and tau (m) and return
## list(P, dc = dP/dc, dtau2 = dP/d(tau^2)), all m x (K-1).

gof_marginal_probit <- function(C,tau) {
    s <- sqrt(1+tau^2)
    U <- C/s
    D <- stats::dnorm(U)
    list(P=stats::pnorm(U),dc=D/s,dtau2=-U*D/(2*s^2))
}

## F, f and f' of the link at a matrix A with shared transcendental work
## (the Z-form's cost is dominated by these evaluations). Arguments of the
## Gumbel links are clamped to +-700, beyond which F, f and f' are already
## exactly 0 or 1 in double precision, so no 0 * Inf products arise.
gof_marginal_fdf <- function(link,A) {
    switch(link,
        logit={e <- exp(-abs(A)); r <- 1/(1+e); f <- e*r*r
               list(F=r*(e+(A>=0)*(1-e)),f=f,df=-sign(A)*f*(1-e)*r)},
        cloglog={A <- pmin(pmax(A,-700),700); ex <- exp(A); f <- ex*exp(-ex)
                 list(F=-expm1(-ex),f=f,df=f*(1-ex))},
        loglog={A <- pmin(pmax(A,-700),700); em <- exp(-A); f <- em*exp(-em)
                list(F=exp(-em),f=f,df=f*(em-1))},
        cauchit={d <- 1/(1+A*A)
                 list(F=stats::pcauchy(A),f=d/pi,df=-2*A*d*d/pi)},
        {lk <- gof_ordinal_link(link)
         list(F=lk$F(A),f=gof_marginal_tail0(lk$f(A)),df=gof_marginal_tail0(lk$df(A)))})
}

## Z-form: E G(c - tau Z) over a normal rule (x, w).
gof_marginal_zform <- function(link,C,tau,rule) {
    out <- list(P=C,dc=C,dtau2=C)
    for (k in seq_len(ncol(C))) {
        v <- gof_marginal_fdf(link,C[,k]-outer(tau,rule$x))
        out$P[,k] <- v$F %*% rule$w
        out$dc[,k] <- v$f %*% rule$w
        out$dtau2[,k] <- (v$df %*% rule$w)/2
    }
    out
}

## e-form: E Phi((c - e)/tau) over a rule (x, w) for the link density.
gof_marginal_eform <- function(C,tau,rule) {
    out <- list(P=C,dc=C,dtau2=C)
    for (k in seq_len(ncol(C))) {
        U <- (C[,k]-matrix(rule$x,nrow(C),length(rule$x),byrow=TRUE))/tau
        D <- stats::dnorm(U)
        out$P[,k] <- stats::pnorm(U) %*% rule$w
        out$dc[,k] <- (D %*% rule$w)/tau
        out$dtau2[,k] <- -((U*D) %*% rule$w)/(2*tau^2)
    }
    out
}

## Trapezoidal rule for the standard normal on [-10, 10] with N points.
gof_marginal_trapezoid <- function(N) {
    x <- seq(-10,10,length.out=N)
    w <- stats::dnorm(x)
    list(x=x,w=w/sum(w))
}

## Non-finite tail products (0 * Inf far in a tail) mapped to their limit 0.
gof_marginal_tail0 <- function(v) {
    v[!is.finite(v)] <- 0
    v
}

## @noRd
## E F(c - tau Z) with its c and tau^2 derivatives for every row of C.
## `nodes` is the minimum number of trapezoid points and of e-form Gauss
## nodes (at least 30 are used for the e-form).
gof_marginal_integrals <- function(link,C,tau,nodes=30L) {
    Q <- gof_count(nodes,"nodes",2L)
    lk <- gof_ordinal_link(link)
    out <- list(P=C,dc=C,dtau2=C)
    put <- function(rows,v) {
        out$P[rows,] <<- v$P; out$dc[rows,] <<- v$dc; out$dtau2[rows,] <<- v$dtau2
    }
    sub <- function(rows) C[rows,,drop=FALSE]
    zero <- tau==0
    if (link=="probit") {
        put(!zero,gof_marginal_probit(sub(!zero),tau[!zero]))
    } else {
        high <- tau>=gof_marginal_switch[[link]]
        if (any(high))
            put(high,gof_marginal_eform(sub(high),tau[high],gof_marginal_density_rule(link,max(Q,30L))))
        ## Trapezoid: h = 20/(N-1) <= min(1/2, 0.18 D/tau), N = 40*2^j+1 >= nodes.
        low <- !zero & !high
        need <- 20/pmin(0.5,0.18*gof_marginal_strip[[link]]/tau[low])
        j <- pmax(ceiling(log2(need/40)),ceiling(log2((Q-1)/40)),0)
        N <- 40*2^j+1
        for (g in unique(N)) {
            rows <- which(low)[N==g]
            put(rows,gof_marginal_zform(link,sub(rows),tau[rows],gof_marginal_trapezoid(g)))
        }
    }
    if (any(zero)) {
        C0 <- sub(zero)
        put(zero,list(P=lk$F(C0),dc=gof_marginal_tail0(lk$f(C0)),
                      dtau2=gof_marginal_tail0(lk$df(C0))/2))
    }
    out$P <- pmin(pmax(out$P,0),1)
    out
}

## @noRd
## Marginal cumulative probabilities and their theta-gradient.
##
## @param problem A `gof_clmm_problem`.
## @param theta Parameter vector c(xi, beta, omega) in the problem layout.
## @param nodes Minimum number of quadrature nodes for the one-dimensional
##   integrals (see gof_marginal_integrals; probit is exact).
## @returns list(cdf = n x (K-1) matrix of P_ik = P(Y_i <= k);
##   gradient = (n(K-1)) x d matrix whose row (k-1)n + i is the derivative of
##   the residual 1{Y_i <= k} - P_ik, i.e. -dP_ik/dtheta (the convention of
##   cumres.ordreg); tau2 = n-vector of tau_i^2). Observation weights are not
##   applied.
gof_clmm_marginal <- function(problem,theta=problem$theta,nodes=30L) {
    if (!inherits(problem,"gof_clmm_problem")) stop("problem must be a gof_clmm_problem",call.=FALSE)
    if (!is.numeric(theta) || any(!is.finite(theta))) stop("theta must be finite",call.=FALSE)
    par <- gof_clmm_unpack(problem,theta)
    n <- problem$n; K1 <- problem$K-1L; d <- length(problem$theta)
    eta <- problem$offset
    if (ncol(problem$X)) eta <- eta+drop(problem$X %*% par$beta)

    ## tau_i^2 = sum_t |L_t' z_it|^2 and d tau_i^2 / d omega = z' dSigma z.
    omega_cols <- unlist(lapply(problem$terms,`[[`,"omega"))
    tau2 <- numeric(n)
    dtau2 <- matrix(0,n,length(omega_cols))
    pos <- 0L
    for (t in seq_along(problem$terms)) {
        Z <- problem$terms[[t]]$Z
        tau2 <- tau2+rowSums((Z %*% par$re[[t]]$L)^2)
        for (dS in par$re[[t]]$dSigma) {
            pos <- pos+1L
            dtau2[,pos] <- rowSums((Z %*% dS)*Z)
        }
    }
    C <- outer(-eta,par$alpha,`+`)
    I <- gof_marginal_integrals(problem$link,C,sqrt(tau2),nodes)

    ## Residual derivative = -dP/dtheta, stacked by threshold.
    xi_cols <- problem$index$xi; beta_cols <- problem$index$beta
    gradient <- matrix(0,n*K1,d,dimnames=list(NULL,names(problem$theta)))
    for (k in seq_len(K1)) {
        rows <- (k-1L)*n+seq_len(n)
        gradient[rows,xi_cols] <- -outer(I$dc[,k],problem$tJac[k,])
        if (length(beta_cols)) gradient[rows,beta_cols] <- problem$X*I$dc[,k]
        if (length(omega_cols)) gradient[rows,omega_cols] <- -I$dtau2[,k]*dtau2
    }
    dimnames(I$P) <- NULL
    list(cdf=I$P,gradient=gradient,tau2=tau2)
}
