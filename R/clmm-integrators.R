## Marginal-likelihood integrators for cumulative link mixed models.
##
## 1.4.0 ships one integrator: nested adaptive Gauss--Hermite quadrature
## (Rabe-Hesketh, Skrondal & Pickles, 2005) for a single grouping factor or a
## nesting chain. Every integrator implements the `gof_integrator()` contract
## from clmm-contracts.R, so a future crossed-effects integrator (e.g. a sparse
## Laplace) can be registered in `gof_integrator_registry` without touching the
## polish step, the marginal layer or the engine.
##
## Estimator definition. The proposal ("state") holds, for every random-effect
## unit, a mode and a Cholesky scale computed at some theta. With the state
## frozen, the per-block log-likelihood
##   l_c(theta; nu) = log sum_q W_q f(y_c | b_q; theta) phi(b_q; Sigma(theta)) / psi_nu(b_q)
## (nested appropriately) is a smooth function of theta, its gradient is the
## exact posterior-weighted complete-data score (Fisher/Louis identity), and
## the polished estimate is the root of the total score with nodes frozen at
## nu(theta.hat).

## ---------------------------------------------------------------------------
## Gauss--Hermite rule for the standard normal measure
## ---------------------------------------------------------------------------

## @noRd
## Gauss--Hermite nodes and weights for the standard normal (probabilists')
## measure: sum(w * g(x)) = E g(Z), Z ~ N(0,1), exactly for polynomials of
## degree <= 2Q-1. Nodes are the eigenvalues of the Jacobi matrix
## (Golub--Welsch), polished by Newton steps on the orthonormal Hermite
## polynomial; weights use the Christoffel function 1/sum_k p_k(x)^2, which
## keeps tail weights accurate in relative terms.
gof_gauss_hermite <- function(Q) {
    Q <- gof_count(Q,"Q",1L)
    if (Q==1L) return(list(nodes=0,weights=1))
    off <- sqrt(seq_len(Q-1L))
    J <- matrix(0,Q,Q)
    J[cbind(seq_len(Q-1L),2:Q)] <- off
    J[cbind(2:Q,seq_len(Q-1L))] <- off
    x <- sort(eigen(J,symmetric=TRUE,only.values=TRUE)$values)
    ## Orthonormal probabilists' Hermite polynomials p_0..p_Q at x.
    orth <- function(x) {
        P <- matrix(0,length(x),Q+1L)
        P[,1L] <- 1
        P[,2L] <- x
        if (Q>=2L) for (k in 2:Q) P[,k+1L] <- (x*P[,k]-sqrt(k-1)*P[,k-1L])/sqrt(k)
        P
    }
    for (it in 1:3) {
        P <- orth(x)
        x <- x-P[,Q+1L]/(sqrt(Q)*P[,Q])          # p_Q' = sqrt(Q) p_{Q-1}
    }
    x <- (x-rev(x))/2                            # exact symmetry
    P <- orth(x)
    w <- 1/rowSums(P[,seq_len(Q),drop=FALSE]^2)
    w <- (w+rev(w))/2
    list(nodes=x,weights=w/sum(w))
}

## Tensor-product grid in q dimensions: nodes (Q^q x q) and log weights.
gof_gh_grid <- function(Q,q) {
    gh <- gof_gauss_hermite(Q)
    idx <- as.matrix(expand.grid(rep(list(seq_len(Q)),q)))
    z <- matrix(gh$nodes[idx],ncol=q)
    lw <- rowSums(matrix(log(gh$weights[idx]),ncol=q))
    list(z=z,logw=lw,Q=Q,q=q,n=nrow(z))
}

## ---------------------------------------------------------------------------
## Row kernel: log category probabilities and eta-derivatives
## ---------------------------------------------------------------------------

## The contract link table extended with log f and d log f / dx, so that
## density/probability ratios are formed on the log scale.
gof_clmm_link <- function(name) {
    lk <- gof_ordinal_link(name)
    lk$logf <- switch(name,
        logit=function(x) stats::dlogis(x,log=TRUE),
        probit=function(x) stats::dnorm(x,log=TRUE),
        cloglog=function(x) x-exp(x),
        loglog=function(x) -x-exp(-x),
        cauchit=function(x) stats::dcauchy(x,log=TRUE))
    lk$dlogf <- switch(name,
        logit=function(x) -tanh(x/2),
        probit=function(x) -x,
        cloglog=function(x) -expm1(x),
        loglog=function(x) expm1(-x),
        cauchit=function(x) -2*x/(1+x^2))
    lk
}

## log(1 - exp(-a)) for a > 0 (Maechler 2012); -Inf for a <= 0.
gof_log1mexp <- function(a) {
    out <- rep(-Inf,length(a))
    big <- which(a>log(2))
    out[big] <- log1p(-exp(-a[big]))
    small <- which(a>0 & a<=log(2))
    out[small] <- log(-expm1(-a[small]))
    out[is.na(a)] <- NaN
    out
}

## @noRd
## For linear predictors eta (P(Y <= k) = F(alpha_k - eta)) and categories y:
## lp = log P(Y = y), rhi = f(hi)/p, rlo = f(lo)/p, d1 = d lp / d eta and,
## for deriv >= 2, d2 = d^2 lp / d eta^2. The interval probability is formed
## on the log scale from the lower tail (midpoint <= 0) or the upper tail.
gof_clmm_kernel <- function(eta,y,alpha,K,lk,deriv=1L) {
    ## alpha_0 = -Inf and alpha_K = Inf fold the end categories into the two
    ## tail formulas: log1mexp(Inf) = 0 gives logF(hi) and logS(lo) exactly.
    top <- y==K
    bot <- y==1L
    hi <- c(alpha,Inf)[y]-eta
    lo <- c(-Inf,alpha)[y]-eta
    lp <- numeric(length(eta))
    lower <- (hi+lo)<=0
    if (any(lower)) {
        a <- lk$logF(hi[lower]); b <- lk$logF(lo[lower])
        lp[lower] <- ifelse(a==-Inf,-Inf,a+gof_log1mexp(a-b))
    }
    upper <- !lower
    if (any(upper)) {
        a <- lk$logS(lo[upper]); b <- lk$logS(hi[upper])
        lp[upper] <- ifelse(a==-Inf,-Inf,a+gof_log1mexp(a-b))
    }
    out <- list(lp=lp)
    if (deriv<1L) return(out)
    rhi <- exp(lk$logf(hi)-lp)
    rlo <- exp(lk$logf(lo)-lp)
    rhi[top] <- 0
    rlo[bot] <- 0
    ## Zero-probability rows carry no posterior weight; keep ratios finite.
    zero <- lp==-Inf
    rhi[zero] <- 0
    rlo[zero] <- 0
    out$rhi <- rhi
    out$rlo <- rlo
    out$d1 <- rlo-rhi
    if (deriv>=2L) {
        ## r * dlogf(x) with r = 0 where the density underflows (0 * Inf -> 0).
        prod0 <- function(r,x) {v <- r*lk$dlogf(x); v[r==0] <- 0; v}
        d2 <- -out$d1^2
        d2 <- d2+prod0(rhi,hi)-prod0(rlo,lo)
        out$d2 <- d2
    }
    out
}

## Expected information for eta: sum_k p_k (d log p_k / d eta)^2.
gof_clmm_expected_eta <- function(eta,alpha,K,lk) {
    E <- numeric(length(eta))
    for (k in seq_len(K)) {
        kk <- gof_clmm_kernel(eta,rep.int(k,length(eta)),alpha,K,lk,1L)
        E <- E+exp(kk$lp)*kk$d1^2
    }
    E
}

## ---------------------------------------------------------------------------
## Batched small-matrix algebra on arrays [n, a, b]
## ---------------------------------------------------------------------------

gof_bchol <- function(A) {
    n <- dim(A)[1L]; s <- dim(A)[2L]
    L <- array(0,dim(A))
    ok <- rep(TRUE,n)
    for (j in seq_len(s)) {
        d <- A[,j,j]
        if (j>1L) for (k in seq_len(j-1L)) d <- d-L[,j,k]^2
        bad <- !is.finite(d) | !(d>1e-14*pmax(abs(A[,j,j]),.Machine$double.xmin))
        ok[bad] <- FALSE
        d[bad] <- 1
        L[,j,j] <- sqrt(d)
        if (j<s) for (i in (j+1L):s) {
            x <- A[,i,j]
            if (j>1L) for (k in seq_len(j-1L)) x <- x-L[,i,k]*L[,j,k]
            L[,i,j] <- x/L[,j,j]
        }
    }
    list(L=L,ok=ok)
}

## Solve (L L') X = B for B [n, s, k].
gof_bcholsolve <- function(L,B) {
    s <- dim(L)[2L]
    Y <- B
    for (i in seq_len(s)) {
        x <- Y[,i,,drop=FALSE]
        if (i>1L) for (j in seq_len(i-1L)) x <- x-L[,i,j]*Y[,j,,drop=FALSE]
        Y[,i,] <- x/L[,i,i]
    }
    X <- Y
    for (i in rev(seq_len(s))) {
        x <- X[,i,,drop=FALSE]
        if (i<s) for (j in (i+1L):s) x <- x-L[,j,i]*X[,j,,drop=FALSE]
        X[,i,] <- x/L[,i,i]
    }
    X
}

## t(A) %*% B for A [n, a, b], B [n, a, c] -> [n, b, c].
gof_bcrossprod <- function(A,B) {
    n <- dim(A)[1L]; a <- dim(A)[2L]; nb <- dim(A)[3L]; nc <- dim(B)[3L]
    out <- array(0,c(n,nb,nc))
    for (i in seq_len(nb)) for (k in seq_len(nc)) {
        s <- 0
        for (j in seq_len(a)) s <- s+A[,j,i]*B[,j,k]
        out[,i,k] <- s
    }
    out
}

## Symmetric inverse square roots C = M^{-1/2} of M = L L' for L [P, q, q];
## returns C [P, q, q] and log det C.
gof_inv_sqrt_batch <- function(L) {
    P <- dim(L)[1L]; q <- dim(L)[2L]
    if (q==1L) return(list(C=array(1/L[,1L,1L],c(P,1L,1L)),logdet=-log(L[,1L,1L])))
    C <- array(0,c(P,q,q))
    logdet <- numeric(P)
    for (p in seq_len(P)) {
        Lp <- matrix(L[p,,],q)
        e <- eigen(tcrossprod(Lp),symmetric=TRUE)
        C[p,,] <- e$vectors %*% (t(e$vectors)/sqrt(e$values))
        logdet[p] <- -sum(log(diag(Lp)))
    }
    list(C=C,logdet=logdet)
}

## ---------------------------------------------------------------------------
## Nesting-chain tree
## ---------------------------------------------------------------------------

## Levels of a nesting chain: per level the stacked design of all terms that
## share the grouping factor (block-diagonal Sigma), unit indices, parent
## units and row lists. Level 1 units are the independent blocks.
gof_agq_tree <- function(problem) {
    if (!identical(problem$structure,"nesting-chain"))
        stop("Nested adaptive quadrature needs a nesting-chain structure",call.=FALSE)
    n <- problem$n
    chain <- problem$chain
    D <- length(chain)
    groups <- vapply(problem$terms,`[[`,"","group")
    lev <- vector("list",D)
    for (j in seq_len(D)) {
        tt <- which(groups==chain[j])
        if (!length(tt)) stop("Nesting chain level '",chain[j],"' has no random-effect term",call.=FALSE)
        idx <- problem$terms[[tt[1L]]]$index
        nu <- problem$terms[[tt[1L]]]$n_levels
        for (t in tt[-1L]) {
            other <- problem$terms[[t]]$index
            if (problem$terms[[t]]$n_levels!=nu || nrow(unique(cbind(idx,other)))!=nu)
                stop("Terms sharing grouping factor '",chain[j],"' must index the same levels",call.=FALSE)
        }
        qs <- vapply(problem$terms[tt],`[[`,0L,"q")
        Z <- do.call(cbind,lapply(problem$terms[tt],`[[`,"Z"))
        lev[[j]] <- list(group=chain[j],terms=tt,q=sum(qs),Z=Z,unit=idx,n_units=nu,
                         cols=split(seq_len(sum(qs)),rep(seq_along(tt),qs)))
    }
    ## Top level <-> independent blocks (one-to-one).
    bid <- problem$blocks$id
    m <- problem$blocks$m
    ub <- unique(cbind(lev[[1L]]$unit,bid))
    if (nrow(ub)!=m || lev[[1L]]$n_units!=m || anyDuplicated(ub[,1L]) || anyDuplicated(ub[,2L]))
        stop("The top-level grouping factor '",chain[1L],"' must coincide with the independent blocks",call.=FALSE)
    lev[[1L]]$unit <- bid
    for (j in seq_len(D)) {
        u <- lev[[j]]$unit
        if (j>1L) {
            parent <- integer(lev[[j]]$n_units)
            parent[u] <- lev[[j-1L]]$unit
            if (any(parent[u]!=lev[[j-1L]]$unit))
                stop("Grouping factor '",chain[j],"' is not nested in '",chain[j-1L],"'",call.=FALSE)
            lev[[j]]$parent <- parent
        }
        rows <- split(seq_len(n),factor(u,levels=seq_len(lev[[j]]$n_units)))
        pos <- integer(n)
        pos[unlist(rows,use.names=FALSE)] <- sequence(lengths(rows))
        lev[[j]]$rows <- unname(rows)
        lev[[j]]$pos <- pos
    }
    for (j in seq_len(D-1L))
        lev[[j]]$children <- unname(split(seq_len(lev[[j+1L]]$n_units),
                                          factor(lev[[j+1L]]$parent,levels=seq_len(lev[[j]]$n_units))))
    list(levels=lev,D=D,n=n,m=m,block=bid)
}

## Level priors at theta: Sigma, inverse and a log-density over effect rows.
gof_agq_priors <- function(tree,up) {
    lapply(tree$levels,function(L) {
        q <- L$q
        Sinv <- matrix(0,q,q)
        chol_t <- list()
        for (a in seq_along(L$terms)) {
            re <- up$re[[L$terms[a]]]
            cc <- L$cols[[a]]
            Sinv[cc,cc] <- chol2inv(t(re$L))
            chol_t[[a]] <- re$L
        }
        logphi <- function(B) {
            out <- 0
            for (a in seq_along(L$terms)) {
                Lt <- chol_t[[a]]
                cc <- L$cols[[a]]
                U <- forwardsolve(Lt,t(B[,cc,drop=FALSE]))
                out <- out-length(cc)/2*log(2*pi)-sum(log(diag(Lt)))-colSums(U^2)/2
            }
            out
        }
        list(Sinv=Sinv,logphi=logphi)
    })
}

## ---------------------------------------------------------------------------
## Vectorised joint-mode search over subtrees
## ---------------------------------------------------------------------------

## @noRd
## Joint mode of the log integrand over the effects of P independent
## subtrees rooted at level j, given fixed ancestor contributions `rc_off`.
## Row copies (rc_row, rc_prob) list the rows of every problem's root unit,
## ordered by problem. Newton steps solve the tree-structured (block-arrow)
## system by eliminating finest levels first (Schur complements); failed
## Cholesky factors (non-PD observed Hessian, e.g. cauchit) switch that
## problem to the expected information. Returns the root modes, the Cholesky
## factor of the root Schur complement, all local modes (for warm starts)
## and iteration counts.
gof_agq_subtree <- function(ctx,j,rc_row,rc_prob,rc_off,P,start=NULL,maxit=100L,tol=1e-20) {
    tree <- ctx$tree
    D <- tree$D
    lev <- j:D
    nl <- length(lev)
    lu <- lu_prob <- lu_unit <- lu_parent <- ids <- vector("list",D)
    for (l in lev) {
        nu <- tree$levels[[l]]$n_units
        key <- (rc_prob-1)*nu+tree$levels[[l]]$unit[rc_row]
        ids[[l]] <- sort(unique(key))
        lu[[l]] <- match(key,ids[[l]])
        lu_prob[[l]] <- as.integer((ids[[l]]-1)%/%nu+1)
        lu_unit[[l]] <- as.integer((ids[[l]]-1)%%nu+1)
        if (l>j) lu_parent[[l]] <- lu[[l-1L]][match(seq_along(ids[[l]]),lu[[l]])]
    }
    q <- vapply(tree$levels,`[[`,0L,"q")
    Zl <- lapply(seq_len(D),function(l) if (l>=j) tree$levels[[l]]$Z[rc_row,,drop=FALSE])
    Zt <- do.call(cbind,Zl[lev])
    s_end <- cumsum(q[lev])
    s_start <- s_end-q[lev]+1L
    sD <- ncol(Zt)
    ij <- which(upper.tri(diag(sD),diag=TRUE),arr.ind=TRUE)
    w <- ctx$w[rc_row]
    y <- ctx$y[rc_row]
    base <- ctx$xb[rc_row]+rc_off
    pri <- ctx$prior
    b <- vector("list",D)
    for (l in lev) {
        b[[l]] <- matrix(0,length(ids[[l]]),q[l])
        if (!is.null(start) && !is.null(start[[l]])) {
            s0 <- start[[l]](lu_prob[[l]],lu_unit[[l]])
            okr <- !is.na(s0[,1L])
            b[[l]][okr,] <- s0[okr,]
        }
    }
    eta_of <- function(b) {
        e <- base
        for (l in lev) e <- e+rowSums(Zl[[l]]*b[[l]][lu[[l]],,drop=FALSE])
        e
    }
    obj_of <- function(b,lp) {
        o <- drop(rowsum(w*lp,rc_prob,reorder=TRUE))
        for (l in lev) o <- o+drop(rowsum(pri[[l]]$logphi(b[[l]]),lu_prob[[l]],reorder=TRUE))
        o
    }
    ## One elimination pass with row curvature cc; returns step, root factor.
    eliminate <- function(b,cc,d1) {
        wc <- w*cc
        wd <- w*d1
        okp <- rep(TRUE,P)
        Dl <- tail(lev,1L)
        nD <- length(ids[[Dl]])
        Mf <- rowsum(Zt[,ij[,1L],drop=FALSE]*Zt[,ij[,2L],drop=FALSE]*wc,lu[[Dl]],reorder=TRUE)
        M <- array(0,c(nD,sD,sD))
        for (r in seq_len(nrow(ij))) {
            M[,ij[r,1L],ij[r,2L]] <- Mf[,r]
            M[,ij[r,2L],ij[r,1L]] <- Mf[,r]
        }
        V <- rowsum(Zt*wd,lu[[Dl]],reorder=TRUE)
        vorig <- vector("list",D)
        for (l in lev) vorig[[l]] <- rowsum(Zl[[l]]*wd,lu[[l]],reorder=TRUE)-b[[l]] %*% pri[[l]]$Sinv
        addprior <- function(M,V,l) {
            self <- s_start[l-j+1L]:s_end[l-j+1L]
            Si <- pri[[l]]$Sinv
            for (a in seq_along(self)) for (c in seq_along(self))
                M[,self[a],self[c]] <- M[,self[a],self[c]]+Si[a,c]
            V[,self] <- V[,self,drop=FALSE]-b[[l]] %*% Si
            list(M=M,V=V)
        }
        mv <- addprior(M,V,Dl); M <- mv$M; V <- mv$V
        store <- vector("list",D)
        if (nl>1L) for (l in rev(lev)[-nl]) {
            k <- l-j+1L
            self <- s_start[k]:s_end[k]
            A <- seq_len(s_start[k]-1L)
            nn <- dim(M)[1L]
            ch <- gof_bchol(M[,self,self,drop=FALSE])
            okp[lu_prob[[l]][!ch$ok]] <- FALSE
            MuA <- M[,self,A,drop=FALSE]
            XA <- gof_bcholsolve(ch$L,MuA)
            Xv <- gof_bcholsolve(ch$L,array(V[,self],c(nn,length(self),1L)))
            R <- M[,A,A,drop=FALSE]-gof_bcrossprod(MuA,XA)
            r <- V[,A,drop=FALSE]-matrix(gof_bcrossprod(MuA,Xv),nn)
            store[[l]] <- list(XA=XA,Xv=Xv)
            na <- length(A)
            np <- length(ids[[l-1L]])
            M <- array(rowsum(matrix(R,nn),lu_parent[[l]],reorder=TRUE),c(np,na,na))
            V <- rowsum(r,lu_parent[[l]],reorder=TRUE)
            mv <- addprior(M,V,l-1L); M <- mv$M; V <- mv$V
        }
        ch <- gof_bchol(M)
        okp[!ch$ok] <- FALSE
        delta <- vector("list",D)
        delta[[j]] <- matrix(gof_bcholsolve(ch$L,array(V,c(P,q[j],1L))),P)
        if (nl>1L) for (l in lev[-1L]) {
            anc <- vector("list",D)
            anc[[l-1L]] <- lu_parent[[l]]
            if (l-2L>=j) for (a in (l-2L):j) anc[[a]] <- lu_parent[[a+1L]][anc[[a+1L]]]
            DA <- do.call(cbind,lapply(j:(l-1L),function(a) delta[[a]][anc[[a]],,drop=FALSE]))
            nn <- nrow(DA)
            st <- store[[l]]
            delta[[l]] <- matrix(st$Xv,nn)-matrix(gof_bcrossprod(aperm(st$XA,c(1L,3L,2L)),
                                                                 array(DA,c(nn,ncol(DA),1L))),nn)
        }
        dec <- numeric(P)
        for (l in lev) dec <- dec+drop(rowsum(rowSums(delta[[l]]*vorig[[l]]),lu_prob[[l]],reorder=TRUE))
        list(delta=delta,ok=okp,dec=dec,L=ch$L)
    }
    newton_dir <- function(b,eta,ker) {
        el <- eliminate(b,-ker$d2,ker$d1)
        fb <- !el$ok
        if (any(fb)) {
            use <- fb[rc_prob]
            cc <- -ker$d2
            cc[use] <- gof_clmm_expected_eta(eta[use],ctx$alpha,ctx$K,ctx$lk)
            el2 <- eliminate(b,cc,ker$d1)
            for (l in lev) {
                r <- fb[lu_prob[[l]]]
                el$delta[[l]][r,] <- el2$delta[[l]][r,]
            }
            el$dec[fb] <- el2$dec[fb]
            el$L[fb,,] <- el2$L[fb,,]
            if (any(!el2$ok[fb])) stop("Expected information is not positive definite in the mode search",call.=FALSE)
        }
        el$fallback <- fb
        el
    }
    conv <- rep(FALSE,P)
    iters <- integer(P)
    halvings <- 0L
    ever <- rep(FALSE,P)
    eta <- eta_of(b)
    ker <- gof_clmm_kernel(eta,y,ctx$alpha,ctx$K,ctx$lk,2L)
    obj <- obj_of(b,ker$lp)
    if (any(!is.finite(obj))) stop("Non-finite log integrand at the start of the mode search",call.=FALSE)
    for (it in 0:maxit) {
        el <- newton_dir(b,eta,ker)
        ever <- ever | el$fallback
        conv <- conv | (is.finite(el$dec) & el$dec<tol)
        if (all(conv) || it==maxit) break
        iters[!conv] <- it+1L
        t <- ifelse(conv,0,1)
        for (h in 0:60) {
            bn <- b
            for (l in lev) bn[[l]] <- b[[l]]+t[lu_prob[[l]]]*el$delta[[l]]
            etan <- eta_of(bn)
            kn <- gof_clmm_kernel(etan,y,ctx$alpha,ctx$K,ctx$lk,2L)
            objn <- obj_of(bn,kn$lp)
            bad <- t>0 & !(is.finite(objn) & objn>=obj-1e-11*(1+abs(obj)))
            if (!any(bad)) break
            halvings <- halvings+1L
            t[bad] <- if (h<60L) t[bad]/2 else 0
        }
        if (any(bad)) {
            ## Give up moving those problems; restore their values exactly.
            for (l in lev) bn[[l]][bad[lu_prob[[l]]],] <- b[[l]][bad[lu_prob[[l]]],]
            etan <- eta_of(bn)
            kn <- gof_clmm_kernel(etan,y,ctx$alpha,ctx$K,ctx$lk,2L)
            objn <- obj_of(bn,kn$lp)
            conv[bad] <- TRUE
        }
        b <- bn; eta <- etan; ker <- kn; obj <- objn
    }
    list(mu=b[[j]],L=el$L,b=b,ids=ids,lu_prob=lu_prob,lu_unit=lu_unit,
         iterations=iters,converged=is.finite(el$dec) & el$dec<tol,
         fallback=el$fallback,fallback_any=ever,halvings=halvings)
}

## ---------------------------------------------------------------------------
## Nested adaptive Gauss--Hermite integrator
## ---------------------------------------------------------------------------

## @noRd
## Nested adaptive Gauss--Hermite quadrature for a nesting chain of any
## depth. `nodes`: Gauss--Hermite nodes per dimension, one value per chain
## level (recycled); NULL uses `gof_clmm_default_nodes()` per level.
gof_integrator_nested_agq <- function(problem,nodes=NULL) {
    if (!inherits(problem,"gof_clmm_problem")) stop("problem must be a gof_clmm_problem",call.=FALSE)
    tree <- gof_agq_tree(problem)
    D <- tree$D
    qv <- vapply(tree$levels,`[[`,0L,"q")
    if (is.null(nodes)) nodes <- vapply(qv,gof_clmm_default_nodes,0L,control=list(nodes=NULL))
    nodes <- unname(vapply(nodes,gof_count,0L,name="nodes",minimum=1L))
    if (length(nodes)==1L) nodes <- rep(nodes,D)
    if (length(nodes)!=D) stop("nodes needs one value per nesting level (",D,")",call.=FALSE)
    grids <- lapply(seq_len(D),function(l) gof_gh_grid(nodes[l],qv[l]))
    lk <- gof_clmm_link(problem$link)
    K <- problem$K
    y <- problem$y
    w <- problem$weights
    X <- problem$X
    off <- problem$offset
    blk <- problem$blocks$id
    m <- problem$blocks$m
    dnames <- names(problem$theta)
    d <- length(problem$theta)

    context <- function(theta) {
        up <- gof_clmm_unpack(problem,theta)
        xb <- off
        if (ncol(X)) xb <- xb+drop(X %*% up$beta)
        list(tree=tree,up=up,alpha=up$alpha,K=K,lk=lk,y=y,w=w,xb=xb,prior=gof_agq_priors(tree,up))
    }

    refresh <- function(theta) {
        theta <- stats::setNames(as.numeric(theta),dnames)
        ctx <- context(theta)
        levels <- vector("list",D)
        rc_row <- unlist(tree$levels[[1L]]$rows,use.names=FALSE)
        rc_prob <- rep.int(seq_len(m),lengths(tree$levels[[1L]]$rows))
        rc_off <- numeric(length(rc_row))
        prob_unit <- seq_len(m)
        prob_parent <- rep(NA_integer_,m)
        prob_block <- seq_len(m)
        start <- NULL
        for (j in seq_len(D)) {
            P <- length(prob_unit)
            sol <- gof_agq_subtree(ctx,j,rc_row,rc_prob,rc_off,P,start=start)
            g <- grids[[j]]
            q <- qv[j]
            Q <- g$n
            ## Nodes b = mu + C z with C the symmetric inverse square root of
            ## the root Schur complement M = L L' (C C' = M^{-1}). The
            ## symmetric root keeps the node set equivariant under reordering
            ## the stacked effects (e.g. term order), unlike a Cholesky root.
            sq <- gof_inv_sqrt_batch(sol$L)
            C <- sq$C
            logdetC <- sq$logdet
            inst_b <- matrix(0,P*Q,q)
            pid <- rep(seq_len(P),each=Q)
            zz <- g$z[rep(seq_len(Q),P),,drop=FALSE]
            for (a in seq_len(q)) {
                v <- sol$mu[pid,a]
                for (k in seq_len(q)) v <- v+C[pid,a,k]*zz[,k]
                inst_b[,a] <- v
            }
            inst_lconst <- g$logw[rep(seq_len(Q),P)]+q/2*log(2*pi)+rowSums(zz^2)/2+logdetC[pid]
            ## Pairs (instance, row) with cumulative random linear predictor.
            r_p <- tabulate(rc_prob,P)
            s_p <- cumsum(c(1L,r_p))[seq_len(P)]
            pos <- sequence(rep(r_p,each=Q),from=rep(s_p,each=Q))
            pr_row <- rc_row[pos]
            pr_inst <- rep.int(seq_len(P*Q),rep(r_p,each=Q))
            pr_zb <- rc_off[pos]+rowSums(tree$levels[[j]]$Z[pr_row,,drop=FALSE]*inst_b[pr_inst,,drop=FALSE])
            inst_block <- rep(prob_block,each=Q)
            levels[[j]] <- list(group=tree$levels[[j]]$group,q=q,nodes=nodes[j],Q=Q,P=P,
                                prob_unit=prob_unit,prob_parent=prob_parent,
                                mu=sol$mu,L=sol$L,inst_b=inst_b,inst_lconst=inst_lconst,
                                inst_block=inst_block,iterations=sol$iterations,
                                converged=sol$converged,fallback=sol$fallback,
                                fallback_any=sol$fallback_any,
                                halvings=sol$halvings)
            if (j==D) {
                pairs <- list(row=pr_row,inst=pr_inst,zb=pr_zb)
                break
            }
            ## Next level: one problem per (instance, child unit).
            inst_unit <- rep(prob_unit,each=Q)
            ch <- tree$levels[[j]]$children[inst_unit]
            nch <- lengths(ch)
            new_parent <- rep.int(seq_len(P*Q),nch)
            new_unit <- unlist(ch,use.names=FALSE)
            rows_new <- tree$levels[[j+1L]]$rows[new_unit]
            new_rc_row <- unlist(rows_new,use.names=FALSE)
            new_rc_prob <- rep.int(seq_along(new_unit),lengths(rows_new))
            pair_start <- cumsum(c(1L,rep(r_p,each=Q)))[seq_len(P*Q)]
            rc_off <- pr_zb[pair_start[new_parent[new_rc_prob]]+tree$levels[[j]]$pos[new_rc_row]-1L]
            ## Warm start from the joint mode of the parent subtree.
            oldp <- ceiling(new_parent/Q)
            start <- vector("list",D)
            for (l in (j+1L):D) local({
                ids_l <- sol$ids[[l]]; b_l <- sol$b[[l]]; nu <- tree$levels[[l]]$n_units
                op <- oldp
                start[[l]] <<- function(lp,lu) b_l[match((op[lp]-1)*nu+lu,ids_l),,drop=FALSE]
            })
            rc_row <- new_rc_row
            rc_prob <- new_rc_prob
            prob_block <- inst_block[new_parent]
            prob_parent <- new_parent
            prob_unit <- new_unit
        }
        structure(list(theta=theta,nodes=nodes,levels=levels,pairs=pairs,
                       pair_block=blk[pairs$row]),class="gof_agq_state")
    }

    ## Forward pass (log-likelihood) and, optionally, the Louis score.
    evaluate <- function(theta,state,score=FALSE) {
        if (!inherits(state,"gof_agq_state")) stop("state must come from refresh()",call.=FALSE)
        ctx <- context(theta)
        pr <- state$pairs
        eta <- ctx$xb[pr$row]+pr$zb
        ker <- gof_clmm_kernel(eta,y[pr$row],ctx$alpha,K,lk,if (score) 1L else 0L)
        S <- drop(rowsum(w[pr$row]*ker$lp,pr$inst,reorder=TRUE))
        pis <- vector("list",D)
        for (l in rev(seq_len(D))) {
            L <- state$levels[[l]]
            val <- S+ctx$prior[[l]]$logphi(L$inst_b)+L$inst_lconst
            Mx <- matrix(val,L$Q)
            mx <- Mx[1L,]
            if (L$Q>1L) for (r in 2:L$Q) mx <- pmax(mx,Mx[r,])
            mx0 <- ifelse(is.finite(mx),mx,0)
            E <- exp(Mx-rep(mx0,each=L$Q))
            g <- mx0+log(colSums(E))
            if (score) pis[[l]] <- as.vector(E/rep(colSums(E),each=L$Q))
            if (l>1L) S <- drop(rowsum(g,L$prob_parent,reorder=TRUE))
        }
        ll <- g
        if (!score) return(ll)
        ## Joint posterior weights of instances: products of conditional ones.
        Pw <- vector("list",D)
        Pw[[1L]] <- pis[[1L]]
        if (D>1L) for (l in 2:D) {
            L <- state$levels[[l]]
            Pw[[l]] <- Pw[[l-1L]][rep(L$prob_parent,each=L$Q)]*pis[[l]]
        }
        pp <- Pw[[D]][pr$inst]*w[pr$row]
        G <- rowsum(cbind(pp*ker$rhi,pp*ker$rlo),pr$row,reorder=TRUE)   # n x 2
        Ghi <- G[,1L]; Glo <- G[,2L]
        A <- matrix(0,length(y),K-1L)
        hi_ok <- y<K
        lo_ok <- y>1L
        A[cbind(which(hi_ok),y[hi_ok])] <- Ghi[hi_ok]
        A[cbind(which(lo_ok),y[lo_ok]-1L)] <- A[cbind(which(lo_ok),y[lo_ok]-1L)]-Glo[lo_ok]
        out <- matrix(0,m,d,dimnames=list(NULL,dnames))
        out[,problem$index$xi] <- rowsum(A,blk,reorder=TRUE) %*% problem$tJac
        if (ncol(X)) out[,problem$index$beta] <- rowsum(X*(Glo-Ghi),blk,reorder=TRUE)
        for (l in seq_len(D)) {
            L <- state$levels[[l]]
            TL <- tree$levels[[l]]
            for (a in seq_along(TL$terms)) {
                re <- ctx$up$re[[TL$terms[a]]]
                Si <- chol2inv(t(re$L))
                B <- L$inst_b[,TL$cols[[a]],drop=FALSE]
                AB <- B %*% Si
                om <- problem$terms[[TL$terms[a]]]$omega
                gom <- vapply(seq_along(om),function(k) {
                    dS <- re$dSigma[[k]]
                    (rowSums((AB %*% dS)*AB)-sum(Si*dS))/2
                },numeric(nrow(B)))
                out[,om] <- out[,om]+rowsum(Pw[[l]]*matrix(gom,nrow(B)),L$inst_block,reorder=TRUE)
            }
        }
        out
    }

    loglik <- function(theta,state) evaluate(theta,state,FALSE)
    score <- function(theta,state) stats::setNames(list(evaluate(theta,state,TRUE)),problem$chain[1L])
    information <- function(theta,state,check=TRUE) {
        if (!requireNamespace("numDeriv",quietly=TRUE))
            stop("The observed information needs the 'numDeriv' package",call.=FALSE)
        f <- function(th) colSums(evaluate(th,state,TRUE))
        ## Two Richardson rounds: truncation error is far below rounding here.
        J <- numDeriv::jacobian(f,as.numeric(theta),method="Richardson",method.args=list(r=2L))
        gof_symmetric_information(J,dnames,check)
    }
    diagnostics <- function(state) {
        data.frame(level=seq_len(D),group=vapply(state$levels,`[[`,"","group"),
                   q=vapply(state$levels,`[[`,0L,"q"),
                   nodes=vapply(state$levels,`[[`,0L,"nodes"),
                   grid=vapply(state$levels,`[[`,0L,"Q"),
                   problems=vapply(state$levels,`[[`,0L,"P"),
                   newton_max=vapply(state$levels,function(L) max(L$iterations),0L),
                   newton_mean=vapply(state$levels,function(L) mean(L$iterations),0),
                   unconverged=vapply(state$levels,function(L) sum(!L$converged),0L),
                   expected_information=vapply(state$levels,function(L) sum(L$fallback),0L),
                   expected_steps=vapply(state$levels,function(L) sum(L$fallback_any),0L),
                   halvings=vapply(state$levels,`[[`,0L,"halvings"),
                   stringsAsFactors=FALSE)
    }
    out <- gof_integrator("nested-agq",refresh=refresh,loglik=loglik,score=score,
                          information=information,diagnostics=diagnostics)
    out$nodes <- nodes
    out
}

## -J symmetrised, with asymmetry and positive-definiteness checks.
gof_symmetric_information <- function(J,names,check=TRUE) {
    asym <- max(abs(J-t(J)))/max(max(abs(J)),.Machine$double.xmin)
    H <- -(J+t(J))/2
    dimnames(H) <- list(names,names)
    attr(H,"asymmetry") <- asym
    if (check) {
        if (!is.finite(asym) || asym>1e-4)
            stop("Observed information is not symmetric (relative asymmetry ",signif(asym,3),
                 "); the score is not the gradient of the log-likelihood",call.=FALSE)
        if (inherits(try(chol(H),silent=TRUE),"try-error"))
            stop("Observed information is not positive definite",call.=FALSE)
    }
    H
}

## ---------------------------------------------------------------------------
## Registry
## ---------------------------------------------------------------------------

gof_integrator_registry <- list(
    "nesting-chain"=list(name="nested-agq",
                         make=function(problem,nodes) gof_integrator_nested_agq(problem,nodes)))

## @noRd
## Select the integrator for a problem's random-effect structure. Returns a
## factory `function(nodes = NULL)` building the integrator; `nodes`
## overrides the per-level nodes (used by the node check with Q + 4).
gof_select_integrator <- function(problem,control=gof_clmm_control()) {
    if (!inherits(problem,"gof_clmm_problem")) stop("problem must be a gof_clmm_problem",call.=FALSE)
    if (identical(problem$structure,"crossed"))
        stop("Crossed random effects are not supported in this version: crossed grouping factors leave ",
             "no independent clusters, which the cluster-multiplier theory requires. ",
             "A future version may add them.",call.=FALSE)
    entry <- gof_integrator_registry[[problem$structure]]
    if (is.null(entry)) stop("No integrator for random-effect structure '",problem$structure,"'",call.=FALSE)
    groups <- vapply(problem$terms,`[[`,"","group")
    q <- vapply(problem$chain,function(g) sum(vapply(problem$terms[groups==g],`[[`,0L,"q")),0L)
    default <- unname(vapply(q,gof_clmm_default_nodes,0L,control=control))
    factory <- function(nodes=NULL) entry$make(problem,if (is.null(nodes)) default else nodes)
    attr(factory,"integrator") <- entry$name
    attr(factory,"nodes") <- unname(default)
    factory
}
