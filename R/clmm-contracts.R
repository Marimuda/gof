## Shared contracts for ordinal multiplier processes and cumulative link
## mixed-model (clmm) diagnostics. Every component codes against these
## constructors/validators; the files implementing the components
## (ordinal-engine.R, clmm-structure.R, clmm-integrators.R, clmm-polish.R,
## clmm-marginal.R, cumres.clmm.R) never change these structures ad hoc.
##
## The structures are deliberately general enough for crossed random
## effects (several multiplier factors, several score matrices), even though
## 1.4.0 only accepts nesting chains with one independent top-level factor.

## ---------------------------------------------------------------------------
## Multiplier design
## ---------------------------------------------------------------------------

## Registered multiplier draw and combination rules. A future crossed design
## registers new rules here instead of changing the engine.
gof_multiplier_rules <- list(
    draw=list(
        ## One iid N(0,1) multiplier per level, factors drawn in order. With a
        ## single factor this is exactly rnorm(L), preserving legacy streams.
        "independent-normal"=function(levels) lapply(levels,stats::rnorm)),
    combine=list(
        ## Row multiplier: sum of factor multipliers (times optional loadings).
        "sum"=function(g,idx,loading) {
            if (length(g)==1L && is.null(loading[[1L]])) return(g[[1L]][idx[[1L]]])
            out <- 0
            for (f in seq_along(g)) {
                term <- g[[f]][idx[[f]]]
                if (!is.null(loading[[f]])) term <- loading[[f]]*term
                out <- out+term
            }
            out
        }))

## @noRd
## Multiplier design for an ordinal cumulative residual process.
##
## @param factors List of multiplier factors. Each is a list with
##   `name` (string), `index` (integer row -> level map with levels 1..L, all
##   used), `influence` (L x d coefficient influence contributions on the
##   sum scale: theta.hat - theta0 ~ sum over factors and levels) and
##   optional `loading` (numeric per-row multiplier loading, default NULL).
## @param n_units Number of independent units (used for reporting).
## @param divisor Positive process scale divisor applied to cumulative sums.
## @param draw,combine Names of registered rules in `gof_multiplier_rules`.
gof_multiplier_design <- function(factors,n_units,divisor,
                                  draw="independent-normal",combine="sum") {
    if (!is.list(factors) || !length(factors)) stop("A multiplier design needs at least one factor",call.=FALSE)
    if (!is.character(draw) || length(draw)!=1L || !draw %in% names(gof_multiplier_rules$draw))
        stop("Unknown multiplier draw rule: ",paste(draw,collapse=","),call.=FALSE)
    if (!is.character(combine) || length(combine)!=1L || !combine %in% names(gof_multiplier_rules$combine))
        stop("Unknown multiplier combination rule: ",paste(combine,collapse=","),call.=FALSE)
    n <- NULL; d <- NULL; nm <- character()
    for (i in seq_along(factors)) {
        f <- factors[[i]]
        if (!is.list(f) || !all(c("name","index","influence") %in% names(f)))
            stop("Each multiplier factor needs name, index and influence",call.=FALSE)
        if (!is.character(f$name) || length(f$name)!=1L || is.na(f$name) || !nzchar(f$name))
            stop("Multiplier factor names must be nonempty strings",call.=FALSE)
        idx <- f$index
        if (!is.numeric(idx) || !is.null(dim(idx)) || !length(idx) || any(!is.finite(idx)) ||
            any(idx!=floor(idx)) || any(idx<1))
            stop("Multiplier factor '",f$name,"' needs a positive integer level index",call.=FALSE)
        if (is.null(n)) n <- length(idx) else if (length(idx)!=n)
            stop("All multiplier factors must index the same rows",call.=FALSE)
        ic <- f$influence
        if (!is.matrix(ic) || !is.numeric(ic) || is.complex(ic) || any(!is.finite(ic)))
            stop("Multiplier factor '",f$name,"' needs a finite numeric influence matrix",call.=FALSE)
        L <- nrow(ic)
        if (max(idx)>L || length(unique(idx))!=L)
            stop("Multiplier factor '",f$name,"': every influence row must be a used level",call.=FALSE)
        if (is.null(d)) d <- ncol(ic) else if (ncol(ic)!=d)
            stop("All multiplier influence matrices need the same parameter columns",call.=FALSE)
        if (!is.null(f$loading) && (!is.numeric(f$loading) || length(f$loading)!=n ||
            any(!is.finite(f$loading))))
            stop("Multiplier loadings must be finite and have one value per row",call.=FALSE)
        factors[[i]] <- list(name=f$name,index=as.integer(idx),influence=ic,loading=f$loading)
        nm <- c(nm,f$name)
    }
    if (anyDuplicated(nm)) stop("Multiplier factor names must be distinct",call.=FALSE)
    if (!is.numeric(n_units) || length(n_units)!=1L || !is.finite(n_units) ||
        n_units<1 || n_units!=floor(n_units))
        stop("n_units must be a positive whole number",call.=FALSE)
    if (!is.numeric(divisor) || length(divisor)!=1L || !is.finite(divisor) || divisor<=0)
        stop("divisor must be a positive finite number",call.=FALSE)
    structure(list(factors=factors,n_units=as.integer(n_units),divisor=divisor,
                   draw=draw,combine=combine,n=n,d=d),class="gof_multiplier_design")
}

## One row per factor, stored on diagnostics as `multiplier`.
gof_multiplier_summary <- function(design) {
    data.frame(factor=vapply(design$factors,`[[`,"","name"),
               levels=vapply(design$factors,function(f) nrow(f$influence),0),
               rows=rep(design$n,length(design$factors)),
               draw=design$draw,combine=design$combine,stringsAsFactors=FALSE)
}

## ---------------------------------------------------------------------------
## Cumulative link functions (shared by integrators and marginal layer)
## ---------------------------------------------------------------------------

gof_link_names <- c("logit","probit","cloglog","loglog","cauchit")

## @noRd
## Link table: F (CDF), f (density), df (density derivative), and stable
## logF = log F(x) and logS = log(1-F(x)). Conventions follow ordinal:
## cloglog F = 1-exp(-exp(x)), loglog F = exp(-exp(-x)).
gof_ordinal_link <- function(name) {
    if (!is.character(name) || length(name)!=1L || !name %in% gof_link_names)
        stop("Unsupported ordinal link: ",paste(name,collapse=","),call.=FALSE)
    log1mexp <- function(a) ifelse(a<=log(2),log(-expm1(-a)),log1p(-exp(-a)))  # log(1-exp(-a)), a>0
    ## Densities and their derivatives vanish in both tails for every link;
    ## guard 0*Inf and Inf/Inf at infinite arguments. The double-exponential
    ## links are evaluated at |x| <= 700, where their densities have already
    ## underflowed to zero, so the clamp is exact.
    tails <- function(fn) {force(fn); function(x) {v <- fn(x); v[!is.finite(x)] <- 0; v}}
    hi <- function(x) pmin(x,700)
    lo <- function(x) pmax(x,-700)
    lk <- switch(name,
        logit=list(name=name,F=stats::plogis,f=stats::dlogis,
                   df=function(x) {p <- stats::plogis(x); stats::dlogis(x)*(1-2*p)},
                   logF=function(x) stats::plogis(x,log.p=TRUE),
                   logS=function(x) stats::plogis(x,lower.tail=FALSE,log.p=TRUE)),
        probit=list(name=name,F=stats::pnorm,f=stats::dnorm,
                    df=function(x) -x*stats::dnorm(x),
                    logF=function(x) stats::pnorm(x,log.p=TRUE),
                    logS=function(x) stats::pnorm(x,lower.tail=FALSE,log.p=TRUE)),
        cloglog=list(name=name,F=function(x) -expm1(-exp(x)),
                     f=function(x) {x <- hi(x); exp(x-exp(x))},
                     df=function(x) {x <- hi(x); exp(x-exp(x))*(1-exp(x))},
                     logF=function(x) log1mexp(exp(x)),
                     logS=function(x) -exp(x)),
        loglog=list(name=name,F=function(x) exp(-exp(-x)),
                    f=function(x) {x <- lo(x); exp(-x-exp(-x))},
                    df=function(x) {x <- lo(x); exp(-x-exp(-x))*(exp(-x)-1)},
                    logF=function(x) -exp(-x),
                    logS=function(x) log1mexp(exp(-x))),
        cauchit=list(name=name,F=stats::pcauchy,f=stats::dcauchy,
                     df=function(x) -2*x/(pi*(1+x^2)^2),
                     logF=function(x) stats::pcauchy(x,log.p=TRUE),
                     logS=function(x) stats::pcauchy(x,lower.tail=FALSE,log.p=TRUE)))
    lk$f <- tails(lk$f)
    lk$df <- tails(lk$df)
    lk
}

## ---------------------------------------------------------------------------
## Random-effect covariance parametrisation (log-Cholesky)
## ---------------------------------------------------------------------------

## omega holds the lower triangle of L (column-major), with log diagonal.
gof_logchol_names <- function(q) {
    ij <- which(lower.tri(diag(q),diag=TRUE),arr.ind=TRUE)
    ij <- ij[order(ij[,2L],ij[,1L]),,drop=FALSE]
    paste0("L",ij[,1L],ij[,2L])
}

gof_logchol_L <- function(omega,q) {
    if (length(omega)!=q*(q+1L)/2L) stop("Log-Cholesky parameter length does not match dimension",call.=FALSE)
    L <- matrix(0,q,q)
    ij <- which(lower.tri(L,diag=TRUE),arr.ind=TRUE)
    ij <- ij[order(ij[,2L],ij[,1L]),,drop=FALSE]
    val <- ifelse(ij[,1L]==ij[,2L],exp(omega),omega)
    L[ij] <- val
    L
}

## Inverse map: a positive-definite Sigma -> omega.
gof_logchol_omega <- function(Sigma) {
    q <- nrow(Sigma)
    L <- t(chol(Sigma))
    ij <- which(lower.tri(L,diag=TRUE),arr.ind=TRUE)
    ij <- ij[order(ij[,2L],ij[,1L]),,drop=FALSE]
    val <- L[ij]
    diagonal <- ij[,1L]==ij[,2L]
    val[diagonal] <- log(val[diagonal])
    stats::setNames(val,gof_logchol_names(q))
}

## Sigma = L L' and its derivatives with respect to each omega entry.
gof_logchol_sigma <- function(omega,q) {
    L <- gof_logchol_L(omega,q)
    ij <- which(lower.tri(L,diag=TRUE),arr.ind=TRUE)
    ij <- ij[order(ij[,2L],ij[,1L]),,drop=FALSE]
    dSigma <- lapply(seq_len(nrow(ij)),function(k) {
        dL <- matrix(0,q,q)
        dL[ij[k,1L],ij[k,2L]] <- if (ij[k,1L]==ij[k,2L]) L[ij[k,1L],ij[k,2L]] else 1
        dL %*% t(L)+L %*% t(dL)
    })
    list(Sigma=L %*% t(L),L=L,dSigma=dSigma)
}

## ---------------------------------------------------------------------------
## clmm problem object
## ---------------------------------------------------------------------------

## @noRd
## Validated cumulative link mixed-model problem: everything the integrator,
## marginal layer and method need, independent of the ordinal package.
##
## Model: P(Y_i <= k | b) = F(alpha_k - x_i'beta - offset_i - sum_t z_it' b_t[level_it]),
## alpha = tJac %*% xi, b_t[l] ~ N(0, Sigma_t) independently over terms and levels.
## theta = c(xi, beta, omega_1, ..., omega_T) with log-Cholesky omega_t.
##
## @param terms List of random-effect terms: label, group (grouping factor
##   name), index (row -> level 1..n_levels), n_levels, Z (n x q, colnames),
##   q. `omega` indices are assigned here.
## @param blocks list(id = row -> independent block 1..m, m).
## @param structure "nesting-chain" or "crossed"; `chain` lists the term
##   groups from coarsest to finest for a nesting chain.
new_gof_clmm_problem <- function(y,K,X,offset,weights,tJac,link,terms,blocks,
                                 structure,chain=NULL,theta,levels=NULL,
                                 relations=NULL,source=list()) {
    n <- length(y)
    if (!is.numeric(K) || length(K)!=1L || K<2 || K!=floor(K)) stop("K must be a whole number >= 2",call.=FALSE)
    if (!n || !is.numeric(y) || any(!is.finite(y)) || any(y!=floor(y)) || any(y<1 | y>K))
        stop("y must contain category indices 1..K",call.=FALSE)
    if (!is.matrix(X) || nrow(X)!=n || !is.numeric(X) || any(!is.finite(X)))
        stop("X must be a finite numeric matrix with one row per observation",call.=FALSE)
    if (ncol(X) && is.null(colnames(X))) stop("X needs column names",call.=FALSE)
    if (!is.numeric(offset) || length(offset)!=n || any(!is.finite(offset)))
        stop("offset must be finite with one value per observation",call.=FALSE)
    if (!is.numeric(weights) || length(weights)!=n || any(!is.finite(weights)) || any(weights<=0))
        stop("weights must be positive and finite (drop zero-weight rows first)",call.=FALSE)
    if (!is.matrix(tJac) || nrow(tJac)!=K-1L || !ncol(tJac) || any(!is.finite(tJac)))
        stop("tJac must be a finite (K-1) x n_threshold matrix",call.=FALSE)
    if (!link %in% gof_link_names) stop("Unsupported ordinal link: ",link,call.=FALSE)
    if (!structure %in% c("nesting-chain","crossed")) stop("structure must be 'nesting-chain' or 'crossed'",call.=FALSE)
    if (!is.list(terms) || !length(terms)) stop("A mixed model needs at least one random-effect term",call.=FALSE)
    nxi <- ncol(tJac); p <- ncol(X)
    pos <- nxi+p
    for (t in seq_along(terms)) {
        tm <- terms[[t]]
        if (!all(c("label","group","index","n_levels","Z") %in% names(tm)))
            stop("Each random-effect term needs label, group, index, n_levels and Z",call.=FALSE)
        if (!is.numeric(tm$index) || length(tm$index)!=n || any(tm$index!=floor(tm$index)) ||
            any(tm$index<1 | tm$index>tm$n_levels) || length(unique(tm$index))!=tm$n_levels)
            stop("Term '",tm$label,"' needs a row -> level index using every level",call.=FALSE)
        if (!is.matrix(tm$Z) || nrow(tm$Z)!=n || !ncol(tm$Z) || any(!is.finite(tm$Z)) || is.null(colnames(tm$Z)))
            stop("Term '",tm$label,"' needs a finite named n x q design Z",call.=FALSE)
        q <- ncol(tm$Z)
        k <- q*(q+1L)/2L
        terms[[t]]$q <- q
        terms[[t]]$index <- as.integer(tm$index)
        terms[[t]]$omega <- pos+seq_len(k)
        pos <- pos+k
    }
    if (!is.list(blocks) || !is.numeric(blocks$id) || length(blocks$id)!=n ||
        any(blocks$id<1 | blocks$id>blocks$m) || length(unique(blocks$id))!=blocks$m)
        stop("blocks must map every row to one of m independent blocks",call.=FALSE)
    blocks$id <- as.integer(blocks$id)
    if (!is.numeric(theta) || length(theta)!=pos || any(!is.finite(theta)) || is.null(names(theta)))
        stop("theta must be a finite named vector of length ",pos,call.=FALSE)
    if (structure=="nesting-chain") {
        groups <- vapply(terms,`[[`,"","group")
        if (is.null(chain)) chain <- unique(groups)
        if (!setequal(chain,groups)) stop("chain must list every term group from coarsest to finest",call.=FALSE)
    }
    structure(list(y=as.integer(y),K=as.integer(K),X=X,offset=as.numeric(offset),
                   weights=as.numeric(weights),n=n,N=sum(weights),tJac=tJac,link=link,
                   terms=terms,blocks=blocks,structure=structure,chain=chain,
                   relations=relations,theta=theta,
                   index=list(xi=seq_len(nxi),beta=nxi+seq_len(p)),
                   levels=levels,source=source),class="gof_clmm_problem")
}

## Unpack theta into model quantities: alpha, beta and per-term Sigma (with
## derivatives dSigma/domega). Every component uses this single mapping.
gof_clmm_unpack <- function(problem,theta=problem$theta) {
    if (length(theta)!=length(problem$theta)) stop("theta has the wrong length",call.=FALSE)
    xi <- theta[problem$index$xi]
    beta <- theta[problem$index$beta]
    re <- lapply(problem$terms,function(tm) gof_logchol_sigma(theta[tm$omega],tm$q))
    list(alpha=drop(problem$tJac %*% xi),xi=xi,beta=beta,re=re)
}

## ---------------------------------------------------------------------------
## Integrator interface
## ---------------------------------------------------------------------------

## @noRd
## A marginal-likelihood integrator for a clmm problem.
##  refresh(theta)            -> state: adaptive nodes etc. (frozen proposal)
##  loglik(theta, state)      -> numeric(m): per-block log-likelihood
##  score(theta, state)       -> named list(<factor> = levels x d) of score
##                               contributions; 1.4.0 returns one element for
##                               the top-level factor, rows = blocks 1..m
##  information(theta, state) -> d x d observed information (frozen nodes)
##  diagnostics(state)        -> list of descriptive diagnostics
gof_integrator <- function(name,refresh,loglik,score,information,diagnostics=function(state) list()) {
    fns <- list(refresh=refresh,loglik=loglik,score=score,information=information,diagnostics=diagnostics)
    if (!is.character(name) || length(name)!=1L || !nzchar(name)) stop("Integrator needs a name",call.=FALSE)
    bad <- names(fns)[!vapply(fns,is.function,logical(1L))]
    if (length(bad)) stop("Integrator components must be functions: ",paste(bad,collapse=", "),call.=FALSE)
    structure(c(list(name=name),fns),class="gof_integrator")
}

## ---------------------------------------------------------------------------
## Control
## ---------------------------------------------------------------------------

##' Control settings for mixed-model ordinal diagnostics
##'
##' Numerical and acceptance settings used by [cumres.clmm()] when it
##' re-evaluates the marginal likelihood with adaptive Gauss--Hermite
##' quadrature and polishes the fitted parameters to a frozen-node score root.
##' Finite quadrature approximates the integrated likelihood; a higher-node
##' calculation checks numerical sensitivity.
##'
##' @param nodes Quadrature nodes per dimension for each random-effect level,
##'   or `NULL` for the defaults 11, 7 and 5 for 1-, 2- and 3-dimensional
##'   effects.
##' @param marginal.nodes Minimum number of quadrature points for the
##'   one-dimensional marginal probabilities and their derivatives.
##' @param score.tol Convergence tolerance for the standardized score
##'   \eqn{\sqrt{S^\top H^{-1}S}}.
##' @param node.tol Maximum standardized score allowed when the polished root
##'   is re-scored with four more nodes per dimension.
##' @param shift.tol Maximum allowed shift from the supplied fit, in standard
##'   errors, for any parameter; larger shifts are an error. The default `Inf`
##'   never refuses: the diagnostic is computed at the polished estimate and
##'   the shift is always reported.
##' @param shift.warn Shift, in standard errors, above which a warning reports
##'   that the diagnosed estimates differ from the supplied fit.
##' @param boundary.tol Relative tolerance defining a covariance at its
##'   boundary: a term whose smallest covariance eigenvalue is at most
##'   `boundary.tol * max(1, largest eigenvalue)`, or whose smallest-eigenvalue
##'   direction is not identified by the information, is refused.
##' @param boundary.z Warn when a term's smallest covariance eigenvalue is
##'   within `boundary.z` standard errors of zero. This flags imprecisely
##'   estimated covariances; it does not refuse the fit.
##' @param max.refresh Maximum number of node refreshes while polishing.
##' @param maxit Maximum Newton iterations per refresh.
##' @param max.dim Largest random-effect dimension per level accepted.
##' @returns A list of class `gof_clmm_control`.
##' @seealso [cumres.clmm()]
##' @export
gof_clmm_control <- function(nodes=NULL,marginal.nodes=30L,score.tol=1e-8,node.tol=0.05,
                             shift.tol=Inf,shift.warn=0.5,boundary.tol=1e-6,boundary.z=1,
                             max.refresh=5L,maxit=50L,max.dim=3L) {
    if (!is.null(nodes)) nodes <- gof_count(nodes,"nodes",1L)
    marginal.nodes <- gof_count(marginal.nodes,"marginal.nodes",2L)
    max.refresh <- gof_count(max.refresh,"max.refresh",1L)
    maxit <- gof_count(maxit,"maxit",1L)
    max.dim <- gof_count(max.dim,"max.dim",1L)
    positive <- function(x,name,inf=FALSE)
        if (!is.numeric(x) || length(x)!=1L || is.na(x) || x<=0 || (!inf && !is.finite(x)))
            stop(name," must be a positive",if (!inf) " finite"," number",call.=FALSE)
    positive(score.tol,"score.tol"); positive(node.tol,"node.tol")
    positive(shift.tol,"shift.tol",TRUE); positive(shift.warn,"shift.warn",TRUE)
    positive(boundary.tol,"boundary.tol"); positive(boundary.z,"boundary.z")
    structure(list(nodes=nodes,marginal.nodes=marginal.nodes,score.tol=score.tol,node.tol=node.tol,
                   shift.tol=shift.tol,shift.warn=shift.warn,boundary.tol=boundary.tol,
                   boundary.z=boundary.z,max.refresh=max.refresh,
                   maxit=maxit,max.dim=max.dim),class="gof_clmm_control")
}

gof_clmm_default_nodes <- function(q,control) {
    if (!is.null(control$nodes)) return(control$nodes)
    c(11L,7L,5L)[min(q,3L)]
}
