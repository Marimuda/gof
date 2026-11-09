## Structure layer for cumulative link mixed models (ordinal::clmm): turns a
## fitted clmm object into a validated `gof_clmm_problem` (see
## clmm-contracts.R). Everything here is base R + stats; the ordinal package is
## never called, only the fields of the fitted object are read.
##
## Verified conventions of ordinal (2026.7.26) that this file relies on:
## * Random-effect bars are found with findbars() and nested grouping
##   `(1|s/c)` is expanded to the terms `c:s` and `s` (and `(1|a/b/c)` to
##   `c:(b:a)`, `b:a`, `a`). Each grouping factor is evaluated on the model
##   frame restricted to rows with positive weight, after droplevels().
## * Terms are reordered by decreasing number of levels with
##   `rev(order(nlev))` whenever any two terms differ in their number of
##   levels; `fit$ST` and `VarCorr(fit)` follow that order and are named by
##   the grouping expression (names repeat when terms share a factor), while
##   `fit$gfList` holds each distinct grouping factor once.
## * `ST[[t]]` is a lower-triangular factor with Sigma_t = ST %*% t(ST)
##   (VarCorr uses tcrossprod); the factor is not unique in sign, so gof
##   always re-factors Sigma with chol().
## * `fit$alpha` holds the threshold *parameters* (gof's xi; names "1|2",
##   "central", "spacing.1", "threshold.1", ...) and the cutpoints are
##   `fit$Theta = tJac %*% alpha` (gof's alpha).
## * Model: P(Y <= k | b) = F(Theta_k - x'beta - offset - z'b).
## * Weights are stored in the model frame as "(weights)"; `fit$model` keeps
##   zero-weight rows. The fixed design is model.matrix(fit$terms, fit$model,
##   fit$contrasts) with columns possibly dropped for rank deficiency; the
##   kept columns are names(fit$beta).
## * An `offset=` argument is swallowed by clmm's `...`, passed to the
##   optimiser control and ignored with only a warning. A formula offset()
##   is read with model.offset() on a frame whose terms come from the
##   fixed-effects formula, so its column position is only right when the
##   offset() term precedes all random-effect terms.

## Tolerance on clmm's reported max |gradient| (log-likelihood scale), scaled
## by sqrt(N) so the implied standardized shift is about 1e-3.
gof_clmm_gradient_tol <- 1e-3

## ---------------------------------------------------------------------------
## Random-effect bar parsing
## ---------------------------------------------------------------------------

gof_deparse <- function(expr)
    paste(deparse(expr,width.cutoff=500L,backtick=mode(expr) %in% c("call","expression","(","function")),
          collapse=" ")

## Bars `(lhs | g)` in a formula right-hand side (same traversal as lme4/ordinal).
gof_find_bars <- function(term) {
    if (is.name(term) || !is.language(term)) return(NULL)
    if (identical(term[[1L]],as.name("("))) return(gof_find_bars(term[[2L]]))
    if (identical(term[[1L]],as.name("|"))) return(list(term))
    if (length(term)==2L) return(gof_find_bars(term[[2L]]))
    c(gof_find_bars(term[[2L]]),gof_find_bars(term[[3L]]))
}

## Remove bars from a formula right-hand side (NULL if nothing remains).
gof_no_bars <- function(term) {
    if (!("|" %in% all.names(term))) return(term)
    if (is.call(term) && identical(term[[1L]],as.name("|"))) return(NULL)
    if (length(term)==2L) {
        if (identical(term[[1L]],as.name("("))) return(gof_no_bars(term[[2L]]))
        nb <- gof_no_bars(term[[2L]])
        if (is.null(nb)) return(NULL)
        term[[2L]] <- nb
        return(term)
    }
    nb2 <- gof_no_bars(term[[2L]])
    nb3 <- gof_no_bars(term[[3L]])
    if (is.null(nb2)) return(nb3)
    if (is.null(nb3)) return(nb2)
    term[[2L]] <- nb2
    term[[3L]] <- nb3
    term
}

## `s/c/d` -> list of grouping expressions named exactly as clmm names them:
## d:(c:s), c:s, s (finest first).
gof_expand_slash <- function(bar) {
    slash <- function(x) {
        if (!("/" %in% all.names(x))) return(x)
        if (!identical(x[[1L]],as.name("/")))
            stop("Cannot parse the grouping factor '",gof_deparse(x),"' in a random-effect term; ",
                 "use nesting of the form (1 | s/c)",call.=FALSE)
        list(slash(x[[2L]]),slash(x[[3L]]))
    }
    interact <- function(x) {
        if (!is.list(x)) return(x)
        inner <- interact(x[[1L]])
        first <- if (is.list(inner)) inner[[1L]] else inner
        list(substitute(a:b,list(a=x[[2L]],b=first)),inner)
    }
    parts <- slash(bar[[3L]])
    if (!is.list(parts)) return(list(bar))
    lapply(unlist(interact(parts)),function(g) substitute(a | b,list(a=bar[[2L]],b=g)))
}

## @noRd
## Parse the random-effect terms of a clmm formula.
##
## @param formula Model formula with bars `(lhs | g)`.
## @param frame Data frame in which the grouping factors and the left-hand
##   side designs are evaluated (the model frame restricted to rows with
##   positive weight, after droplevels()).
## @returns list(terms, fixed): `terms` in clmm's order (that of `fit$ST`),
##   each list(label, group, index, n_levels, levels, Z, q); `fixed` is the
##   fixed-effects formula with the bars removed.
gof_re_terms <- function(formula,frame) {
    if (!inherits(formula,"formula") || length(formula)!=3L)
        stop("A two-sided model formula is needed",call.=FALSE)
    if ("||" %in% all.names(formula[[3L]]))
        stop("Uncorrelated random-effect terms '(x || g)' are not supported (clmm does not support them ",
             "either); write them as separate terms, e.g. (1 | g) + (0 + x | g)",call.=FALSE)
    env <- environment(formula)
    if (is.null(env)) env <- parent.frame()
    bars <- gof_find_bars(formula[[3L]])
    if (!length(bars)) stop("The formula has no random-effect terms of the form (1 | g)",call.=FALSE)
    bars <- unlist(lapply(bars,gof_expand_slash))
    terms <- lapply(bars,function(bar) {
        label <- gof_deparse(bar)
        group <- gof_deparse(bar[[3L]])
        gf <- tryCatch(eval(bar[[3L]],frame,env),error=function(e)
            stop("Cannot evaluate the grouping factor '",group,"' of random-effect term (",label,"): ",
                 conditionMessage(e),call.=FALSE))
        if (length(gf)!=nrow(frame))
            stop("Grouping factor '",group,"' has length ",length(gf)," but the model frame has ",
                 nrow(frame)," rows",call.=FALSE)
        if (anyNA(gf)) stop("Grouping factor '",group,"' has missing values",call.=FALSE)
        gf <- droplevels(as.factor(gf))
        lhs <- stats::as.formula(call("~",bar[[2L]]),env=env)
        Z <- tryCatch(stats::model.matrix(lhs,frame),error=function(e)
            stop("Cannot build the design of random-effect term (",label,"): ",conditionMessage(e),call.=FALSE))
        if (!ncol(Z)) stop("Random-effect term (",label,") has an empty design",call.=FALSE)
        Z <- matrix(as.numeric(Z),nrow(Z),ncol(Z),dimnames=list(NULL,colnames(Z)))
        list(label=label,group=group,index=as.integer(gf),n_levels=nlevels(gf),
             levels=levels(gf),Z=Z,q=ncol(Z))
    })
    ## clmm's term order (getREterms): decreasing number of levels via
    ## rev(order(.)) whenever the numbers of levels are not all equal.
    nlev <- vapply(terms,`[[`,0L,"n_levels")
    if (any(diff(nlev)!=0L)) terms <- terms[rev(order(nlev))]
    fixed <- gof_no_bars(formula[[3L]])
    if (is.null(fixed)) fixed <- 1
    fixed <- stats::as.formula(call("~",formula[[2L]],fixed),env=env)
    list(terms=terms,fixed=fixed)
}

## ---------------------------------------------------------------------------
## Factor relations and independent blocks
## ---------------------------------------------------------------------------

## Does every level of a map to exactly one level of b (a nested in b)?
gof_nested_in <- function(a,b) {
    pairs <- unique(cbind(a,b))
    !anyDuplicated(pairs[,1L])
}

## @noRd
## Pairwise relations between the distinct grouping factors of `terms`:
## matrix with entry [A, B] one of "identical" (same partition of the rows),
## "nested_in" (each level of A lies within exactly one level of B),
## "contains" (B nested in A) or "crossed".
gof_factor_relations <- function(terms) {
    groups <- vapply(terms,`[[`,"","group")
    keep <- !duplicated(groups)
    groups <- groups[keep]
    idx <- lapply(terms[keep],`[[`,"index")
    G <- length(groups)
    rel <- matrix("identical",G,G,dimnames=list(groups,groups))
    for (i in seq_len(G)) for (j in seq_len(G)) if (i!=j) {
        ab <- gof_nested_in(idx[[i]],idx[[j]])
        ba <- gof_nested_in(idx[[j]],idx[[i]])
        rel[i,j] <- if (ab && ba) "identical" else if (ab) "nested_in" else if (ba) "contains" else "crossed"
    }
    rel
}

## @noRd
## Independent blocks: connected components of the bipartite graph between
## rows and (term, level) nodes, computed by label propagation (each pass
## replaces a row's label by the smallest label among rows sharing a level).
## Blocks are numbered 1..m by first appearance in the rows.
gof_independent_blocks <- function(n,terms) {
    lab <- seq_len(n)
    idx <- lapply(terms,`[[`,"index")
    if (any(vapply(idx,length,0L)!=n)) stop("Every term index needs one entry per row",call.=FALSE)
    repeat {
        old <- lab
        for (ix in idx) {
            o <- order(ix,lab)
            first <- !duplicated(ix[o])
            low <- integer(max(ix))
            low[ix[o][first]] <- lab[o][first]
            lab <- low[ix]
        }
        if (identical(lab,old)) break
    }
    id <- match(lab,unique(lab))
    list(id=id,m=max(id))
}

## "nesting-chain" when no two distinct grouping factors are crossed; the
## chain lists the groups from coarsest (fewest levels) to finest.
gof_clmm_chain <- function(terms,relations) {
    if (any(relations=="crossed")) return(list(structure="crossed",chain=NULL))
    groups <- rownames(relations)
    nlev <- vapply(terms[match(groups,vapply(terms,`[[`,"","group"))],`[[`,0L,"n_levels")
    list(structure="nesting-chain",chain=groups[order(nlev)])
}

## ---------------------------------------------------------------------------
## Parameter mapping
## ---------------------------------------------------------------------------

## @noRd
## Map clmm estimates to gof's theta = c(xi, beta, omega_1, ..., omega_T).
##
## @param model A clmm fit.
## @param problem_parts A list with `terms` as returned by gof_re_terms()
##   (clmm's term order).
## @returns Named theta. xi = fit$alpha (threshold parameters, cutpoints
##   tJac %*% xi), beta = fit$beta, omega_t = log-Cholesky of
##   tcrossprod(ST[[t]]), named "<term label>.L11", ... Attribute
##   `boundary` lists terms whose clmm covariance was singular and was
##   regularised by a tiny ridge so that omega is finite.
gof_clmm_theta <- function(model,problem_parts) {
    terms <- problem_parts$terms
    ST <- model$ST
    if (!is.list(ST) || length(ST)!=length(terms))
        stop("The clmm fit has ",length(ST)," random-effect covariance factors but the formula gives ",
             length(terms)," random-effect terms",call.=FALSE)
    xi <- model$alpha
    if (is.null(names(xi)) || !length(xi))
        stop("The clmm fit has no named threshold parameters (alpha)",call.=FALSE)
    beta <- model$beta
    if (is.null(beta)) beta <- numeric(0)
    if (length(beta) && is.null(names(beta))) stop("The clmm fixed effects (beta) are unnamed",call.=FALSE)
    omega <- list(); boundary <- character()
    for (t in seq_along(terms)) {
        tm <- terms[[t]]; st <- ST[[t]]
        if (!is.null(names(ST)) && names(ST)[t]!=tm$group)
            stop("Cannot match random-effect term (",tm$label,") to clmm's covariance factor '",
                 names(ST)[t],"'",call.=FALSE)
        if (!is.matrix(st) || nrow(st)!=tm$q || ncol(st)!=tm$q ||
            (!is.null(colnames(st)) && !identical(colnames(st),colnames(tm$Z))))
            stop("clmm's covariance factor for term (",tm$label,") does not match its design columns (",
                 paste(colnames(tm$Z),collapse=", "),")",call.=FALSE)
        if (any(!is.finite(st))) stop("clmm's covariance factor for term (",tm$label,") is not finite",call.=FALSE)
        Sigma <- tcrossprod(unname(st))
        ok <- tryCatch({chol(Sigma); min(eigen(Sigma,symmetric=TRUE,only.values=TRUE)$values) >
                        1e-12*max(1,max(diag(Sigma)))},error=function(e) FALSE)
        if (!ok) {
            Sigma <- Sigma+diag(1e-10*max(1,max(diag(Sigma))),tm$q)
            boundary <- c(boundary,tm$label)
        }
        om <- gof_logchol_omega(Sigma)
        omega[[t]] <- stats::setNames(om,paste0(tm$label,".",gof_logchol_names(tm$q)))
    }
    theta <- c(xi,beta,unlist(omega))
    attr(theta,"boundary") <- boundary
    theta
}

## ---------------------------------------------------------------------------
## Problem construction
## ---------------------------------------------------------------------------

## Offsets from offset() terms of the fixed-effects terms, looked up by name.
gof_formula_offset <- function(tt,mf) {
    off <- numeric(nrow(mf))
    pos <- attr(tt,"offset")
    if (!length(pos)) return(off)
    vars <- vapply(as.list(attr(tt,"variables"))[-1L],gof_deparse,"")
    for (v in vars[pos]) {
        if (!v %in% names(mf)) stop("Offset term '",v,"' is missing from the model frame",call.=FALSE)
        off <- off+as.numeric(mf[[v]])
    }
    off
}

## @noRd
## Validate a clmm fit and build the corresponding gof_clmm_problem.
gof_clmm_problem <- function(model,control=gof_clmm_control()) {
    if (!inherits(model,"clmm")) stop("model must be a fitted ordinal::clmm object",call.=FALSE)
    if (!inherits(control,"gof_clmm_control")) stop("control must come from gof_clmm_control()",call.=FALSE)
    call <- model$call
    if (!is.null(call) && "offset" %in% names(as.list(call)))
        stop("The clmm fit was called with an 'offset' argument, which clmm silently ignores ",
             "(it is passed on to the optimiser control); put the offset in the formula as ",
             "offset(<variable>) before the random-effect terms and refit",call.=FALSE)
    mf <- model$model
    if (!is.data.frame(mf) || !nrow(mf))
        stop("The clmm fit has no model frame; refit with model = TRUE",call.=FALSE)
    link <- model$link
    if (!is.character(link) || length(link)!=1L || !link %in% gof_link_names)
        stop("Unsupported clmm link '",paste(link,collapse=","),"'; supported links are ",
             paste(gof_link_names,collapse=", "),call.=FALSE)

    ## Response, weights and categories.
    yf <- stats::model.response(mf)
    if (!is.factor(yf)) stop("The clmm response is not a factor",call.=FALSE)
    lev <- levels(yf)
    K <- length(lev)
    if (!is.null(model$y.levels) && !identical(lev,model$y.levels))
        stop("The response levels in the model frame do not match the fit",call.=FALSE)
    w <- stats::model.weights(mf)
    if (is.null(w)) w <- rep(1,nrow(mf))
    w <- as.numeric(w)
    if (any(!is.finite(w)) || any(w<0)) stop("clmm weights must be finite and nonnegative",call.=FALSE)
    keep <- w>0
    if (!any(keep)) stop("Every observation has zero weight",call.=FALSE)
    y <- as.integer(yf)
    tot <- vapply(seq_len(K),function(k) sum(w[keep][y[keep]==k]),0)
    if (any(tot<=0))
        stop("Response categor",if (sum(tot<=0)>1L) "ies " else "y ",
             paste0("'",lev[tot<=0],"'",collapse=", ")," ",if (sum(tot<=0)>1L) "have" else "has",
             " no observations with positive weight; its thresholds are not identified. ",
             "Merge or drop empty categories and refit",call.=FALSE)

    ## Convergence of the supplied fit.
    conv <- model$optRes$convergence
    if (is.null(conv) || !identical(as.integer(conv),0L))
        stop("The clmm fit did not converge (optimizer convergence code ",
             if (is.null(conv)) "missing" else conv,
             if (!is.null(model$optRes$message)) paste0(": ",model$optRes$message),
             "); refit with different starting values or control settings first",call.=FALSE)
    N <- sum(w[keep])
    if (!is.null(model$gradient)) {
        mg <- max(abs(model$gradient))
        tol <- gof_clmm_gradient_tol*max(1,sqrt(N))
        if (!is.finite(mg) || mg>tol)
            stop("The clmm fit has not converged: max |gradient| = ",format(mg,digits=3),
                 " exceeds ",format(tol,digits=3),"; refit (e.g. with other starting values) first",
                 call.=FALSE)
    }

    ## Fixed design and offset.
    tt <- model$terms
    if (is.null(tt)) tt <- attr(mf,"terms")
    beta <- model$beta
    if (is.null(beta)) beta <- numeric(0)
    Xfull <- stats::model.matrix(tt,mf,contrasts.arg=model$contrasts)
    miss <- setdiff(names(beta),colnames(Xfull))
    if (length(miss))
        stop("Cannot rebuild the fixed-effects design: columns ",paste(miss,collapse=", "),
             " of the fit are missing",call.=FALSE)
    X <- Xfull[keep,names(beta),drop=FALSE]
    X <- matrix(as.numeric(X),nrow(X),ncol(X),dimnames=list(NULL,names(beta)))
    off <- gof_formula_offset(tt,mf)
    used <- stats::model.offset(mf)
    if (is.null(used)) used <- numeric(nrow(mf))
    if (length(used)!=nrow(mf) || !is.numeric(used) || !isTRUE(all.equal(as.numeric(used),off)))
        stop("clmm read the formula offset from the wrong model-frame column (this happens when ",
             "offset() follows a random-effect term in the formula); move offset() before the ",
             "random-effect terms, e.g. y ~ x + offset(o) + (1 | g), and refit",call.=FALSE)
    if (any(!is.finite(off[keep]))) stop("The formula offset has non-finite values",call.=FALSE)

    ## Thresholds.
    tJac <- model$tJac
    if (!is.matrix(tJac) || nrow(tJac)!=K-1L || ncol(tJac)!=length(model$alpha))
        stop("The clmm threshold Jacobian (tJac) does not match the ",K," response categories",call.=FALSE)

    ## Random effects on the positive-weight rows (as clmm builds them).
    rf <- droplevels(mf[keep,,drop=FALSE])
    re <- gof_re_terms(model$formula,rf)
    terms <- re$terms
    gfl <- model$gfList
    for (tm in terms) {
        if (tm$n_levels<3L)
            stop("Grouping factor '",tm$group,"' has ",tm$n_levels," levels; at least 3 are needed",
                 call.=FALSE)
        if (tm$q>control$max.dim)
            stop("Random-effect term (",tm$label,") has dimension ",tm$q," per level, more than ",
                 "max.dim = ",control$max.dim,"; simplify the term or raise gof_clmm_control(max.dim)",
                 call.=FALSE)
        if (!is.null(gfl)) {
            f <- gfl[[tm$group]]
            if (is.null(f) || length(f)!=length(tm$index) ||
                !identical(as.integer(f),tm$index) || !identical(levels(f),tm$levels))
                stop("The parsed grouping factor '",tm$group,"' does not match the clmm fit; ",
                     "please report this with a reproducible example",call.=FALSE)
        }
    }
    theta <- gof_clmm_theta(model,list(terms=terms))
    boundary <- attr(theta,"boundary")
    attr(theta,"boundary") <- NULL

    relations <- gof_factor_relations(terms)
    n <- sum(keep)
    blocks <- gof_independent_blocks(n,terms)
    st <- gof_clmm_chain(terms,relations)
    if (st$structure=="nesting-chain") {
        if (blocks$m<=length(theta))
            stop("Only ",blocks$m," independent blocks (levels of '",st$chain[1L],"') for ",
                 length(theta)," parameters; the multiplier correction needs more blocks than ",
                 "parameters. Simplify the model or use more groups",call.=FALSE)
        if (blocks$m<30L)
            warning("Only ",blocks$m," independent blocks (levels of '",st$chain[1L],"'); the ",
                    "cluster multiplier approximation may be inaccurate with fewer than 30",call.=FALSE)
    }
    src <- list(call=call,dropped=sum(!keep),theta_clmm=theta,data=mf,fixed=re$fixed,
                threshold=model$threshold,nAGQ=model$nAGQ,boundary=boundary,rows=which(keep))
    new_gof_clmm_problem(y=y[keep],K=K,X=X,offset=off[keep],weights=w[keep],tJac=tJac,link=link,
                         terms=terms,blocks=blocks,structure=st$structure,chain=st$chain,
                         theta=theta,levels=lev,relations=relations,source=src)
}
