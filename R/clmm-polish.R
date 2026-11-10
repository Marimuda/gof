## Polishing clmm estimates to a frozen-node adaptive-quadrature score root,
## score, acceptance checks (nodes, shift, boundary) and the cluster
## multiplier design built from the polished per-block influence.

## Central-difference Jacobian of the total score with frozen state, used for
## the Newton direction (the final information uses Richardson extrapolation).
gof_polish_hessian <- function(integ,theta,state) {
    d <- length(theta)
    J <- matrix(0,d,d)
    for (k in seq_len(d)) {
        h <- 1e-5*max(1,abs(theta[k]))
        tp <- theta; tp[k] <- tp[k]+h
        tm <- theta; tm[k] <- tm[k]-h
        J[,k] <- (colSums(integ$score(tp,state)[[1L]])-colSums(integ$score(tm,state)[[1L]]))/(2*h)
    }
    H <- -(J+t(J))/2
    dimnames(H) <- list(names(theta),names(theta))
    H
}

## Cholesky factor of H, with Levenberg damping when H is not PD.
gof_damped_chol <- function(H) {
    R <- tryCatch(chol(H),error=function(e) NULL)
    lambda <- 0
    if (is.null(R)) {
        lambda <- 1e-8*max(1,max(abs(diag(H))))
        repeat {
            R <- tryCatch(chol(H+diag(lambda,nrow(H))),error=function(e) NULL)
            if (!is.null(R)) break
            lambda <- lambda*10
        }
    }
    list(R=R,lambda=lambda)
}

gof_standardized_score <- function(S,H) {
    ch <- gof_damped_chol(H)
    sqrt(max(0,sum(backsolve(ch$R,S,transpose=TRUE)^2)))
}

## Newton with Armijo step halving on the frozen-state total log-likelihood.
## The central-difference Hessian H is reused across iterations (and passed
## between refreshes): it is recomputed only when the line search has to
## shorten a step taken with a stale H or progress stalls. Returns early
## (refresh=TRUE) once theta has moved more than `max_move` standardized
## units from the refresh point: frozen nodes placed for a distant theta
## integrate poorly, so the proposal is refreshed first.
gof_polish_newton <- function(integ,theta,state,control,H=NULL,max_move=Inf) {
    theta_start <- theta
    R0 <- NULL
    fresh <- FALSE
    dec_old <- Inf
    ## Trial points can leave the parameter space numerically (e.g. a
    ## variance underflowing to zero); they count as -Inf in the line search.
    total <- function(th) tryCatch(sum(integ$loglik(th,state)),error=function(e) -Inf)
    ll <- total(theta)
    if (!is.finite(ll)) stop("Non-finite marginal log-likelihood at the starting values",call.=FALSE)
    S <- colSums(integ$score(theta,state)[[1L]])
    it <- 0L
    repeat {
        if (is.null(H)) {
            H <- gof_polish_hessian(integ,theta,state)
            fresh <- TRUE
        }
        ch <- gof_damped_chol(H)
        if (is.null(R0)) R0 <- ch$R
        step <- drop(chol2inv(ch$R) %*% S)
        dec2 <- sum(S*step)
        dec <- sqrt(max(dec2,0))
        if (dec<control$score.tol)
            return(list(theta=theta,iterations=it,converged=TRUE,H=H,loglik=ll))
        if (sqrt(sum((R0 %*% (theta-theta_start))^2))>max_move)
            return(list(theta=theta,iterations=it,converged=FALSE,refresh=TRUE,H=H,loglik=ll))
        if (!fresh && dec>.25*dec_old) {
            H <- NULL                                     # slow progress: new Hessian
            next
        }
        if (it>=control$maxit) break
        t <- 1
        repeat {
            tn <- theta+t*step
            lln <- total(tn)
            if (is.finite(lln) && lln>=ll+1e-4*t*dec2-1e-12*(1+abs(ll))) break
            t <- t/2
            if (t<1e-12) break
        }
        if (t<1 && !fresh) {
            H <- NULL                                     # stale Hessian: retry fresh
            next
        }
        if (t<1e-12) return(list(theta=theta,iterations=it,converged=FALSE,H=H,loglik=ll,
                                 reason="line search failed"))
        theta <- tn
        ll <- lln
        S <- colSums(integ$score(theta,state)[[1L]])
        dec_old <- dec
        fresh <- FALSE
        it <- it+1L
    }
    list(theta=theta,iterations=it,converged=FALSE,H=H,loglik=ll,reason="maxit reached")
}

## Smallest eigenvalue of each term's Sigma, delta-method SE and z = lambda/SE.
## The SE uses an eigen-decomposition of H so a singular or indefinite
## information along a direction that moves lambda gives SE = Inf (z = 0).
gof_clmm_boundary <- function(problem,theta,H,boundary.tol=1e-6) {
    up <- gof_clmm_unpack(problem,theta)
    e <- eigen((H+t(H))/2,symmetric=TRUE)
    tol <- max(abs(e$values))*1e-12
    rows <- lapply(seq_along(problem$terms),function(t) {
        tm <- problem$terms[[t]]
        re <- up$re[[t]]
        es <- eigen(re$Sigma,symmetric=TRUE)
        lam <- min(es$values)
        v <- es$vectors[,which.min(es$values)]
        g <- numeric(length(theta))
        g[tm$omega] <- vapply(re$dSigma,function(D) drop(crossprod(v,D %*% v)),0)
        proj <- drop(crossprod(e$vectors,g))
        ## SE = Inf only if a material part (> 1%) of the gradient lies in the
        ## (numerically) singular or indefinite directions of H.
        bad <- e$values<=tol & proj^2>1e-4*max(sum(g^2),.Machine$double.xmin)
        se <- if (any(bad)) Inf else sqrt(sum(proj[e$values>tol]^2/e$values[e$values>tol]))
        ## Keep a small covariance eigenvalue separate from failure to
        ## identify its uncertainty using the observed information.
        data.frame(term=tm$label,group=tm$group,lambda_min=lam,se=se,
                   z=if (is.finite(se) && se>0) lam/se else 0,
                   at_boundary=lam<=boundary.tol*max(1,max(es$values)),
                   identified=is.finite(se),
                   stringsAsFactors=FALSE)
    })
    do.call(rbind,rows)
}

gof_boundary_stop <- function(labels,detail="",stage="polish",diagnostics=list()) {
    gof_clmm_abort(paste0("Random-effect covariance at or near its boundary (variance ~ 0, |correlation| ~ 1 or ",
         "singular covariance) for term(s): ",paste(labels,collapse=", "),detail,
         ". The cluster-multiplier approximation is not valid at the boundary; ",
         "simplify the random-effect structure"),"covariance_boundary",stage,diagnostics)
}

## @noRd
## Newton-polish clmm estimates to the root of the adaptive-quadrature score
## with nodes frozen at nu(theta.hat), then run the node, shift and boundary
## acceptance checks. `integrator_factory(nodes)` builds the integrator
## (NULL = default nodes); see `gof_select_integrator()`.
gof_polish <- function(problem,integrator_factory,control=gof_clmm_control()) {
    result <- tryCatch(gof_polish_attempt(problem,integrator_factory,control),gof_clmm_error=identity)
    if (!inherits(result,"gof_clmm_error")) return(result)
    recoverable <- c("covariance_boundary","polish_nonconvergence","quadrature_unstable",
                     "information_unidentified","information_indefinite")
    if (!isTRUE(control$refit) || !result$code %in% recoverable) stop(result)
    gof_retry_polish(problem,integrator_factory,control,result)
}

gof_polish_attempt <- function(problem,integrator_factory,control) {
    snapshot <- new.env(parent=emptyenv())
    tryCatch(gof_polish_impl(problem,integrator_factory,control,snapshot),error=function(e) {
        diagnostics <- snapshot$diagnostics
        if (is.null(diagnostics)) diagnostics <- list()
        if (inherits(e,"gof_clmm_error")) {
            e$diagnostics <- utils::modifyList(diagnostics,e$diagnostics)
            stop(e)
        }
        gof_clmm_abort(conditionMessage(e),"numerical_failure","polish",diagnostics,e)
    })
}

gof_polish_impl <- function(problem,integrator_factory,control,snapshot) {
    on.exit(if (!isTRUE(snapshot$successful))
        snapshot$diagnostics <- gof_polish_snapshot(environment()),add=TRUE)
    if (!inherits(problem,"gof_clmm_problem")) stop("problem must be a gof_clmm_problem",call.=FALSE)
    if (!inherits(control,"gof_clmm_control")) stop("control must come from gof_clmm_control()",call.=FALSE)
    theta0 <- problem$theta
    theta <- theta0
    if (length(problem$source$boundary))
        gof_boundary_stop(problem$source$boundary," (singular covariance in the supplied fit)",
                          stage="input_validation",diagnostics=list(input_covariance=problem$source$covariance))
    integ <- integrator_factory(NULL)
    ## The estimator freezes the nodes: with Q = 1 the log-determinant of the
    ## Laplace approximation is frozen too and the fixed point is a joint-mode
    ## estimate whose variance collapses to zero; Q = 2 does not settle.
    if (!is.null(integ$nodes) && any(integ$nodes<3L))
        stop("Polishing needs at least 3 quadrature nodes per dimension (got ",
             paste(integ$nodes,collapse="/"),"): with frozen nodes, Laplace (1) and 2-node ",
             "rules do not define a usable estimator; increase nodes via gof_clmm_control(nodes=)",
             call.=FALSE)
    refreshes <- 0L
    iterations <- integer()
    delta <- NULL
    converged <- FALSE
    failure <- NULL
    trace <- data.frame(refresh=integer(),newton=integer(),step=numeric(),score=numeric())
    early <- 0L
    H <- NULL
    while (refreshes-early<control$max.refresh) {
        state <- integ$refresh(theta)
        refreshes <- refreshes+1L
        if (!is.null(delta)) {
            S <- colSums(integ$score(theta,state)[[1L]])
            step <- sqrt(sum((gof_damped_chol(H)$R %*% delta)^2))
            std <- gof_standardized_score(S,H)
            trace[nrow(trace),c("step","score")] <- c(step,std)
            if (step<1e-6 && std<1e-6) {
                converged <- TRUE
                break
            }
        }
        nw <- gof_polish_newton(integ,theta,state,control,H=H,max_move=.5)
        H <- nw$H
        iterations <- c(iterations,nw$iterations)
        trace[nrow(trace)+1L,] <- list(refreshes,nw$iterations,NA_real_,NA_real_)
        if (isTRUE(nw$refresh) && early<control$maxit) {
            ## Early refresh after a large move; not counted against max.refresh.
            early <- early+1L
            delta <- nw$theta-theta
            theta <- nw$theta
            next
        }
        if (!nw$converged) {
            failure <- nw$reason
            theta <- nw$theta
            break
        }
        delta <- nw$theta-theta
        theta <- nw$theta
        if (nw$iterations==0L) {
            ## Already a root at freshly refreshed nodes: a fixed point.
            converged <- TRUE
            trace[nrow(trace),c("step","score")] <- c(0,NA_real_)
            break
        }
    }
    if (!converged) {
        ## A run-away towards a covariance boundary is the usual cause.
        H <- tryCatch(gof_polish_hessian(integ,theta,integ$refresh(theta)),error=function(e) NULL)
        if (!is.null(H)) {
            bd <- gof_clmm_boundary(problem,theta,H,control$boundary.tol)
            boundary <- bd
        }
        gof_clmm_abort(paste0("Polishing the clmm estimates did not converge",
             if (!is.null(failure)) paste0(" (",failure,")"),
             " after ",refreshes," node refreshes; increase max.refresh or maxit in gof_clmm_control()"),
             "polish_nonconvergence","polish")
    }
    score <- integ$score(theta,state)
    ## Validation happens below, after the boundary check (at the boundary
    ## the information is singular and the boundary message is the useful one).
    info <- if ("check" %in% names(formals(integ$information)))
        integ$information(theta,state,check=FALSE) else integ$information(theta,state)
    if (is.null(attr(info,"asymmetry"))) attr(info,"asymmetry") <- 0
    loglik <- integ$loglik(theta,state)

    ## Boundary first: at the boundary the information is (near) singular.
    ## Only a real boundary refuses; imprecision (small z) is reported.
    boundary <- gof_clmm_boundary(problem,theta,info,control$boundary.tol)
    if (any(boundary$at_boundary)) gof_boundary_stop(boundary$term[boundary$at_boundary],
        paste0(" (smallest eigenvalue ",paste(signif(boundary$lambda_min[boundary$at_boundary],3),collapse=", "),
               "; boundary.tol = ",control$boundary.tol,")"))
    if (any(!boundary$identified))
        gof_clmm_abort("Random-effect covariance is not identified by the observed information; ordinary inference is unavailable",
                       "information_unidentified","polish")
    weak <- boundary$z<control$boundary.z
    if (any(weak))
        warning("Random-effect covariance imprecisely estimated for term(s) ",
                paste(boundary$term[weak],collapse=", ")," (smallest eigenvalue within ",
                paste(signif(boundary$z[weak],3),collapse=", ")," SE of zero); ",
                "the cluster-multiplier approximation may be less accurate",call.=FALSE)
    asym <- attr(info,"asymmetry")
    if (!is.finite(asym) || asym>1e-4)
        gof_clmm_abort(paste0("Observed information is not symmetric (relative asymmetry ",signif(asym,3),")"),
                       "information_asymmetric","polish")
    Rinfo <- tryCatch(chol(info),error=function(e) NULL)
    if (is.null(Rinfo)) gof_clmm_abort("Observed information at the polished estimate is not positive definite",
                                     "information_indefinite","polish")
    Hinv <- chol2inv(Rinfo)
    dimnames(Hinv) <- dimnames(info)
    S <- colSums(score[[1L]])
    std_score <- sqrt(sum(backsolve(Rinfo,S,transpose=TRUE)^2))

    ## Node check: re-score at Q + 4 nodes per dimension.
    integ4 <- integrator_factory(integ$nodes+4L)
    state4 <- integ4$refresh(theta)
    S4 <- colSums(integ4$score(theta,state4)[[1L]])
    std4 <- sqrt(sum(backsolve(Rinfo,S4,transpose=TRUE)^2))
    node_check <- list(nodes=integ$nodes,check_nodes=integ4$nodes,standardized_score=std4,
                       tolerance=control$node.tol,loglik=sum(integ4$loglik(theta,state4)))
    if (std4>control$node.tol)
        gof_clmm_abort(paste0("Quadrature node check failed: the score at ",paste(integ4$nodes,collapse="/"),
             " nodes per dimension has standardized norm ",signif(std4,3)," > node.tol = ",
             control$node.tol,"; increase nodes via gof_clmm_control(nodes=)"),"quadrature_unstable","polish",
             diagnostics=list(node_check=node_check))

    ## Shift from the supplied (clmm) estimates.
    se <- sqrt(diag(Hinv))
    reference <- problem$source$shift_reference
    if (is.null(reference)) reference <- theta0
    if (anyNA(reference) && is.finite(control$shift.tol))
        gof_clmm_abort("A finite shift.tol cannot be checked for a singular input covariance on the log-Cholesky scale",
                       "parameter_shift_unavailable","polish")
    dl <- theta-reference
    shift <- data.frame(parameter=names(theta),estimate_clmm=unname(reference),estimate=unname(theta),
                        se=unname(se),shift_se=unname(dl/se),stringsAsFactors=FALSE)
    if (!is.null(problem$source$shift_reference)) shift$estimate_start <- unname(theta0)
    ms <- max(abs(shift$shift_se),na.rm=TRUE)
    if (ms>control$shift.tol)
        gof_clmm_abort(paste0("Polished estimates moved ",signif(ms,3)," standard errors from the supplied fit (> shift.tol = ",
             control$shift.tol,"); refit more accurately (e.g. clmm(nAGQ=)) or use ",
             "gof_clmm_control(shift.tol=Inf) to accept"),"parameter_shift","polish",diagnostics=list(shift=shift))
    if (ms>control$shift.warn)
        warning("Polished estimates moved up to ",signif(ms,3)," standard errors from the supplied fit; ",
                "the diagnostic uses the polished estimates (see fit_check$shift)",call.=FALSE)

    dg <- integ$diagnostics(state)
    if (any(dg$expected_information>0))
        warning("The adaptive quadrature scale used the expected information for ",
                sum(dg$expected_information)," random-effect unit(s) whose observed curvature is not ",
                "positive definite at the stationary point (possibly a multimodal posterior); ",
                "quadrature may be inaccurate for them",call.=FALSE)
    snapshot$successful <- TRUE
    structure(list(theta=theta,loglik=sum(loglik),loglik_blocks=loglik,score=score,
                   standardized_score=std_score,information=info,vcov=Hinv,state=state,
                   iterations=iterations,refreshes=refreshes,early_refreshes=early,trace=trace,
                   shift=shift,node_check=node_check,
                   boundary=boundary,integrator=integ$name,nodes=integ$nodes,
                   diagnostics=dg),class="gof_clmm_fit")
}

## @noRd
## Cluster multiplier design for a polished nesting-chain clmm fit: one
## N(0,1) multiplier per top-level block with influence s_c H^{-1}.
gof_clmm_multiplier <- function(problem,fit) {
    if (!inherits(problem,"gof_clmm_problem")) stop("problem must be a gof_clmm_problem",call.=FALSE)
    if (!inherits(fit,"gof_clmm_fit")) stop("fit must come from gof_polish()",call.=FALSE)
    if (!identical(problem$structure,"nesting-chain"))
        stop("Multiplier designs are only available for a single grouping factor or a nesting chain; ",
             "structure '",problem$structure,"' is not supported",call.=FALSE)
    if (length(fit$score)!=1L || !identical(names(fit$score),problem$chain[1L]))
        stop("The fit must carry one score matrix for the top-level factor '",problem$chain[1L],"'",call.=FALSE)
    sc <- fit$score[[1L]]
    if (nrow(sc)!=problem$blocks$m || ncol(sc)!=length(problem$theta))
        stop("Score matrix does not match the problem's blocks and parameters",call.=FALSE)
    gof_multiplier_design(list(list(name=problem$chain[1L],index=problem$blocks$id,
                                    influence=sc %*% solve(fit$information))),
                          n_units=problem$blocks$m,divisor=problem$N/sqrt(problem$blocks$m))
}
