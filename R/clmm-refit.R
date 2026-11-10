## Bounded, deterministic initialization retries for the same integrated
## estimating equation. Every returned root passes the ordinary checks.
## This is optional numerical recovery, not boundary-valid inference.
gof_retry_polish <- function(problem,factory,control,original) {
    if (length(problem$source$boundary) && is.finite(control$shift.tol))
        gof_clmm_abort("A finite shift.tol cannot be checked for a singular input covariance on the log-Cholesky scale",
                       "parameter_shift_unavailable","polish",diagnostics=list(input_covariance=problem$source$covariance))
    base <- attr(factory,"nodes")
    if (is.null(base)) stop(original)
    groups <- vapply(problem$terms,`[[`,"","group")
    q <- vapply(problem$chain,function(g) sum(vapply(problem$terms[groups==g],`[[`,0L,"q")),0L)
    initial <- tryCatch(gof_clmm_unpack(problem,problem$theta),error=function(e) NULL)
    if (is.null(initial)) stop(original)
    attempts <- list(list(strategy="input",condition=original))
    for (retry in 1:2) {
        working <- problem
        working$source$boundary <- character()
        reference <- problem$theta
        if (length(problem$source$boundary)) for (tm in problem$terms)
            if (tm$label %in% problem$source$boundary) reference[tm$omega] <- NA_real_
        working$source$shift_reference <- reference
        if (retry==2L) {
            frequency <- vapply(seq_len(problem$K),function(k) sum(problem$weights[problem$y==k]),0)
            probability <- cumsum(frequency)[seq_len(problem$K-1L)]/sum(frequency)
            cutpoints <- switch(problem$link,logit=stats::qlogis(probability),probit=stats::qnorm(probability),
                cauchit=stats::qcauchy(probability),cloglog=log(-log1p(-probability)),loglog=-log(-log(probability)))
            xi <- tryCatch(qr.solve(problem$tJac,cutpoints),error=function(e) NULL)
            if (!is.null(xi) && all(is.finite(xi)) && all(diff(drop(problem$tJac %*% xi))>0))
                working$theta[problem$index$xi] <- xi
            working$theta[problem$index$beta] <- 0
        }
        valid <- TRUE
        for (i in seq_along(problem$terms)) {
            tm <- problem$terms[[i]]
            scale <- colSums(tm$Z^2*(problem$weights/problem$N))
            if (any(!is.finite(scale)) || any(scale<=.Machine$double.eps)) {valid<-FALSE;break}
            S <- diag(.25/scale,tm$q)
            if (retry==1L) S <- S+initial$re[[i]]$Sigma
            working$theta[tm$omega] <- gof_logchol_omega(S)
        }
        if (!valid) break
        nodes <- if (!is.null(control$nodes)) base else pmin(base+4L*retry,ifelse(q>=3L,9L,21L))
        # Bound the largest row/node expansion, including the Q+4 check.
        budget <- if (is.null(control$max.refit.pairs)) 5e6 else control$max.refit.pairs
        log_pairs <- log(problem$n)+sum(q*log(nodes+4L))
        if (!is.finite(log_pairs) || log_pairs>log(budget)) {
            attempts[[length(attempts)+1L]] <- list(strategy="budget",nodes=nodes,
                condition="Refit quadrature row/node expansion exceeds max.refit.pairs")
            next
        }
        refactory <- function(override=NULL) factory(if (is.null(override)) nodes else override)
        result <- tryCatch(gof_polish_attempt(working,refactory,control),gof_clmm_error=identity)
        strategy <- if (retry==1L) "interior ridge" else "marginal cutpoints and diagonal covariance"
        entry <- list(strategy=strategy,nodes=nodes,theta_start=working$theta)
        if (inherits(result,"gof_clmm_error")) {
            entry$condition <- result
            attempts[[length(attempts)+1L]] <- entry
        } else {
            entry$loglik <- result$loglik
            attempts[[length(attempts)+1L]] <- entry
            result$initialization <- list(strategy=strategy,attempts=attempts,
                input_covariance=problem$source$covariance,shift_reference=reference)
            warning("An interior score root was obtained after numerical refitting; inspect fit_check$initialization and fit_check$shift",
                    call.=FALSE)
            return(result)
        }
    }
    gof_clmm_abort("Numerical refitting did not produce a checked interior score root; ordinary inference is unavailable",
                   "refit_failed","polish",diagnostics=list(attempts=attempts),parent=original)
}
