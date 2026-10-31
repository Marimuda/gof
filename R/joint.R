#' Joint cumulative residual tests over prespecified ordering variables
#'
#' Generate one multiplier per independent residual unit and reuse it for
#' every ordering. This supplies a joint null law, family tests and a
#' simultaneous null reference envelope on the evaluated grids.
#'
#' @param model A supported `lm` or `glm` fit, or `NULL` for explicit residual,
#'   derivative and influence inputs as described in [cumres.default()].
#' @param ... Arguments passed to [cumres()], including `variable`, `R`,
#'   `plots`, `b`, `subset.max` and, for `model=NULL`, `r`, `dr` and `ic`.
#'   `simulation` is selected internally and must not also be supplied.
#' @param level Reference-envelope probability, strictly between zero and one.
#' @details
#' The observed processes and nuisance corrections are those of Lin, Wei and
#' Ying (2002), equations (2.4), (2.6), and (2.8). In each realization the same
#' multiplier belongs to the same subject under every ordering. Taking a
#' maximum over a fixed, finite set of orderings is a continuous functional of
#' their joint process. The family Sup test uses
#' \eqn{\max_j\sup_t |W_j(t)|}. The family `L2.range` test uses the maximum
#' over orderings of \eqn{\{\int W_j(t)^2dt / (\max t-\min t)\}^{1/2}}.
#' This range normalization makes the family L2 test invariant to positive
#' changes of coordinate units and translations. A one-point grid contributes
#' zero to L2. Marginal `L2` results retain the historical unnormalized norm.
#'
#' Monte Carlo tail estimates are `(1 + number at least as extreme)/(R + 1)`.
#' This avoids zero estimates but does not make the underlying asymptotic
#' approximation exact. The common envelope is `[-critical, critical]` for
#' every selected process and every evaluated point, where `critical` is the
#' order statistic at `ceiling(level * (R + 1))` of the simulated family Sup
#' statistics. It is infinite when the requested probability exceeds the
#' simulation resolution. The envelope describes the joint null residual
#' processes; it is not a confidence interval for fitted means or parameters.
#'
#' The selected orderings, bandwidth and grid rule must be prespecified for
#' formal family inference. A finite grid limits both the tests and envelope
#' to that grid. This family test is not the multivariate lower-orthant omnibus
#' test in equation (2.5) and cannot detect every model departure. Marginal
#' `cumres()` simulation columns from separate orderings cannot be combined
#' retrospectively to obtain this joint distribution. SEM collections and
#' ordinal adapters are not accepted by this entry point.
#'
#' Structurally zero-information processes contribute zero. If any selected
#' process is numerically indeterminate, the whole family's p-values and
#' envelope are withheld; the family is never silently reduced.
#' This also applies to an unresolved fitted ordering under the linear/GLM
#' numerical check. Its marginal curve is descriptive; other resolved marginal
#' tails remain available. Named `resolution` and `on_unresolved` arguments are
#' passed to [cumres.glm()]; declared fixed-coordinate coarsening precedes every
#' family statistic and simulation.
#' The maximum-statistic adjusted marginal tails use the common full-model
#' null. They do not assert strong familywise control for separately defined
#' component hypotheses under arbitrary partial alternatives.
#' @returns An object with classes `cumres_joint` and `cumres`, retaining the
#'   ordinary marginal curves, grids and metadata. Its `joint` component
#'   contains the observed and simulated family statistics, family p-values,
#'   maximum-statistic adjusted marginal p-values, the requested `level`,
#'   `critical`, grid-specific `envelope` data frames, and `status`.
#'   `as.data.frame()` includes marginal and family rows distinguished by
#'   `scope`. `summary()` and optional `tidy()` retain these rows. `plot()`
#'   displays each process with the common Sup reference envelope.
#'   `coef()` retains the marginal p-value table; family p-values are in
#'   `joint$p.value` and the rows whose `scope` is `"family"`.
#' @references D. Y. Lin, L. J. Wei and Z. Ying (2002).
#'   Model-Checking Techniques Based on Cumulative Residuals.
#'   Biometrics 58, 1--12. \doi{10.1111/j.0006-341X.2002.00001.x}.
#'   The finite-family maximum and range-normalized L2 aggregation here are
#'   derived extensions of that joint multiplier process.
#' @seealso [cumres()], [cumres.default()], [cumres.glm()]
#' @examples
#' set.seed(51)
#' d <- data.frame(x = seq(-2, 2, length.out = 120), z = rnorm(120))
#' d$y <- 1 + d$x + 0.3*d$z + rnorm(120)
#' diagnostic <- lm(y ~ x + z, data = d) |>
#'   cumres_joint(variable = c("x", "z"), R = 99, plots = 3)
#' as.data.frame(diagnostic)
#' diagnostic$joint$critical
#' plot(diagnostic, idx = "x")
#' @export
cumres_joint <- function(model, ..., level=0.95) {
    if (!is.null(model) && (!inherits(model,"lm") || inherits(model,c("ordreg","lvmfit"))))
        stop("Joint ordering inference supports lm/glm fits or explicit inputs with model=NULL",call.=FALSE)
    dots <- list(...)
    if ("simulation" %in% names(dots))
        stop("cumres_joint selects simulation='joint'; do not also supply simulation",call.=FALSE)
    gof_joint_level(level)
    result <- cumres(model,...,simulation="joint")
    gof_joint_finalize(result,level)
}

gof_joint_level <- function(level) {
    if (!is.numeric(level) || is.complex(level) || length(level)!=1L ||
        !is.finite(level) || level<=0 || level>=1)
        stop("level must be a finite probability strictly between zero and one",call.=FALSE)
    invisible(level)
}

gof_joint_engine <- function(model,variable,r,dr,ic,R,b,plots,subset,subset.max,vnames) {
    n <- length(r)
    p <- ncol(variable)
    q <- ncol(dr)
    orders <- lapply(seq_len(p),function(j) order(variable[,j]))
    if (b==0) {
        idx <- lapply(seq_len(p),function(j) gof_grid(variable[,j],subset.max,subset))
        x <- lapply(seq_len(p),function(j) variable[orders[[j]][idx[[j]]],j])
        project <- function(w,j) w[idx[[j]]]
    } else {
        x <- lapply(seq_len(p),function(j) gof_window_grid(variable[,j],b,subset.max,subset))
        upper <- lapply(seq_len(p),function(j) findInterval(x[[j]],variable[orders[[j]],j]))
        lower <- lapply(seq_len(p),function(j) findInterval(x[[j]],variable[orders[[j]],j]+b))
        project <- function(w,j) {
            prefix <- c(0,w)
            prefix[upper[[j]]+1L]-prefix[lower[[j]]+1L]
        }
    }
    derivatives <- lapply(orders,function(o) {
        if (!q) return(matrix(numeric(),n,0L))
        matrix(apply(dr[o,,drop=FALSE],2L,cumsum),nrow=n,ncol=q)
    })
    W <- lapply(seq_len(p),function(j) project(cumsum(r[orders[[j]]])/sqrt(n),j))
    realization <- function() {
        g <- rnorm(n)
        influence <- if (q) drop(crossprod(ic,g)) else numeric()
        lapply(seq_len(p),function(j) {
            # Crucially, g is reordered with the SUBJECT, not assigned to its rank.
            prefix <- cumsum((r*g)[orders[[j]]])
            if (q) prefix <- prefix+drop(derivatives[[j]]%*%influence)
            project(prefix/sqrt(n),j)
        })
    }
    statistics <- function(w) unlist(lapply(seq_len(p),function(j)
        c(SupTest(w[[j]]),L2Test(w[[j]],x[[j]]))),use.names=FALSE)
    null <- t(vapply(seq_len(R),function(k) statistics(realization()),numeric(2L*p)))
    What <- lapply(x,function(t) matrix(0,length(t),plots))
    for (k in seq_len(plots)) {
        draw <- realization()
        for (j in seq_len(p)) What[[j]][,k] <- draw[[j]]
    }
    result <- structure(list(W=W,What=What,x=x,KS=rep(NA_real_,p),CvM=rep(NA_real_,p),
        R=R,n=n,sd=NULL,cvalues=NULL,variable=vnames,type=rep("residual",p),
        untie=rep(0,n),model=class(model)[1],joint=list(statistics=null),
        design_info=gof_design_info(correction_dimension=q)),class=c("cumres_joint","cumres"))
    gof_joint_finalize(result)
}

gof_joint_finalize <- function(res,level=0.95) {
    gof_joint_level(level)
    p <- length(res$W)
    sup <- seq.int(1L,2L*p,2L)
    l2 <- sup+1L
    null <- res$joint$statistics
    if (!is.matrix(null) || nrow(null)!=res$R || ncol(null)!=2L*p)
        stop("Joint simulation statistics are missing or inconsistent",call.=FALSE)
    zero <- which(res$informative %in% FALSE)
    for (j in zero) null[,c(sup[j],l2[j])] <- 0
    observed <- unlist(lapply(seq_len(p),function(j)
        c(SupTest(res$W[[j]]),L2Test(res$W[[j]],res$x[[j]]))),use.names=FALSE)
    # range normalization is defined as zero on one-point grids, whose L2 is zero.
    ranges <- vapply(res$x,function(t) tail(t,1L)-t[1L],numeric(1))
    if (any(!is.finite(ranges)))
        stop("Ordering-variable ranges must be finite; rescale the coordinates",call.=FALSE)
    normalizer <- sqrt(ranges)
    normalizer[normalizer==0] <- 1
    normalized <- sweep(null[,l2,drop=FALSE],2L,normalizer,"/")
    family <- c(Sup=max(observed[sup]),L2.range=max(observed[l2]/normalizer))
    simulated <- cbind(Sup=apply(null[,sup,drop=FALSE],1L,max),
                       L2.range=apply(normalized,1L,max))
    tail_estimate <- function(draws,threshold) (1+sum(draws>=threshold))/(nrow(null)+1)
    marginal <- vapply(seq_len(2L*p),function(j) tail_estimate(null[,j],observed[j]),numeric(1))
    res$KS <- marginal[sup]
    res$CvM <- marginal[l2]
    unavailable <- if (is.null(res$informative)) rep(FALSE,p) else is.na(res$informative)
    unavailable <- unavailable | gof_ordering_unresolved(res)
    indeterminate <- any(unavailable)
    res$KS[unavailable] <- res$CvM[unavailable] <- NA_real_
    probability <- vapply(seq_len(2L),function(j) tail_estimate(simulated[,j],family[j]),numeric(1))
    names(probability) <- names(family)
    adjusted <- cbind(Sup=vapply(seq_len(p),function(j)
        tail_estimate(simulated[,1L],observed[sup[j]]),numeric(1)),
        L2.range=vapply(seq_len(p),function(j)
        tail_estimate(simulated[,2L],observed[l2[j]]/normalizer[j]),numeric(1)))
    rownames(adjusted) <- res$variable
    k <- ceiling(level*(res$R+1))
    critical <- if (k>res$R) Inf else sort(simulated[,1L],partial=k)[k]
    if (indeterminate) {
        probability[] <- NA_real_
        adjusted[] <- NA_real_
        critical <- NA_real_
    }
    envelope <- lapply(res$x,function(t) data.frame(x=t,lower=rep(-critical,length(t)),
                                                   upper=rep(critical,length(t))))
    names(envelope) <- res$variable
    res$joint <- list(statistics=null,observed=family,simulated=simulated,p.value=probability,
        adjusted.p.value=adjusted,level=level,critical=critical,envelope=envelope,
        status=if (indeterminate) "indeterminate" else if (length(zero)==p) "uninformative" else "ok",
        scope="prespecified ordering variables and evaluated grids",
        L2.normalization="coordinate range")
    class(res) <- unique(c("cumres_joint",class(res)))
    res
}

#' @rdname cumres_joint
#' @param x,object A `cumres_joint` diagnostic or its summary.
#' @param row.names,optional Arguments accepted by [base::as.data.frame()].
#' @export
as.data.frame.cumres_joint <- function(x,row.names=NULL,optional=FALSE,...) {
    marginal <- gof_result_table(x)
    marginal$scope <- "marginal"
    family <- marginal[rep(1L,2L),,drop=FALSE]
    family$variable <- "<family>"
    family$test <- names(x$joint$observed)
    family$statistic <- unname(x$joint$observed)
    family$p.value <- unname(x$joint$p.value)
    family$informative <- if (is.null(x$informative) || identical(x$joint$status,"indeterminate")) NA else
        !identical(x$joint$status,"uninformative")
    family$aggregation <- "maximum across orderings"
    family$scope <- "family"
    if ("ordering_status" %in% names(family)) {
        unresolved <- gof_ordering_unresolved(x)
        family$ordering_source <- "family"
        family$ordering_status <- if (any(unresolved)) "unresolved" else "family"
        family$ordering_reason <- if (any(unresolved)) paste0("Unresolved ordering: ",
            paste(x$variable[unresolved],collapse=", ")) else NA_character_
        for (field in c("resolution","distinct_values","tie_blocks","tied_units","changed_units"))
            family[[field]] <- NA_real_
    }
    result <- rbind(marginal,family)
    rownames(result) <- if (is.null(row.names)) NULL else row.names
    result
}

#' @rdname cumres_joint
#' @export
summary.cumres_joint <- function(object,...) {
    result <- summary.cumres(object,...)
    result$joint <- object$joint
    class(result) <- c("summary.cumres_joint","summary.cumres")
    result
}

#' @rdname cumres_joint
#' @export
print.cumres_joint <- function(x,...) {
    print(as.data.frame(x),row.names=FALSE,...)
    cat("\nJoint subject-level multiplier simulation; ",x$R," realizations.\n",sep="")
    cat("Family scope: ",x$joint$scope,".\n",sep="")
    cat("Simultaneous null reference envelope level: ",x$joint$level,".\n",sep="")
    if (identical(x$joint$status,"indeterminate"))
        cat("Family inference withheld: a selected process is numerically indeterminate.\n")
    invisible(x)
}

#' @rdname cumres_joint
#' @export
print.summary.cumres_joint <- function(x,...) {
    print(x$table,row.names=FALSE,...)
    cat("\nJoint subject-level multiplier simulation; ",x$R," realizations.\n",sep="")
    if (identical(x$joint$status,"indeterminate"))
        cat("Family inference withheld: a selected process is numerically indeterminate.\n")
    invisible(x)
}

#' @rdname cumres_joint
#' @param idx Numeric process indices or ordering-variable names.
#' @param envelope Whether to draw the simultaneous Sup null reference envelope.
#' @param envelope.col,envelope.lty Colour and line type for the envelope.
#' @export
plot.cumres_joint <- function(x,idx=seq_along(x$W),envelope=TRUE,
                              envelope.col="steelblue",envelope.lty=2,...) {
    if (!is.logical(envelope) || length(envelope)!=1L || is.na(envelope))
        stop("envelope must be TRUE or FALSE",call.=FALSE)
    if (is.character(idx)) {
        idx <- match(idx,x$variable)
        if (anyNA(idx)) stop("Unknown process name in idx",call.=FALSE)
    }
    if (!is.numeric(idx) || any(!is.finite(idx)) ||
        any(idx!=floor(idx) | idx<1 | idx>length(x$W))) stop("Invalid process index",call.=FALSE)
    dots <- list(...)
    if (envelope && !is.finite(x$joint$critical))
        warning("No finite simultaneous envelope is available at this level",call.=FALSE)
    for (j in idx) {
        arguments <- c(list(x=x,idx=j),dots)
        if (is.null(dots$ylim) && envelope && is.finite(x$joint$critical)) {
            peak <- max(abs(x$W[[j]]),x$joint$critical)
            arguments$ylim <- c(-1,1)*if (peak>0) 1.1*peak else 1
        }
        do.call(plot.cumres,arguments)
        if (envelope && is.finite(x$joint$critical))
            graphics::abline(h=c(-x$joint$critical,x$joint$critical),
                             col=envelope.col,lty=envelope.lty)
    }
    invisible(x)
}
