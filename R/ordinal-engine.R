## Shared multiplier engine for ordinal cumulative residual processes.
##
## Extracted from cumres.ordreg (gof 1.3.0) so that fixed-effect ordinal fits
## (one multiplier per row) and mixed-model fits (one multiplier per
## independent block) use the same sort, tie grid, cumulative gradients,
## aggregation, statistics and result layout. With a single "rows" factor the
## engine performs exactly the floating-point operations and random draws of
## the 1.3.0 method.

## Map the user's aggregate argument to the aggregation function and mode.
gof_ordinal_aggregation <- function(aggregate) {
    if (is.character(aggregate)) {
        aggregation <- match.arg(aggregate,c("omnibus","score","cutpoints"))
        aggregate <- switch(aggregation,omnibus=gof_ordinal_norm,score=sum,cutpoints=identity)
    } else {
        if (!is.function(aggregate)) stop("aggregate must be a function or a named aggregation mode")
        aggregation <- if (identical(aggregate,sum)) "score" else
            if (identical(aggregate,identity)) "cutpoints" else "custom"
    }
    list(fn=aggregate,aggregation=aggregation)
}

## Reject linear/GLM simulation controls passed through `...`. Errors are
## reported against the calling method, as before the extraction.
gof_ordinal_check_dots <- function(...) {
    call <- sys.call(-1L)
    if (any(c("resolution","on_unresolved") %in% names(list(...))))
        stop(simpleError("resolution and on_unresolved are not implemented for ordinal diagnostics"))
    if ("simulation" %in% names(list(...)))
        stop(simpleError("Ordinal diagnostics use shared cutpoint multipliers for one ordering; simulation is a default/linear/GLM option",call))
    invisible(NULL)
}

gof_ordinal_norm <- function(x) {
    scale <- max(abs(x))
    if (scale==0) return(0)
    scale*sqrt(sum((x/scale)^2))
}

## @noRd
## Multiplier simulation of an ordinal cumulative residual process.
##
## @param residual n x k matrix of cumulative-cutpoint residuals.
## @param gradient (n*k) x d matrix; row (j-1)*n+i is the derivative of
##   residual[i,j] with respect to the parameters.
## @param multiplier A gof_multiplier_design over the same n rows and d
##   parameter columns. Cumulative sums are divided by multiplier$divisor.
## @param ordering Finite ordering values, one per row.
## @param aggregate,aggregation Function and mode from gof_ordinal_aggregation().
## @param informative TRUE, FALSE (processes set to zero, p-values one) or NA
##   (curves retained, missing p-values).
## @param R,plots Validated simulation and retained-curve counts.
## @param variable,model Labels stored on the result.
## @param design_info One-row gof_design_info() data frame.
## @param metadata Named list of further result fields.
## @returns A cumres object. The caller plots it if requested.
gof_ordinal_engine <- function(residual,gradient,multiplier,ordering,aggregate,aggregation,
                               informative,R,plots,variable,model,design_info,metadata=list()) {
    if (!inherits(multiplier,"gof_multiplier_design"))
        stop("multiplier must be a gof_multiplier_design",call.=FALSE)
    if (!is.matrix(residual) || !is.matrix(gradient)) stop("Ordinal residuals and gradients must be matrices",call.=FALSE)
    n <- nrow(residual)
    k <- ncol(residual)
    if (!n || !k || length(ordering)!=n || nrow(gradient)!=n*k || multiplier$n!=n ||
        ncol(gradient)!=multiplier$d)
        stop("Incompatible ordinal residual, gradient, ordering or multiplier dimensions",call.=FALSE)
    if (!is.list(metadata) || (length(metadata) && (is.null(names(metadata)) || any(!nzchar(names(metadata))))))
        stop("metadata must be a named list",call.=FALSE)
    divisor <- multiplier$divisor
    x <- ordering
    o <- order(x)
    x <- x[o]
    ends <- cumsum(rle(x)$lengths)
    x <- x[ends]
    residual <- residual[o,,drop=FALSE]
    eta <- lapply(seq_len(k),function(j) {
        z <- gradient[seq_len(n)+(j-1L)*n,,drop=FALSE][o,,drop=FALSE]
        apply(z,2,cumsum)
    })
    ## Levels are renumbered by first appearance in the sorted rows, so draws
    ## do not depend on the labelling of levels. For a "rows" factor with
    ## index 1:n this is the identity on sorted rows and IC = influence[o,].
    factors <- lapply(multiplier$factors,function(f) {
        index <- f$index[o]
        first <- unique(index)
        list(idx=match(index,first),influence=f$influence[first,,drop=FALSE],
             loading=if (!is.null(f$loading)) f$loading[o])
    })
    if (any(!is.finite(residual)) || any(!is.finite(gradient)) ||
        any(!vapply(factors,function(f) all(is.finite(f$influence)),logical(1))))
        stop("Ordinal residuals, gradients and influences must be finite")
    levels <- vapply(factors,function(f) nrow(f$influence),integer(1))
    idx <- lapply(factors,`[[`,"idx")
    loading <- lapply(factors,`[[`,"loading")
    ic <- lapply(factors,`[[`,"influence")
    draw <- gof_multiplier_rules$draw[[multiplier$draw]]
    rowmult <- gof_multiplier_rules$combine[[multiplier$combine]]
    combine <- function(w) {
        rows <- lapply(seq_len(nrow(w)),function(i) aggregate(w[i,]))
        sizes <- lengths(rows)
        if (!length(sizes) || !sizes[1] || any(sizes!=sizes[1]) ||
            any(!vapply(rows,function(z) is.numeric(z) && !is.complex(z) && is.null(dim(z)),logical(1))) ||
            any(!is.finite(unlist(rows))))
            stop("aggregate must return a fixed-length finite numeric vector")
        do.call(rbind,rows)
    }
    observed_cdf <- apply(residual,2,cumsum)[ends,,drop=FALSE]/divisor
    if (identical(informative,FALSE)) observed_cdf[] <- 0
    observed <- combine(observed_cdf)
    q <- ncol(observed)
    What <- lapply(seq_len(q),function(i) matrix(0,length(ends),plots))
    obs.ks <- apply(abs(observed),2,max)
    obs.l2 <- apply(observed,2,function(w) L2Test(w,x))
    ks <- l2 <- matrix(0,R,q)
    for (j in seq_len(max(R,plots))) {
        g <- draw(levels)
        u <- rowmult(g,idx,loading)
        B <- colSums(ic[[1L]]*g[[1L]])
        for (f in seq_along(g)[-1L]) B <- B+colSums(ic[[f]]*g[[f]])
        w <- vapply(seq_len(k),function(i)
            (cumsum(residual[,i]*u)+drop(eta[[i]] %*% B))[ends]/divisor,
            numeric(length(ends)))
        w <- matrix(w,ncol=k)
        if (identical(informative,FALSE)) w[] <- 0
        w <- combine(w)
        if (ncol(w)!=q) stop("aggregate returned a different number of processes")
        if (j<=R) {
            ks[j,] <- apply(abs(w),2,max)
            l2[j,] <- apply(w,2,function(z) L2Test(z,x))
        }
        if (j<=plots) for (i in seq_len(q)) What[[i]][,j] <- w[,i]
    }
    labels <- if (q==1L) variable else paste(variable,seq_len(q),sep=".")
    result <- structure(list(W=lapply(seq_len(q),function(i) observed[,i]),
                             What=What,x=rep(list(x),q),
                             KS=colMeans(sweep(ks,2,obs.ks,">=")),
                             CvM=colMeans(sweep(l2,2,obs.l2,">=")),
                             n=multiplier$n_units,R=R,
                             type=rep(if (aggregation=="omnibus") "ordinal-norm" else "normal",q),
                             aggregation=aggregation,informative=rep(informative,q),sd=NULL,cvalues=NULL,
                             variable=labels,model=model,design_info=design_info),class="cumres")
    if (is.na(informative)) result$KS[] <- result$CvM[] <- NA_real_
    result$multiplier <- gof_multiplier_summary(multiplier)
    if (length(metadata)) result[names(metadata)] <- metadata
    result
}
