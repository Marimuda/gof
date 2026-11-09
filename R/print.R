##' @rdname cumres-object
##' @export
print.cumres <- function(x,...) {
    cat("\n")
    print(coef(x), ...)
    cat("\nBased on ", x$R, " realizations.\n",sep="")
    gof_print_design(x)
    if (!is.null(x$aggregation)) cat("Ordinal aggregation: ",x$aggregation,".\n",sep="")
    if (anyNA(x$informative))
        cat("Numerically indeterminate score process: ",
            paste(x$variable[is.na(x$informative)],collapse=", "),".\n",sep="")
    if (any(gof_ordering_unresolved(x)))
        cat("Fitted ordering unresolved: ",paste(x$variable[gof_ordering_unresolved(x)],collapse=", "),
            ". Curves are descriptive; inferential p-values are unavailable.\n",sep="")
    invisible(x)
}

##' @rdname cumres-object
##' @export
summary.cumres <- function(object,...) {
    res <- with(object, list(KS=KS, CvM=CvM, R=R, n=n, type=type, model=model, variable=variable))
    res$nobs <- object$nobs
    res$response <- object$response
    res$informative <- object$informative
    res$aggregation <- object$aggregation
    res$ordering <- object$ordering
    res$ordering_check <- object$ordering_check
    res$design_info <- object$design_info
    res$multiplier <- object$multiplier
    res$random_effects <- object$random_effects
    res$fit_check <- object$fit_check
    res$table <- as.data.frame(object)
    return(structure(res, class="summary.cumres"))
}

##' @rdname cumres-object
##' @export
print.summary.cumres <- function(x,...) {
    cat("\n");
    gof_print_design(x)
    for (i in seq_along(x$variable)) {
        if (!is.null(x$KS)) cat("Sup-test: p-value=", x$KS[i], "\n", sep="")
        if (!is.null(x$CvM)) cat("L2-test: p-value=", x$CvM[i], "\n", sep="")
        if (!is.null(x$aggregation)) cat("Ordinal aggregation: ",x$aggregation,".\n",sep="")
        if (!is.null(x$informative) && is.na(x$informative[i]))
            cat("Numerically indeterminate score process.\n")
        if (gof_ordering_unresolved(x)[i])
            cat("Fitted ordering unresolved; curve is descriptive and inferential p-values are unavailable.\n")
        cat("Based on ", x$R, " realizations. Cumulated residuals ordered by '", x$variable[i], "'.\n", sep="")
        cat("---\n");
    }
    invisible(x)
}


