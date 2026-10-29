##' @rdname cumres-object
##' @export
coef.cumres <- function(object,...) {
    columns <- list("p-value(Sup)"=object$KS,"p-value(L2)"=object$CvM)
    if (all(vapply(columns,is.null,logical(1)))) return(data.frame(row.names=object$variable))
    res <- as.data.frame(columns[!vapply(columns,is.null,logical(1))],check.names=FALSE)
    rownames(res) <- object$variable
    res
}
