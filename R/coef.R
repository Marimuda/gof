##' @export
coef.cumres <- function(object,...) {
    columns <- list("p-value(Sup)"=object$KS,"p-value(L2)"=object$CvM)
    res <- as.data.frame(columns[!vapply(columns,is.null,logical(1))],check.names=FALSE)
    rownames(res) <- object$variable
    res
}
