#' Work with cumulative residual diagnostics
#'
#' Extract test results for reporting or further analysis. Base R pipelines
#' work directly with `cumres()`, `summary()`, `as.data.frame()` and `plot()`.
#'
#' @name cumres-object
#' @aliases cumres-object
#' @param object,x A `cumres` diagnostic or, for the summary methods, an object
#'   returned by `summary()` for that diagnostic.
#' @param row.names Optional row names for the returned data frame.
#' @param optional Logical argument accepted for compatibility with
#'   [base::as.data.frame()]. Column names are fixed by the result contract.
#' @param ... Further arguments. The data-frame and tidying methods reserve
#'   these for future use; print methods pass applicable arguments to printing.
#' @returns `as.data.frame()` returns a data frame with two rows per ordering
#'   variable, in variable order and then Sup/L2 order, and these columns:
#'   \describe{
#'     \item{variable}{The ordering variable or process label.}
#'     \item{test}{`"Sup"` or `"L2"`.}
#'     \item{statistic}{Observed test statistic, calculated from the returned
#'       process and evaluation grid. Missing for saved summaries without curves
#'       or an extracted result table.}
#'     \item{p.value}{Monte Carlo tail estimate. Missing when unavailable or
#'       when the model adapter reports numerical indeterminacy.}
#'     \item{n}{Number of independent residual units. Binomial trial expansion
#'       can make this larger than the number of fitted rows.}
#'     \item{nobs}{Positive-weight fitted rows when reported by the adapter;
#'       otherwise missing.}
#'     \item{R}{Number of null simulations used for the p-values.}
#'     \item{informative}{The adapter's logical information flag. `FALSE`
#'       denotes structurally recognized zero-information processes. `NA`
#'       denotes either unavailable metadata or numerical indeterminacy; only
#'       the latter also withholds otherwise available p-values.}
#'     \item{aggregation}{Ordinal aggregation mode (`"omnibus"`, `"score"`,
#'       `"cutpoints"` or `"custom"`); missing for other model adapters.}
#'   }
#'   `summary()` returns a `summary.cumres` object with the existing p-values
#'   and metadata plus a `table` containing this data frame. `coef()` retains
#'   the historical wide p-value table, with process names as row names; these
#'   entries are diagnostic p-values, not fitted regression coefficients.
#'   `print()` returns its input invisibly. Optional `generics::tidy()` and
#'   `broom::tidy()` methods return the same result columns as a tibble.
#' @details Extraction consumes no random draws and does not rerun simulations.
#'   The additional endpoints retain their interpretation in extracted tables:
#'   [cumres_joint()] adds marginal/family `scope`; [cumres_normality()] adds
#'   calibration, diagnostic and residual labels; [cumres_bootstrap()] adds
#'   calibration/null labels, unresolved-refit counts and simulated tail bounds.
#'   Their method-specific help defines these additional columns.
#'   The data-frame methods also accept saved matrix-based diagnostic objects.
#'   A saved summary made before result tables were added retains its p-values
#'   but cannot reconstruct observed statistics without the original curves.
#'
#'   A large p-value is not proof of model adequacy. The Sup and L2 tests use
#'   different summaries of the same process; neither supplies an independent
#'   replication. See [cumres()] and the model-specific help for assumptions.
#' @seealso [cumres()], [plot.cumres()]
#' @examples
#' set.seed(42)
#' d <- data.frame(x = seq(-2, 2, length.out = 120))
#' d$y <- 1 + d$x + rnorm(nrow(d))
#' g <- lm(y ~ x, data = d) |> cumres(variable = "x", R = 99, plots = 5)
#' g |> as.data.frame()
#' g |> summary() |> as.data.frame()
#' coef(g)
#' if (requireNamespace("generics", quietly = TRUE) &&
#'     requireNamespace("tibble", quietly = TRUE)) {
#'   generics::tidy(g)
#' }
NULL

gof_result_table <- function(object) {
    variables <- object$variable
    p <- length(variables)
    result <- data.frame(variable=rep(as.character(variables),each=2L),
                         test=rep(c("Sup","L2"),p),
                         statistic=rep(NA_real_,2L*p),
                         p.value=rep(NA_real_,2L*p),
                         n=rep(if (is.null(object$n)) NA_real_ else object$n,2L*p),
                         nobs=rep(if (is.null(object$nobs)) NA_real_ else object$nobs,2L*p),
                         R=rep(if (is.null(object$R)) NA_real_ else object$R,2L*p),
                         informative=rep(NA,p*2L),
                         aggregation=rep(if (is.null(object$aggregation)) NA_character_ else object$aggregation,2L*p),
                         stringsAsFactors=FALSE)
    if (!p) return(result)
    if (!is.null(object$KS)) result$p.value[seq.int(1L,2L*p,2L)] <- object$KS
    if (!is.null(object$CvM)) result$p.value[seq.int(2L,2L*p,2L)] <- object$CvM
    if (!is.null(object$informative)) result$informative <- rep(object$informative,each=2L)
    columns <- function(z) if (is.matrix(z)) lapply(seq_len(ncol(z)),function(i) z[,i]) else z
    observed <- columns(object$W)
    grid <- columns(object$x)
    if (!is.null(observed) && !is.null(grid)) {
        for (i in seq_len(p)) {
            result$statistic[2L*i-1L] <- SupTest(observed[[i]])
            result$statistic[2L*i] <- L2Test(observed[[i]],grid[[i]])
        }
    }
    result
}

#' @rdname cumres-object
#' @export
as.data.frame.cumres <- function(x,row.names=NULL,optional=FALSE,...) {
    result <- gof_result_table(x)
    if (!is.null(row.names)) rownames(result) <- row.names
    result
}

#' @rdname cumres-object
#' @export
as.data.frame.summary.cumres <- function(x,row.names=NULL,optional=FALSE,...) {
    result <- if (is.null(x$table)) gof_result_table(x) else x$table
    if (!is.null(row.names)) rownames(result) <- row.names
    result
}

#' @rdname cumres-object
#' @exportS3Method generics::tidy
tidy.cumres <- function(x,...) {
    if (!requireNamespace("tibble",quietly=TRUE))
        stop("Install the optional 'tibble' package to use tidy()",call.=FALSE)
    tibble::as_tibble(as.data.frame(x))
}

#' @rdname cumres-object
#' @exportS3Method generics::tidy
tidy.summary.cumres <- function(x,...) {
    tidy.cumres(x,...)
}
