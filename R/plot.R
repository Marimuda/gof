##' Plot observed and simulated cumulative residual processes
##'
##' Display the observed residual process together with saved realizations
##' under the null. These curves provide a visual comparison and are not
##' simultaneous confidence bands.
##'
##' @param x An object produced by \code{cumres}.
##' @param idx Numeric process indices or character labels from
##' \code{x$variable}. An empty vector draws no processes.
##' @param col Colour or vector of colours for saved null processes. Use
##' \code{"none"}, \code{NULL} or an empty vector to suppress simulated lines.
##' @param col.alpha Transparency of saved null-process lines, between zero
##' and one. The observed curve uses the current graphical foreground colour.
##' @param legend Legend style: \code{"type1"} displays available p-values
##' and information-status annotations; \code{"type2"} identifies observed
##' and simulated lines; \code{"none"} suppresses the legend.
##' @param xlab Optional horizontal-axis label.
##' @param ylab Optional vertical-axis label.
##' @param vs Logical; use the ordering-variable coordinates when \code{TRUE}
##' (the default), or sequential evaluation-point indices when \code{FALSE}.
##' @param ylim Optional numeric vertical-axis limits.
##' @param title Optional plot title, overriding a supplied \code{main}.
##' @param ... Named graphical arguments passed to \code{graphics::plot},
##' such as \code{xlim}, \code{axes}, \code{cex.axis} or \code{main}.
##' @returns The original \code{cumres} object, invisibly, so plotting can
##' appear within a pipeline without replacing its diagnostic result.
##' @seealso \code{\link{cumres}}.
##' @author Klaus K. Holst
##' @keywords hplot regression
##' @examples
##' set.seed(16)
##' d <- data.frame(x = seq(-2, 2, length.out = 100))
##' d$y <- 1 + d$x + rnorm(nrow(d))
##' fit <- lm(y ~ x, data = d)
##' diagnostic <- cumres(fit, variable = "x", R = 99, plots = 5)
##' plot(diagnostic, idx = "x", legend = "type2")
##' plot(diagnostic, col = NULL, legend = "none", vs = FALSE)
##'
##' # Plotting returns the diagnostic invisibly.
##' invisible(diagnostic |> plot(legend = "type1"))
##' @method plot cumres
##' @export
plot.cumres <- function(x,idx=seq_along(x$W),col="purple",col.alpha=0.3,
                         legend=c("type1","type2","none"),xlab,ylab,vs=TRUE,
                         ylim=NULL,title,...) {
    original <- x
    # Accept saved pre-1.0.4 matrix-based objects as well as current lists.
    columns <- function(z) if (is.matrix(z)) lapply(seq_len(ncol(z)),function(i) z[,i]) else z
    x$W <- columns(x$W)
    x$x <- columns(x$x)
    if (is.character(idx)) {
        idx <- match(idx,x$variable)
        if (anyNA(idx)) stop("Unknown process name in idx")
    }
    if (!is.numeric(idx) || any(!is.finite(idx)) ||
        any(idx!=floor(idx) | idx<1 | idx>length(x$W))) stop("Invalid process index")
    if (!is.logical(vs) || length(vs)!=1L || is.na(vs)) stop("vs must be TRUE or FALSE")
    if (!is.numeric(col.alpha) || length(col.alpha)!=1L || !is.finite(col.alpha) ||
        col.alpha<0 || col.alpha>1) stop("col.alpha must be between zero and one")
    if (is.null(legend) || identical(legend,FALSE)) legend <- "none" else
        if (identical(legend,TRUE)) legend <- "type1" else
            legend <- match.arg(legend,c("type1","type2","none"))
    if (length(col)>1L && any(col=="none",na.rm=TRUE)) stop("Use col='none' alone to disable simulated lines")
    colors <- if (is.null(col) || !length(col) || identical(col,"none")) NULL else
        grDevices::adjustcolor(col,alpha.f=col.alpha)
    extra <- list(...)
    if (length(extra) && (is.null(names(extra)) || any(!nzchar(names(extra)))))
        stop("Additional plotting arguments must be named")
    for (i in idx) {
        observed <- x$W[[i]]; locations <- x$x[[i]]; simulated <- x$What[[i]]
        if (!is.numeric(observed) || !length(observed) || any(!is.finite(observed)) ||
            !is.numeric(locations) || length(locations)!=length(observed) || any(!is.finite(locations)) ||
            !is.matrix(simulated) || nrow(simulated)!=length(observed) || any(!is.finite(simulated)))
            stop("Incompatible or nonfinite plotted process data")
        score <- length(x$type)>=i && x$type[i]=="score"
        ordinal_norm <- length(x$type)>=i && x$type[i]=="ordinal-norm"
        xlabel <- if (!missing(xlab)) xlab else if (!vs) "Observation" else if (score) "Time" else x$variable[i]
        ylabel <- if (!missing(ylab)) ylab else if (ordinal_norm)
            expression(sqrt(sum(W[k](x)^2,k==1,K))) else as.expression(list(bquote(W[.(x$variable[i])](x))))
        main <- if (!missing(title)) title else if (!is.null(extra$main)) extra$main else
            if (!is.null(x$response)) x$response else if (score) x$variable[i] else ""
        x0 <- if (vs) locations else seq_along(locations)
        limits <- if (!is.null(ylim)) ylim else if (ordinal_norm) {
            peak <- max(c(observed,simulated))
            c(0,if (peak==0) 1 else 1.25*peak)
        } else max(abs(observed))*2*c(-1,1)
        arguments <- list(x=x0,y=observed,type="n",lwd=2,xlab=xlabel,ylab=ylabel,ylim=limits,main=main)
        for (name in names(extra)) arguments[[name]] <- extra[[name]]
        if (!missing(title)) arguments$main <- title
        do.call(graphics::plot,arguments)
        has_samples <- !is.null(colors) && ncol(simulated)>0L
        if (has_samples)
            graphics::matlines(x0,simulated,type="s",col=colors,lwd=1,lty=1)
        observed_color <- graphics::par("col")
        graphics::lines(x0,observed,type="s",lwd=2,col=observed_color)
        if (legend=="type1") {
            labels <- character()
            if (!is.null(x$KS)) labels <- c(labels,paste0("Sup-test p=",format(x$KS[i],digits=3)))
            if (!is.null(x$CvM)) labels <- c(labels,paste0("L2-test p=",format(x$CvM[i],digits=3)))
            if (!is.null(x$informative)) {
                if (identical(x$informative[i],FALSE)) labels <- c(labels,"No lack-of-fit information")
                if (is.na(x$informative[i])) labels <- c(labels,"Numerically indeterminate score process")
            }
            if (length(labels)) graphics::legend("topright",labels,bg="white")
        } else if (legend=="type2") {
            labels <- if (ordinal_norm) "Observed cutpoint norm" else "Observed"
            line_colors <- observed_color; line_widths <- 2
            if (has_samples) {
                labels <- c(labels,if (ordinal_norm) "Simulated cutpoint norm" else "Sample process under null")
                line_colors <- c(line_colors,colors[1])
                line_widths <- c(line_widths,1)
            }
            graphics::legend("topright",labels,col=line_colors,lty=1,lwd=line_widths,bg="white")
        }
    }
    invisible(original)
}



