##' Plot cumulative residuals from a 'cumres' object
##'
##' \code{plot} displays the observed cumulative residual process with
##' realizations under the null.
##'
##'
##' @param x Object produced by the function \code{cumres}.
##' @param idx vector of numbers (or variable names) indicating which processes
##' from the \code{x} to plot.
##' @param col Color of the sample processes. By setting this parameter to "none"
##' or \code{NULL} no realizations will be drawn. The number of realizations is
##' determined by the \code{cumres}-object.
##' @param col.alpha Transparency-level of plotted sample processes
##' @param legend Type of legend where "type1" gives p-values of GoF-tests and
##' "type2" gives usual type of legends.
##' @param xlab Optional label of x-axis
##' @param ylab Optional label of y-axis
##' @param vs Label of predictor
##' @param ylim Range of y axis
##' @param title Main title
##' @param ... Additional arguments passed to the plot-routine.
##' @author Klaus K. Holst
##' @keywords hplot regression
##' @examples
##'
##' n <- 500; x <- abs(rnorm(n,sd=0.2))+0.01; y <- sqrt(x) + rnorm(n,sd=0.2)
##' l <- lm(y ~ x)
##' g <- cumres(l, R=500)
##' plot(g, idx=1, legend="type2")
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
        xlabel <- if (!missing(xlab)) xlab else if (!vs) "Observation" else if (score) "Time" else x$variable[i]
        ylabel <- if (!missing(ylab)) ylab else as.expression(list(bquote(W[.(x$variable[i])](x))))
        main <- if (!missing(title)) title else if (!is.null(extra$main)) extra$main else
            if (!is.null(x$response)) x$response else if (score) x$variable[i] else ""
        x0 <- if (vs) locations else seq_along(locations)
        limits <- if (is.null(ylim)) max(abs(observed))*2*c(-1,1) else ylim
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
            labels <- "Observed"; line_colors <- observed_color; line_widths <- 2
            if (has_samples) {
                labels <- c(labels,"Sample process under null")
                line_colors <- c(line_colors,colors[1])
                line_widths <- c(line_widths,1)
            }
            graphics::legend("topright",labels,col=line_colors,lty=1,lwd=line_widths,bg="white")
        }
    }
    invisible(original)
}



