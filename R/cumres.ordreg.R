##' Cumulative residual diagnostics for an ordinal regression
##'
##' @param model A fitted lava ordinal regression.
##' @param variable Covariate name, or "predicted" for the linear predictor.
##' @param R Number of null processes to simulate.
##' @param plot Display the resulting process.
##' @param aggregate Function combining cumulative residuals across categories.
##' @param plots Number of simulated processes to retain.
##' @param ... Arguments passed to plot when plot=TRUE.
##' @return A cumres object with Sup and L2 tests.
##' @export
cumres.ordreg <- function(model,variable="predicted",R=100,plot=FALSE,
                          aggregate=sum,plots=min(R,50L),...) {
    R <- gof_count(R,"R",1L)
    plots <- gof_count(plots,"plots")
    if (!requireNamespace("numDeriv",quietly=TRUE)) stop("package 'numDeriv' not available")
    if (!is.function(aggregate)) stop("aggregate must be a function")
    env <- model$up
    if (length(variable)!=1L || !is.character(variable)) stop("Select one ordinal covariate")
    if (variable=="predicted") x <- env$lp[,1] else {
        if (!variable %in% colnames(env$X)) stop("Variable not found: ",variable)
        x <- env$X[,variable]
    }
    n <- length(x)
    k <- length(env$alpha)
    if (!n || !k || any(!is.finite(x))) stop("Ordinal covariates must be finite")
    cpr <- as.matrix(predict(model,type="cumulative"))
    Y <- outer(env$y,seq_len(k),"<=")*1
    residual <- Y-cpr
    gradient <- -numDeriv::jacobian(function(p) predict(model,p=p,type="cumulative"),coef(model))
    ic <- as.matrix(lava::iid(model))
    o <- order(x)
    x <- x[o]
    ends <- cumsum(rle(x)$lengths)
    x <- x[ends]
    residual <- residual[o,,drop=FALSE]
    eta <- lapply(seq_len(k),function(j) {
        z <- gradient[seq_len(n)+(j-1L)*n,,drop=FALSE][o,,drop=FALSE]
        apply(z,2,cumsum)
    })
    ic <- ic[o,,drop=FALSE]
    if (any(!is.finite(residual)) || any(!is.finite(gradient)) || any(!is.finite(ic)))
        stop("Ordinal residuals, gradients and influences must be finite")
    combine <- function(w) {
        rows <- lapply(seq_len(nrow(w)),function(i) aggregate(w[i,]))
        sizes <- lengths(rows)
        if (!length(sizes) || !sizes[1] || any(sizes!=sizes[1]) ||
            any(!vapply(rows,is.numeric,logical(1))) || any(!is.finite(unlist(rows))))
            stop("aggregate must return a fixed-length finite numeric vector")
        do.call(rbind,rows)
    }
    observed <- combine(apply(residual,2,cumsum)[ends,,drop=FALSE]/sqrt(n))
    q <- ncol(observed)
    What <- lapply(seq_len(q),function(i) matrix(0,length(ends),plots))
    obs.ks <- apply(abs(observed),2,max)
    obs.l2 <- apply(observed,2,function(w) L2Test(w,x))
    ks <- l2 <- matrix(0,R,q)
    for (j in seq_len(max(R,plots))) {
        g <- rnorm(n)
        B <- colSums(ic*g)
        w <- vapply(seq_len(k),function(i)
            (cumsum(residual[,i]*g)+drop(eta[[i]] %*% B))[ends]/sqrt(n),
            numeric(length(ends)))
        w <- combine(matrix(w,ncol=k))
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
                             n=n,R=R,type=rep("normal",q),sd=NULL,cvalues=NULL,
                             variable=labels,model="ordreg"),class="cumres")
    if (isTRUE(plot)) graphics::plot(result,...)
    result
}
