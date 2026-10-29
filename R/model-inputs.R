gof_information_inverse <- function(model,p) {
    if (!p) return(matrix(0,0L,0L))
    if (is.null(model$qr)) stop("The fitted model must retain its QR decomposition")
    R <- qr.R(model$qr)[seq_len(p),seq_len(p),drop=FALSE]
    inverse <- chol2inv(R)
    result <- matrix(0,p,p)
    pivot <- model$qr$pivot[seq_len(p)]
    result[pivot,pivot] <- inverse
    result
}

gof_fitted_data <- function(data,mf) {
    data <- as.data.frame(data)
    positions <- match(rownames(mf),rownames(data))
    if (!anyNA(positions)) return(data[positions,,drop=FALSE])
    if (nrow(data)==nrow(mf)) return(data)
    stop("data must contain the fitted rows, identified by row names, or be in fitted-row order")
}

gof_informative_processes <- function(X,weights,ordering,grid=NULL,b=0) {
    # The fitted score equations have the form X' diag(h) r=0. If each indicator
    # is a union of saturated design cells with constant precision within
    # cells, its residual sum is constrained to zero. This is a structural
    # lack of test information, not evidence inferred from small p-values.
    p <- ncol(X)
    informative <- rep(TRUE,ncol(ordering))
    if (!p) return(informative)
    for (j in seq_len(ncol(ordering))) {
        levels <- unique(ordering[,j])
        if (length(levels)>p) next
        constant <- rep(TRUE,p)
        cells <- matrix(0,length(levels),p)
        precision_constant <- TRUE
        for (k in seq_along(levels)) {
            use <- ordering[,j]==levels[k]
            if (length(unique(weights[use]))!=1L || !is.finite(weights[use][1]) || weights[use][1]==0) {
                precision_constant <- FALSE
                break
            }
            cells[k,] <- X[which(use)[1],]
            constant <- constant & vapply(seq_len(p),function(column)
                all(X[use,column]==cells[k,column]),logical(1))
        }
        if (precision_constant && any(constant) &&
            qr(cells[,constant,drop=FALSE])$rank==length(levels)) informative[j] <- FALSE
    }
    # A general reparameterization can hide cell indicators in combinations of
    # nonconstant columns. This QR check is a numerical diagnostic, never an
    # exact-saturation certificate: leave the curves intact and report NA when
    # all ACTUAL evaluated indicators are indistinguishable from score space.
    candidates <- which(informative)
    if (!length(candidates)) return(informative)
    if (is.null(grid)) grid <- lapply(seq_len(ncol(ordering)),function(j) sort(unique(ordering[,j])))
    if (length(grid)!=ncol(ordering)) stop("Evaluation grids must match ordering variables")
    column_scale <- apply(abs(X),2,max)
    column_scale[column_scale==0] <- 1
    score <- sweep(X,2,column_scale,"/")*weights
    if (any(!is.finite(score))) return(informative)
    column_scale <- apply(abs(score),2,max)
    column_scale[column_scale==0] <- 1
    score <- sweep(score,2,column_scale,"/")
    factor <- qr(score,LAPACK=TRUE)
    diagonal <- abs(diag(qr.R(factor)))
    # Do not use the arbitrary orthogonal completion of a rank-deficient
    # factorization as score information. Full-rank fitted designs are the
    # normal case; ambiguous score rank is outside this sufficient check.
    rank_bound <- .Machine$double.eps*max(dim(score))*max(diagonal)
    if (length(diagonal)!=p || any(diagonal<=rank_bound)) return(informative)
    Q <- qr.Q(factor)
    # gamma_k is the standard accumulated-roundoff form. The operation budget
    # covers the Householder basis and two projection products conservatively.
    # It depends only on dimensions and machine precision, never residuals,
    # simulated variances or observed p-values. Near-zero positive departures
    # above this backward-error scale remain informative.
    operations <- 8*as.double(nrow(score))*max(1,p)
    accumulated <- operations*.Machine$double.eps
    if (accumulated>=1) return(informative)
    gamma <- accumulated/(1-accumulated)
    for (j in candidates) {
        unresolved <- length(grid[[j]])>0L
        for (t in grid[[j]]) {
            indicator <- as.numeric(ordering[,j]<=t)
            if (b>0) indicator <- indicator*as.numeric(ordering[,j]>t-b)
            projection <- drop(Q %*% crossprod(Q,indicator))
            error <- sqrt(sum((indicator-projection)^2))
            bound <- gamma*(sqrt(sum(indicator^2))+sqrt(sum(projection^2)))
            if (error>bound) {
                unresolved <- FALSE
                break
            }
        }
        if (unresolved) informative[j] <- NA
    }
    informative
}

