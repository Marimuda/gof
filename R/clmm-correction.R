## Regular-interior experimental corrections. If A is the matrix of cluster
## residual contributions and Q spans the cluster likelihood scores, the
## projection reference is (I-QQ')A. HC2 scales its cluster rows by
## 1/sqrt(1-diag(QQ')). The observed residual process is never rescaled.
## Encode this using the shared engine's existing row loadings and correction
## basis; no new random stream or fixed-effect ordinal behavior is introduced.
gof_clmm_correction <- function(problem,fit,residual,gradient,method="derivative") {
    multiplier <- gof_clmm_multiplier(problem,fit)
    if (is.null(method)) method <- "derivative"
    if (method=="derivative") return(list(gradient=gradient,multiplier=multiplier,
        metadata=list(method=method)))
    if (!method %in% c("projection","projection-HC2")) stop("Unknown CLMM correction",call.=FALSE)
    S <- fit$score[[1L]]
    scale <- sqrt(colSums(S^2))
    if (any(!is.finite(scale)) || any(scale==0))
        gof_clmm_abort("Cluster scores do not identify the full correction space",
                       "score_rank_deficient","correction",list(score_scale=scale))
    decomposition <- qr(sweep(S,2L,scale,"/"),tol=1e-10)
    if (decomposition$rank!=ncol(S))
        gof_clmm_abort("Cluster scores do not identify the full correction space",
                       "score_rank_deficient","correction",list(rank=decomposition$rank,parameters=ncol(S)))
    Q <- qr.Q(decomposition)
    h <- rowSums(Q^2)
    loading <- rep(1,nrow(S))
    if (method=="projection-HC2") {
        if (any(1-h<=sqrt(.Machine$double.eps)))
            gof_clmm_abort("A cluster has unit score leverage; HC2 inference is unavailable",
                           "score_leverage_singular","correction",list(leverage=h))
        loading <- 1/sqrt(1-h)
    }
    index <- problem$blocks$id
    gradient <- do.call(rbind,lapply(seq_len(ncol(residual)),function(k)
        -Q[index,,drop=FALSE]*residual[,k]))
    multiplier <- gof_multiplier_design(list(list(name=problem$chain[1L],index=index,
        influence=Q*loading,loading=loading[index])),n_units=problem$blocks$m,
        divisor=problem$N/sqrt(problem$blocks$m))
    list(gradient=gradient,multiplier=multiplier,metadata=list(method=method,
        rank=decomposition$rank,leverage=h,max_leverage=max(h),
        scope="Experimental regular-interior correction; full calibration qualification pending"))
}
