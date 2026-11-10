## Machine-readable refusals. A failed optimization is not evidence that a
## covariance maximum lies on the boundary; retain the numerical state so
## callers can distinguish those cases without parsing an error message.
gof_clmm_abort <- function(message,code,stage,diagnostics=list(),parent=NULL) {
    stop(structure(list(message=message,call=NULL,code=code,stage=stage,
                        diagnostics=diagnostics,parent=parent),
                   class=c(paste0("gof_clmm_",code),"gof_clmm_error","error","condition")))
}

gof_polish_snapshot <- function(frame) {
    keep <- c("theta","theta0","trace","H","S","info","boundary","refreshes",
              "iterations","early","failure","node_check","std_score")
    out <- lapply(keep,function(name)
        if (exists(name,envir=frame,inherits=FALSE)) get(name,envir=frame,inherits=FALSE) else NULL)
    names(out) <- keep
    problem <- get("problem",envir=frame)
    if (!is.null(out$theta)) out$covariance <- tryCatch(
        lapply(gof_clmm_unpack(problem,out$theta)$re,`[[`,"Sigma"),error=function(e) NULL)
    if (exists("integ",envir=frame,inherits=FALSE)) {
        integ <- get("integ",envir=frame)
        out$nodes <- integ$nodes
        if (exists("state",envir=frame,inherits=FALSE))
            out$integration <- tryCatch(integ$diagnostics(get("state",envir=frame)),error=function(e) NULL)
        if (!is.null(out$theta) && exists("state",envir=frame,inherits=FALSE))
            out$last_frozen_score <- tryCatch(colSums(integ$score(out$theta,get("state",envir=frame))[[1L]]),
                                              error=function(e) NULL)
    }
    out
}
