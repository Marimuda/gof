## Keep the historical R names while native clients use valid C++ identifiers.
.fastapprox <- function(time,newtime,equal=FALSE,type=0L) {
    type <- gof_count(type,"type")
    if (type>2L) stop("type must be 0, 1 or 2")
    if (!is.logical(equal) || length(equal)!=1L || is.na(equal)) stop("equal must be TRUE or FALSE")
    if (!is.numeric(time) || is.complex(time) || !is.null(dim(time)) ||
        !is.numeric(newtime) || is.complex(newtime) || !is.null(dim(newtime)))
        stop("time and newtime must be real numeric vectors")
    FastApprox(time,newtime,equal,type)
}
.clusterindex <- function(t) {
    if (!is.numeric(t) || is.complex(t) || !is.null(dim(t)) || any(!is.finite(t)) ||
        any(t<0 | t!=floor(t) | t>.Machine$integer.max))
        stop("Cluster identifiers must be nonnegative R integers")
    ClusterIndex(t)
}
