# Ordering policies do not alter residual units, the estimator, or calibration.
# Exact ties enter one step; measurement coarsening must be declared explicitly.
gof_ordering_coordinates <- function(values,resolution=NULL,units=NULL) {
    values <- as.matrix(values)
    if (is.null(units)) units <- rep.int(1,nrow(values))
    if (!is.numeric(units) || is.complex(units) || length(units)!=nrow(values) ||
        any(!is.finite(units) | units<=0 | units!=floor(units)))
        stop("Ordering unit counts must be positive finite integers per coordinate",call.=FALSE)
    variables <- colnames(values)
    p <- ncol(values)
    steps <- rep(NA_real_,p)
    if (!is.null(resolution)) {
        if (!is.numeric(resolution) || is.complex(resolution) || !is.null(dim(resolution)) ||
            !length(resolution) || any(!is.finite(resolution) | resolution<=0) ||
            is.null(names(resolution)) || anyNA(names(resolution)) ||
            any(!nzchar(names(resolution))) || anyDuplicated(names(resolution)))
            stop("resolution must be a named numeric vector of positive finite measurement steps",call.=FALSE)
        if ("predicted" %in% names(resolution))
            stop("resolution cannot coarsen fitted predictions; use a prespecified observed ordering",call.=FALSE)
        if (length(setdiff(names(resolution),variables)))
            stop("resolution must name selected observed ordering variables",call.=FALSE)
        steps[match(names(resolution),variables)] <- resolution
    }
    changed <- numeric(p)
    for (j in which(!is.na(steps))) {
        quotient <- values[,j]/steps[j]
        if (any(!is.finite(quotient) | abs(quotient)>2^53))
            stop("resolution is too fine to represent its grid indices exactly; rescale the ordering",call.=FALSE)
        rounded <- round(quotient)*steps[j]
        if (any(!is.finite(rounded))) stop("resolution produces nonfinite coordinates",call.=FALSE)
        changed[j] <- sum(units[rounded!=values[,j]])
        values[,j] <- rounded
    }
    counts <- lapply(seq_len(p),function(j) {
        z <- values[,j]
        # Binomial response bootstrap stores one residual total per group,
        # but its underlying units are trials. Count ties without expanding.
        as.numeric(rowsum(units,match(z,unique(z)),reorder=FALSE)[,1L])
    })
    metadata <- data.frame(variable=variables,source=rep("fixed",p),ties=rep("group",p),
        resolution=steps,distinct_values=vapply(counts,length,integer(1)),
        tie_blocks=vapply(counts,function(n) sum(n>1L),integer(1)),
        tied_units=vapply(counts,function(n) sum(n[n>1L]),numeric(1)),
        changed_units=changed,status=rep("fixed",p),reason=rep(NA_character_,p),
        stringsAsFactors=FALSE)
    list(values=values,metadata=metadata)
}

gof_ordering_unresolved <- function(object) {
    if (is.null(object$ordering)) return(rep(FALSE,length(object$variable)))
    object$ordering$status=="unresolved"
}

gof_fitted_ordering_check <- function(X,beta,offset,mu,success,exposure,family,b=0) {
    # A local numerical diagnostic, NOT an exact-root certificate. A Fisher
    # score correction estimates fitting error; gamma_k allowances screen for
    # floating-point cancellation. No coordinates or fitted means are changed.
    # Validate this diagnostic against independent high-precision references.
    fail <- function(reason,details=list()) list(status="unresolved",reason=reason,details=details)
    pass <- function(status,details=list()) list(status=status,reason=NA_character_,details=details)
    n <- nrow(X); p <- ncol(X)
    if (!p) return(pass("specified"))
    if (any(!is.finite(c(X,beta,offset,mu,success,exposure))))
        return(fail("Nonfinite inputs to the fitted-order numerical check"))
    functional <- data.frame(X,offset=offset,check.names=FALSE)
    representatives <- which(!duplicated(functional))
    k <- length(representatives)
    if (k==1L) return(pass("structural"))
    # With a full-rank fit and exactly p natural-parameter cells, every union
    # of those cells is saturated. Binomial expansion gives unit trial weights;
    # other mean-response processes require constant precision within a cell.
    if (family$family=="binomial") constant_precision <- TRUE else {
        augmented <- data.frame(functional,exposure=exposure,check.names=FALSE)
        constant_precision <- sum(!duplicated(augmented))==k
    }
    if (k==p && constant_precision) return(pass("saturated"))
    scale <- apply(abs(X),2L,max)
    if (any(!is.finite(scale) | scale==0)) return(fail("Singular fitted design in the ordering check"))
    Z <- sweep(X,2L,scale,"/")
    eta <- drop(X%*%beta)+offset
    derivative <- as.numeric(family$mu.eta(eta))
    variance <- as.numeric(family$variance(mu))
    if (any(!is.finite(c(derivative,variance))) || any(variance<=0) || any(derivative==0))
        return(fail("Degenerate link derivative or variance in the ordering check"))
    ratio <- derivative/variance
    score <- (success-exposure*mu)*ratio
    magnitude <- (abs(success)+exposure*abs(mu))*abs(ratio)
    information_weight <- exposure*derivative^2/variance
    H <- crossprod(Z,Z*information_weight)
    inverse <- tryCatch(chol2inv(chol(H)),error=function(e) NULL)
    if (is.null(inverse) || any(!is.finite(inverse)))
        return(fail("Ill-conditioned information in the fitted-order numerical check"))
    gamma <- function(operations) {
        x <- as.double(operations)*.Machine$double.eps
        if (!is.finite(x) || x>=1) Inf else x/(1-x)
    }
    # Summation plus predictor/link/weight evaluation allowance. This is a
    # conservative engineering screen; libm and nonlinear root errors are not
    # enclosed by a formal outward-rounded interval proof.
    allowance <- gamma(n+2*as.double(p)+32)
    absolute_Z <- abs(Z); absolute_inverse <- abs(inverse)
    H_error <- allowance*crossprod(absolute_Z,absolute_Z*abs(information_weight))
    product_error <- abs(diag(p)-inverse%*%H)+absolute_inverse%*%H_error+
        gamma(p+2)*absolute_inverse%*%abs(H)
    amplification <- max(rowSums(product_error))
    if (!is.finite(amplification) || amplification>=1)
        return(fail("Information roundoff is unresolved in the fitted-order numerical check"))
    gradient <- drop(crossprod(Z,score))
    correction <- drop(inverse%*%gradient)
    score_error <- allowance*drop(crossprod(absolute_Z,magnitude))
    correction_error <- (drop(absolute_inverse%*%score_error)+
        gamma(p+2)*drop(absolute_inverse%*%abs(gradient)))/(1-amplification)
    radius <- abs(correction)+correction_error
    if (any(!is.finite(radius))) return(fail("Nonfinite fitting-error allowance for fitted ordering"))
    eta_allowance <- drop(abs(Z)%*%radius)+gamma(2*as.double(p)+8)*
        (drop(abs(X)%*%abs(beta))+abs(offset))
    if ((family$family=="poisson" && family$link %in% c("identity","sqrt") &&
         any(eta<=eta_allowance)) ||
        (family$family=="binomial" && family$link=="log" && any(eta+eta_allowance>=0)))
        return(fail("The local fitting allowance crosses the link's mean domain"))
    o <- representatives[order(mu[representatives])]
    left <- head(o,-1L); right <- tail(o,-1L)
    difference_X <- X[right,,drop=FALSE]-X[left,,drop=FALSE]
    difference_Z <- sweep(difference_X,2L,scale,"/")
    offset_difference <- offset[right]-offset[left]
    contrast <- drop(difference_X%*%beta)+offset_difference
    arithmetic <- gamma(2*as.double(p)+8)*
        (drop(abs(difference_X)%*%abs(beta))+abs(offset_difference))+
        (2*as.double(p)+8)*.Machine$double.xmin
    error <- drop(abs(difference_Z)%*%radius)+arithmetic
    increasing <- derivative[left]>0 & derivative[right]>0
    decreasing <- derivative[left]<0 & derivative[right]<0
    stable <- (increasing & contrast>error) | (decreasing & -contrast>error)
    # The Gaussian inverse link is decreasing on each side of its pole.
    # Across the pole, signs determine order only when neither eta can cross 0.
    if (family$family=="gaussian" && family$link=="inverse") {
        same_side <- sign(eta[left])==sign(eta[right])
        stable <- (same_side & stable) | (!same_side & mu[right]>mu[left])
        stable <- stable & abs(eta[left])>eta_allowance[left] & abs(eta[right])>eta_allowance[right]
    }
    stable <- stable & mu[right]>mu[left] & is.finite(error) & is.finite(contrast)
    details <- list(method="local score correction and accumulated-roundoff screen",
        roundoff_factor=allowance,information_amplification=amplification,
        max_scaled_score_correction=max(abs(correction)),
        comparisons=length(stable),unresolved_comparisons=sum(!stable),bandwidth=b)
    if (any(!stable)) return(fail("Fitted comparisons are sensitive to local fitting error or roundoff",details))
    if (b>0) {
        # Window entry/exit order also matters: stable ordering of mu alone
        # does not settle a comparison of mu_i with mu_j+b. Screen the cross-
        # event comparisons using the same local fitting allowance. Within a
        # single event type, adding a constant preserves the checked order.
        lower <- family$linkinv(eta-eta_allowance)
        upper <- family$linkinv(eta+eta_allowance)
        mean_error <- pmax(abs(lower-mu),abs(upper-mu))+gamma(8)*abs(mu)
        if (any(!is.finite(mean_error)))
            return(fail("Moving-window link evaluation is numerically unresolved",details))
        ids <- rep(representatives,2L)
        type <- rep(0:1,each=k)
        event <- mu[ids]+type*b
        error_event <- mean_error[ids]+gamma(2)*(abs(mu[ids])+type*b)
        if (any(!is.finite(c(event,error_event))))
            return(fail("Moving-window endpoint arithmetic is nonfinite",details))
        o <- order(event)
        a <- head(o,-1L); z <- tail(o,-1L)
        stable_event <- event[z]>event[a]
        different <- type[a]!=type[z] & ids[a]!=ids[z]
        stable_event[different] <- (event[z]-event[a])[different]>
            (error_event[a]+error_event[z])[different]
        details$window_comparisons <- length(stable_event)
        details$unresolved_window_comparisons <- sum(!stable_event)
        if (any(!stable_event))
            return(fail("Fitted moving-window entry/exit comparisons are sensitive to fitting error or roundoff",details))
    }
    pass("checked",details)
}

gof_ordering_report <- function(metadata,on_unresolved) {
    variables <- metadata$variable[metadata$status=="unresolved"]
    if (!length(variables)) return(invisible(NULL))
    message <- paste0("Fitted ordering unresolved for ",paste(variables,collapse=", "),
        "; coordinates are retained and inferential p-values withheld. Use a prespecified fixed ordering.")
    if (on_unresolved=="error") stop(message,call.=FALSE)
    warning(message,call.=FALSE)
}
