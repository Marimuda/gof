# Internal regularity checks for the MLE-based residual multiplier expansion.
# A finite number returned by IRLS is not evidence that an interior MLE exists.

gof_mle_direction <- function(X,nonnegative,equal=integer()) {
    # Find d with C d >= 0, E d = 0 and sum(C d) >= 1. Homogeneity makes
    # the last condition equivalent to at least one strictly improving row;
    # there is no coefficient-size, p-value or solver-objective threshold.
    if (!length(nonnegative)) return(FALSE)
    scale <- apply(abs(X),2L,max)
    if (any(!is.finite(scale)) || any(scale==0))
        stop("Cannot qualify maximum-likelihood regularity for a singular design",call.=FALSE)
    Z <- sweep(X,2L,scale,"/")
    C <- Z[abs(nonnegative),,drop=FALSE]*sign(nonnegative)
    E <- Z[equal,,drop=FALSE]
    row_scale <- apply(abs(C),1L,max)
    C <- C[row_scale>0,,drop=FALSE]
    row_scale <- row_scale[row_scale>0]
    if (!nrow(C)) return(FALSE)
    C <- C/row_scale
    if (nrow(E)) {
        row_scale <- apply(abs(E),1L,max)
        E <- E[row_scale>0,,drop=FALSE]
        if (nrow(E)) E <- E/row_scale[row_scale>0]
    }
    A <- rbind(C,E,colSums(C))
    directions <- c(rep(">=",nrow(C)),rep("=",nrow(E)),">=")
    rhs <- c(rep(0,nrow(C)+nrow(E)),1)
    # lpSolve variables are nonnegative: d=d.plus-d.minus is unrestricted.
    solution <- lpSolve::lp("min",rep(0,2L*ncol(Z)),cbind(A,-A),directions,rhs)
    # Solver infeasibility is a numerical overlap diagnostic, not an exact
    # rational-arithmetic Farkas certificate for every ill-conditioned input.
    if (solution$status==2L) return(FALSE)
    if (solution$status!=0L)
        stop("Cannot qualify maximum-likelihood regularity: linear-programming solver status ",
             solution$status,call.=FALSE)
    d <- solution$solution[seq_len(ncol(Z))]-solution$solution[ncol(Z)+seq_len(ncol(Z))]
    margin <- drop(A%*%d)-rhs
    # Independently verify a reported feasible direction against accumulated
    # matrix-product roundoff. A poor numerical solution is an unresolved
    # regularity check, never an ordinary test p-value.
    error <- 64*.Machine$double.eps*(ncol(Z)+1)*
        (drop(abs(A)%*%abs(d))+abs(rhs)+.Machine$double.xmin)
    if (any(!is.finite(margin)) || tail(margin,1L)+1<=tail(error,1L) ||
        any(margin[directions==">="]< -error[directions==">="]) ||
        any(abs(margin[directions=="="])>error[directions=="="]))
        stop("Cannot qualify maximum-likelihood regularity: imprecise separation solution",call.=FALSE)
    TRUE
}

gof_mle_interior <- function(X,beta,offset,kind,success,failure,weights) {
    # Krawczyk-style numerical qualification of an interior likelihood root.
    # For the constrained non-Gaussian links, negative log likelihood is convex:
    #   binomial/log: -s*eta-f*log(1-exp(eta)), eta<0;
    #   Poisson/identity: -w*y*log(eta)+w*eta, eta>0;
    #   Poisson/sqrt: -2*w*y*log(eta)+w*eta^2, eta>0.
    # Gaussian log/inverse links use the unscaled mean score (dispersion
    # cancels). Their observed Hessian need not be globally positive; the
    # enclosure qualifies a local regular mean root, not global optimality.
    # A positive observed Hessian plus small score ALONE is insufficient:
    # a boundary optimum can have both. We enclose the actual design domain
    # and Hessian over a parameter box and require an interior fixed point.
    #
    # In exact interval arithmetic K=beta-R*g+(I-R*H(B))*(B-beta) inside B
    # proves an interior root. Here conservative gamma_k roundoff envelopes
    # implement a numerical qualification, not a formal outward-rounded
    # interval-arithmetic proof for arbitrary platforms/libm implementations.
    # Any unresolved enclosure fails closed. Bounds depend only on arithmetic,
    # dimensions and the actual likelihood derivatives, never fitted p-values.
    n <- nrow(X); p <- ncol(X)
    unit <- .Machine$double.eps
    gamma <- function(k) {
        error <- as.double(k)*unit
        if (!is.finite(error) || error>=1) return(Inf)
        error/(1-error)
    }
    rounded_up <- function(value,k) value/(1-gamma(k))+as.double(k)*.Machine$double.xmin
    fail <- function() stop("Cannot qualify a regular interior maximum-likelihood mean fit for this link; ",
                            "the optimum may be on the mean-domain boundary or numerically unresolved",call.=FALSE)
    column_scale <- apply(abs(X),2L,max)
    if (any(column_scale==0 | !is.finite(column_scale))) fail()
    Z <- sweep(X,2L,column_scale,"/")
    centre <- beta*column_scale
    eta <- drop(X%*%beta)+offset
    scaled_eta <- drop(Z%*%centre)+offset
    # The budget includes column division, centre multiplication, the dot
    # product and offset addition, not just the final matrix multiplication.
    eta_error <- rounded_up(gamma(p+8L)*(drop(abs(X)%*%abs(beta))+abs(offset))+
                            abs(scaled_eta-eta),p+10L)
    domain <- function(lo,hi) {
        if (kind=="binomial.log") return(all(hi<0))
        if (kind=="gaussian.log") return(all(is.finite(lo)) && all(is.finite(hi)))
        if (kind=="gaussian.inverse") return(all((lo>0 & hi>0) | (lo<0 & hi<0)))
        all(lo>0)
    }
    eta_interval <- function(radius=rep(0,n)) {
        endpoint_error <- rounded_up(gamma(4L)*(abs(eta)+eta_error+radius),8L)
        list(lo=eta-eta_error-radius-endpoint_error,
             hi=eta+eta_error+radius+endpoint_error)
    }
    initial_eta <- eta_interval()
    if (!domain(initial_eta$lo,initial_eta$hi)) fail()
    derivatives <- function(z) {
        if (kind=="binomial.log") {
            ratio <- exp(z)/(-expm1(z))
            first <- success
            second <- failure*ratio
            h <- second/(-expm1(z))
        } else if (kind=="poisson.identity") {
            first <- weights*success/z
            second <- weights
            h <- first/z
        } else if (kind=="poisson.sqrt") {
            first <- 2*weights*success/z
            second <- 2*weights*z
            h <- first/z+2*weights
        } else if (kind=="gaussian.log") {
            mu <- exp(z)
            first <- weights*success*mu
            second <- weights*mu^2
            h <- weights*(2*mu^2-success*mu)
        } else {
            mu <- 1/z
            first <- weights*mu^3
            second <- weights*success*mu^2
            h <- weights*(3*mu^4-2*success*mu^3)
        }
        h_magnitude <- if (kind=="gaussian.log") weights*(2*mu^2+abs(success*mu)) else
            if (kind=="gaussian.inverse") weights*(3*mu^4+2*abs(success*mu^3)) else h
        list(score=first-second,magnitude=abs(first)+abs(second),h=h,h_magnitude=h_magnitude)
    }
    curvature <- function(lo,hi) {
        if (kind=="binomial.log") {
            high <- derivatives(hi)$h
            return(list(lo=derivatives(lo)$h,hi=high,magnitude=high))
        }
        if (kind %in% c("poisson.identity","poisson.sqrt"))
            return(list(lo=derivatives(hi)$h,hi=derivatives(lo)$h,magnitude=derivatives(lo)$h))
        # Gaussian Hessian weights are polynomial in mu. Include every
        # stationary point in the interval, not only its endpoints.
        if (kind=="gaussian.log") {
            lower <- exp(lo); upper <- exp(hi)
            hmu <- function(mu) weights*(2*mu^2-success*mu)
            stationary <- success/4
            magnitude <- weights*(2*upper^2+abs(success)*upper)
        } else {
            lower <- 1/hi; upper <- 1/lo
            hmu <- function(mu) weights*(3*mu^4-2*success*mu^3)
            stationary <- success/2
            maximum_mu <- pmax(abs(lower),abs(upper))
            magnitude <- weights*(3*maximum_mu^4+2*abs(success)*maximum_mu^3)
        }
        at_lower <- hmu(lower); at_upper <- hmu(upper)
        at_stationary <- hmu(stationary)
        inside <- stationary>=lower & stationary<=upper
        low <- pmin(at_lower,at_upper,ifelse(inside,at_stationary,Inf))
        high <- pmax(at_lower,at_upper,ifelse(inside,at_stationary,-Inf))
        if (kind=="gaussian.inverse") {
            inside_zero <- lower<=0 & upper>=0
            low <- pmin(low,ifelse(inside_zero,0,Inf))
            high <- pmax(high,ifelse(inside_zero,0,-Inf))
        }
        list(lo=low,hi=high,magnitude=magnitude)
    }
    at_centre <- derivatives(eta)
    initial <- curvature(initial_eta$lo,initial_eta$hi)
    if (any(!is.finite(unlist(c(at_centre,initial))))) fail()
    g <- -drop(crossprod(Z,at_centre$score))
    H <- crossprod(Z,Z*at_centre$h)
    inverse <- tryCatch(chol2inv(chol(H)),error=function(e) NULL)
    if (is.null(inverse) || any(!is.finite(inverse))) fail()
    absolute_Z <- abs(Z); absolute_inverse <- abs(inverse)
    # Include summation, derivative evaluation, and uncertain eta values.
    score_error <- rounded_up(gamma(n+p+64L)*drop(crossprod(absolute_Z,at_centre$magnitude))+
                        drop(crossprod(absolute_Z,(pmax(abs(initial$lo),abs(initial$hi))+
                             gamma(64L)*initial$magnitude)*eta_error)),n+p+64L)
    correction <- drop(inverse%*%g)
    correction_error <- rounded_up(drop(absolute_inverse%*%score_error)+
                            gamma(p+2L)*drop(absolute_inverse%*%abs(g)),p+4L)
    radius <- rounded_up(2*(abs(correction)+correction_error+
                            unit*pmax(1,abs(centre))),p+8L)
    if (any(!is.finite(radius)) || any(radius<=0)) fail()
    eta_radius <- rounded_up(drop(absolute_Z%*%radius),p+4L)
    enclosed_eta <- eta_interval(eta_radius)
    lo <- enclosed_eta$lo
    hi <- enclosed_eta$hi
    if (!domain(lo,hi)) fail()
    interval <- curvature(lo,hi)
    if (any(!is.finite(unlist(interval)))) fail()
    row_difference <- pmax(abs(interval$lo-at_centre$h),abs(interval$hi-at_centre$h))
    H_error <- rounded_up(crossprod(absolute_Z,absolute_Z*row_difference)+
                         gamma(n+64L)*crossprod(absolute_Z,absolute_Z*interval$magnitude),n+64L)
    contraction <- rounded_up(abs(diag(p)-inverse%*%H)+absolute_inverse%*%H_error+
                             gamma(p+4L)*(absolute_inverse%*%abs(H)),p+8L)
    image_radius <- rounded_up(abs(correction)+correction_error+
                               drop(contraction%*%radius),p+8L)
    if (any(!is.finite(image_radius)) || any(image_radius>=radius)) fail()
    invisible(TRUE)
}

gof_validate_glm_mle <- function(model,X,y,weights,eta=NULL) {
    f <- stats::family(model)
    if (!f$family %in% c("binomial","poisson","gaussian")) return(invisible(TRUE))
    active <- weights>0
    if (f$family=="binomial" && is.matrix(y)) active <- active & rowSums(y)>0
    X <- X[active,,drop=FALSE]
    y <- if (is.matrix(y)) y[active,,drop=FALSE] else y[active]
    weights <- weights[active]
    offset <- model$offset
    if (is.null(offset)) offset <- rep(0,length(active))
    offset <- offset[active]
    beta <- stats::coef(model)
    if (is.null(eta)) eta <- drop(X%*%beta)+offset else eta <- eta[active]
    if (any(!is.finite(eta)))
        stop("Maximum-likelihood regularity requires a finite linear predictor",call.=FALSE)
    mu <- f$linkinv(eta)
    if (any(!is.finite(mu)) ||
        (f$family %in% c("binomial","poisson") && any(mu<=0)) ||
        (f$family=="binomial" && any(mu>=1)))
        stop("Diagnostics require fitted means strictly inside the response mean domain",call.=FALSE)
    # Known means have no estimated MLE to run into an estimation boundary.
    # In particular glm.fit intentionally labels empty models boundary=TRUE.
    if (!length(beta)) return(invisible(TRUE))
    if (f$family=="binomial") {
        if (is.factor(y)) y <- as.numeric(y!=levels(y)[1L])
        success <- if (is.matrix(y)) y[,1L]*weights else y*weights
        failure <- if (is.matrix(y)) y[,2L]*weights else (1-y)*weights
        improvement <- c(which(success>0),-which(failure>0))
        if (gof_mle_direction(X,improvement))
            stop("Binomial separation or a numerically unresolved separating direction: ",
                 "diagnostics require a regular interior maximum-likelihood fit",call.=FALSE)
        if (f$link=="log")
            gof_mle_interior(X,beta,offset,"binomial.log",success,failure,weights)
    } else if (f$family=="poisson") {
        success <- as.numeric(y)
        # For a Poisson log link this is the finite-MLE recession condition.
        # For identity/sqrt it is sufficient to rule out an interior optimum,
        # but is not a complete boundary test; the enclosure below is needed.
        if (gof_mle_direction(X,-which(success==0),which(success>0)))
            stop("A Poisson mean-domain boundary or numerically unresolved recession direction ",
                 "prevents a regular interior maximum-likelihood fit",call.=FALSE)
        if (f$link %in% c("identity","sqrt"))
            gof_mle_interior(X,beta,offset,paste0("poisson.",f$link),success,NULL,weights)
    } else if (f$link %in% c("log","inverse")) {
        gof_mle_interior(X,beta,offset,paste0("gaussian.",f$link),as.numeric(y),NULL,weights)
    }
    if (isTRUE(model$boundary))
        stop("The model fitter stopped at a boundary; a regular interior maximum-likelihood fit is required",call.=FALSE)
    invisible(TRUE)
}
