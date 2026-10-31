#' Residual distribution diagnostics for Gaussian structural equation models
#'
#' Compare standardized marginal or projected structural residuals with the
#' normal distribution, accounting for estimation of the model parameters.
#'
#' @param model A regular Gaussian maximum-likelihood `lava::lvmfit` object,
#'   fitted to complete independent-subject data without weights or clusters.
#' @param variable `NULL` pools the symmetrically whitened endogenous manifest
#'   residuals. A single endogenous or latent variable name selects its
#'   standardized projected structural error, conditional on all observed
#'   endogenous variables and the fixed exogenous design.
#' @param R Positive whole number of null simulations.
#' @param plots Nonnegative whole number of simulated processes to retain.
#' @param grid Strictly increasing standard-normal probabilities between zero
#'   and one. The tests use this prespecified finite grid. The Sup statistic
#'   is the maximum absolute ECDF discrepancy on the grid; L2 is the square
#'   root of its left-step integral over this probability grid, approximating
#'   integration against the standard-normal distribution.
#' @param score.tol Positive finite numerical tolerance for the norm of the
#'   total fitted Gaussian score standardized by expected Fisher information.
#'   This checks local optimizer adequacy in standard-error units, not global
#'   maximization or finite-sample calibration.
#' @details
#' This implements the complete-data Gaussian setting of the approximate
#' parametric simulation in Sanchez, Houseman and Ryan (2009), Section 3.1,
#' equations (10)--(11). Each simulation generates the entire endogenous
#' observed vector conditional on the original exogenous covariates. The ECDF
#' process and Gaussian score are calculated from the same generated data.
#' A first-order correction based on expected Fisher information and the
#' derivative of the residual ECDF accounts for nuisance-parameter estimation.
#' Models are not refitted for each simulation.
#' Monte Carlo p-values use `(1 + exceedances)/(R + 1)`, including ties in
#' the upper tail. This avoids reporting zero and has minimum `1/(R + 1)`;
#' it does not make the fitted one-step approximation finite-sample exact.
#'
#' For a residual map `A(theta)` and observed vector with Gaussian mean and
#' covariance, derivatives of its exact univariate Gaussian CDF give the
#' correction. Both the changing mean and changing standardization/projection
#' are included. This is not an ordinary fitted-residual KS test.
#'
#' The observed process is `sqrt(n) * (Fhat - grid)`. The returned `x` is the
#' probability grid, so the usual L2 result extractor integrates against
#' standard-normal probability, not Lebesgue measure on normal quantiles.
#' The grid omits the two infinite tails and approximates the continuous
#' statistics in the paper. Increase its density and tail range for a
#' sensitivity analysis specified independently of the observed discrepancy.
#'
#' The marginal pool checks the one-dimensional standardized marginal
#' distributions; it is not a universally consistent test of multivariate
#' normality or every form of dependence. A selected structural projection is
#' a mixture of observed error sources: rejection cannot uniquely identify
#' which source is nonnormal, and a large p-value does not establish that an
#' isolated latent error is normal. Correct mean/covariance specification,
#' regular interior identifiable Gaussian estimation and sufficient effective
#' information remain assumptions. Missing, weighted and clustered fits, active
#' boundary error covariances and singular expected information are unsupported.
#' First-order calibration can depart from the nominal level in finite samples;
#' numerical convergence and passing a bounded simulation screen do not supply
#' exact size control. The installed `joint-normality-bootstrap` vignette
#' records the qualification studies and their uncertainty.
#' @returns A `cumres_normality` object inheriting from `cumres`, with one
#'   process, compatible with `plot()`,
#'   `summary()`, `as.data.frame()` and optional tidiers. The `calibration`
#'   field is `"gaussian-one-step"`; `normality` records the projection,
#'   normal quantiles, fitted standardized residuals, ECDF derivative,
#'   expected information and standardized fitted score norm.
#' @references B. N. Sanchez, E. A. Houseman and L. M. Ryan (2009).
#'   Residual-Based Diagnostics for Structural Equation Models. Biometrics
#'   65, 104--115. doi:10.1111/j.1541-0420.2008.01022.x.
#' @seealso [cumres.lvmfit()], [cumres()], [plot.cumres()]
#' @examples
#' if (requireNamespace("numDeriv", quietly = TRUE)) {
#'   set.seed(83)
#'   structure <- lava::lvm(y ~ x)
#'   fit <- lava::estimate(structure, lava::sim(structure, 150))
#'   diagnostic <- fit |> cumres_normality(R = 99, plots = 3)
#'   as.data.frame(diagnostic)
#'   plot(diagnostic, xlab = "Standard-normal probability")
#' }
#' @export
cumres_normality <- function(model, variable=NULL, R=999L, plots=min(R,50L),
                            grid=seq(.001,.999,length.out=401L), score.tol=1e-3) {
    R <- gof_count(R,"R",1L)
    plots <- gof_count(plots,"plots")
    if (!is.numeric(grid) || is.complex(grid) || !is.null(dim(grid)) ||
        length(grid)<2L || any(!is.finite(grid)) || any(grid<=0 | grid>=1) ||
        any(diff(grid)<=0))
        stop("grid must contain at least two strictly increasing probabilities between zero and one",call.=FALSE)
    if (!is.numeric(score.tol) || is.complex(score.tol) || length(score.tol)!=1L || !is.finite(score.tol) ||
        score.tol<=0)
        stop("score.tol must be a positive finite numerical tolerance",call.=FALSE)
    setup <- gof_normality_setup(model,variable,grid,score.tol)
    n <- setup$n
    quantiles <- stats::qnorm(grid)
    ecdf_process <- function(z) {
        # Pool components only after generating a complete subject vector.
        z <- as.numeric(z)
        sqrt(n)*(findInterval(quantiles,sort(z))/length(z)-grid)
    }
    statistic <- function(w) c(SupTest(w),L2Test(w,grid))
    observed <- ecdf_process(setup$standardized)
    null_process <- function() {
        error <- matrix(stats::rnorm(n*setup$d),n,setup$d) %*% setup$chol
        w <- ecdf_process(error %*% t(setup$map))
        if (length(setup$parameters)) {
            score <- gof_normality_score(error,setup)
            w <- w+sqrt(n)*drop(setup$delta %*% setup$information.inverse %*% score)
        }
        w
    }
    simulated <- t(vapply(seq_len(R),function(i) statistic(null_process()),numeric(2L)))
    curves <- matrix(0,length(grid),plots)
    for (i in seq_len(plots)) curves[,i] <- null_process()
    observed.stat <- statistic(observed)
    p <- (1+colSums(sweep(simulated,2L,observed.stat,">=")))/(R+1)
    label <- if (is.null(variable)) "marginal normality" else paste0(variable," normality")
    structure(list(W=list(observed),What=list(curves),x=list(grid),KS=p[1L],CvM=p[2L],
        variable=label,R=R,n=n,nobs=n,type="normality",model="lvmfit",
        design_info=gof_design_info(unit="subject",n_rows=n,
            fit_parameters=length(setup$parameters),estimated_parameters=length(setup$parameters),
            correction_dimension=length(setup$parameters),n_components=ncol(setup$standardized)),
        calibration="gaussian-one-step",normality=list(variable=variable,
            quantiles=quantiles,residuals=setup$standardized,projection=setup$map,
            delta=setup$delta,information=setup$information,score.norm=setup$score.norm,
            score.tol=score.tol)),class=c("cumres_normality","cumres"))
}

#' @rdname cumres_normality
#' @param x,object A normality diagnostic or its summary.
#' @param row.names,optional Arguments passed to [base::as.data.frame()].
#' @param ... Additional arguments passed to the corresponding result or plot
#'   method.
#' @export
as.data.frame.cumres_normality <- function(x,row.names=NULL,optional=FALSE,...) {
    result <- as.data.frame.cumres(x,row.names=row.names,optional=optional,...)
    result$calibration <- rep("gaussian-one-step",nrow(result))
    result$diagnostic <- rep("normality",nrow(result))
    result$residual <- rep(if (is.null(x$normality$variable)) "marginal" else
        x$normality$variable,nrow(result))
    result
}

#' @rdname cumres_normality
#' @export
summary.cumres_normality <- function(object,...) {
    result <- summary.cumres(object,...)
    result$normality <- object$normality
    result$calibration <- object$calibration
    class(result) <- c("summary.cumres_normality","summary.cumres")
    result
}

#' @rdname cumres_normality
#' @export
print.cumres_normality <- function(x,...) {
    cat("\nGaussian SEM residual normality diagnostic\n")
    print(coef.cumres(x),...)
    cat("\nOne-step Gaussian calibration based on ",x$R," simulations.\n",sep="")
    cat("Statistics use a prespecified finite standard-normal probability grid.\n")
    invisible(x)
}

#' @rdname cumres_normality
#' @export
print.summary.cumres_normality <- function(x,...) {
    print.cumres_normality(x,...)
}

#' @rdname cumres_normality
#' @param xlab,ylab Optional plot axis labels.
#' @export
plot.cumres_normality <- function(x,xlab="Standard-normal probability (u)",
    ylab=expression(sqrt(n)*(hat(F)(Phi^{-1}(u))-u)),...) {
    plot.cumres(x,xlab=xlab,ylab=ylab,...)
}

gof_normality_moments <- function(model,parameters,data,variable) {
    ii <- lava::index(model)
    moment <- lava::moments(model,p=parameters,data=data,conditional=TRUE)
    covariance <- as.matrix(moment$C)
    mean <- as.matrix(moment$xi)
    n <- nrow(data)
    if (nrow(mean)==1L && n>1L) mean <- mean[rep(1L,n),,drop=FALSE]
    if (nrow(mean)!=n || ncol(mean)!=nrow(covariance) ||
        any(!is.finite(mean)) || any(!is.finite(covariance)) ||
        !isSymmetric(covariance))
        stop("Invalid conditional Gaussian mean or covariance",call.=FALSE)
    covariance.chol <- tryCatch(chol(covariance),error=function(e) NULL)
    if (is.null(covariance.chol) || isTRUE(attr(moment$IAi,"pseudo")))
        stop("Normality diagnostics require an invertible structural system and positive-definite observed covariance",call.=FALSE)
    if (is.null(variable)) {
        spectral <- eigen(covariance,symmetric=TRUE)
        if (any(spectral$values<=0)) stop("Observed covariance must be positive definite",call.=FALSE)
        map <- tcrossprod(sweep(spectral$vectors,2L,sqrt(spectral$values),"/"),spectral$vectors)
    } else {
        selected <- match(variable,ii$vars)
        observed <- match(ii$endogenous,ii$vars)
        # Cov(error, Y | exogenous) = P (I-A')^-T, restricted to Y.
        error.covariance <- moment$P
        error.covariance[ii$exo.idx,] <- 0
        error.covariance[,ii$exo.idx] <- 0
        cross <- error.covariance[selected,,drop=FALSE] %*%
            t(moment$IAi[observed,,drop=FALSE])
        map <- cross %*% chol2inv(covariance.chol)
        variance <- drop(map %*% covariance %*% t(map))
        if (!is.finite(variance) || variance<=0)
            stop("The selected structural error has no positive-variance observable projection",call.=FALSE)
        map <- map/sqrt(variance)
    }
    list(mean=mean,covariance=covariance,map=map,chol=covariance.chol,moment=moment)
}

gof_normality_score <- function(error,setup) {
    whitened <- error %*% setup$covariance.inverse
    vapply(seq_along(setup$parameters),function(j) {
        mean.term <- sum(whitened*setup$mean.derivative[[j]])
        covariance.term <- sum((whitened %*% setup$covariance.derivative[[j]])*whitened)
        mean.term+.5*(covariance.term-setup$n*sum(diag(
            setup$covariance.inverse %*% setup$covariance.derivative[[j]])))
    },numeric(1L))
}

gof_normality_setup <- function(model,variable,grid,score.tol) {
    if (!inherits(model,"lvmfit")) stop("Supply a fitted lava Gaussian structural-equation model",call.=FALSE)
    if (!is.null(model$opt$convergence) && !isTRUE(model$opt$convergence==0L))
        stop("Model has not converged",call.=FALSE)
    if (!is.null(model$cluster)) stop("Normality diagnostics require independent subjects, not clustered fits",call.=FALSE)
    if (!is.null(model$weights)) stop("Weighted normality diagnostics are not supported",call.=FALSE)
    if (!is.null(model$data$data2)) stop("Multi-sample Gaussian fits require a separate normality construction",call.=FALSE)
    if (length(model$estimator)!=1L || !is.character(model$estimator) ||
        is.na(model$estimator) || !model$estimator %in% c("gaussian","gaussian1","gaussian2","normal"))
        stop("Normality diagnostics require a Gaussian maximum-likelihood fit",call.=FALSE)
    data <- as.data.frame(model.frame(model))
    ii <- lava::index(model)
    required <- ii$manifest
    if (!length(ii$endogenous) || !nrow(data) || !all(required %in% names(data)) ||
        any(!vapply(data[required],function(x) is.numeric(x) && !is.complex(x),logical(1L))) ||
        any(!is.finite(as.matrix(data[required]))))
        stop("Normality diagnostics require complete finite numeric manifest data",call.=FALSE)
    if (!is.null(variable) && (!is.character(variable) || length(variable)!=1L ||
        is.na(variable) || !variable %in% c(ii$endogenous,ii$latent)))
        stop("variable must be NULL or one endogenous or latent structural-error name",call.=FALSE)
    other <- setdiff(seq_along(ii$vars),ii$exo.idx)
    if (length(ii$exo.idx) && any(ii$P[ii$exo.idx,other,drop=FALSE]!=0))
        stop("Exogenous variables must have structurally zero error covariance with endogenous and latent variables",call.=FALSE)
    parameters <- lava::pars(model)
    if (!is.numeric(parameters) || any(!is.finite(parameters))) stop("Fitted parameters must be finite",call.=FALSE)
    base <- gof_normality_moments(model,parameters,data,variable)
    structural <- as.matrix(base$moment$P[other,other,drop=FALSE])
    active <- which(diag(structural)!=0)
    free.variance <- vapply(other,function(i) ii$P[i,i]!=0,logical(1L))
    inactive <- setdiff(seq_len(nrow(structural)),active)
    if (any(!is.finite(structural)) || !isSymmetric(structural) ||
        any(diag(structural)[free.variance]<=0) || any(diag(structural)<0) ||
        (length(inactive) && any(structural[inactive,,drop=FALSE]!=0)) ||
        (length(active) && is.null(tryCatch(chol(structural[active,active,drop=FALSE]),error=function(e) NULL))))
        stop("Normality diagnostics require interior positive-definite active structural error covariance",call.=FALSE)
    n <- nrow(data); d <- ncol(base$mean); k <- nrow(base$map); q <- length(parameters)
    setup <- list(n=n,d=d,parameters=parameters,map=base$map,chol=base$chol,
        covariance.inverse=chol2inv(base$chol),mean.derivative=list(),
        covariance.derivative=list(),information=matrix(0,q,q),
        information.inverse=matrix(0,q,q),delta=matrix(0,length(grid),q),score.norm=0,
        standardized=(as.matrix(data[ii$endogenous])-base$mean) %*% t(base$map))
    if (!q) return(setup)
    if (!requireNamespace("numDeriv",quietly=TRUE)) stop("Install 'numDeriv' for normality nuisance corrections",call.=FALSE)
    derivative <- numDeriv::jacobian(function(p) {
        names(p) <- names(parameters)
        moment <- gof_normality_moments(model,p,data,variable)
        c(moment$mean,moment$covariance,moment$map)
    },parameters,method="Richardson")
    if (any(!is.finite(derivative))) stop("Nonfinite Gaussian moment derivatives",call.=FALSE)
    for (j in seq_len(q)) {
        setup$mean.derivative[[j]] <- matrix(derivative[seq_len(n*d),j],n,d)
        setup$covariance.derivative[[j]] <- matrix(derivative[n*d+seq_len(d*d),j],d,d)
    }
    for (j in seq_len(q)) for (l in seq_len(j)) {
        entry <- sum((setup$mean.derivative[[j]] %*% setup$covariance.inverse)*setup$mean.derivative[[l]])+
            n/2*sum((setup$covariance.inverse %*% setup$covariance.derivative[[j]])*
                t(setup$covariance.inverse %*% setup$covariance.derivative[[l]]))
        setup$information[j,l] <- setup$information[l,j] <- entry
    }
    information.chol <- tryCatch(chol(setup$information),error=function(e) NULL)
    if (is.null(information.chol)) stop("Expected Gaussian information must be positive definite; unidentified normality corrections are unsupported",call.=FALSE)
    setup$information.inverse <- chol2inv(information.chol)
    observed.score <- gof_normality_score(as.matrix(data[ii$endogenous])-base$mean,setup)
    setup$score.norm <- sqrt(drop(crossprod(observed.score,setup$information.inverse %*% observed.score)))
    if (!is.finite(setup$score.norm) || setup$score.norm>score.tol)
        stop("Fitted Gaussian score is not sufficiently close to an interior root; refit the supplied model",call.=FALSE)
    z <- stats::qnorm(grid)
    for (j in seq_len(q)) {
        map.derivative <- matrix(derivative[n*d+d*d+seq_len(k*d),j],k,d)
        # Average over subjects and pooled projection components. Both factors
        # in A Sigma A' move, hence dVar/2 = dA Sigma A'.
        location <- mean(setup$mean.derivative[[j]] %*% t(base$map))
        scale <- mean(rowSums((map.derivative %*% base$covariance)*base$map))
        setup$delta[,j] <- stats::dnorm(z)*(location-z*scale)
    }
    setup
}
