##' @rdname cumres.glm
##' @export
`cumres.lm` <- function(model,...) {
  cumres.glm(model,...)
}

##' Cumulative residual diagnostics for linear and generalized linear models
##' 
##' Given the generalized linear models model \deqn{g(E(Y_i|X_{i1},...,X_{ik})) =
##' \sum_{j=1}^k \beta_jX_{ij}} the \code{cumres}-function calculates the
##' observed cumulative sum of residual process, cumulating the residuals,
##' \eqn{e_i}, by the jth covariate: \deqn{W_j(t) = n^{-1/2}\sum_{i=1}^n
##' 1_{\{X_{ij}\le t\}}e_i} and Sup and L2 test
##' statistics are calculated via simulation from the asymptotic distribution of
##' the cumulative residual process under the null (Lin et al., 2002).
##' 
##' 
##' @aliases cumres.glm cumres.lm
##' @param model Model object (\code{lm} or \code{glm})
##' @param variable Names of ordering variables. The default selects fitted
##' predictions and design columns with more than two distinct values. Explicit
##' selections may also use binary or constant numeric variables, or ordered
##' factors supplied through \code{data}. Unordered factors are rejected.
##' @param data Data containing additional ordering variables. Rows are aligned
##' to the fitted model by row names, or must already be in fitted-row order.
##' @param R Positive whole number of multiplier simulations.
##' @param b Nonnegative moving-window bandwidth. Zero gives cumulative
##' residuals; positive values give residual sums over the window \eqn{(t-b,t]}.
##' @param plots Nonnegative whole number of null processes to retain for plotting.
##' @param ... additional arguments passed to \code{cumres.default}. In
##' particular \code{subset.max} (default 500): each process is evaluated at
##' at most this many distinct covariate values, retaining every observation
##' in the fitted model. Use \code{subset.max=Inf} or \code{NULL} for all
##' distinct values. Tied observations are accumulated together. The deprecated
##' \code{breakties} argument is accepted with a warning when nonzero; no jitter
##' is added. A custom \code{subset} function receives sorted distinct values
##' and \code{n.qt}, \code{n.eq}, and returns their 1-based evaluation positions.
##' @returns A \code{cumres} object. See \code{\link{cumres}} for its observed
##' and simulated curves, grids, Monte Carlo p-values and fitted-row metadata.
##' @note Normal, binomial and Poisson models with their standard links are
##' supported. Normal and Poisson prior weights describe inverse dispersion of
##' independently observed response means, as in \code{stats::glm}; their scores
##' and coefficient influences account for those weights. Zero-weight rows are
##' excluded. Binomial weights specify trial counts; weighted proportions and
##' grouped responses must reconstruct nonnegative integer successes/failures
##' within floating-point roundoff. Binomial trials are expanded without refitting.
##' Diagnostics retain the fitted coefficients, covariance and offsets, and omit
##' the same incomplete rows as the fit. Normal models require positive residual
##' degrees of freedom. A known mean with no estimated coefficients is supported.
##' Estimated mean parameters must have a regular interior likelihood root.
##' Binomial separation and Poisson recession checks use linear programming;
##' constrained/noncanonical links additionally require numerical local-root
##' qualification. A finite coefficient and \code{converged=TRUE} alone do not
##' establish these conditions. Unresolved/boundary cases fail explicitly.
##' These checks are numerical qualifications, not exact-arithmetic interval
##' proofs, global optimization certificates, or proof that a loose/custom
##' fitter estimated the coefficients accurately. Appropriate estimator
##' consistency, precision and finite-information conditions remain required.
##' Predicted coordinates are reconstructed from the fitted coefficients so
##' identical design rows are not separated by subtraction roundoff.
##' Structurally recognized saturated ordering cells have no lack-of-fit
##' information: \code{informative=FALSE}, zero curves and p-values of one.
##' Other processes whose selected indicators are numerically indistinguishable
##' from the fitted score space are reported as \code{informative=NA}, with
##' missing p-values and the computed curves retained. This avoids interpreting
##' numerical cancellation as evidence of lack of fit. These checks use the
##' design/score structure, not a threshold on observed p-values or residuals.
##' Monte Carlo comparisons count simulated statistics greater than or equal to
##' the observed statistic. \code{n} is the number of independent residual units
##' (trials for expanded binomial data); \code{nobs} records positive-weight fitted rows.
##' @author Klaus K. Holst
##' @seealso \code{\link{cumres}}, \code{\link{cumres.default}} and
##' \code{\link{plot.cumres}}. The \href{https://CRAN.R-project.org/package=timereg}{timereg}
##' package provides similar GoF methods for survival data.
##' @references D.Y. Lin and L.J. Wei and Z. Ying (2002) \emph{Model-Checking
##' Techniques Based on Cumulative Residuals}. Biometrics, Volume 58, pp 1-12.
##' 
##' John Q. Su and L.J. Wei (1991) \emph{A lack-of-fit test for the mean function
##' in a generalized linear model}. Journal. Amer. Statist. Assoc., Volume 86,
##' Number 414, pp 420-426.
##' @keywords models regression
##' @export
##' @examples
##' set.seed(13)
##' d <- data.frame(x = seq(-2, 2, length.out = 120))
##' d$y <- rbinom(nrow(d), 1, plogis(-0.5 + d$x))
##' fit <- glm(y ~ x, data = d, family = binomial())
##' diagnostic <- cumres(fit, variable = "x", R = 99, plots = 5)
##' coef(diagnostic)
##' plot(diagnostic, idx = "x")
##'
##' # The limit applies to evaluation points, not observations.
##' all_points <- cumres(fit, variable = "x", R = 99, plots = 0,
##'   subset.max = Inf)
##' all_points$nobs
##'
##' # A positive bandwidth requests residual sums over (t - b, t].
##' window <- cumres(fit, variable = "x", R = 99, plots = 0, b = 0.5)
##' coef(window)
`cumres.glm` <- function(model,
                  variable=c("predicted",colnames(model.matrix(model))),
                  data=model.frame(model),
                  R=1000,b=0,plots=min(R,50),...) {
    if (!inherits(model,"lm")) stop("Unsupported model!")
    is_glm <- inherits(model,"glm")
    if (is_glm && !family(model)$family %in% c("gaussian","binomial","poisson"))
        stop("Unsupported family: ",family(model)$family)
    beta <- coef(model)
    if (anyNA(beta)) stop("Over-parametrized model")
    if (!is.numeric(beta) || any(!is.finite(beta))) stop("Coefficients must be finite numeric values")
    if (is_glm && !isTRUE(model$converged)) stop("Model has not converged")
    mf <- tryCatch(model.frame(model),error=function(e)
        stop("Cannot recover fitted rows; fit with model=TRUE or retain the model data: ",conditionMessage(e)))
    Xfit <- model.matrix(model)
    if (nrow(Xfit)!=nrow(mf) || ncol(Xfit)!=length(beta) ||
        length(model$fitted.values)!=nrow(mf)) stop("Inconsistent fitted model dimensions")
    lhs <- formula(model)[[2]]
    response <- paste(deparse(lhs,width.cutoff=500),collapse="")
    # Keep the historical first-response-variable alias. Other LHS components
    # can be legitimate covariates, e.g. trial size in successes/size.
    response_aliases <- unique(c(response,names(mf)[1],head(all.vars(lhs),1L)))
    automatic_variables <- missing(variable)
    force(variable)
    if (!is.character(variable) || !length(variable) || anyNA(variable))
        stop("variable must contain one or more names")
    variable[variable %in% response_aliases] <- "predicted"
    variable <- unique(variable)
    w <- model.weights(mf)
    if (is.null(w)) w <- rep(1,nrow(mf))
    if (!is.numeric(w) || length(w)!=nrow(mf) || any(!is.finite(w)) || any(w<0))
        stop("Model weights must be finite and nonnegative")
    active <- w>0
    if (!any(active)) stop("The fitted model has no positive-weight observations")
    if ((!is_glm || family(model)$family=="gaussian") && sum(active)<=length(beta))
        stop("Normal-model diagnostics require positive residual degrees of freedom")
    supplied_weights <- !is.null(model.weights(mf))
    y <- model.response(mf)
    # lm stores fitted = y-residual, which can manufacture tiny differences
    # between identical design rows. Reconstruct from the actual coefficients
    # so mathematically tied predictions remain tied without rounding data.
    eta_fit <- drop(Xfit%*%beta)
    if (!is.null(model$offset)) eta_fit <- eta_fit+model$offset
    yhat_fit <- if (is_glm) family(model)$linkinv(eta_fit) else eta_fit
    is_binomial <- is_glm && family(model)$family=="binomial"
    grouped <- is_binomial && (is.matrix(y) || supplied_weights)
    if (is.matrix(y) && !grouped) stop("Only two-column binomial matrix responses are supported")
    if (is_binomial && is.factor(y)) y <- as.numeric(y!=levels(y)[1])
    if (grouped) {
        if (is.matrix(y)) {
            if (ncol(y)!=2L) stop("Only two-column binomial matrix responses are supported")
            counts <- y*w
        } else counts <- cbind(y,1-y)*w
        # Allow only floating-point roundoff when reconstructing integer trials.
        tolerance <- 64*.Machine$double.eps*pmax(1,abs(counts))
        if (any(!is.finite(counts)) || any(counts<0 | abs(counts-round(counts))>tolerance))
            stop("Binomial responses and weights must define nonnegative integer success/failure counts")
        counts <- round(counts)
        active <- rowSums(counts)>0
        total <- sum(counts)
        if (!total || total>.Machine$integer.max) stop("The expanded trial count must be between 1 and .Machine$integer.max")
        gof_validate_glm_mle(model,Xfit,y,w,eta_fit)
        rows <- rep(seq_len(nrow(mf)),rowSums(counts))
        y0 <- unlist(lapply(seq_len(nrow(counts)),function(i) rep(c(1,0),counts[i,])))
        X <- Xfit[rows,,drop=FALSE]
        yhat <- yhat_fit[rows]
        eta <- eta_fit[rows]
        f <- family(model)
        dmu <- f$mu.eta(eta)
        if (any(yhat<=0 | yhat>=1)) stop("Binomial fitted probabilities must be strictly between zero and one")
        r <- y0-yhat
        dr <- -X*as.numeric(dmu)
        ic <- if (length(beta)) (X*as.numeric(r*dmu/f$variance(yhat))) %*% vcov(model) else matrix(0,length(r),0L)
    } else {
        if (is_binomial && any(!y %in% c(0,1)))
            stop("Binomial proportions require trial counts supplied as weights")
        if (is_glm) gof_validate_glm_mle(model,Xfit,y,w,eta_fit)
        rows <- which(active)
        X <- Xfit[rows,,drop=FALSE]
        yhat <- yhat_fit[rows]
        if (is_glm) {
            f <- family(model)
            dmu <- as.numeric(f$mu.eta(eta_fit[rows]))
            r <- as.numeric(y)[rows]-yhat
            dr <- -X*dmu
        } else {
            r <- as.numeric(model$residuals)[rows]
            dr <- -X
        }
        if (!length(beta) || all(r==0)) {
            ic <- matrix(0,length(r),length(beta))
        } else if (supplied_weights) {
            # R prior weights describe response-mean precision. Each active
            # row is one independent mean; account for weights in the score
            # and coefficient influence, without reinterpreting them as rows.
            score <- if (is_glm) r*w[rows]*dmu/f$variance(yhat) else r*w[rows]
            ic <- (X*as.numeric(score)) %*% gof_information_inverse(model,length(beta))
        } else {
            fit <- model
            fit$na.action <- NULL
            ic <- lava::iid(fit)[rows,,drop=FALSE]
        }
    }

    needed <- setdiff(variable,c("predicted",colnames(Xfit)))
    extra <- if (length(needed)) gof_fitted_data(data,mf) else NULL
    unknown <- setdiff(needed,names(extra))
    if (length(unknown)) stop("Variable not found: ",paste(unknown,collapse=", "))
    values <- list()
    for (v in variable) {
        z <- if (v=="predicted") yhat else if (v %in% colnames(Xfit)) X[,v] else extra[[v]][rows]
        if (is.ordered(z)) z <- as.numeric(z)
        if (!is.numeric(z) || is.complex(z) || !is.null(dim(z)) || any(!is.finite(z)))
            stop("Ordering variable '",v,"' must be a finite numeric vector or ordered factor")
        # Preserve automatic exclusion of intercepts and binary design columns.
        # An explicitly requested numeric variable can include any tie count.
        if (!automatic_variables || v=="predicted" || length(unique(z))>2L) values[[v]] <- z
    }
    if (!length(values)) stop("No selected ordering variables")
    ordering <- do.call(cbind,values)
    colnames(ordering) <- names(values)
    res <- cumres.default(NULL,variable=ordering,r=r,dr=dr,ic=ic,coef=NULL,
                          R=R,plots=plots,b=b,...)
    res$model <- class(model)[1]
    res$response <- response
    res$nobs <- sum(active)
    canonical <- !is_glm || (family(model)$family=="gaussian" && family(model)$link=="identity") ||
        (family(model)$family=="binomial" && family(model)$link=="logit") ||
        (family(model)$family=="poisson" && family(model)$link=="log")
    precision <- if (grouped) rep(1,length(rows)) else w[rows]
    if (!canonical) precision <- precision*as.numeric(dmu)/f$variance(yhat)
    res$informative <- gof_informative_processes(X,precision,ordering,res$x,b)
    for (i in which(res$informative %in% FALSE)) {
        res$W[[i]][] <- 0
        res$What[[i]][] <- 0
        res$KS[i] <- res$CvM[i] <- 1
    }
    indeterminate <- which(is.na(res$informative))
    res$KS[indeterminate] <- res$CvM[indeterminate] <- NA_real_
    res
}

