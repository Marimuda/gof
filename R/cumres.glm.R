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
##' @param resolution Optional named numeric vector of positive finite measurement
##' steps for selected observed ordering variables, for example
##' \code{c(x=0.25)}. These coordinates become \code{round(x / step) * step},
##' with origin zero and R's rounding of represented values to the nearest even
##' integer at exact half steps. The model and residual units are unchanged.
##' Choose the resolution from the recording convention before inspecting test
##' results. Fitted predictions cannot be coarsened by this argument.
##' @param on_unresolved Policy for a fitted ordering flagged by the numerical
##' check. \code{"report"} (default) warns, retains the numerical curves and
##' statistics, and returns missing inferential p-values for that ordering.
##' \code{"error"} stops before simulation. Other selected marginal orderings
##' remain available in report mode; a joint family containing an unresolved
##' member has no inferential p-values or simultaneous envelope. These controls
##' are implemented for linear/GLM methods; \code{resolution} also supports
##' \code{cumres_bootstrap()}.
##' @param ... additional arguments passed to \code{cumres.default}. In
##' particular \code{subset.max} (default 500): each process is evaluated at
##' at most this many distinct covariate values, retaining every observation
##' in the fitted model. Use \code{subset.max=Inf} or \code{NULL} for all
##' distinct values. Tied observations are accumulated together. The deprecated
##' \code{breakties} argument is accepted with a warning when nonzero; no jitter
##' is added. A custom \code{subset} function receives sorted distinct values
##' and \code{n.qt}, \code{n.eq}, and returns their 1-based evaluation positions.
##' Named additional controls must unambiguously match \code{breakties},
##' \code{subset}, \code{subset.max} or \code{simulation}; unused names are
##' rejected so misspelled ordering controls cannot be silently ignored.
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
##' This does not certify ties between different design rows at the exact
##' likelihood solution. R's GLM convergence criterion concerns deviance,
##' not the accuracy of every fitted-value comparison. With a flat or nearly
##' flat fitted mean, numerical fitting error can split an exact tie, merge
##' distinct values or reverse their order; the Sup statistic can change
##' substantially. A local Fisher-score correction and accumulated-roundoff
##' screen flag comparisons whose sign is not numerically resolved. For moving
##' windows the screen also compares entry/exit events, including a fitted
##' coordinate against another coordinate plus the bandwidth. The
##' coordinates and supplied coefficients are retained. A passed check is an
##' ordinary numerical diagnostic, not an exact-root or interval-arithmetic
##' certificate, nor a guarantee of calibration. Prespecified
##' fixed covariate orderings avoid this estimated-coordinate ambiguity.
##' Arbitrary rounding of predictions or tighter fitting tolerances alone
##' do not establish correct fitted ordering or statistical calibration.
##' The \code{ordering} component records one row per ordering: source, grouped
##' tie rule, declared resolution, counts of distinct values/tie blocks/tied
##' residual units, changed coordinates, check status and reason. Status is
##' \code{"fixed"} for observed coordinates, \code{"specified"} for a known
##' mean with no estimated coefficients, \code{"structural"} for identical
##' design/offset rows, \code{"saturated"} when ordering cannot add lack-of-fit
##' information across saturated natural-parameter cells, \code{"checked"} for
##' other fitted orderings that pass the local screen, or \code{"unresolved"}.
##' This status is separate from the \code{informative} score-space flag.
##' Diagnostic arithmetic details are retained in \code{ordering_check}.
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
##' @section Coordinate convention:
##' The \code{"predicted"} coordinate is the fitted response mean. Lin, Wei
##' and Ying (2002), equation (2.10), instead use the fitted linear predictor.
##' A strictly increasing change of coordinate preserves the all-point
##' cumulative Sup statistic, but generally changes L2's integration measure
##' and a fixed-bandwidth moving-window statistic. The coordinate scale is
##' therefore part of the diagnostic specification. The half-open window
##' \eqn{(t-b,t]} agrees with their equation (2.7).
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
##'
##' # Declare a recording resolution for the diagnostic coordinate only.
##' recorded <- cumres(fit, variable = "x", resolution = c(x = 0.25),
##'   R = 39, plots = 0)
##' recorded$ordering
##' # A strict workflow can stop on unresolved fitted ordering.
##' cumres(fit, variable = "predicted", on_unresolved = "error", R = 39,
##'   plots = 0)$ordering
`cumres.glm` <- function(model,
                  variable=c("predicted",colnames(model.matrix(model))),
                  data=model.frame(model),
                  R=1000,b=0,plots=min(R,50),...,
                  resolution=NULL,on_unresolved=c("report","error")) {
    on_unresolved <- tryCatch(match.arg(on_unresolved,c("report","error")),
        error=function(e) stop("on_unresolved must be 'report' or 'error'",call.=FALSE))
    if (!is.numeric(b) || is.complex(b) || !is.null(dim(b)) || length(b)!=1L || !is.finite(b) || b<0)
        stop("b must be a finite nonnegative bandwidth",call.=FALSE)
    extra_names <- names(substitute(list(...)))[-1L]
    extra_names <- extra_names[nzchar(extra_names)]
    unknown <- extra_names[is.na(pmatch(extra_names,c("breakties","subset","subset.max","simulation"),duplicates.ok=TRUE))]
    if (length(unknown)) stop("Unused diagnostic argument(s): ",paste(unknown,collapse=", "),call.=FALSE)
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
    # so identical design rows remain tied without rounding data. This does
    # not establish exact-MLE ties between different design rows.
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
    coordinates <- gof_ordering_coordinates(ordering,resolution)
    ordering <- coordinates$values
    predicted <- match("predicted",colnames(ordering))
    ordering_check <- NULL
    if (!is.na(predicted)) {
        offset <- model$offset
        if (is.null(offset)) offset <- rep(0,nrow(Xfit))
        success <- if (grouped) counts[,1L] else as.numeric(y)*w
        exposure <- if (grouped) rowSums(counts) else w
        ordering_check <- gof_fitted_ordering_check(Xfit[active,,drop=FALSE],beta,
            offset[active],yhat_fit[active],success[active],exposure[active],
            if (is_glm) family(model) else stats::gaussian(),b=b)
        coordinates$metadata$source[predicted] <- "fitted"
        coordinates$metadata$status[predicted] <- ordering_check$status
        coordinates$metadata$reason[predicted] <- ordering_check$reason
    }
    gof_ordering_report(coordinates$metadata,on_unresolved)
    res <- cumres.default(NULL,variable=ordering,r=r,dr=dr,ic=ic,coef=NULL,
                          R=R,plots=plots,b=b,...)
    res$model <- class(model)[1]
    res$response <- response
    res$nobs <- sum(active)
    res$ordering <- coordinates$metadata
    res$ordering_check <- ordering_check$details
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
    unresolved <- gof_ordering_unresolved(res)
    res$KS[unresolved] <- res$CvM[unresolved] <- NA_real_
    if (inherits(res,"cumres_joint")) res <- gof_joint_finalize(res)
    res
}

