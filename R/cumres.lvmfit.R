cumresMaxSet <- function(m,var,...) {
  ii <- lava::index(m)
  if (length(var)!=1L || !is.character(var) || is.na(var) || !var %in% ii$vars)
      stop("Select one named structural residual for full prediction")
  ## In V=(I-A')^{-1} epsilon, Cov(epsilon,V)=P(I-A)^{-1}.
  ## Covariance paths determine zeros for every admissible parameter value.
  ## Boolean closure works even when a matrix of 0/1 edges is singular.
  reach <- ii$A!=0
  diag(reach) <- TRUE
  for (k in seq_len(nrow(reach)))
      reach <- reach | outer(reach[,k],reach[k,],"&")
  related <- (ii$P!=0)[var,,drop=FALSE] %*% (reach*1)>0
  candidates <- intersect(ii$manifest,colnames(reach))
  candidates[!as.logical(related[1,candidates])]
}


##' Cumulative residual processes for structural equation models
##' 
##' Calculates GoF statistics based on cumulative residual processes for
##' structural equation models fitted with the \code{lava} package.
##' 
##' With \code{y} and \code{x} given as functions the user can decide which
##' variables to use in the prediction of the outcome and predictor (use the
##' \code{predict} method as below).
##' 
##' @param model A fitted \code{lvmfit} object from \code{lava::estimate}.
##' @param y A formula specifying the association to be checked. Alternatively
##' the outcome specified as a function or a string with the name of the outcome
##' in the model. Lists of formulas retain a list of results. Multiple named
##' outcomes or lists of residual functions are crossed with the requested
##' predictors to produce separately named diagnostics. A residual function
##' may attach a matrix or function-valued \code{grad} attribute.
##' @param x Predictor name, numeric vector, or function of \code{p}. Character
##' vectors and lists request separate diagnostics for each predictor. Explicit
##' numeric/function predictors must satisfy conditional centering of the
##' selected residual, as well as the estimated-ordering conditions in the
##' details; pointwise conditional centering alone is insufficient.
##' @param full If FALSE, named latent predictors use their conditional means
##' given the exogenous covariates; observed predictors use their data values.
##' If TRUE, named predictors are projected onto all manifest variables with no
##' structural covariance path from the selected residual's error. These
##' variables are guaranteed uncorrelated with that error for every admissible
##' parameter value. The selection is conservative in models with exact path
##' cancellations. Supplied numeric/function predictors are used directly.
##' @details This method implements mean/linearity residual-process diagnostics;
##' use \code{\link{cumres_normality}} for residual-distribution diagnostics
##' from the SEM reference. For a parameter-dependent callback \eqn{x_\theta}, its indicator
##' derivative contributes unless
##' \eqn{E[r_{\theta_0}\partial_\theta x_\theta\mid x_{\theta_0}=t]=0}.
##' A sufficient construction uses the same fixed admissible conditioning data
##' \eqn{Z} for the ordering and its parameter derivatives, with
##' \eqn{E[r_{\theta_0}\mid Z]=0}. Named predictions use structural safe sets;
##' arbitrary callbacks are not automatically certified. Smoothness, moment
##' and stochastic-equicontinuity conditions for estimated orderings are also
##' required. The adapter differentiates residuals, not arbitrary ordering
##' callbacks. The caller must justify these stronger conditions or supply a
##' different, fully derived influence process through the default interface.
##' @param data Complete finite numeric model data (default \code{model.frame(model)}).
##' Residuals, predictions, scores and information all use these same rows.
##' @param p Optional parameter vector, with names matching the fitted
##' coefficients. It is used consistently in all diagnostic quantities; no
##' refitting is performed. Calibrated inference assumes these parameters were
##' estimated from the supplied data using the model's estimation procedure.
##' @param R Positive whole number of multiplier simulations.
##' @param b Nonnegative moving-window bandwidth. Zero gives cumulative
##' residuals; positive values give residual sums over \eqn{(t-b,t]}.
##' @param plots Nonnegative whole number of null processes to retain for plotting.
##' @param ... Additional arguments passed to \code{cumres.default}, including
##' \code{subset.max}, which limits grid values without discarding data rows.
##' @method cumres lvmfit
##' @export
##' @returns A \code{cumres} object for one residual/predictor pair, a named
##' list for multiple pairs, or a list of results when \code{y} is a list of
##' formulas. See \code{\link{cumres}} for each diagnostic's common result
##' components; normal/GLM-specific metadata are not added by this method.
##' @author Klaus K. Holst
##' @note The residual process assumes a Gaussian SEM with independent subjects
##' and complete data, a positive semidefinite structural error covariance and
##' an invertible structural system. With \code{full=FALSE}, named observed predictors with a
##' covariance path from the selected structural error are rejected; use
##' \code{full=TRUE} to construct an orthogonal prediction instead.
##' Declared exogenous variables must have structurally zero error covariance
##' with endogenous/latent variables; contradictory forced declarations are
##' rejected because all lava predictions condition on the exogenous set.
##' Clustered, weighted and incomplete-data SEM fits need different influence
##' processes and are rejected explicitly. Collections are evaluated directly
##' from the supplied objects, including models fitted inside functions.
##' @references B.N. Sanchez and E. A. Houseman and L. M. Ryan (2009)
##' \emph{Residual-Based Diagnostics for Structural Equation Models}. Biometrics
##' Volume 65 (1), pp 104-115.
##' @keywords models regression
##' @seealso \code{\link{cumres}}, \code{\link{cumres.default}},
##' \code{\link{plot.cumres}}.
##' @examples
##' if (requireNamespace("numDeriv", quietly = TRUE)) {
##'   set.seed(14)
##'   model <- lava::lvm(y ~ x)
##'   d <- lava::sim(model, 100)
##'   fit <- lava::estimate(model, d)
##'   diagnostic <- cumres(fit, y ~ x, R = 49, plots = 3)
##'   coef(diagnostic)
##'
##'   # A named collection retains one diagnostic for each formula.
##'   diagnostics <- cumres(fit, list(check_x = y ~ x), R = 49, plots = 0)
##'   names(diagnostics)
##' }
cumres.lvmfit <- function(model,y,x,full=FALSE,
                   data=model.frame(model),p,
                   R=1000, b=0, plots=min(R,50),
                   ...) {
  if (!inherits(model,"lvmfit")) stop("Supply a fitted lava structural-equation model")
  if (!requireNamespace("numDeriv",quietly=TRUE)) stop("package 'numDeriv' not available")
  if (!is.logical(full) || length(full)!=1L || is.na(full)) stop("full must be TRUE or FALSE")
  if (!is.null(model$opt$convergence) && !isTRUE(model$opt$convergence==0L))
      stop("Model has not converged")
  if (!is.null(model$cluster)) stop("Structural-equation diagnostics require independent subjects, not clustered fits")
  if (!is.null(model$weights)) stop("Weighted structural-equation fits require a different residual process and are not supported")
  if (missing(y)) stop("Supply a residual y or a formula y ~ x")
  fitted.p <- lava::pars(model)
  if (length(fitted.p) && (length(model$estimator)!=1L || !is.character(model$estimator) ||
      is.na(model$estimator) || !model$estimator %in% c("gaussian","gaussian1","gaussian2","normal")))
      stop("Structural-equation diagnostics require a Gaussian fitted model")
  if (missing(p)) p <- fitted.p
  if (!is.numeric(p) || is.complex(p) || !is.null(dim(p)) ||
      length(p)!=length(fitted.p) || any(!is.finite(p)))
      stop("p must contain one finite parameter per fitted coefficient")
  if (!is.null(names(p))) {
      if (anyDuplicated(names(p)) || !setequal(names(p),names(fitted.p)))
          stop("p names must match the fitted coefficients")
      p <- p[names(fitted.p)]
  } else names(p) <- names(fitted.p)
  data <- as.data.frame(data)
  required <- names(model.frame(model))
  if (!all(required %in% names(data))) stop("data must contain all manifest model variables")
  if (!nrow(data) || any(!vapply(data[required],function(z) is.numeric(z) && !is.complex(z),logical(1))) ||
      any(!is.finite(as.matrix(data[required]))))
      stop("Structural-equation diagnostics require complete finite numeric model data")
  ii <- lava::index(model)
  other <- setdiff(seq_along(ii$vars),ii$exo.idx)
  if (length(ii$exo.idx) && length(other) &&
      any(ii$P[ii$exo.idx,other,drop=FALSE]!=0))
      stop("Declared exogenous variables have structural error covariance with other variables; ",
           "fit without a forced exogeneity override")
  moment <- tryCatch(lava::moments(model,p=p,data=data),error=function(e)
      stop("Invalid Gaussian structural moments: ",conditionMessage(e),call.=FALSE))
  covariance <- as.matrix(moment$P)
  if (!length(covariance) || any(!is.finite(covariance)) || any(diag(covariance)<0) ||
      !isSymmetric(covariance)) stop("Gaussian structural error covariance must be finite and positive semidefinite")
  ## Only allow negative eigenvalues attributable to floating-point backward
  ## error; do not choose a data/seed-dependent near-zero covariance threshold.
  roundoff <- 100*.Machine$double.eps*nrow(covariance)*max(abs(covariance))
  if (min(eigen(covariance,symmetric=TRUE,only.values=TRUE)$values)< -roundoff)
      stop("Gaussian structural error covariance must be positive semidefinite")
  if (isTRUE(attr(moment$IAi,"pseudo")))
      stop("Gaussian structural system must be invertible; a pseudoinverse is not a valid generative SEM")
  dots <- list(...)
  if ("simulation" %in% names(dots) && !identical(dots$simulation,"marginal"))
      stop("Joint ordering simulation supports the default and linear/GLM methods; SEM collections remain separate diagnostics")
  ## IC.default evaluates information at pars(model). Update a copy so explicit
  ## p/data options affect scores, information and residuals consistently.
  fit <- model
  fit$opt$estimate <- fit$opt$par <- p
  fit$data$model.frame <- data
  fit$data$n <- nrow(data)
  influence <- if (!length(p)) matrix(numeric(),nrow(data),0L) else
      if (identical(p,fitted.p) && identical(data,model.frame(model)))
      lava::iid(model) else lava::iid(fit,data=data)
  influence <- as.matrix(influence)
  one <- function(yy,xx) {
      if (!is.function(yy) && (!is.character(yy) || length(yy)!=1L || is.na(yy) ||
          !yy %in% c(lava::endogenous(model),lava::latent(model))))
          stop("Select one endogenous or latent residual, or a residual function")
      if (is.character(xx)) {
          if (length(xx)!=1L || is.na(xx) || !xx %in% lava::vars(model))
              stop("Predictor not found in the structural model")
          if (full) {
              if (!is.character(yy)) stop("full prediction requires a named structural residual")
              by <- cumresMaxSet(model,yy)
              prediction <- if (length(by)) predict(model,x=by,y=xx,p=p,data=data) else
                  predict(model,x=~1,p=p,data=data)
              x0 <- prediction[,xx]
          } else if (xx %in% lava::latent(model)) {
              x0 <- predict(model,x=~1,p=p,data=data)[,xx]
          } else {
              if (is.character(yy) && !xx %in% cumresMaxSet(model,yy))
                  stop("Observed predictor may correlate with the structural residual; use full=TRUE ",
                       "or supply a predictor with a justified conditional-centering property")
              x0 <- data[[xx]]
          }
      } else x0 <- if (is.function(xx)) xx(p) else xx
      if (!is.numeric(x0) || is.complex(x0) || length(x0)!=nrow(data) ||
          (!is.null(dim(x0)) && NCOL(x0)!=1L) || any(!is.finite(x0)))
          stop("Each predictor must return one finite numeric value per data row")
      residual <- if (is.function(yy)) yy else function(q) {
          rr <- predict(model,x=lava::vars(model),residual=TRUE,p=q,data=data)
          rr[,yy]
      }
      r <- residual(p)
      if (!is.numeric(r) || is.complex(r) || length(r)!=nrow(data) ||
          (!is.null(dim(r)) && NCOL(r)!=1L) || any(!is.finite(r)))
          stop("Each residual must return one finite numeric value per data row")
      grad <- attr(r,"grad")
      if (is.null(grad)) grad <- if (!length(p)) matrix(numeric(),nrow(data),0L) else
          numDeriv::jacobian(residual,p,method=lava::lava.options()$Dmethod)
      if (is.function(grad)) grad <- grad(p)
      variable <- matrix(as.numeric(x0),ncol=1L,dimnames=list(NULL,if (is.character(xx)) xx else "x"))
      result <- do.call(cumres.default,c(list(model=NULL,variable=variable,r=r,dr=grad,
                          ic=influence,coef=NULL,R=R,plots=plots,b=b),dots))
      result$type <- rep("sem",length(result$W))
      result$model <- "lvmfit"
      result
  }
  collection <- function(yy,xx) {
      if (inherits(yy,"formula")) {
          outcome <- lava::getoutcome(yy)
          if (!length(outcome)) stop("Formula must name one or more residuals on the left")
          raw.predictor <- attr(outcome,"x")
          predictor <- if (!length(raw.predictor)) rep(1,nrow(data)) else
              lava::decomp.specials(raw.predictor)
          return(collection(lava::decomp.specials(outcome),
                            predictor))
      }
      if (is.list(yy) && length(yy) && all(vapply(yy,inherits,logical(1),"formula")))
          return(lapply(yy,function(f) collection(f,NULL)))
      ys <- if (is.list(yy) || (is.character(yy) && length(yy)>1L)) as.list(yy) else list(yy)
      xs <- if (is.list(xx) || (is.character(xx) && length(xx)>1L)) as.list(xx) else list(xx)
      if (!length(ys) || !length(xs) || is.null(xx)) stop("Supply nonempty residuals and predictors")
      if (length(ys)==1L && length(xs)==1L) return(one(ys[[1]],xs[[1]]))
      out <- list()
      for (i in seq_along(ys)) for (j in seq_along(xs)) {
          label.y <- if (is.character(ys[[i]])) ys[[i]] else paste0("y",i)
          label.x <- if (is.character(xs[[j]])) xs[[j]] else paste0("x",j)
          out[[length(out)+1L]] <- one(ys[[i]],xs[[j]])
          names(out)[length(out)] <- paste(label.y,label.x,sep="<-")
      }
      out
  }
  collection(y,if (missing(x)) NULL else x)
}
