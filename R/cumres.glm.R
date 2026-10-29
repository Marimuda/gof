##' @export
`cumres.lm` <- function(model,...) {
  cumres.glm(model,...)
}

##' Calculates GoF statistics based on cumulative residual processes
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
##' @aliases cumres cumres.glm cumres.lm
##' @param model Model object (\code{lm} or \code{glm})
##' @param variable List of variable to order the residuals after
##' @param data data.frame used to fit model (complete cases)
##' @param R Number of samples used in simulation
##' @param b Nonnegative moving-window bandwidth. Zero gives cumulative
##' residuals; positive values give residual sums over the window \eqn{(t-b,t]}.
##' @param plots Number of realizations to save for use in the plot-routine
##' @param ... additional arguments passed to \code{cumres.default}. In
##' particular \code{subset.max} (default 500): each process is evaluated at
##' at most this many distinct covariate values, retaining every observation
##' in the fitted model. Use \code{subset.max=Inf} or \code{NULL} for all
##' distinct values. Tied observations are accumulated together. The deprecated
##' \code{breakties} argument is accepted with a warning when nonzero; no jitter
##' is added. A custom \code{subset} function receives sorted distinct values
##' and \code{n.qt}, \code{n.eq}, and returns their 1-based evaluation positions.
##' @return Returns an object of class 'cumres'.
##' @note Unweighted normal, binomial and Poisson models are supported;
##' noncanonical binomial links are also tested. Grouped binomial responses
##' must contain integer success/failure counts. Diagnostics use the original
##' maximum-likelihood fit and do not refit grouped data. Rows omitted during
##' fitting are also omitted by the diagnostic. Monte Carlo p-values count
##' simulated statistics greater than or equal to the observed statistic.
##' @author Klaus K. Holst
##' @seealso The \href{https://CRAN.R-project.org/package=timereg}{timereg}
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
##' 
##' sim1 <- function(n=100, f=function(x1,x2) {10+x1+x2^2}, sd=1, seed=1) {
##'   if (!is.null(seed))
##'     set.seed(seed)
##'   x1 <- rnorm(n);
##'   x2 <- rnorm(n)
##'   X <- cbind(1,x1,x2)
##'   y <- f(x1,x2) + rnorm(n,sd=sd)
##'   d <- data.frame(y,x1,x2)
##'   return(d)
##' }
##' d <- sim1(100); l <- lm(y ~ x1 + x2,d)
##' system.time(g <- cumres(l, R=100, plots=50))
##' g
##' \donttest{plot(g)}
##' g1 <- cumres(l, c("y"), R=100, plots=50)
##' g1
##' g2 <- cumres(l, c("y"), R=100, plots=50, b=0.5)
##' g2
##' 
`cumres.glm` <- function(model,
                  variable=c("predicted",colnames(model.matrix(model))),
                  data=data.frame(model.matrix(model)), 
                  R=1000, b=0, plots=min(R,50),
                  ...) {
  if (!inherits(model,"lm")) stop("Unsupported model!")
  if (inherits(model,"glm") && !family(model)$family %in% c("gaussian","binomial","poisson"))
      stop("Unsupported family: ",family(model)$family)
  if (anyNA(coef(model))) stop("Over-parametrized model")
  if (inherits(model,"glm") && !isTRUE(model$converged)) stop("Model has not converged")
  mf <- model.frame(model)
  if (!is.null(model.weights(mf))) stop("Weight argument not supported")
  X <- model.matrix(model)
  response <- all.vars(formula(model))[1]
  force(variable)
  y <- model.response(mf)
  yhat <- as.numeric(model$fitted.values)
  if (is.matrix(y)) {
      if (!inherits(model,"glm") || family(model)$family!="binomial" || ncol(y)!=2L)
          stop("Only two-column binomial matrix responses are supported")
      if (any(!is.finite(y)) || any(y<0 | y!=floor(y)))
          stop("Grouped binomial responses must be nonnegative integer counts")
      rows <- rep(seq_len(nrow(y)),rowSums(y))
      y0 <- unlist(lapply(seq_len(nrow(y)),function(i) rep(c(1,0),y[i,])))
      ## Keep the actual fit, including its coefficients, controls and method.
      ## Expand scores per Bernoulli trial; the coefficient covariance from
      ## the grouped likelihood is the same as for the expanded likelihood.
      X <- X[rows,,drop=FALSE]
      yhat <- yhat[rows]
      eta <- model$linear.predictors[rows]
      f <- family(model)
      dmu <- f$mu.eta(eta)
      r <- y0-yhat
      dr <- -X*as.numeric(dmu)
      ic <- (X*as.numeric(r*dmu/f$variance(yhat))) %*% vcov(model)
  } else {
      ## Stored fit vectors contain complete cases, unlike residuals()/predict()
      ## which pad omitted observations with NA for na.exclude fits.
      if (inherits(model,"glm")) {
          r <- if (is.factor(y)) as.numeric(y!=levels(y)[1])-yhat else as.numeric(y)-yhat
          dr <- -X*as.numeric(family(model)$mu.eta(model$linear.predictors))
      } else {
          r <- as.numeric(model$residuals)
          dr <- -X
      }
      fit <- model
      fit$na.action <- NULL
      ic <- lava::iid(fit)
  }

  if (!is.na(match(response, variable))) variable[match(response, variable)] <- "predicted"
  variable <- unique(variable)
  if (any(!variable %in% c("predicted",colnames(X))))
      stop("Variable not found: ",paste(setdiff(variable,c("predicted",colnames(X))),collapse=", "))
  UsedData <- X[,na.omit(match(variable, colnames(X))),drop=FALSE]
  myvars <- colnames(UsedData)[apply(UsedData,2,function(x) length(unique(x))>2)] ## Only consider variables with more than two levels
  if ("predicted"%in%variable) myvars <- c("predicted",myvars)  
  if (!length(myvars)) stop("No selected variable has more than two distinct values")

  variable <- c()
  for (v in myvars) {
    x <- NULL
    if (v=="predicted") {
        x <- yhat
    } else  if (v %in% colnames(X)) {
        x <- X[,v]    
    } else warning("Variable '", v , "' not found.\n",sep="")
    if (!is.null(x)) {
        variable <- cbind(variable, x)
        colnames(variable)[ncol(variable)] <- v
    }
}
  
  res <- cumres.default(NULL, variable=variable,
                        r=r, dr=dr, ic=ic, coef=NULL,
                        R=R, plots=plots, b=b, ...)
  
  return(res)
}

