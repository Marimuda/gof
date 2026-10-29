##' Model diagnostics based on cumulative residuals
##'
##' The \pkg{gof} package checks fitted mean structures using cumulative
##' residual processes, multiplier simulations and Sup/L2 tests. Start with
##' \code{\link{cumres}} for a model-to-diagnostic workflow, then use the
##' model-specific help for its assumptions and available options.
##'
##' Supported adapters cover normal/binomial/Poisson regression, Gaussian
##' structural equation models and ordinal regression fitted with \pkg{lava}.
##' The low-level default method accepts custom residual and influence
##' contributions when their statistical assumptions can be justified.
##' \code{\link{cumres_joint}} provides simultaneous ordering-variable inference;
##' \code{\link{cumres_normality}} checks Gaussian SEM residual distributions;
##' \code{\link{cumres_bootstrap}} supplies parametric calibration for canonical
##' binomial and Poisson residual diagnostics.
##' @author Klaus K. Holst <klaus@@holst.it>
##' @seealso \code{\link{cumres}}, \code{\link{cumres.glm}},
##' \code{\link{cumres.lvmfit}}, \code{\link{cumres.ordreg}},
##' \code{\link{cumres_joint}}, \code{\link{cumres_normality}},
##' \code{\link{cumres_bootstrap}},
##' \code{\link{plot.cumres}}, \code{\link{surgunit}}.
##' @references D. Y. Lin, L. J. Wei and Z. Ying (2002).
##' \emph{Model-Checking Techniques Based on Cumulative Residuals}.
##' Biometrics 58, 1--12. \doi{10.1111/j.0006-341X.2002.00001.x}.
##'
##' John Q. Su and L. J. Wei (1991). \emph{A lack-of-fit test for the mean
##' function in a generalized linear model}. Journal of the American
##' Statistical Association 86, 420--426.
##' @keywords internal
##' @importFrom stats quantile runif coef family formula model.frame model.matrix
##' @importFrom stats model.response model.weights vcov predict rnorm
##' @importFrom utils head tail
##' @useDynLib gof, .registration = TRUE
##' @importFrom Rcpp loadModule
##' @importFrom methods new
"_PACKAGE"

loadModule("gofmod", TRUE)


##' Surgical unit data
##'
##' Survival time and four covariates for 54 patients undergoing liver surgery.
##' This is the surgical unit example used by Lin, Wei and Ying (2002).
##' @name surgunit
##' @docType data
##' @format A data frame with 54 rows and five variables:
##' \describe{
##'   \item{\code{bloodclot}}{Numeric blood-clotting measure.}
##'   \item{\code{prognostic}}{Integer prognostic measure.}
##'   \item{\code{enzyme}}{Integer enzyme measure.}
##'   \item{\code{liver function}}{Numeric liver-function measure. Its name
##'   contains a space; use backticks in formulas or \code{[[} for extraction.}
##'   \item{\code{survival}}{Integer survival-time measure.}
##' }
##' @details Measurement units and fuller variable definitions are not stored
##' as metadata in the distributed data. Consult the cited original sources
##' before interpreting these measurements in substantive analyses.
##' @source Neter, J., Kutner, M. H., Nachtsheim, C. J., and Wasserman, W.
##' (1996), \emph{Applied Linear Statistical Models}, 4th edition. Chicago: Irwin.
##' @references D. Y. Lin, L. J. Wei and Z. Ying (2002).
##' \emph{Model-Checking Techniques Based on Cumulative Residuals}.
##' Biometrics 58, 1--12. \doi{10.1111/j.0006-341X.2002.00001.x}.
##' @examples
##' data(surgunit)
##' names(surgunit)
##' head(surgunit[["liver function"]])
NULL

