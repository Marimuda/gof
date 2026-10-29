\name{cumres.default}
\alias{cumres.default}
\title{Cumulative residual processes from supplied residuals and influences}
\usage{
\method{cumres}{default}(model, variable, r = stats::residuals(model), dr,
  R = 1e3L, b = 0, plots = min(R, 50L), coef = stats::coef(model),
  ic = lava::iid(model), breakties = 0, subset = subset_t,
  subset.max = 500L, ...)
}
\arguments{
\item{model}{Model providing coefficients, residuals and influences, or \code{NULL}
when these quantities are supplied explicitly.}
\item{variable}{Finite real numeric ordering vector or matrix with one row per
independent residual unit. Tied values are evaluated as complete groups.}
\item{r}{Real residual vector, or function of \code{coef} and \code{...} returning
that vector. A \code{grad} attribute can supply the residual gradient.}
\item{dr}{Residual gradient matrix, with one column per fitted parameter, or a
function of \code{coef} and \code{...}. If omitted, use the residual's gradient
attribute or differentiate the residual function with \code{numDeriv}.}
\item{R}{Positive whole number of multiplier simulations.}
\item{b}{Nonnegative bandwidth; zero gives cumulative residuals and positive
values give sums over the window \eqn{(t-b,t]}.}
\item{plots}{Nonnegative whole number of additional null processes to retain.}
\item{coef}{Finite real coefficient vector. An empty vector describes a model
with no estimated coefficients and no projection adjustment.}
\item{ic}{Finite real coefficient influence matrix with one row per residual
unit and the same parameter columns as \code{dr}. Matched zero-column matrices
are supported; when no coefficients/gradient columns exist they default to zero.}
\item{breakties}{Deprecated jitter width. Nonzero values warn; ties are retained.}
\item{subset}{Function receiving sorted distinct values and \code{n.qt},
\code{n.eq}, and returning unique valid 1-based evaluation positions.}
\item{subset.max}{Positive whole-number cap on distinct evaluation values per
process. \code{Inf} or \code{NULL} evaluates all distinct values.}
\item{...}{Additional arguments passed to the residual and gradient callbacks.}
}
\description{Constructs cumulative residual diagnostics using supplied residual,
gradient and influence contributions. The caller must justify conditional
centering and the independent-unit influence representation for custom models.}
\value{A \code{cumres} object with observed/simulated processes, grids and Sup/L2
Monte Carlo p-values. Comparisons include equality.}
\seealso{\code{\link{cumres.glm}}, \code{\link{cumres.lvmfit}}}
