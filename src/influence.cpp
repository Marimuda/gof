/*!
  @file influence.cpp
  @author Klaus K. Holst
  @copyright 2018-2020, Klaus Kähler Holst

  @brief Utility functions for Generalized Linear Models
*/
#include <target/influence.hpp>
#include <cmath>
#include <stdexcept>

namespace target {
namespace {
void validate_iid(const arma::vec &y, const arma::vec &p, const arma::mat &x, const arma::vec &w) {
    if (y.is_empty() || x.n_cols == 0 || y.n_elem != p.n_elem || y.n_elem != x.n_rows ||
        y.n_elem != w.n_elem || !y.is_finite() || !p.is_finite() || !x.is_finite() ||
        !w.is_finite() || arma::any(w < 0) || !arma::any(w > 0))
        throw std::invalid_argument(
            "Supply matching finite outcomes, predictions, design and nonnegative informative weights");
}

arma::mat inverse_information(const arma::mat &information) {
    arma::mat inverse;
    if (!information.is_finite() || !arma::inv_sympd(inverse, information) || !inverse.is_finite())
        throw std::domain_error("Information matrix must be finite and positive definite");
    return inverse;
}
} // namespace

IID::IID(arma::mat score, arma::mat v) {
    if (v.n_rows != v.n_cols || score.n_cols != v.n_rows || !score.is_finite() || !v.is_finite())
        throw std::invalid_argument(
            "Score and covariance dimensions must match and contain finite values");
    if (v.n_rows == 0)
        iid.set_size(score.n_rows, 0);
    else
        iid = score * v;
    vcov = v;
    if (!iid.is_finite())
        throw std::domain_error("Influence computation overflowed");
}

IID logistic_iid(const arma::vec &y, const arma::vec &p, const arma::mat &x, const arma::vec &w) {
    validate_iid(y, p, x, w);
    if (arma::any(y < 0) || arma::any(y > 1) || arma::any(p < 0) || arma::any(p > 1))
        throw std::invalid_argument("Logistic outcomes and probabilities must be in [0,1]");
    const arma::vec s = (y - p) % w;
    arma::mat U = x;
    for (arma::uword i = 0; i < x.n_cols; ++i)
        U.col(i) %= s;
    const arma::vec v = p % (1 - p) % w;
    arma::mat H = arma::zeros(x.n_cols, x.n_cols);
    for (arma::uword i = 0; i < x.n_rows; ++i)
        H += v(i) * (x.row(i).t() * x.row(i));
    return IID(U, inverse_information(H));
}

IID linear_iid(const arma::vec &y, const arma::vec &p, const arma::mat &x, const arma::vec &w) {
    validate_iid(y, p, x, w);
    const arma::uword active = arma::accu(w > 0);
    if (active <= x.n_cols)
        throw std::invalid_argument("Linear influence needs positive residual degrees of freedom");
    const double df = static_cast<double>(active) - static_cast<double>(x.n_cols);
    const arma::vec r = y - p;
    const arma::vec s = r % w;
    const double sigma2 = arma::dot(r, s) / df;
    if (!std::isfinite(sigma2))
        throw std::domain_error("Residual variance overflowed");
    arma::mat U = x;
    for (arma::uword i = 0; i < x.n_cols; ++i)
        U.col(i) %= s;
    arma::mat H = arma::zeros(x.n_cols, x.n_cols);
    for (arma::uword i = 0; i < x.n_rows; ++i)
        H += w(i) * (x.row(i).t() * x.row(i));
    // sigma2 cancels in score * inverse information. Avoid 0/0 for exact fits.
    IID result(U, inverse_information(H));
    result.vcov *= sigma2;
    return result;
}
} // namespace target
