/*!
  @file probability.cpp
  @author Klaus K. Holst
  @copyright 2018-2020, Klaus Kähler Holst

  @brief Focused component extracted without changing its algorithms.
*/
#include <target/probability.hpp>
#include <cmath>
#include <complex>
#include <stdexcept>

namespace target {
// Subtract a row maximum before exponentiating to avoid overflow.
arma::mat softmax(arma::mat lp, bool ref, bool log) {
    if (!lp.is_finite())
        throw std::invalid_argument("Softmax logits must be finite");
    if (ref)
        lp.insert_cols(0, arma::zeros(lp.n_rows));
    if (lp.n_cols == 0)
        throw std::invalid_argument("Softmax needs at least one category");
    if (lp.n_rows == 0)
        return lp;
    const arma::colvec lpmax = arma::max(lp, 1);
    lp.each_col() -= lpmax;
    const arma::colvec denom = sum(exp(lp), 1);
    lp.each_col() -= arma::log(denom);
    return log ? lp : arma::exp(lp);
}

arma::vec softmax(arma::vec u) {
    if (u.is_empty() || !u.is_finite())
        throw std::invalid_argument("Softmax logits must be nonempty and finite");
    u -= u.max();
    arma::vec result = arma::exp(u);
    result /= arma::accu(result);
    return result;
}

arma::mat expit(arma::mat x) {
    if (x.has_nan())
        throw std::invalid_argument("Expit arguments must not be NaN");
    for (arma::uword i = 0; i < x.n_elem; ++i) {
        const double z = x(i);
        if (z >= 0)
            x(i) = 1.0 / (1.0 + std::exp(-z));
        else {
            const double e = std::exp(z);
            x(i) = e / (1.0 + e);
        }
    }
    return x;
}

arma::cx_mat expit(arma::cx_mat x) {
    if (!x.is_finite())
        throw std::invalid_argument("Complex expit arguments must be finite");
    for (arma::uword i = 0; i < x.n_elem; ++i) {
        const std::complex<double> z = x(i);
        if (z.real() >= 0)
            x(i) = 1.0 / (1.0 + std::exp(-z));
        else {
            const std::complex<double> e = std::exp(z);
            x(i) = e / (1.0 + e);
        }
    }
    if (!x.is_finite())
        throw std::domain_error("Complex expit is singular or overflowed at this argument");
    return x;
}
} // namespace target
