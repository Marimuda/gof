/*!
  @file differentiation.cpp
  @author Klaus K. Holst
  @copyright 2020, Klaus Kähler Holst

  @brief Focused component extracted without changing its algorithms.
*/
#include <target/differentiation.hpp>
#include <stdexcept>

namespace target {
// Complex step derivative
arma::mat deriv(cx_func f, arma::vec theta) {
    if (!f || !theta.is_finite())
        throw std::invalid_argument(
            "Supply a function and finite parameters for complex-step differentiation");
    arma::cx_vec thetac = arma::conv_to<arma::cx_vec>::from(theta);
    arma::cx_mat val0 = f(thetac);
    if (!val0.is_finite() || arma::any(arma::vectorise(arma::imag(val0)) != 0))
        throw std::invalid_argument(
            "Complex-step function must be finite and real at real parameters");
    arma::uword n = val0.n_elem;
    arma::uword p = theta.n_elem;
    arma::mat res(n, p);
    // DBL_MIN underflows intermediate products even for modest derivatives.
    const double h = 1e-20;
    cx_dbl h0 = cx_dbl(0, h);
    for (arma::uword i = 0; i < p; i++) {
        arma::cx_vec theta0 = thetac;
        theta0[i] += h0;
        const arma::cx_mat perturbed = f(theta0);
        if (perturbed.n_rows != val0.n_rows || perturbed.n_cols != val0.n_cols ||
            !perturbed.is_finite())
            throw std::invalid_argument(
                "Complex-step function must retain its output shape and finite values");
        arma::mat val = imag(perturbed) / h;
        if (!val.is_finite())
            throw std::domain_error("Complex-step derivative overflowed");
        for (arma::uword j = 0; j < n; j++)
            res(j, i) = val[j];
    }
    return (res);
}
} // namespace target
