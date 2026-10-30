/*!
  @file statistics.cpp
  @author Klaus K. Holst
  @copyright 2020, Klaus Kähler Holst

  @brief Focused component extracted without changing its algorithms.
*/
#include <target/statistics.hpp>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace target {
double SupTest(const arma::vec &D) {
    if (D.is_empty() || !D.is_finite())
        throw std::invalid_argument("Process must be nonempty and finite");
    // Avoid a lazy Armadillo abs expression: its generic element accessor
    // passes an uninitialized, unused auxiliary scalar to the abs operation.
    // A scalar maximum is exact for finite inputs and needs no temporary.
    double maximum = 0;
    for (const double value : D)
        maximum = std::max(maximum, std::fabs(value));
    return maximum;
}

double L2Test(const arma::vec &D, const arma::vec &t) {
    if (D.is_empty() || D.n_elem != t.n_elem || !D.is_finite() || !t.is_finite() || !t.is_sorted())
        throw std::invalid_argument(
            "Process and sorted times must have equal nonzero lengths and finite values");
    double norm = 0.0;
    for (arma::uword i = 0; i + 1 < t.n_elem; ++i) {
        // sqrt(b-a) without overflowing an interval straddling zero, then
        // hypot instead of squaring tiny/large residuals prematurely.
        const double root_delta = t(i) < 0 && t(i + 1) > 0
                                      ? std::hypot(std::sqrt(-t(i)), std::sqrt(t(i + 1)))
                                      : std::sqrt(t(i + 1) - t(i));
        norm = std::hypot(norm, root_delta * D(i));
    }
    if (!std::isfinite(norm))
        throw std::domain_error("L2 statistic overflowed");
    return norm;
}

// Comparison of ECDF of (x1,...,xn) with null CDF
// evaluated in G = (G(x1),...,G(xn))
double CramerVonMises(const arma::vec &x, arma::vec G) {
    if (x.is_empty() || x.n_elem != G.n_elem || !x.is_finite() || !G.is_finite() ||
        arma::any(G < 0) || arma::any(G > 1))
        throw std::invalid_argument(
            "Supply nonempty matching observations and CDF values in [0,1]");
    // back to original order of input data:
    arma::uvec ord = arma::stable_sort_index(x);
    G = G.elem(ord);
    const double n = static_cast<double>(G.n_elem);
    double res = 1.0 / (12.0 * n);
    for (arma::uword i = 0; i < G.n_elem; i++) {
        double val = (2.0 * static_cast<double>(i) + 1.0) / (2.0 * n) - G(i);
        res += val * val;
    }
    return res;
}
} // namespace target
