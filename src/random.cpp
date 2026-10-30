// RNG boundary extracted from Klaus K. Holst's 2020 cumres implementation.
#include <gof/detail/random.hpp>

namespace gof {
namespace detail {

arma::vec normal_draw(arma::uword count) {
#ifdef ARMA_R
    const int r_count = checked_r_normal_count(count);
    Rcpp::RNGScope scope;
    return Rcpp::as<arma::vec>(Rcpp::rnorm(r_count));
#else
    return arma::randn<arma::vec>(count);
#endif
}

} // namespace detail
} // namespace gof
