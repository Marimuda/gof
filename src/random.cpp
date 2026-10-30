// RNG boundary extracted from Klaus K. Holst's 2020 cumres implementation.
#include <gof/detail/random.hpp>

namespace gof {
namespace detail {

arma::vec normal_draw(arma::uword count) {
#ifdef ARMA_R
    Rcpp::RNGScope scope;
    return Rcpp::as<arma::vec>(Rcpp::rnorm(count));
#else
    return arma::randn<arma::vec>(count);
#endif
}

} // namespace detail
} // namespace gof
