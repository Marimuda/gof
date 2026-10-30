// Internal implementation detail, not a stable downstream API.
#pragma once

#include <gof/detail/armadillo.hpp>
#include <limits>
#include <stdexcept>

namespace gof {
namespace detail {

// Rcpp::rnorm takes an int even when the matrix backend has wider indices.
// Validate before casting, allocating or entering the R RNG scope. Kept
// separate so the boundary can be tested without a multi-gigabyte allocation.
inline int checked_r_normal_count(arma::uword count) {
    if (count > static_cast<arma::uword>(std::numeric_limits<int>::max()))
        throw std::invalid_argument("Normal draw count exceeds the Rcpp integer limit");
    return static_cast<int>(count);
}

// Uses the caller's existing R RNG stream in the R build, or Armadillo's
// existing stream in the standalone build. Never seeds or stores RNG state.
arma::vec normal_draw(arma::uword count);

} // namespace detail
} // namespace gof
