// Internal implementation detail, not a stable downstream API.
#pragma once

#include <gof/detail/armadillo.hpp>

namespace gof {
namespace detail {

// Uses the caller's existing R RNG stream in the R build, or Armadillo's
// existing stream in the standalone build. Never seeds or stores RNG state.
arma::vec normal_draw(arma::uword count);

} // namespace detail
} // namespace gof
