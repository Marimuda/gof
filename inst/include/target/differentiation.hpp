// Klaus K. Holst's complex-step differentiation interface, separated from utils.
#pragma once

#include <gof/detail/types.hpp>

namespace target {

using cx_dbl = gof::detail::complex_scalar;
using cx_func = gof::detail::complex_function;

// Real Jacobian of a holomorphic extension real-valued at theta. Output rows
// follow the column-major vectorization of the function's matrix result.
arma::mat deriv(cx_func f, arma::vec theta);

} // namespace target
