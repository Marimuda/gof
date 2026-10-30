// Internal shared aliases. Legacy headers retain their original public names.
#pragma once

#include <gof/detail/armadillo.hpp>
#include <complex>
#include <functional>
#include <vector>

namespace gof {
namespace detail {

using complex_scalar = std::complex<double>;
using complex_function = std::function<arma::cx_mat(arma::cx_vec)>;
using matrix_list = std::vector<arma::mat>;

} // namespace detail
} // namespace gof
