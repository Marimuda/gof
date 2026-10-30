/*!
  @file glm.hpp
  @author Klaus K. Holst
  @copyright 2018-2020, Klaus Kähler Holst

  @brief Utility functions for Generalized Linear Models

*/

#pragma once

// Preserve the historical standalone macro and global aliases for clients.
#ifndef ARMA_R
#define MATHLIB_STANDALONE
#endif
#include <gof/detail/armadillo.hpp>
#include <gof/detail/types.hpp>
#include <cmath>
#include <complex>
#include <cfloat>
#include <functional>
#include <vector>

#include <target/probability.hpp>
#include <target/influence.hpp>

using cx_dbl = gof::detail::complex_scalar;
using cx_func = gof::detail::complex_function;
using matlist = gof::detail::matrix_list;
