// Klaus K. Holst's process statistic interfaces, separated from utils.
#pragma once

#include <gof/detail/armadillo.hpp>

namespace target {

double SupTest(const arma::vec &process);
double L2Test(const arma::vec &process, const arma::vec &time);
double CramerVonMises(const arma::vec &x, arma::vec cdf);

} // namespace target
