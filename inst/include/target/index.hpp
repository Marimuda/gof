// Klaus K. Holst's index lookup interface, separated from utils.
#pragma once

#include <gof/detail/armadillo.hpp>

namespace target {

// type: 0 nearest (right on an exact distance tie), 1 right, 2 left.
// The nonconst time reference is retained for source/binary compatibility;
// the implementation does not modify it. The first column is zero-based;
// the optional equality column retains a 1-based match or zero for no match.
arma::umat fastapprox(arma::vec &time, const arma::vec &newtime, bool equal = false,
                      unsigned type = 0);

} // namespace target
