// Klaus K. Holst's grouping interfaces, separated from utils.
#pragma once

#include <gof/detail/armadillo.hpp>

namespace target {

arma::umat clusterid(const arma::uvec &id);

// Cluster starts are zero-based, strictly increasing and begin at zero.
arma::mat groupsum(const arma::mat &x, const arma::uvec &cluster, bool reduce = true);

// Distinct rows in first-seen order, with zero-based group membership.
// categories must be positive; the labels themselves may be arbitrary.
void fastpattern(const arma::umat &y, arma::umat &pattern, arma::uvec &group,
                 unsigned categories = 2);

} // namespace target
