// Klaus K. Holst's influence interfaces, separated from probability transforms.
#pragma once

#include <gof/detail/armadillo.hpp>

namespace target {

class IID {
  public:
    arma::mat iid;
    arma::mat vcov;
    IID() : iid(arma::zeros(1, 1)), vcov(arma::zeros(1, 1)) {}
    IID(arma::mat score, arma::mat v);
};

IID logistic_iid(const arma::vec &y, const arma::vec &p, const arma::mat &x, const arma::vec &w);
IID linear_iid(const arma::vec &y, const arma::vec &p, const arma::mat &x, const arma::vec &w);

} // namespace target
