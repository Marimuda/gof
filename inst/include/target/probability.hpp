// Klaus K. Holst's probability transforms, separated from glm helpers.
#pragma once

#include <gof/detail/armadillo.hpp>

namespace target {

// Finite row-wise logits. The matrix overload optionally adds a reference
// category or returns log probabilities. Inputs are owned work copies.
arma::vec softmax(arma::vec logits);
arma::mat softmax(arma::mat logits, bool ref = true, bool log = false);

// Real expit includes the limits at +/-infinity. The complex overload supports
// finite analytic arguments, including complex-step perturbations.
arma::mat expit(arma::mat x);
arma::cx_mat expit(arma::cx_mat x);

} // namespace target
