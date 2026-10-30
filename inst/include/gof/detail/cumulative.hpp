// Internal implementation detail, not a stable downstream API.
// Original cumulative-process implementation: Klaus K. Holst, 2020.
#pragma once

#include <gof/detail/armadillo.hpp>

namespace gof {
namespace detail {

// Borrowed inputs: callers retain ownership and must validate dimensions,
// finiteness and ordering permutations before calling the kernels. No input
// references are stored beyond the call. The legacy facade owns validation.
struct cumulative_process_view {
    const arma::vec &residuals;
    const arma::mat &derivatives;
    const arma::mat &influences;
    const arma::umat &ordering;
};

arma::mat observed_cumulative_process(const arma::vec &residuals, const arma::umat &ordering);

// One normal multiplier per SORTED RANK, deliberately preserving the legacy
// marginal coupling and arithmetic order. This is not the subject-coupled
// engine used by cumres_joint(). Indices are validated zero-based positions;
// a 0-by-0 matrix requests every position.
//
// derivative_prefix is explicit reusable scratch/output storage. It is
// overwritten on every ordering and retains the final ordering's prefix,
// matching the legacy eta field. Its old contents never affect the result.
// The draw is supplied by the caller; this function never accesses an RNG.
arma::mat legacy_multiplier_process(cumulative_process_view data, const arma::vec &rank_multipliers,
                                    const arma::umat &indices, arma::mat &derivative_prefix);

} // namespace detail
} // namespace gof
