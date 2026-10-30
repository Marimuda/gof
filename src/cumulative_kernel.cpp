/*
 * Deterministic cumulative-process kernels.
 * Extracted from Klaus K. Holst's 2020 cumres implementation without changing
 * summation order, projection arithmetic or legacy multiplier coupling.
 */
#include <gof/detail/cumulative.hpp>
#include <cmath>
#include <stdexcept>

namespace gof {
namespace detail {

arma::mat observed_cumulative_process(const arma::vec& residuals,
                                      const arma::umat& ordering) {
    const arma::uword n = residuals.n_elem;
    arma::mat result(n, ordering.n_cols);
    for (arma::uword j = 0; j < ordering.n_cols; ++j) {
        const arma::vec prefix = arma::cumsum(residuals.elem(ordering.col(j))) /
            std::sqrt(static_cast<double>(n));
        result.col(j) = prefix;
    }
    return result;
}

arma::mat legacy_multiplier_process(cumulative_process_view data,
                                    const arma::vec& rank_multipliers,
                                    const arma::umat& indices,
                                    arma::mat& derivative_prefix) {
    const arma::uword n = data.residuals.n_elem;
    if (rank_multipliers.n_elem != n || !rank_multipliers.is_finite())
        throw std::invalid_argument("Supply one finite multiplier per sorted rank");
    const arma::uword rows = indices.is_empty() ? n : indices.n_rows;
    arma::mat result(rows, data.ordering.n_cols);
    for (arma::uword j = 0; j < data.ordering.n_cols; ++j) {
        const arma::uvec order = data.ordering.col(j);
        arma::vec residual_prefix = arma::cumsum(data.residuals.elem(order) % rank_multipliers);
        if (data.influences.n_cols == 0) {
            derivative_prefix.set_size(n, 0);
            if (!indices.is_empty()) residual_prefix = residual_prefix.elem(indices.col(j));
            result.col(j) = residual_prefix / std::sqrt(static_cast<double>(n));
            continue;
        }
        derivative_prefix = arma::cumsum(data.derivatives.rows(order), 0);
        arma::rowvec coefficient_draw(data.influences.n_cols);
        for (arma::uword k = 0; k < data.influences.n_cols; ++k) {
            const arma::vec contribution = data.influences.col(k);
            coefficient_draw(k) = arma::sum(contribution.elem(order) % rank_multipliers);
        }
        arma::vec correction(rows);
        if (!indices.is_empty()) {
            residual_prefix = residual_prefix.elem(indices.col(j));
            for (arma::uword i = 0; i < rows; ++i)
                correction(i) = arma::as_scalar(coefficient_draw * derivative_prefix.row(indices(i, j)).t());
        } else {
            for (arma::uword i = 0; i < n; ++i)
                correction(i) = arma::as_scalar(coefficient_draw * derivative_prefix.row(i).t());
        }
        result.col(j) = (residual_prefix + correction) / std::sqrt(static_cast<double>(n));
    }
    return result;
}

} // namespace detail
} // namespace gof
