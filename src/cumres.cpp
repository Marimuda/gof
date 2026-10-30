/*!
  @file cumres.cpp
  @author Klaus K. Holst
  @copyright 2020, Klaus Kähler Holst

  @brief Generic function for calculating cumulative residuals

  Test statistics

*/
#include <target/cumres.hpp>
#include <gof/detail/cumulative.hpp>
#include <gof/detail/random.hpp>
#include <stdexcept>
#include <limits>

namespace target {

cumres::cumres(const arma::vec &r, const arma::mat &dr, const arma::mat &ic)
    : r(r), dr(dr), ic(ic) {
    if (r.n_elem > std::numeric_limits<unsigned>::max())
        throw std::invalid_argument("Too many observations for the cumulative-process index type");
    n = static_cast<unsigned>(r.n_elem);
    validate_inputs();
    arma::vec inp(n, 1);
    for (unsigned i = 0; i < n; i++)
        inp(i) = i;
    this->ord = arma::conv_to<arma::uvec>::from(inp);
    this->inp = inp;
    this->order(inp);
}

void cumres::validate_inputs() const {
    if (n == 0 || r.n_elem != n || dr.n_rows != n || ic.n_rows != n || dr.n_cols != ic.n_cols ||
        !r.is_finite() || !dr.is_finite() || !ic.is_finite())
        throw std::invalid_argument(
            "Incompatible or nonfinite residual, gradient or influence inputs");
}

void cumres::validate_order() const {
    validate_inputs();
    if (inp.n_rows != n || inp.n_cols == 0 || !inp.is_finite() || ord.n_rows != n ||
        ord.n_cols != inp.n_cols || ord.max() >= n)
        throw std::invalid_argument("Invalid cumulative process state");
    for (arma::uword j = 0; j < ord.n_cols; ++j) {
        for (arma::uword i = 1; i < n; ++i)
            if (inp(i, j) < inp(i - 1, j))
                throw std::invalid_argument("Process times must be sorted");
        arma::uvec seen(n, arma::fill::zeros);
        for (arma::uword i = 0; i < n; ++i)
            if (seen(ord(i, j))++)
                throw std::invalid_argument("Process order must be a permutation");
    }
}

void cumres::validate_indices(const arma::umat &idx) const {
    if (idx.n_rows == 0 && idx.n_cols == 0)
        return;
    if (idx.n_rows == 0 || idx.n_cols != ord.n_cols || idx.max() >= n)
        throw std::invalid_argument("Evaluation indices have invalid dimensions or positions");
}

void cumres::order(const arma::mat &inp) {
    validate_inputs();
    if (inp.n_rows != n || inp.n_cols == 0 || !inp.is_finite())
        throw std::invalid_argument("Ordering variables must have one finite row per observation");
    arma::uword p = inp.n_cols;
    arma::umat ord(inp.n_rows, p);
    arma::mat sorted(inp.n_rows, p);
    for (arma::uword j = 0; j < p; j++) {
        arma::vec inpj = inp.col(j);
        ord.col(j) = arma::stable_sort_index(inpj);
        sorted.col(j) = inpj.elem(ord.col(j));
    }
    this->inp = sorted;
    this->ord = ord;
    if (p == 1) {
        if (dr.n_cols == 0)
            eta.set_size(n, 0);
        else
            eta = arma::cumsum(dr.rows(ord), 0);
    }
}

arma::vec cumres::rnorm() {
    validate_inputs();
    return gof::detail::normal_draw(n);
}

arma::mat cumres::obs() {
    validate_order();
    return gof::detail::observed_cumulative_process(r, ord);
}

// Sample single process
arma::mat cumres::sample(const arma::umat &idx) {
    validate_order();
    validate_indices(idx);
    const arma::vec draw = rnorm();
    return gof::detail::legacy_multiplier_process({r, dr, ic, ord}, draw, idx, eta);
}

// Sample 'r' processes
arma::mat cumres::sample(unsigned R,            // Number of samples
                         const arma::umat &idx, // subset of 'time-points' of process to sample
                         bool quantiles) {
    if (R == 0 || R > static_cast<unsigned>(std::numeric_limits<int>::max()))
        throw std::invalid_argument("R must be positive and fit in an R integer");
    if (quantiles)
        throw std::invalid_argument("Quantile output is not supported; use quantiles=false");
    validate_order();
    validate_indices(idx);
    for (arma::uword j = 0; j < idx.n_cols; ++j)
        for (arma::uword i = 1; i < idx.n_rows; ++i)
            if (idx(i, j) < idx(i - 1, j))
                throw std::invalid_argument("Evaluation indices must be sorted for statistics");
    arma::uword p = ord.n_cols;
    arma::mat res(R, 2 * p);
    // qt.fill(0);
    // unsigned n = this->n;
    // arma::mat qt(n, std::ceil(R*0.05));
    for (unsigned i = 0; i < R; i++) {
        arma::mat wi = this->sample(idx);
        for (arma::uword j = 0; j < p; j++) {
            arma::vec t0 = inp.col(j);
            if (!idx.is_empty()) {
                t0 = t0.elem(idx.col(j));
            }
            res(i, j * 2) = SupTest(wi.col(j));
            res(i, j * 2 + 1) = L2Test(wi.col(j), t0);
        }
        /* if (quantiles) { // TODO: Disable for now. Capture quantiles
	 for (unsigned j=0; j<n; j++) {
	 wi = abs(wi);
	 arma::rowvec qtj = qt.row(j);
	 unsigned i = qtj.index_min();
	 if (wi(j)>qt(j,i)) qt(j,i) = wi(j);
	 }
	 } */
    }
    // this->qt = qt;
    return res;
}

} // namespace target
