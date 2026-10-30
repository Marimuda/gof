/*!
  @file index.cpp
  @author Klaus K. Holst
  @copyright 2020, Klaus Kähler Holst

  @brief Focused component extracted without changing its algorithms.
*/
#include <target/index.hpp>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace target {
arma::umat fastapprox(arma::vec &time, // sorted times
                      const arma::vec &newtime, bool equal,
                      // type: (0: nearedst, 1: right, 2: left)
                      unsigned type) {
    if (time.is_empty() || !time.is_finite() || !newtime.is_finite())
        throw std::invalid_argument("Approximation times must be finite and time must be nonempty");
    if (!time.is_sorted() || type > 2)
        throw std::invalid_argument("time must be sorted and type must be 0, 1 or 2");
    arma::uvec idx(newtime.n_elem);
    arma::uvec eq(newtime.n_elem);

    for (arma::uword i = 0; i < newtime.n_elem; i++) {
        auto it = std::lower_bound(time.begin(), time.end(), newtime(i));
        arma::uword pos;
        if (type == 2) {
            auto after = std::upper_bound(time.begin(), time.end(), newtime(i));
            pos = after == time.begin() ? 0 : static_cast<arma::uword>(after - time.begin() - 1);
        } else if (it == time.end()) {
            pos = time.n_elem - 1;
        } else {
            pos = static_cast<arma::uword>(it - time.begin());
            if (type == 0 && pos > 0 &&
                std::fabs(newtime(i) - time(pos - 1)) < std::fabs(newtime(i) - time(pos)))
                --pos;
        }
        idx(i) = pos;
        eq(i) = time(pos) == newtime(i) ? pos + 1 : 0;
    }
    if (equal) {
        return (arma::join_horiz(idx, eq));
    }
    return (idx);
}
} // namespace target
