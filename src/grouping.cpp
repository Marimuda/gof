/*!
  @file grouping.cpp
  @author Klaus K. Holst
  @copyright 2020, Klaus Kähler Holst

  @brief Focused component extracted without changing its algorithms.
*/
#include <target/grouping.hpp>
#include <stdexcept>

namespace target {
// Find index of clusters/ids in sorted integer vector
arma::umat clusterid(const arma::uvec &id) {
    if (id.is_empty())
        return arma::umat(0, 2);
    arma::uword ncl = 1;
    arma::uword n = id.size();
    arma::uword id0 = id(0);
    for (arma::uword i = 0; i < n; i++) {
        if (id(i) != id0) {
            ncl++;
            id0 = id(i);
        }
    }
    arma::umat res(ncl, 2); // column 1: index, column 2: size
    res.fill(0);
    arma::uword cl = 0;
    id0 = id(0);
    for (arma::uword i = 0; i < n; i++) {
        if (id(i) != id0) {
            cl++;
            res(cl, 0) = i; // index
            id0 = id(i);
        }
        res(cl, 1)++; // size
    }
    return res;
}

arma::mat groupsum(const arma::mat &x, const arma::uvec &cluster, bool reduce) {
    if (!x.is_finite())
        throw std::invalid_argument("Group values must be finite");
    if (x.n_rows == 0) {
        if (!cluster.is_empty())
            throw std::invalid_argument("Empty data cannot have cluster starts");
        return arma::mat(0, x.n_cols);
    }
    if (cluster.is_empty() || cluster(0) != 0 || cluster.max() >= x.n_rows)
        throw std::invalid_argument("Cluster starts must begin at zero and lie within the data");
    for (arma::uword i = 1; i < cluster.n_elem; ++i)
        if (cluster(i) <= cluster(i - 1))
            throw std::invalid_argument("Cluster starts must be strictly increasing");
    if (x.n_cols == 0)
        return arma::mat(reduce ? cluster.n_elem : x.n_rows, 0);
    arma::uword ncl = cluster.n_elem;
    arma::uword n = ncl;
    if (!reduce) {
        n = x.n_rows;
    }
    arma::uword stop;
    arma::mat res(n, x.n_cols);
    arma::rowvec tmp(x.n_cols);
    for (arma::uword i = 0; i < ncl; i++) { // Iterate over individuals
        arma::uword start = cluster(i);
        if (i == (ncl - 1)) {
            stop = x.n_rows;
        } else {
            stop = cluster(i + 1);
        }
        tmp.fill(0);
        for (arma::uword j = start; j < stop; j++) {
            tmp += x.row(j);
        }
        if (reduce) {
            res.row(i) = tmp;
        } else {
            for (arma::uword j = start; j < stop; j++) {
                res.row(j) = tmp;
            }
        }
    }
    if (!res.is_finite())
        throw std::domain_error("Group sum overflowed");
    return (res);
}

void fastpattern(const arma::umat &y, arma::umat &pattern, arma::uvec &group, unsigned categories) {
    if (categories == 0)
        throw std::invalid_argument("categories must be positive");
    arma::uword n = y.n_rows;
    arma::uword k = y.n_cols;
    if (n == 0) {
        pattern.set_size(0, k);
        group.set_size(0);
        return;
    }
    if (k == 0) {
        pattern.set_size(1, 0);
        group.zeros(n);
        return;
    }
    arma::uvec mygroup(n);
    // At most n distinct rows exist. The old categories^k allocation could
    // overflow or be too small for arbitrary labels and floating rounding.
    arma::umat mypattern(n, k);
    arma::uword K = 0;

    for (arma::uword i = 0; i < n; i++) {
        arma::urowvec ri = y.row(i);
        bool found = false;
        for (arma::uword j = 0; j < K; j++) {
            if (sum(ri != mypattern.row(j)) == 0) {
                found = true;
                mygroup(i) = j;
                break;
            }
        }
        if (!found) {
            K++;
            mypattern.row(K - 1) = ri;
            mygroup(i) = K - 1;
        }
    }
    pattern = mypattern.rows(0, K - 1);
    group = mygroup;
}
} // namespace target
