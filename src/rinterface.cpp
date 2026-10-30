/*!
  @file rinterface.cpp
  @author Klaus K. Holst
  @copyright 2020, Klaus Kähler Holst

  @brief R interface for the cumres class. 

  The relevant bindings are created in \c RcppExports.cpp, \c RcppExports.h
  after running \c Rcpp::compileAttributes()
*/

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::plugins(cpp17)]]

#include <target/cumres.hpp>
#include <RcppArmadillo.h>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace {
bool real_numeric(SEXP x) { return TYPEOF(x) == REALSXP || TYPEOF(x) == INTSXP; }

arma::umat checked_indices(SEXP x, arma::uword n, arma::uword p) {
    if (!Rf_isMatrix(x) || !real_numeric(x))
        throw std::invalid_argument("Evaluation indices must be a real numeric matrix");
    arma::mat value = Rcpp::as<arma::mat>(x);
    if (!value.is_empty()) {
        if (value.n_cols != p)
            throw std::invalid_argument(
                "Evaluation indices must be whole valid zero-based positions");
        // Check scalar values directly. Lazy floor expressions in Armadillo's
        // generic accessor pass an uninitialized, unused auxiliary argument.
        for (const double position : value)
            if (!std::isfinite(position) || position < 0 || position >= n ||
                position != std::floor(position))
                throw std::invalid_argument(
                    "Evaluation indices must be whole valid zero-based positions");
    }
    return arma::conv_to<arma::umat>::from(value);
}

arma::mat sample_checked(target::cumres *engine, SEXP idx) {
    return engine->sample(checked_indices(idx, engine->n, engine->ord.n_cols));
}

arma::mat statistics_checked(target::cumres *engine, SEXP count, SEXP idx, SEXP quantiles) {
    if (!real_numeric(count) || Rf_xlength(count) != 1)
        throw std::invalid_argument("R must be a positive whole number");
    double value = Rcpp::as<double>(count);
    if (!std::isfinite(value) || value < 1 || value != std::floor(value) ||
        value > std::numeric_limits<int>::max())
        throw std::invalid_argument("R must be a positive whole number within the R matrix limit");
    if (TYPEOF(quantiles) != LGLSXP || Rf_xlength(quantiles) != 1 ||
        LOGICAL(quantiles)[0] == NA_LOGICAL)
        throw std::invalid_argument("quantiles must be TRUE or FALSE");
    return engine->sample(static_cast<unsigned>(value),
                          checked_indices(idx, engine->n, engine->ord.n_cols),
                          LOGICAL(quantiles)[0] != 0);
}

void order_checked(target::cumres *engine, SEXP x) {
    if (!real_numeric(x))
        throw std::invalid_argument("Ordering variables must be real numeric values");
    arma::mat value = Rf_isMatrix(x) ? Rcpp::as<arma::mat>(x) : arma::mat(Rcpp::as<arma::vec>(x));
    engine->order(value);
}
} // namespace

// [[Rcpp::export]]
double SupTest(const arma::vec &x) { return target::SupTest(x); }

// [[Rcpp::export]]
double L2Test(const arma::vec &x, const arma::vec &t) { return target::L2Test(x, t); }

// [[Rcpp::export]]
arma::umat FastApprox(arma::vec &time, arma::vec &newtime, bool equal = false, unsigned type = 0) {
    arma::umat result = target::fastapprox(time, newtime, equal, type);
    result.col(0) += 1;
    return result;
}

// [[Rcpp::export]]
arma::umat ClusterIndex(const arma::uvec &t) {
    arma::umat res = target::clusterid(t);
    res.col(0) += 1;
    return res;
}

RCPP_MODULE(gofmod) {
    using namespace Rcpp;
    class_<target::cumres>("CumRes")
        // expose the constructor
        .constructor<arma::vec, arma::mat, arma::mat>("Constructor")
        .field("inp", &target::cumres::inp)
        .field("ord", &target::cumres::ord)
        .field("r", &target::cumres::r)
        .field("qt", &target::cumres::qt)

        .method("samplestat", &statistics_checked, "sample process and return Sup and L2 statistic")
        .method("sample1", &sample_checked, "sample process")
        .method("obs", &target::cumres::obs, "Return observed process")
        .method("rnorm", &target::cumres::rnorm, "Sample from Gaussian")
        .method("order", &order_checked, "Order observations after input variable");
}
