#!/usr/bin/env Rscript
# Load gof from the intended library before sourcing this file, or configure
# R_LIBS when running it with Rscript. This exercises the INSTALLED C++ API.
if (!requireNamespace("gof",quietly=TRUE)) stop("Install gof in an isolated library first")
if (!requireNamespace("Rcpp",quietly=TRUE)) stop("Rcpp is required")
probe <- new.env(parent=globalenv())
requested <- Sys.getenv("GOF_CXX_STD","CXX17")
if (!nzchar(requested)) requested <- "CXX17"
plugins <- c(CXX14="cpp14",CXX17="cpp17",CXX20="cpp20",CXX23="cpp23")
if (!requested %in% names(plugins)) stop("Unsupported GOF_CXX_STD for the C++ client")
Rcpp::sourceCpp(env=probe,code=paste0("// [[Rcpp::plugins(",plugins[[requested]],")]]\n",'
// [[Rcpp::depends(gof,RcppArmadillo)]]
#include <gof.h>
#include <cmath>
#include <stdexcept>

void client_require(bool value) {
    if (!value) throw std::runtime_error("downstream C++ API regression");
}

// [[Rcpp::export]]
bool gof_cpp_client() {
    arma::vec x = {1,-2,3}, t = {1,2,4}, q = {.5,2.5,9};
    arma::uvec ids = {1,1,3};
    // Repeat calls to check both signature validation and cached callable paths.
    for (unsigned i=0;i<2;++i) {
        client_require(gof::SupTest(x)==3);
        client_require(std::fabs(gof::L2Test(x,t)-3)<1e-14);
        arma::umat index = gof::FastApprox(t,q);
        client_require(index.n_rows==3 && index.n_cols==1);
        client_require(index(0)==1 && index(1)==2 && index(2)==3);
        arma::umat clusters = gof::ClusterIndex(ids);
        client_require(clusters.n_rows==2 && clusters.n_cols==2);
        client_require(clusters(0,0)==1 && clusters(0,1)==2);
        client_require(clusters(1,0)==3 && clusters(1,1)==1);
    }
    bool rejected=false;
    try { gof::SupTest(arma::vec()); }
    catch (const std::exception&) { rejected=true; }
    client_require(rejected);
    rejected=false;
    try { gof::L2Test(x,arma::vec{1,2}); }
    catch (const std::exception&) { rejected=true; }
    client_require(rejected);
    arma::vec unsorted = {2,1};
    rejected=false;
    try { gof::FastApprox(unsorted,q); }
    catch (const std::exception&) { rejected=true; }
    client_require(rejected);
    client_require(gof::ClusterIndex(arma::uvec()).n_rows==0);
    return true;
}
// [[Rcpp::export]]
int gof_cpp_language_standard() { return __cplusplus; }
'))
stopifnot(isTRUE(probe$gof_cpp_client()))
actual <- probe$gof_cpp_language_standard()
minimum <- c(CXX14=201402L,CXX17=201703L,CXX20=202002L,CXX23=202100L)
stopifnot(actual>=minimum[[requested]])
cat("Requested client plugin:",plugins[[requested]],"; actual __cplusplus:",actual,"\n")
if (requested=="CXX14" && actual>201402L)
    cat("Host R selected a newer standard; this is not an actual C++14 client build.\n")
cat("Installed downstream C++ header/callable checks passed for gof ",
    as.character(utils::packageVersion("gof"))," at ",find.package("gof"),"\n",sep="")
