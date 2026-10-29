// Standalone regression tests for native helpers not exposed by the R API.
// Compile this file with src/utils.cpp and src/cumres.cpp, -I inst/include,
// -std=c++20 and -larmadillo, then run the resulting executable.
#include <target/utils.hpp>
#include <target/cumres.hpp>
#include <cmath>
#include <iostream>
#include <stdexcept>

void require(bool ok) {
    if (!ok) throw std::runtime_error("native regression failure");
}

int main() {
    arma::vec x = {3,1,2};
    arma::vec G = {.8,.2,.4};
    double expected = 1.0/36.0;
    arma::vec sorted = {.2,.4,.8};
    for (arma::uword i=0;i<3;++i) {
        const double error = (2.0*i+1.0)/6.0-sorted(i);
        expected += error*error;
    }
    require(std::fabs(target::CramerVonMises(x,G)-expected)<1e-14);
    require(std::fabs(target::CramerVonMises(arma::vec{2},arma::vec{.5})-1.0/12.0)<1e-14);
    bool rejected = false;
    try { target::CramerVonMises(arma::vec(),arma::vec()); }
    catch (const std::invalid_argument&) { rejected=true; }
    require(rejected);
    arma::vec t = {1,2,2,4}, queries = {0,2,9};
    arma::umat left = target::fastapprox(t,queries,false,2);
    require(left(0)==0 && left(1)==2 && left(2)==3);
    require(target::clusterid(arma::uvec()).n_rows==0);
    arma::vec r = {1,-1,2};
    arma::mat dr(3,1,arma::fill::ones), ic(3,1,arma::fill::zeros);
    target::cumres engine(r,dr,ic);
    engine.order(arma::mat(t.head(3)));
    require(engine.obs().n_rows==3);
    rejected=false;
    try { engine.sample(arma::umat(1,2,arma::fill::zeros)); }
    catch (const std::invalid_argument&) { rejected=true; }
    require(rejected);
    std::cout << "Native helper regressions passed\n";
}
