// Separate consumer translation unit for the historical umbrella interfaces.
#ifdef GOF_REVERSE_HEADER_ORDER
#include <target/utils.hpp>
#include <target/glm.hpp>
#else
#include <target/glm.hpp>
#include <target/utils.hpp>
#endif
#include <target/cumres.hpp>
#include <type_traits>
#include <vector>

static_assert(std::is_same<cx_dbl, target::cx_dbl>::value, "complex scalar alias changed");
static_assert(std::is_same<cx_func, target::cx_func>::value, "complex function alias changed");
static_assert(std::is_same<matlist, std::vector<arma::mat>>::value, "matrix list alias changed");
static_assert(
    std::is_same<decltype(&target::deriv), arma::mat (*)(target::cx_func, arma::vec)>::value,
    "derivative signature changed");
static_assert(std::is_same<decltype(&target::fastapprox),
                           arma::umat (*)(arma::vec &, const arma::vec &, bool, unsigned)>::value,
              "index signature changed");

bool legacy_header_client() {
    auto softmax = static_cast<arma::vec (*)(arma::vec)>(&target::softmax);
    auto expit = static_cast<arma::mat (*)(arma::mat)>(&target::expit);
    const arma::vec p = softmax(arma::vec{0, 0});
    const arma::mat q = expit(arma::mat(1, 1, arma::fill::zeros));
    target::IID influence(arma::mat(2, 1, arma::fill::ones), arma::mat(1, 1, arma::fill::ones));
    target::cumres process(arma::vec{1, -1}, arma::mat(2, 0), arma::mat(2, 0));
    return p(0) == .5 && p(1) == .5 && q(0) == .5 && influence.iid(1) == 1 &&
           process.obs()(1) == 0 && target::EmptyMat.is_empty() && target::COL_RESET != nullptr;
}
