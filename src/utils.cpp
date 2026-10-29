
/*!
  @file utils.cpp
  @author Klaus K. Holst
  @copyright 2020, Klaus Kähler Holst

  @brief Various utility functions and constants

*/
#include <target/utils.hpp>
#include <algorithm>
#include <stdexcept>

namespace target {

  // Complex step derivative
  arma::mat deriv(cx_func f, arma::vec theta) {
    if (!f || !theta.is_finite())
      throw std::invalid_argument("Supply a function and finite parameters for complex-step differentiation");
    arma::cx_vec thetac = arma::conv_to<arma::cx_vec>::from(theta);
    arma::cx_mat val0 = f(thetac);
    if (!val0.is_finite() || arma::any(arma::vectorise(arma::imag(val0)) != 0))
      throw std::invalid_argument("Complex-step function must be finite and real at real parameters");
    arma::uword n = val0.n_elem;
    arma::uword p = theta.n_elem;
    arma::mat res(n, p);
    // DBL_MIN underflows intermediate products even for modest derivatives.
    const double h = 1e-20;
    cx_dbl h0 = cx_dbl(0, h);
    for (arma::uword i=0; i < p; i++) {
      arma::cx_vec theta0 = thetac;
      theta0[i] += h0;
      const arma::cx_mat perturbed = f(theta0);
      if (perturbed.n_rows != val0.n_rows || perturbed.n_cols != val0.n_cols || !perturbed.is_finite())
        throw std::invalid_argument("Complex-step function must retain its output shape and finite values");
      arma::mat val = imag(perturbed)/h;
      if (!val.is_finite()) throw std::domain_error("Complex-step derivative overflowed");
      for (arma::uword j=0; j < n; j++)
	res(j, i) = val[j];
    }
    return(res);
  }


  // Find index of clusters/ids in sorted integer vector
  arma::umat clusterid(const arma::uvec &id) {
    if (id.is_empty()) return arma::umat(0, 2);
    arma::uword ncl = 1;
    arma::uword n = id.size();
    arma::uword id0 = id(0);
    for (arma::uword i=0; i < n; i++) {
      if (id(i) != id0) {
	ncl++;
	id0 = id(i);
      }
    }
    arma::umat res(ncl, 2);  // column 1: index, column 2: size
    res.fill(0);
    arma::uword cl = 0;
    id0 = id(0);
    for (arma::uword i=0; i < n; i++) {
      if (id(i) != id0) {
	cl++;
	res(cl, 0) = i;  // index
	id0 = id(i);
      }
      res(cl, 1)++;  // size
    }
    return res;
  }


  arma::mat groupsum(const arma::mat &x,
		     const arma::uvec &cluster,
		     bool reduce) {
    if (!x.is_finite()) throw std::invalid_argument("Group values must be finite");
    if (x.n_rows == 0) {
      if (!cluster.is_empty()) throw std::invalid_argument("Empty data cannot have cluster starts");
      return arma::mat(0,x.n_cols);
    }
    if (cluster.is_empty() || cluster(0) != 0 || cluster.max() >= x.n_rows)
      throw std::invalid_argument("Cluster starts must begin at zero and lie within the data");
    for (arma::uword i=1; i < cluster.n_elem; ++i)
      if (cluster(i) <= cluster(i-1))
        throw std::invalid_argument("Cluster starts must be strictly increasing");
    if (x.n_cols == 0) return arma::mat(reduce ? cluster.n_elem : x.n_rows,0);
    arma::uword ncl = cluster.n_elem;
    arma::uword n = ncl;
    if (!reduce) {
      n = x.n_rows;
    }
    arma::uword stop;
    arma::mat res(n, x.n_cols);
    arma::rowvec tmp(x.n_cols);
    for (arma::uword i=0; i < ncl; i++) {  // Iterate over individuals
      arma::uword start = cluster(i);
      if (i == (ncl-1)) {
	stop = x.n_rows;
      } else {
	stop = cluster(i+1);
      }
      tmp.fill(0);
      for (arma::uword j=start; j < stop; j++) {
	tmp += x.row(j);
      }
      if (reduce) {
	res.row(i) = tmp;
      } else {
        for (arma::uword j=start; j < stop; j++) {
	  res.row(j) = tmp;
	}
      }
    }
    if (!res.is_finite()) throw std::domain_error("Group sum overflowed");
    return(res);
  }


  void fastpattern(const arma::umat &y,
		   arma::umat &pattern,
		   arma::uvec &group,
		   unsigned categories) {
    if (categories == 0) throw std::invalid_argument("categories must be positive");
    arma::uword n = y.n_rows;
    arma::uword k = y.n_cols;
    if (n == 0) {
      pattern.set_size(0,k);
      group.set_size(0);
      return;
    }
    if (k == 0) {
      pattern.set_size(1,0);
      group.zeros(n);
      return;
    }
    arma::uvec mygroup(n);
    // At most n distinct rows exist. The old categories^k allocation could
    // overflow or be too small for arbitrary labels and floating rounding.
    arma::umat mypattern(n, k);
    arma::uword K = 0;

    for (arma::uword i=0; i < n; i++) {
      arma::urowvec ri = y.row(i);
      bool found = false;
      for (arma::uword j=0; j < K; j++) {
	if (sum(ri != mypattern.row(j)) == 0) {
	  found = true;
	  mygroup(i) = j;
	  break;
	}
      }
      if (!found) {
	K++;
	mypattern.row(K-1) = ri;
	mygroup(i) = K-1;
      }
    }
    pattern = mypattern.rows(0, K-1);
    group = mygroup;
  }


  arma::umat fastapprox(arma::vec &time,  // sorted times
			const arma::vec &newtime,
			bool equal,
			// type: (0: nearedst, 1: right, 2: left)
			unsigned type) {
    if (time.is_empty() || !time.is_finite() || !newtime.is_finite())
      throw std::invalid_argument("Approximation times must be finite and time must be nonempty");
    if (!time.is_sorted() || type > 2)
      throw std::invalid_argument("time must be sorted and type must be 0, 1 or 2");
    arma::uvec idx(newtime.n_elem);
    arma::uvec eq(newtime.n_elem);

    for (arma::uword i=0; i < newtime.n_elem; i++) {
      auto it = std::lower_bound(time.begin(), time.end(), newtime(i));
      arma::uword pos;
      if (type == 2) {
        auto after = std::upper_bound(time.begin(), time.end(), newtime(i));
        pos = after == time.begin() ? 0 : static_cast<arma::uword>(after-time.begin()-1);
      } else if (it == time.end()) {
        pos = time.n_elem-1;
      } else {
        pos = static_cast<arma::uword>(it-time.begin());
        if (type == 0 && pos > 0 && std::fabs(newtime(i)-time(pos-1)) <
            std::fabs(newtime(i)-time(pos))) --pos;
      }
      idx(i) = pos;
      eq(i) = time(pos) == newtime(i) ? pos+1 : 0;
    }
    if (equal) {
      return( arma::join_horiz(idx, eq) );
    }
    return( idx );
  }


  double SupTest(const arma::vec &D) {
    if (D.is_empty() || !D.is_finite())
      throw std::invalid_argument("Process must be nonempty and finite");
    return arma::max(arma::abs(D));
  }

  double L2Test(const arma::vec &D, const arma::vec &t) {
    if (D.is_empty() || D.n_elem != t.n_elem || !D.is_finite() || !t.is_finite() || !t.is_sorted())
      throw std::invalid_argument("Process and sorted times must have equal nonzero lengths and finite values");
    double norm = 0.0;
    for (arma::uword i=0; i+1 < t.n_elem; ++i) {
      // sqrt(b-a) without overflowing an interval straddling zero, then
      // hypot instead of squaring tiny/large residuals prematurely.
      const double root_delta = t(i) < 0 && t(i+1) > 0
        ? std::hypot(std::sqrt(-t(i)),std::sqrt(t(i+1)))
        : std::sqrt(t(i+1)-t(i));
      norm = std::hypot(norm,root_delta*D(i));
    }
    if (!std::isfinite(norm)) throw std::domain_error("L2 statistic overflowed");
    return norm;
  }

  // Comparison of ECDF of (x1,...,xn) with null CDF
  // evaluated in G = (G(x1),...,G(xn))
  double CramerVonMises(const arma::vec &x, arma::vec G) {
    if (x.is_empty() || x.n_elem != G.n_elem || !x.is_finite() || !G.is_finite() ||
        arma::any(G < 0) || arma::any(G > 1))
      throw std::invalid_argument("Supply nonempty matching observations and CDF values in [0,1]");
    // back to original order of input data:
    arma::uvec ord = arma::stable_sort_index(x);
    G = G.elem(ord);
    const double n = static_cast<double>(G.n_elem);
    double res = 1.0/(12.0*n);
    for (arma::uword i=0; i < G.n_elem; i++) {
      double val = (2.0*static_cast<double>(i)+1.0)/(2.0*n)-G(i);
      res += val*val;
    }
    return res;
  }


  arma::mat const EmptyMat = arma::mat();
  arma::vec const EmptyVec = arma::vec();


  // Foreground colors are in form of 3x, bacground are 4x
  const char* COL_RESET = "\x1b[0m";
  const char* COL_DEF   = "\x1b[39m";
  const char* BLACK     = "\x1b[30m";
  const char* RED       = "\x1b[31m";
  const char* MAGENTA   = "\x1b[35m";
  const char* YELLOW    = "\x1b[33m";
  const char* GREEN     = "\x1b[32m";
  const char* BLUE      = "\x1b[34m";
  const char* CYAN      = "\x1b[36m";
  const char* WHITE     = "\x1b[37m";
  const char* GRAY      = "\x1b[90m";
  const char* LRED      = "\x1b[91m";
  const char* LGREEN    = "\x1b[92m";
  const char* LYELLOW   = "\x1b[93m";
  const char* LBLUE     = "\x1b[94m";
  const char* LMAGENTA  = "\x1b[95m";
  const char* LCYAN     = "\x1b[96m";
  const char* LWHITE    = "\x1b[97m";

}  // namespace target
