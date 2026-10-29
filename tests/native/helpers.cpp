// Standalone independent algebra/property tests for the entire native surface.
// Build with src/utils.cpp, src/glm.cpp, src/cumres.cpp and -larmadillo.
#include <target/utils.hpp>
#include <target/glm.hpp>
#include <target/cumres.hpp>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {
  unsigned cases=0, failures=0, assertions=0;
  const double inf=std::numeric_limits<double>::infinity();
  const double not_a_number=std::numeric_limits<double>::quiet_NaN();
  void require(bool condition, const char *message="native assertion failed") {
    ++assertions;
    if (!condition) throw std::runtime_error(message);
  }
  void near(double actual, double expected, double tolerance=1e-12) {
    require(std::isfinite(actual) && std::fabs(actual-expected)<=tolerance*(1+std::fabs(expected)),
            "numeric values differ");
  }
  template<class A,class B> void same(const A &a,const B &b,double tolerance=1e-12) {
    const auto left=a.eval();
    const auto right=b.eval();
    require(left.n_rows==right.n_rows && left.n_cols==right.n_cols,"matrix dimensions differ");
    for (arma::uword i=0;i<left.n_elem;++i) near(left(i),right(i),tolerance);
  }
  template<class E=std::invalid_argument,class F> void rejects(F function) {
    bool caught=false;
    try { function(); } catch (const E&) { caught=true; }
    require(caught,"expected explicit error was not raised");
  }
  void test(const char *name,const std::function<void()> &body) {
    ++cases;
    try {body();std::cout<<"PASS "<<name<<"\n";}
    catch(const std::exception &e) {++failures;std::cout<<"FAIL "<<name<<": "<<e.what()<<"\n";}
  }
  arma::mat manual_process(const target::cumres &engine,const arma::vec &g) {
    arma::mat expected(engine.n,engine.ord.n_cols,arma::fill::zeros);
    for(arma::uword j=0;j<engine.ord.n_cols;++j) {
      arma::vec B(engine.ic.n_cols,arma::fill::zeros),derivative(engine.dr.n_cols,arma::fill::zeros);
      for(arma::uword i=0;i<engine.n;++i)
        for(arma::uword k=0;k<B.n_elem;++k) B(k)+=engine.ic(engine.ord(i,j),k)*g(i);
      double residual=0;
      for(arma::uword i=0;i<engine.n;++i) {
        residual+=engine.r(engine.ord(i,j))*g(i);
        for(arma::uword k=0;k<derivative.n_elem;++k) derivative(k)+=engine.dr(engine.ord(i,j),k);
        double projection=0;
        for(arma::uword k=0;k<B.n_elem;++k) projection+=B(k)*derivative(k);
        expected(i,j)=(residual+projection)/std::sqrt(static_cast<double>(engine.n));
      }
    }
    return expected;
  }
  target::cumres fixture() {
    arma::vec r={1,-2,.5,3};
    arma::mat dr={{1,.5},{-.5,2},{2,-1},{.5,.25}};
    arma::mat ic={{.2,-.1},{-.3,.5},{.1,.2},{.4,-.2}};
    return target::cumres(r,dr,ic);
  }
}

int main() {
  std::cout.setf(std::ios::unitbuf);
  test("softmax vector probabilities and translation invariance",[] {
    const arma::vec v={-2,0,1};
    arma::vec expected(3);
    const double denominator=std::exp(-2)+1+std::exp(1);
    for(arma::uword i=0;i<3;++i) expected(i)=std::exp(v(i))/denominator;
    same(target::softmax(v),expected);
    same(target::softmax(arma::vec(v+1000)),expected);
    near(arma::accu(target::softmax(arma::vec{10000,-10000,10000})),1);
    same(target::softmax(arma::vec{3}),arma::vec{1});
  });
  test("softmax matrix reference and log probabilities",[] {
    arma::mat logits={{-2,0,1},{1000,1001,999}};
    const arma::mat probability=target::softmax(logits,false,false);
    for(arma::uword row=0;row<logits.n_rows;++row) {
      arma::vec v=logits.row(row).t();
      same(arma::vec(probability.row(row).t()),target::softmax(v));
      near(arma::accu(probability.row(row)),1);
    }
    same(arma::exp(target::softmax(logits,false,true)).eval(),probability);
    const arma::mat reference=target::softmax(arma::mat{{0,0}},true,false);
    same(reference,arma::mat{{1.0/3,1.0/3,1.0/3}});
    near(target::softmax(arma::mat(1,0),true,false)(0,0),1);
    const arma::mat empty=target::softmax(arma::mat(0,2),true,false);
    require(empty.n_rows==0 && empty.n_cols==3);
  });
  test("softmax invalid inputs",[] {
    rejects([&]{target::softmax(arma::vec());});
    rejects([&]{target::softmax(arma::vec{inf,0});});
    rejects([&]{target::softmax(arma::mat(1,1,arma::fill::value(not_a_number)),false,false);});
    rejects([&]{target::softmax(arma::mat(2,0),false,false);});
  });
  test("real expit stable extremes and limits",[] {
    const arma::mat p=target::expit(arma::mat{{-inf,-1000,-1,0,1,1000,inf}});
    near(p(0),0);near(p(1),0);near(p(2),1/(1+std::exp(1.0)));
    near(p(3),.5);near(p(4),1/(1+std::exp(-1.0)));near(p(5),1);near(p(6),1);
    require(target::expit(arma::mat(0,2)).n_rows==0);
    rejects([&]{target::expit(arma::mat(1,1,arma::fill::value(not_a_number)));});
  });
  test("complex expit stable analytic extension",[] {
    arma::cx_mat z(1,4);
    z(0)={-1000,1e-20};z(1)={1000,1e-20};z(2)={0,1e-20};z(3)={1,.2};
    const arma::cx_mat p=target::expit(z);
    require(p.is_finite());
    near(p(0).real(),0);near(p(1).real(),1);
    near(p(2).imag()/1e-20,.25);
    const auto reference=1.0/(1.0+std::exp(-z(3)));
    near(p(3).real(),reference.real());near(p(3).imag(),reference.imag());
    rejects([&]{target::expit(arma::cx_mat(1,1,arma::fill::value(std::complex<double>(not_a_number,0))));});
  });
  test("IID constructor multiplication and validation",[] {
    target::IID identity;
    require(identity.iid.n_elem==1 && identity.vcov.n_elem==1);
    target::IID fit(arma::mat{{1,2},{3,4}},arma::mat{{2,0},{0,3}});
    same(fit.iid,arma::mat{{2,6},{6,12}});
    rejects([&]{target::IID bad(arma::mat(3,2),arma::mat(3,3));});
    rejects([&]{target::IID bad(arma::mat(3,2),arma::mat(2,3));});
    rejects([&]{target::IID bad(arma::mat(1,1,arma::fill::value(not_a_number)),arma::mat(1,1,arma::fill::ones));});
  });
  test("logistic influence matches explicit weighted 2x2 algebra",[] {
    arma::vec y={0,1,0,1,1},p={.2,.4,.5,.6,.8},w={1,2,0,.5,3};
    arma::mat x={{1,-2},{1,-1},{1,0},{1,1},{1,2}};
    double h00=0,h01=0,h11=0;
    for(arma::uword i=0;i<5;++i) {
      const double v=p(i)*(1-p(i))*w(i);
      h00+=v;h01+=v*x(i,1);h11+=v*x(i,1)*x(i,1);
    }
    const double det=h00*h11-h01*h01;
    arma::mat inverse={{h11/det,-h01/det},{-h01/det,h00/det}};
    const auto fit=target::logistic_iid(y,p,x,w);
    same(fit.vcov,inverse);
    for(arma::uword i=0;i<5;++i) {
      const double score=(y(i)-p(i))*w(i);
      near(fit.iid(i,0),score*(inverse(0,0)+x(i,1)*inverse(1,0)));
      near(fit.iid(i,1),score*(inverse(0,1)+x(i,1)*inverse(1,1)));
    }
  });
  test("linear influence matches explicit precision-weight algebra and zero weights",[] {
    arma::vec y={-1,0,3,1,4},p={-1.5,-.5,.5,1.5,2.5},w={1,2,0,.5,3};
    arma::mat x={{1,-2},{1,-1},{1,0},{1,1},{1,2}};
    double h00=0,h01=0,h11=0,sse=0;
    for(arma::uword i=0;i<5;++i) {
      h00+=w(i);h01+=w(i)*x(i,1);h11+=w(i)*x(i,1)*x(i,1);
      sse+=w(i)*(y(i)-p(i))*(y(i)-p(i));
    }
    const double det=h00*h11-h01*h01,sigma2=sse/(4-2);
    arma::mat inverse={{h11/det,-h01/det},{-h01/det,h00/det}};
    const auto fit=target::linear_iid(y,p,x,w);
    same(fit.vcov,arma::mat(sigma2*inverse));
    for(arma::uword i=0;i<5;++i) {
      const double score=(y(i)-p(i))*w(i);
      near(fit.iid(i,0),score*(inverse(0,0)+x(i,1)*inverse(1,0)));
      near(fit.iid(i,1),score*(inverse(0,1)+x(i,1)*inverse(1,1)));
    }
    same(fit.iid.row(2).eval(),arma::rowvec{0,0});
    const auto exact=target::linear_iid(p,p,x,w);
    require(exact.iid.is_finite() && arma::accu(arma::abs(exact.iid))==0);
    require(exact.vcov.is_finite() && arma::accu(arma::abs(exact.vcov))==0);
  });
  test("influence input contracts and singular information",[] {
    const arma::vec y={0,1,0},p={.25,.5,.75},w={1,1,1};
    arma::mat x={{1,-1},{1,0},{1,1}};
    rejects([&]{target::logistic_iid(y,p.head(2),x,w);});
    rejects([&]{target::linear_iid(y,p,x,arma::vec{1,-1,1});});
    rejects([&]{target::logistic_iid(y,p,x,arma::zeros(3));});
    rejects([&]{target::logistic_iid(arma::vec{0,2,1},p,x,w);});
    rejects([&]{target::logistic_iid(y,arma::vec{0,.5,1.1},x,w);});
    rejects([&]{target::linear_iid(y,p,x,arma::vec{1,0,1});});
    rejects([&]{target::linear_iid(arma::vec{1},arma::vec{0},arma::mat{{1,2}},arma::vec{1});});
    arma::mat singular(3,2,arma::fill::ones);
    rejects<std::domain_error>([&]{target::logistic_iid(y,p,singular,w);});
    rejects<std::domain_error>([&]{target::linear_iid(y,p,singular,w);});
    rejects([&]{target::linear_iid(arma::vec{0,not_a_number,1},p,x,w);});
  });
  test("complex-step Jacobian analytic matrix oracle",[] {
    target::cx_func function=[](arma::cx_vec z) {
      arma::cx_mat result(2,2);
      result(0,0)=z(0)*z(0);result(1,0)=std::sin(z(0));
      result(0,1)=std::exp(z(1));result(1,1)=z(0)*z(1);
      return result;
    };
    same(target::deriv(function,arma::vec{2,3}),arma::mat{{4,0},{std::cos(2.0),0},{0,std::exp(3.0)},{3,2}});
    const arma::mat tiny=target::deriv([](arma::cx_vec z) {
      return arma::cx_mat(1,1,arma::fill::value(z(0)*1e-100));
    },arma::vec{2});
    near(tiny(0)/1e-100,1);
  });
  test("complex-step empty dimensions and invalid analytic contracts",[] {
    target::cx_func constant=[](arma::cx_vec) {return arma::cx_mat(2,1,arma::fill::ones);};
    const auto empty=target::deriv(constant,arma::vec());
    require(empty.n_rows==2 && empty.n_cols==0);
    rejects([&]{target::deriv(target::cx_func(),arma::vec{1});});
    rejects([&]{target::deriv(constant,arma::vec{not_a_number});});
    rejects([&]{target::deriv([](arma::cx_vec) {
      return arma::cx_mat(1,1,arma::fill::value(std::complex<double>(1,1)));
    },arma::vec{1});});
    rejects([&]{target::deriv([](arma::cx_vec z) {
      return z(0).imag()==0 ? arma::cx_mat(1,2,arma::fill::ones) : arma::cx_mat(2,1,arma::fill::ones);
    },arma::vec{1});});
  });
  test("cluster starts and sizes reconstruct all contiguous runs",[] {
    const arma::uvec id={5,5,2,2,2,9};
    const arma::umat clusters=target::clusterid(id);
    same(clusters,arma::umat{{0,2},{2,3},{5,1}});
    const auto empty=target::clusterid(arma::uvec());
    require(empty.n_rows==0 && empty.n_cols==2);
    const arma::uword large=std::numeric_limits<arma::uword>::max();
    same(target::clusterid(arma::uvec{large,large,0}),arma::umat{{0,2},{2,1}});
  });
  test("groupsum reduced and expanded algebra",[] {
    const arma::mat x={{1,10},{2,20},{3,30},{4,40},{5,50}};
    const arma::uvec start={0,2,3};
    same(target::groupsum(x,start,true),arma::mat{{3,30},{3,30},{9,90}});
    same(target::groupsum(x,start,false),arma::mat{{3,30},{3,30},{3,30},{9,90},{9,90}});
    require(target::groupsum(arma::mat(0,2),arma::uvec(),false).n_cols==2);
    require(target::groupsum(arma::mat(3,0),arma::uvec{0,2},true).n_rows==2);
  });
  test("groupsum rejects incomplete or overlapping partitions",[] {
    const arma::mat x={{1,10},{2,20},{3,30},{4,40}};
    for(const arma::uvec &starts: {arma::uvec(),arma::uvec{1},arma::uvec{0,0},arma::uvec{0,3,1},arma::uvec{0,4}})
      rejects([&]{target::groupsum(x,starts,false);});
    rejects([&]{target::groupsum(arma::mat(0,2),arma::uvec{0},true);});
    rejects([&]{target::groupsum(arma::mat(1,1,arma::fill::value(not_a_number)),arma::uvec{0},true);});
  });
  test("fastpattern first-seen rows reconstruct arbitrary category labels",[] {
    const arma::umat y={{10,0},{20,1},{10,0},{30,2},{20,1}};
    arma::umat pattern;arma::uvec group;
    target::fastpattern(y,pattern,group,1);
    same(pattern,arma::umat{{10,0},{20,1},{30,2}});
    same(group,arma::uvec{0,1,0,2,1});
    for(arma::uword i=0;i<y.n_rows;++i) same(pattern.row(group(i)).eval(),y.row(i).eval());
  });
  test("fastpattern empty shapes and reserved category hint",[] {
    arma::umat pattern;arma::uvec group;
    target::fastpattern(arma::umat(0,3),pattern,group,2);
    require(pattern.n_rows==0 && pattern.n_cols==3 && group.n_elem==0);
    target::fastpattern(arma::umat(4,0),pattern,group,2);
    require(pattern.n_rows==1 && pattern.n_cols==0 && group.n_elem==4 && arma::accu(group)==0);
    rejects([&]{target::fastpattern(arma::umat(1,1,arma::fill::zeros),pattern,group,0);});
  });
  test("fastapprox nearest left right boundaries and equality flags",[] {
    arma::vec times={1,2,2,4},query={0,2,2.5,9};
    same(target::fastapprox(times,query,false,0),arma::uvec{0,1,2,3});
    same(target::fastapprox(times,query,false,1),arma::uvec{0,1,3,3});
    same(target::fastapprox(times,query,false,2),arma::uvec{0,2,2,3});
    same(target::fastapprox(times,query,true,1),arma::umat{{0,0},{1,2},{3,0},{3,0}});
    arma::vec extreme={-1e308,1e308},center={0};
    require(target::fastapprox(extreme,center,false,0)(0)==1);
    const auto empty=target::fastapprox(times,arma::vec(),true,0);
    require(empty.n_rows==0 && empty.n_cols==2);
  });
  test("fastapprox exhaustive nearest-neighbour property",[] {
    arma::vec times={-5,-2,-1,0,3,7,9};
    arma::vec query=arma::linspace(-8,12,201);
    const auto found=target::fastapprox(times,query,false,0);
    for(arma::uword i=0;i<query.n_elem;++i) {
      arma::uword best=0;
      for(arma::uword j=1;j<times.n_elem;++j)
        if(std::fabs(query(i)-times(j))<=std::fabs(query(i)-times(best))) best=j;
      require(found(i)==best);
    }
    rejects([&]{target::fastapprox(times,arma::vec{not_a_number},false,0);});
    arma::vec empty,unsorted={2,1};
    rejects([&]{target::fastapprox(empty,query,false,0);});
    rejects([&]{target::fastapprox(unsorted,query,false,0);});
    rejects([&]{target::fastapprox(times,query,false,3);});
  });
  test("Sup and L2 statistics match direct weighted norms",[] {
    const arma::vec d={-2,3,1},time={0,2,5};
    near(target::SupTest(d),3);
    near(target::L2Test(d,time),std::sqrt(2*4+3*9));
    near(target::L2Test(arma::vec{3},arma::vec{7}),0);
    near(target::L2Test(arma::vec{1e200,0},arma::vec{0,1})/1e200,1);
    near(target::L2Test(arma::vec{1e-100,0},arma::vec{0,1e-300})/1e-250,1);
    near(target::L2Test(arma::vec{1e-154,0},arma::vec{-1e308,1e308}),std::sqrt(2.0));
  });
  test("Sup and L2 input contracts",[] {
    rejects([&]{target::SupTest(arma::vec());});
    rejects([&]{target::SupTest(arma::vec{not_a_number});});
    rejects([&]{target::L2Test(arma::vec{1},arma::vec{0,1});});
    rejects([&]{target::L2Test(arma::vec{1,2},arma::vec{1,0});});
    rejects([&]{target::L2Test(arma::vec{1,inf},arma::vec{0,1});});
  });
  test("CramerVonMises rank formula ordering and singleton",[] {
    const arma::vec x={3,1,2},cdf={.8,.2,.4},ordered={.2,.4,.8};
    double expected=1.0/36;
    for(arma::uword i=0;i<3;++i) expected+=std::pow((2.0*i+1)/6-ordered(i),2);
    near(target::CramerVonMises(x,cdf),expected);
    near(target::CramerVonMises(arma::vec{2},arma::vec{.5}),1.0/12);
    rejects([&]{target::CramerVonMises(arma::vec(),arma::vec());});
    rejects([&]{target::CramerVonMises(arma::vec{1},arma::vec{1.1});});
    rejects([&]{target::CramerVonMises(arma::vec{1,2},arma::vec{.5});});
  });
  test("cumres observed processes and alias-safe ordering",[] {
    auto engine=fixture();
    engine.order(arma::mat{{3,0},{1,1},{2,1},{0,2}});
    same(engine.obs(),arma::mat{{3,1},{1,-1},{1.5,-.5},{2.5,2.5}}/2);
    auto single=fixture();
    single.order(arma::vec{3,1,2,0});
    const arma::mat sorted=single.inp;
    single.order(single.inp);
    same(single.inp,sorted);
  });
  test("cumres full and selected simulations match an independent multiplier oracle",[] {
    auto engine=fixture();
    engine.order(arma::mat{{3,0},{1,1},{2,1},{0,2}});
    arma::arma_rng::set_seed(81);const arma::vec g=arma::randn<arma::vec>(engine.n);
    const arma::mat expected=manual_process(engine,g);
    arma::arma_rng::set_seed(81);same(engine.sample(),expected);
    const arma::umat idx={{0,1},{2,3}};
    arma::mat selected(2,2);
    for(arma::uword j=0;j<2;++j) for(arma::uword i=0;i<2;++i) selected(i,j)=expected(idx(i,j),j);
    arma::arma_rng::set_seed(81);same(engine.sample(idx),selected);
    arma::arma_rng::set_seed(82);const arma::vec normal=arma::randn<arma::vec>(engine.n);
    arma::arma_rng::set_seed(82);same(engine.rnorm(),normal);
  });
  test("cumres batch statistics consume the expected random stream",[] {
    auto engine=fixture();engine.order(arma::vec{3,1,2,0});
    const arma::umat idx=arma::uvec{0,2,3};
    arma::mat expected(4,2);
    arma::arma_rng::set_seed(83);
    for(arma::uword draw=0;draw<4;++draw) {
      const arma::vec g=arma::randn<arma::vec>(engine.n);
      const arma::mat full=manual_process(engine,g);
      arma::vec values(3),time(3);
      for(arma::uword i=0;i<3;++i) {values(i)=full(idx(i),0);time(i)=engine.inp(idx(i),0);}
      expected(draw,0)=arma::max(arma::abs(values));
      double integral=0;
      for(arma::uword i=0;i+1<3;++i) integral+=(time(i+1)-time(i))*values(i)*values(i);
      expected(draw,1)=std::sqrt(integral);
    }
    arma::arma_rng::set_seed(83);same(engine.sample(4,idx,false),expected);
  });
  test("cumres construction preserves the caller's seed and shared RNG stream",[] {
    arma::arma_rng::set_seed(42);
    const arma::vec first=arma::randn<arma::vec>(4);
    const arma::vec second=arma::randn<arma::vec>(4);
    arma::arma_rng::set_seed(42);
    auto engine=fixture();
    auto another=fixture();
    same(engine.rnorm(),first,0);
    same(another.rnorm(),second,0);
    arma::arma_rng::set_seed(42);
    auto repeated=fixture();
    same(repeated.rnorm(),first,0);
    arma::arma_rng::set_seed(97);
    const arma::vec expected_next=arma::randn<arma::vec>(4);
    arma::arma_rng::set_seed(97);
    auto unused=fixture();
    (void)unused;
    same(arma::randn<arma::vec>(4),expected_next,0);
  });
  test("cumres derivative mutation invalidates the old prefix cache",[] {
    auto engine=fixture();engine.order(arma::vec{3,1,2,0});
    engine.dr*=2;
    engine.eta.fill(not_a_number); // derived cache is recomputed, not treated as input
    arma::arma_rng::set_seed(84);const arma::vec g=arma::randn<arma::vec>(engine.n);
    const arma::mat expected=manual_process(engine,g);
    arma::arma_rng::set_seed(84);same(engine.sample(),expected);
  });
  test("cumres supports zero-coefficient models without a projection",[] {
    target::cumres engine(arma::vec{1,-1,2},arma::mat(3,0),arma::mat(3,0));
    engine.order(arma::vec{3,1,2});
    same(engine.obs(),arma::vec{-1,1,2}/std::sqrt(3.0));
    arma::arma_rng::set_seed(85);const arma::vec g=arma::randn<arma::vec>(engine.n);
    const arma::mat expected=manual_process(engine,g);
    arma::arma_rng::set_seed(85);same(engine.sample(),expected);
  });
  test("cumres guards mutated fields before access or allocation",[] {
    auto base=fixture();
    auto engine=base;engine.r(0)=not_a_number;rejects([&]{engine.obs();});
    engine=base;engine.dr(0,0)=not_a_number;rejects([&]{engine.sample();});
    engine=base;engine.ic(0,0)=inf;rejects([&]{engine.sample();});
    engine=base;engine.ord(1,0)=engine.ord(0,0);rejects([&]{engine.obs();});
    engine=base;engine.ord(0,0)=engine.n;rejects([&]{engine.obs();});
    engine=base;engine.inp(0,0)=inf;rejects([&]{engine.sample();});
    engine=base;engine.inp(0,0)=100;rejects([&]{engine.obs();});
    engine=base;engine.n=1;rejects([&]{engine.rnorm();});
    engine=base;engine.dr.set_size(4,0);rejects([&]{engine.sample();});
    rejects([&]{target::cumres bad(arma::vec(),arma::mat(0,1),arma::mat(0,1));});
  });
  test("cumres index and simulation-count contracts",[] {
    auto engine=fixture();
    rejects([&]{engine.sample(arma::umat(1,2,arma::fill::zeros));});
    rejects([&]{engine.sample(arma::umat(0,1));});
    rejects([&]{engine.sample(arma::umat(2,0));});
    rejects([&]{engine.sample(arma::umat(1,1,arma::fill::value(engine.n)));});
    rejects([&]{engine.sample(0,arma::umat(),false);});
    rejects([&]{engine.sample(std::numeric_limits<unsigned>::max(),arma::umat(),false);});
    rejects([&]{engine.sample(1,arma::umat(),true);});
    rejects([&]{engine.sample(1,arma::umat(arma::uvec{2,0}),false);});
  });
  test("empty constants and terminal colour constants",[] {
    require(target::EmptyMat.is_empty() && target::EmptyVec.is_empty());
    const char *colours[]={target::COL_RESET,target::COL_DEF,target::BLACK,target::RED,target::GREEN,
      target::YELLOW,target::BLUE,target::MAGENTA,target::CYAN,target::WHITE,target::GRAY,
      target::LRED,target::LGREEN,target::LYELLOW,target::LBLUE,target::LMAGENTA,target::LCYAN,target::LWHITE};
    for(const char *colour:colours) require(colour && std::string(colour).find("\x1b[")==0);
  });
  std::cout<<cases<<" native cases, "<<assertions<<" assertions, "<<failures<<" failures\n";
  return failures ? 1 : 0;
}
