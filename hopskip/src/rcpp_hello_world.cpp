#include "boost/any.hpp"
#include "Rcpp.h"
#include "simple_hazard.hpp"
#include "bugs.hpp"
#include "excess_growth.hpp"
#include "rng.hpp"
#include "smv.hpp"


using namespace Rcpp;
using namespace hsb;


namespace hsb {
class SIRObserver : public hsb::TrajectoryObserver {
 public:
  virtual void event(double when, int what, int64_t who, int64_t who2) {
    when_.push_back(when);
    what_.push_back(what);
    who_.push_back(who);
    who2_.push_back(who2);
  }

  std::vector<double> when_;
  std::vector<int> what_;
  std::vector<int64_t> who_;
  std::vector<int64_t> who2_;
};
}



// [[Rcpp::export]]
DataFrame simple_hazard(SEXP pairwise_distanceS, SEXP parametersS) {
  afidd::LogInit("error");

  NumericVector pairwise(pairwise_distanceS);
  std::vector<double> distance(pairwise.begin(), pairwise.end());

  List parameters(parametersS);
  std::map<std::string, boost::any> params;
  params["individual_cnt"]=int64_t{as<int>(parameters["individual_cnt"])};
  params["beta0"]=double{as<double>(parameters["beta0"])};
  params["beta1"]=double{as<double>(parameters["beta1"])};
  params["beta2"]=double{as<double>(parameters["beta2"])};
  params["gamma"]=double{as<double>(parameters["gamma"])};
  params["cutoff"]=double{as<double>(parameters["cutoff"])};
  int64_t rand_seed=int64_t{as<int>(parameters["seed"])};
  params["seed"]=rand_seed;

  auto observer=std::make_shared<SIRObserver>();

  RandGen rng(rand_seed);
  hsb::simple_hazard::SIR_run(params, distance, observer, rng);

  NumericVector when(observer->when_.begin(), observer->when_.end());
  IntegerVector what(observer->what_.begin(), observer->what_.end());
  IntegerVector who(observer->who_.begin(), observer->who_.end());
  IntegerVector who2(observer->who2_.begin(), observer->who2_.end());

  return DataFrame::create(Named("times")=when,
      Named("event")=what, Named("who")=who, Named("actor")=who2);
}



// [[Rcpp::export]]
DataFrame bugs(SEXP pairwise_distanceS, SEXP parametersS) {
  afidd::LogInit("error");

  NumericVector pairwise(pairwise_distanceS);
  std::vector<double> distance(pairwise.begin(), pairwise.end());

  List parameters(parametersS);
  std::map<std::string, boost::any> params;
  for (auto iname : std::vector<std::string>{"individual_cnt", "seed",
      "initial_bug_cnt"}) {
    try {
      params[iname]=int64_t{as<int>(parameters[iname])};
    } catch (std::exception& e) {
      BOOST_LOG_TRIVIAL(error) << "Could not cast "<<iname<<
        " to an integer";
      throw;
    }
  }
  for (auto dname : std::vector<std::string>{"birth", "death", "carrying",
      "move0", "move1", "gamma", "cutoff"}) {
    try {
      params[dname]=double{as<double>(parameters[dname])};
    } catch (std::exception& e) {
      BOOST_LOG_TRIVIAL(error) << "Could not cast "<<dname<<" to a double";
      throw;
    }
  }

  auto observer=std::make_shared<SIRObserver>();

  int64_t rand_seed=boost::any_cast<int64_t>(params["seed"]);
  RandGen rng(rand_seed);
  hsb::bugs::SIR_run(params, distance, observer, rng);

  NumericVector when(observer->when_.begin(), observer->when_.end());
  IntegerVector what(observer->what_.begin(), observer->what_.end());
  IntegerVector who(observer->who_.begin(), observer->who_.end());
  IntegerVector who2(observer->who2_.begin(), observer->who2_.end());

  return DataFrame::create(Named("times")=when,
      Named("event")=what, Named("who")=who, Named("actor")=who2);
}

// [[Rcpp::export]]
SEXP TestExcessGrowthDistribution() {
  return Rcpp::wrap(hsb::TestExcessGrowthDistribution());
}
