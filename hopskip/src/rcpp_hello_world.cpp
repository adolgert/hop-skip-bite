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
  virtual void final() {}

  std::vector<double> when_;
  std::vector<int> what_;
  std::vector<int64_t> who_;
  std::vector<int64_t> who2_;
};



class CallbackEventObserver : public hsb::TrajectoryObserver {
  Rcpp::Function callback_;
  std::vector<double> when_;
  std::vector<int> what_;
  std::vector<int64_t> who_;
  std::vector<int64_t> who2_;
 public:
  CallbackEventObserver(Rcpp::Function callback) : callback_(callback) {}
  virtual ~CallbackEventObserver() {}
  virtual void event(double when, int what, int64_t who, int64_t who2) {
    when_.push_back(when);
    what_.push_back(what);
    who_.push_back(who);
    who2_.push_back(who2);
  }

  virtual void final() {
    NumericVector when(when_.begin(), when_.end());
    IntegerVector what(what_.begin(), what_.end());
    IntegerVector who(who_.begin(), who_.end());
    IntegerVector who2(who2_.begin(), who2_.end());

    auto df=DataFrame::create(Named("when")=when,
      Named("event")=what, Named("whom")=who, Named("who")=who2);
    BOOST_LOG_TRIVIAL(debug)<<"Writing results to callback";
    callback_(df);
    when_.clear();
    what_.clear();
    who_.clear();
    who2_.clear();
  }
};
}


// [[Rcpp::export]]
SEXP simple_hazard(SEXP pairwise_distanceS, SEXP parametersS,
    Rcpp::Function callback) {
  afidd::LogInit("info");

  NumericVector pairwise(pairwise_distanceS);
  std::vector<double> distance(pairwise.begin(), pairwise.end());

  List parameters(parametersS);
  std::map<std::string, boost::any> params;
  params["individual_cnt"]=int64_t{as<int>(parameters["individual_cnt"])};
  params["initial"]=int64_t{as<int>(parameters["initial"])};
  params["N0"]=double{as<double>(parameters["N0"])};
  params["beta0"]=double{as<double>(parameters["beta0"])};
  params["beta1"]=double{as<double>(parameters["beta1"])};
  params["beta2"]=double{as<double>(parameters["beta2"])};
  params["gamma"]=double{as<double>(parameters["gamma"])};
  params["cutoff"]=double{as<double>(parameters["cutoff"])};
  params["growthrate"]=double{as<double>(parameters["growthrate"])};
  params["carrying"]=double{as<double>(parameters["carrying"])};
  int64_t rand_seed=int64_t{as<int>(parameters["seed"])};
  params["seed"]=rand_seed;

  auto observer=std::make_shared<CallbackEventObserver>(callback);

  RandGen rng(rand_seed);
  auto gspn=hsb::simple_hazard::SimpleHazardGSPN(params, distance, rng);

  int run_cnt=1;
  if (parameters.containsElementNamed("runs")) {
    run_cnt=as<int>(parameters["runs"]);
    BOOST_LOG_TRIVIAL(debug)<<"has run cnt";
  } else {
    BOOST_LOG_TRIVIAL(debug)<<"doesn't have cnt";
  }

  for (int run_idx=0; run_idx<run_cnt; ++run_idx) {
    hsb::simple_hazard::SIR_run(params, gspn, observer, rng);
  }
  return Rcpp::wrap(0);
}



// [[Rcpp::export]]
DataFrame bugs(SEXP pairwise_distanceS, SEXP parametersS) {
  afidd::LogInit("info");

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


template<typename B>
void fillit(B& callback, int idx) {
  std::vector<int> events { 1, 3, 5, 7, 9};
  Rcpp::IntegerVector sevents(events.size());
  for (size_t i=0; i<events.size(); ++i) {
    sevents[i]=events[i];
  }
  sevents[1]=idx;
  callback(sevents);
}

// [[Rcpp::export]]
SEXP TestCallback(Function callback) {
  fillit(callback, 4);
  fillit(callback, 7);
  return Rcpp::wrap(0);
}
