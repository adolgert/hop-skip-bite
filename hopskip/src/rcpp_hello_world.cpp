#include "boost/any.hpp"
#include "Rcpp.h"
#include "simple_hazard.hpp"
#include "bugs.hpp"
#include "rng.hpp"

using namespace Rcpp;
using namespace hsb;

// [[Rcpp::export]]
List rcpp_hello_world() {
  boost::any holdit;
  CharacterVector x = CharacterVector::create( "foo", "bar" )  ;
  NumericVector y   = NumericVector::create( 0.0, 1.0 ) ;
  List z            = List::create( x, y ) ;

  std::map<std::string, boost::any> params;
  params["individual_cnt"]=int64_t{10};
  params["beta0"]=double{0.1};
  params["beta1"]=double{0.1};
  params["beta2"]=double{0.1};
  params["gamma"]=double{0.1};

  int64_t rand_seed=33333;
  RandGen rng(rand_seed);
  hsb::simple_hazard::SIR_run(params, rng);

  return z ;
}
