#include <string>
#include <sstream>
#include "boost/program_options.hpp"



int main(int argc, char *argv[]) {
  namespace po=boost::program_options;
  po::options_description desc("Hop-skip for infection");

  std::map<std::string, boost::any> params;
  params["individual_cnt"]=int64_t{10};
  params["beta0"]=double{0.1};
  params["beta1"]=double{0.1};
  params["beta2"]=double{0.1};
  params["gamma"]=double{0.1};

  int64_t rand_seed=33333;
  RandGen rng(rand_seed);
  simple_hazard::SIR_run(params, rng);
  return 0;
}
