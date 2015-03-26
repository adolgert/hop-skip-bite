#ifndef _BUGS_HPP_
#define _BUGS_HPP_ 1

#include <map>
#include <string>
#include <memory>
#include "boost/any.hpp"
#include "rng.hpp"
#include "trajectory_observer.hpp"



namespace hsb {
namespace bugs {

int64_t SIR_run(std::map<std::string, boost::any> params,
  const std::vector<double>& pairwise_distance,
  std::shared_ptr<TrajectoryObserver> observer, RandGen& rng);

}
}

#endif // _BUGS_HPP_
