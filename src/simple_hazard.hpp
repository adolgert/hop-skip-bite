#ifndef _SIMPLE_HAZARD_HPP_
#define _SIMPLE_HAZARD_HPP_ 1

#include <map>
#include <string>
#include "boost/any.hpp"
#include "rng.hpp"



namespace hsb {
namespace simple_hazard {

int64_t SIR_run(std::map<std::string, boost::any> params, RandGen& rng);

}
}

#endif // _SIMPLE_HAZARD_HPP_
