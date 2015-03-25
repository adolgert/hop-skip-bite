#ifndef _BUGS_HPP_
#define _BUGS_HPP_ 1

#include <map>
#include <string>
#include "boost/any.hpp"
#include "rng.hpp"



namespace hsb {
namespace bugs {

int64_t SIR_run(std::map<std::string, boost::any> params, RandGen& rng);

}
}

#endif // _BUGS_HPP_
