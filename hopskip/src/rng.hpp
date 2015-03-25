#ifndef _RNG_HPP_
#define _RNG_HPP_ 1

#include "boost/random/mersenne_twister.hpp"

namespace hsb {
// Defines the random number generator for the simulation.
using RandGen=boost::mt19937;
}
#endif
