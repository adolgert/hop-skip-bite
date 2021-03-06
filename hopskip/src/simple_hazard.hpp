#ifndef _SIMPLE_HAZARD_HPP_
#define _SIMPLE_HAZARD_HPP_ 1

#include <map>
#include <string>
#include <memory>
#include "boost/any.hpp"
#include "smv.hpp"
#include "rng.hpp"
#include "trajectory_observer.hpp"



namespace hsb {
namespace simple_hazard {

namespace smv=afidd::smv;
using namespace smv;

enum class Parameter : int { none, beta0, beta1, beta2, gamma,
  N0, carrying, growthrate, growthscale, streetfactor };

// A token is an instance of this class.
struct AnonymousToken {
    AnonymousToken()=default;
    inline friend
    std::ostream& operator<<(std::ostream& os, const AnonymousToken& at) {
        return os << "T";
    }
};


// Identifies a place uniquely.
struct SIRPlace
{
  int64_t individual;
  int64_t disease;
  SIRPlace()=default;
  SIRPlace(int64_t i, int64_t d) : individual(i), disease(d) {}
  friend inline
  bool operator<(const SIRPlace& a, const SIRPlace& b) {
    return afidd::smv::LazyLess(a.individual, b.individual,
        a.disease, b.disease);
  }

  friend inline
  bool operator==(const SIRPlace& a, const SIRPlace& b) {
    return (a.individual == b.individual) && (a.disease==b.disease);
  }

  friend inline
  std::ostream& operator<<(std::ostream& os, const SIRPlace& cp) {
    return os << '(' << cp.individual << ", " << cp.disease << ')';
  }
};



// This identifies a transition uniquely.
struct SIRTKey {
    int64_t i;
    int64_t j;
    int64_t kind;
    SIRTKey()=default;
    SIRTKey(int64_t i, int64_t j, int64_t kind) : i(i), j(j), kind(kind) {}

    friend inline
    bool operator<(const SIRTKey& a, const SIRTKey& b) {
        return afidd::smv::LazyLess(a.i, b.i, a.j, b.j, a.kind, b.kind);
    }

    friend inline
    bool operator==(const SIRTKey& a, const SIRTKey& b) {
        return (a.i==b.i) && (a.j==b.j) && (a.kind==b.kind);
    }

    friend inline
    std::ostream& operator<<(std::ostream& os, const SIRTKey& cp) {
        return os << '(' << cp.i << ',' << cp.j << ',' << cp.kind << ')';
    }
};


// This is as much of the marking as the transition will see.
using Local=LocalMarking<Uncolored<AnonymousToken>>;
// Extra state to add to the system state. Will be passed to transitions.
struct WithParams {
  // Put our parameters here.
  std::map<Parameter,double> params;
};


// The transition needs to know the local marking and any extra state.
using SIRTransition=ExplicitTransition<Local,RandGen,WithParams>;

using Dist=TransitionDistribution<RandGen>;
using ExpDist=ExponentialDistribution<RandGen>;



using SIRGSPN=
    ExplicitTransitions<SIRPlace, SIRTKey, Local, RandGen, WithParams>;


void BuildSystem(SIRGSPN& bg, const std::vector<double>& pairwise,
  const std::vector<int>& streets, 
  double cutoff, double street_factor, RandGen& rng);

void BuildGrowthSystem(SIRGSPN& bg, const std::vector<double>& pairwise,
  const std::vector<int>& streets, 
  double cutoff, double street_factor, RandGen& rng);

SIRGSPN SimpleHazardGSPN(std::map<std::string, boost::any> params,
  const std::vector<double>& pairwise_distance,
  const std::vector<int>& streets, RandGen& rng);

int64_t SIR_run(std::map<std::string, boost::any> params,
  SIRGSPN& gspn,
  std::shared_ptr<TrajectoryObserver> observer, RandGen& rng);

}
}

#endif // _SIMPLE_HAZARD_HPP_
