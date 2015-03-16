#include <iostream>
#include <map>
#include <string>
#include "smv.hpp"
#include "rng.hpp"

namespace smv=afidd::smv;
using namespace smv;

namespace hsb {
namespace simple_hazard {

// A token is an instance of this class.
struct AnonymousToken {
    AnonymousToken()=default;
    inline friend
    std::ostream& operator<<(std::ostream& os, const AnonymousToken& at) {
        return os << "T";
    }    I
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
using Local=LocalMarking<Uncolored<IndividualToken>>;
// Extra state to add to the system state. Will be passed to transitions.
struct WithParams {
  // Put our parameters here.
  std::map<std::string,double> params;
};


// The transition needs to know the local marking and any extra state.
using SIRTransition=ExplicitTransition<Local,RandGen,WithParams>;

using Dist=TransitionDistribution<RandGen>;
using ExpDist=ExponentialDistribution<RandGen>;


// Infect a hop-distance away with one hazard.
class Infect0 : public SIRTransition {
  public:
    Infect0() {}

    virtual std::pair<bool, std::unique_ptr<Dist>>
    Enabled(const UserState& s, const Local& lm,
        double te, double t0, RandGen& rng) override {
        int64_t I=lm.template Length<0>(0);
        int64_t S=lm.template Length<0>(1);
        if (S>0 && I>0) {
            double rate=s.params["beta0"];
            return {true, std::unique_ptr<ExpDist>(new ExpDist(rate, te))};
        } else {
            return {false, std::unique_ptr<Dist>(nullptr)};
        }
    }

    virtual void Fire(UserState& s, Local& lm, double t0,
        RandGen& rng) override {
        lm.template Move<0,0>(1, 2, 1);
    }
};



// Infect at other distances with this hazard.
class Infect1 : public SIRTransition {
  public:
    Infect1() {}

    virtual std::pair<bool, std::unique_ptr<Dist>>
    Enabled(const UserState& s, const Local& lm,
        double te, double t0, RandGen& rng) override {
        int64_t I=lm.template Length<0>(0);
        int64_t S=lm.template Length<0>(1);
        if (S>0 && I>0) {
            double rate=s.params["beta1"];
            return {true, std::unique_ptr<ExpDist>(new ExpDist(rate, te))};
        } else {
            return {false, std::unique_ptr<Dist>(nullptr)};
        }
    }

    virtual void Fire(UserState& s, Local& lm, double t0,
        RandGen& rng) override {
        lm.template Move<0,0>(1, 2, 1);
    }
};


using SIRGSPN=
    ExplicitTransitions<SIRPlace, SIRTKey, Local, RandGen, WithParams>;


void BuildSystem(SIRGSPN& bg, std::vector<std::tuple<double,2>>& point,
    RandGen& rng) {
    using Edge=BuildGraph<SIRGSPN>::PlaceEdge;

    for (int64_t pcreate_idx=0; pcreate_idx<point.size(); ++pcreate_idx) {
        for (int64_t dcreate_idx=0; dcreate_idx<3; ++dcreate_idx) {
            bg.AddPlace( {pcreate_idx, dcreate_idx}, 0);
        }
    }

    for (int64_t source_idx=0; source_idx<point.size(); ++source_idx) {
        auto source=SIRPlace{source_idx, 1};
        double source_x=std::get<0>(point[source_idx]);
        double source_y=std::get<1>(point[source_idx]);
        for (int64_t target_idx=0; target_idx<point.size(); ++target_idx) {
            auto susceptible=SIRPlace{target_idx, 0};
            auto infected=SIRPlace{target_idx, 1};
            double target_x=std::get<0>(point[target_idx]);
            double target_y=std::get<1>(point[target_idx]);
            if ((source_x-target_x)**2 + (source_y-target_y)**2<0.09) {
                bg.AddTransition({source_idx, target_idx, 1},
                    {Edge{source, -1}, Edge{susceptible, -1},
                    Edge{infected, 1}},
                    std::unique_ptr<SIRTransition>(new Infect0()));
            } else {
                bg.AddTransition({source_idx, target_idx, 1},
                    {Edge{source, -1}, Edge{susceptible, -1},
                    Edge{infected, 1}},
                    std::unique_ptr<SIRTransition>(new Infect1()));
            }
        }
    }
}







} // namespace simple_hazard
} // namespace hsb