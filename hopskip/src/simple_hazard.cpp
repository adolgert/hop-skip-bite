#include <iostream>
#include <map>
#include <cmath>
#include <string>
#include <vector>
#include "boost/any.hpp"
#include "smv.hpp"
#include "spatial_process.hpp"
#include "simple_hazard.hpp"
#include "excess_growth.hpp"

namespace smv=afidd::smv;
using namespace smv;

namespace hsb {
namespace simple_hazard {

enum class Parameter : int { none, beta0, beta1, beta2, gamma,
  N0, carrying, growthrate };

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
            double rate=s.params.at(Parameter::beta0);
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
            double rate=s.params.at(Parameter::beta1);
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


// Infect a hop-distance away with one hazard.
class ExcessInfect : public SIRTransition {
  public:
    ExcessInfect() {}

    virtual std::pair<bool, std::unique_ptr<Dist>>
    Enabled(const UserState& s, const Local& lm,
        double te, double t0, RandGen& rng) override {
        int64_t I=lm.template Length<0>(0);
        int64_t S=lm.template Length<0>(1);
        if (S>0 && I>0) {
            return {true, std::unique_ptr<ExcessGrowth<RandGen>>(
              new ExcessGrowth<RandGen>(
                s.params.at(Parameter::N0),
                s.params.at(Parameter::carrying),
                s.params.at(Parameter::growthrate),
                te))};
        } else {
            return {false, std::unique_ptr<Dist>(nullptr)};
        }
    }

    virtual void Fire(UserState& s, Local& lm, double t0,
        RandGen& rng) override {
        lm.template Move<0,0>(1, 2, 1);
    }
};


class Notify : public SIRTransition {
  public:
    Notify() {}

    virtual std::pair<bool, std::unique_ptr<Dist>>
    Enabled(const UserState& s, const Local& lm,
        double te, double t0, RandGen& rng) override {
        int64_t I=lm.template Length<0>(0);
        if (I>0) {
            double rate=s.params.at(Parameter::gamma);
            return {true, std::unique_ptr<ExpDist>(new ExpDist(rate, te))};
        } else {
            return {false, std::unique_ptr<Dist>(nullptr)};
        }
    }

    virtual void Fire(UserState& s, Local& lm, double t0,
        RandGen& rng) override {
        lm.template Move<0,0>(0, 1, 1);
    }
};


using SIRGSPN=
    ExplicitTransitions<SIRPlace, SIRTKey, Local, RandGen, WithParams>;


void BuildSystem(SIRGSPN& bg, const std::vector<double>& pairwise,
      double cutoff, RandGen& rng) {
    using Edge=BuildGraph<SIRGSPN>::PlaceEdge;
    int64_t cnt=std::lround(std::sqrt(pairwise.size()));

    for (int64_t pcreate_idx=0; pcreate_idx<cnt; ++pcreate_idx) {
        for (int64_t dcreate_idx=0; dcreate_idx<3; ++dcreate_idx) {
            bg.AddPlace( {pcreate_idx, dcreate_idx}, 0);
        }
    }

    for (int64_t source_idx=0; source_idx<cnt; ++source_idx) {
        auto source=SIRPlace{source_idx, 1};
        for (int64_t target_idx=0; target_idx<cnt; ++target_idx) {
            auto susceptible=SIRPlace{target_idx, 0};
            auto infected=SIRPlace{target_idx, 1};
            if (pairwise[source_idx+cnt*target_idx]<cutoff) {
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

        bg.AddTransition({source_idx, source_idx, 3},
          {Edge{source, -1}, Edge{SIRPlace{source_idx, 2}, 1}},
          std::unique_ptr<SIRTransition>(new Notify()));

        auto nsource=SIRPlace{source_idx, 2};
        for (int64_t target_idx=0; target_idx<cnt; ++target_idx) {
            auto susceptible=SIRPlace{target_idx, 0};
            auto infected=SIRPlace{target_idx, 1};
            if (pairwise[source_idx+cnt*target_idx]<cutoff) {
                bg.AddTransition({source_idx, target_idx, 2},
                    {Edge{nsource, -1}, Edge{susceptible, -1},
                    Edge{infected, 1}},
                    std::unique_ptr<SIRTransition>(new Infect0()));
            } else {
                bg.AddTransition({source_idx, target_idx, 2},
                    {Edge{nsource, -1}, Edge{susceptible, -1},
                    Edge{infected, 1}},
                    std::unique_ptr<SIRTransition>(new Infect1()));
            }
        }
    }
}



template<typename GSPN, typename SIRState>
struct SIROutput {
 public:
  SIROutput(const GSPN& gspn, int64_t individual_cnt,
      std::shared_ptr<TrajectoryObserver> observer)
  : gspn_(gspn), individual_cnt_(individual_cnt), observer_(observer),
  sir_(individual_cnt, std::array<int64_t,3>{0, -1, -1}) {}

  bool operator()(const SIRState& state) {
    times_.push_back(state.CurrentTime());
    size_t time_idx=times_.size();

    auto transition=gspn_.VertexTransition(state.last_transition);
    observer_->event(state.CurrentTime(), transition.kind, transition.j,
      transition.i);
    switch (transition.kind) {
      case 0:
        std::get<1>(sir_[transition.j])=time_idx;
        break;

      case 1:
        std::get<1>(sir_[transition.j])=time_idx;
        break;
        
      case 2:
        std::get<1>(sir_[transition.j])=time_idx;
        break;
        
      case 3:
        std::get<2>(sir_[transition.i])=time_idx;
        break;

      default:
        assert(false);
        break;
    }
    return true;
  }

  void initial(const SIRState& state) {
    for (int64_t ind_idx=0; ind_idx<individual_cnt_; ++ind_idx) {
      int64_t inf_place=gspn_.PlaceVertex({ind_idx, 1});
      if (Length<0>(state.marking, inf_place)>0) {
        observer_->event(0, 1, ind_idx, 0);
        std::get<0>(sir_[ind_idx])=-1;
        std::get<1>(sir_[ind_idx])=0;
      }
    }
  }

  void final(const SIRState& state) {}

 private:
  const GSPN& gspn_;
  int64_t individual_cnt_;
  std::vector<double> times_;
  std::vector<std::array<int64_t,3>> sir_;
  std::shared_ptr<TrajectoryObserver> observer_;
};



int64_t SIR_run(std::map<std::string, boost::any> params,
    const std::vector<double>& pairwise_distance,
    std::shared_ptr<TrajectoryObserver> observer, RandGen& rng) {
  assert(params["individual_cnt"].type()==typeid(int64_t));
  int64_t individual_cnt=boost::any_cast<int64_t>(params["individual_cnt"]);

  int64_t place_cnt=3*individual_cnt;
  int64_t transition_cnt=(individual_cnt*individual_cnt)*2 + individual_cnt;
  SIRGSPN gspn(place_cnt+transition_cnt);
  double cutoff=boost::any_cast<double>(params["cutoff"]);
  BuildSystem(gspn, pairwise_distance, cutoff, rng);

  using Mark=Marking<SIRGSPN::PlaceKey, Uncolored<AnonymousToken>>;
  using SIRState=GSPNState<Mark, SIRGSPN::TransitionKey,WithParams>;

  SIRState state;

  state.user.params[Parameter::beta0]=boost::any_cast<double>(params["beta0"]);
  state.user.params[Parameter::beta1]=boost::any_cast<double>(params["beta1"]);
  state.user.params[Parameter::beta2]=boost::any_cast<double>(params["beta2"]);
  state.user.params[Parameter::gamma]=boost::any_cast<double>(params["gamma"]);

  std::uniform_int_distribution<int64_t> random_individual(0, individual_cnt-1);
  int64_t infected_start=random_individual(rng);
  for (int64_t init_idx=0; init_idx<individual_cnt; ++init_idx) {
    if (init_idx!=infected_start) {
      Add<0>(state.marking, gspn.PlaceVertex({init_idx, 0}), AnonymousToken{});
    } else {
      Add<0>(state.marking, gspn.PlaceVertex({init_idx, 1}), AnonymousToken{});
    }
  }

  //using Propagator=PropagateCompetingProcesses<int64_t,RandGen>;
  using Propagator=NonHomogeneousPoissonProcesses<int64_t,RandGen>;
  Propagator competing;
  using Dynamics=StochasticDynamics<SIRGSPN,SIRState,RandGen>;
  Dynamics dynamics(gspn, {&competing});

  SIROutput<SIRGSPN,SIRState> output_function(gspn, individual_cnt, observer);

  dynamics.Initialize(&state, &rng);

  BOOST_LOG_TRIVIAL(info)<<"Starting main loop";
  bool running=true;
  auto nothing=[](SIRState&)->void {};
  double last_time=state.CurrentTime();
  while (running) {
    running=dynamics(state);
    if (running) {
      double new_time=state.CurrentTime();
      if (new_time-last_time<-1e-12) {
        BOOST_LOG_TRIVIAL(warning) << "last time "<<last_time <<" "
          << " new_time "<<new_time;
      }
      last_time=new_time;
      running=output_function(state);
    } else {
      BOOST_LOG_TRIVIAL(info)<<"No transitions left to fire "
          <<state.CurrentTime();
    }
    // auto v=competing.content_size();
    // SMVLOG(BOOST_LOG_TRIVIAL(debug)<<"Competing Processes: size "
    //     <<v.first<<" infinities "<<v.second);
  }
  BOOST_LOG_TRIVIAL(info)<<"Reached end time "<<state.CurrentTime();
  output_function.final(state);
  return 0;
}



} // namespace simple_hazard
} // namespace hsb