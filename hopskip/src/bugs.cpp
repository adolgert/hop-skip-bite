#include <iostream>
#include <map>
#include <string>
#include <vector>
#include "boost/any.hpp"
#include "smv.hpp"
#include "spatial_process.hpp"
#include "bugs.hpp"

namespace smv=afidd::smv;
using namespace smv;

namespace hsb {
namespace bugs {

enum class Parameter : int { none, birth, death, carrying, move0,
    move1, gamma };

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
class Birth : public SIRTransition {
  public:
    Birth() {}

    virtual std::pair<bool, std::unique_ptr<Dist>>
    Enabled(const UserState& s, const Local& lm,
        double te, double t0, RandGen& rng) override {
        int64_t S=lm.template Length<0>(0);
        if (S>0) {
            double rate=S*s.params.at(Parameter::birth)*
              (1-S/s.params.at(Parameter::carrying));
            return {true, std::unique_ptr<ExpDist>(new ExpDist(rate, te))};
        } else {
            return {false, std::unique_ptr<Dist>(nullptr)};
        }
    }

    virtual void Fire(UserState& s, Local& lm, double t0,
        RandGen& rng) override {
        lm.template Add<0>(0, AnonymousToken{});
    }
};



// Infect a hop-distance away with one hazard.
class Death : public SIRTransition {
  public:
    Death() {}

    virtual std::pair<bool, std::unique_ptr<Dist>>
    Enabled(const UserState& s, const Local& lm,
        double te, double t0, RandGen& rng) override {
        int64_t S=lm.template Length<0>(0);
        if (S>0) {
            double rate=S*s.params.at(Parameter::death);
            return {true, std::unique_ptr<ExpDist>(new ExpDist(rate, te))};
        } else {
            return {false, std::unique_ptr<Dist>(nullptr)};
        }
    }

    virtual void Fire(UserState& s, Local& lm, double t0,
        RandGen& rng) override {
        lm.template Remove<0>(0, 1, rng);
    }
};



// Infect a hop-distance away with one hazard.
class Move0 : public SIRTransition {
  public:
    Move0() {}

    virtual std::pair<bool, std::unique_ptr<Dist>>
    Enabled(const UserState& s, const Local& lm,
        double te, double t0, RandGen& rng) override {
        int64_t I=lm.template Length<0>(0);
        if (I>0) {
            double rate=I*s.params.at(Parameter::move0);
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


// Infect a hop-distance away with one hazard.
class Move1 : public SIRTransition {
  public:
    Move1() {}

    virtual std::pair<bool, std::unique_ptr<Dist>>
    Enabled(const UserState& s, const Local& lm,
        double te, double t0, RandGen& rng) override {
        int64_t I=lm.template Length<0>(0);
        if (I>0) {
            double rate=I*s.params.at(Parameter::move1);
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




class Notify : public SIRTransition {
  public:
    Notify() {}

    virtual std::pair<bool, std::unique_ptr<Dist>>
    Enabled(const UserState& s, const Local& lm,
        double te, double t0, RandGen& rng) override {
        int64_t infested=lm.template Length<0>(0);
        int64_t cryptic=lm.template Length<0>(1);
        if (infested>0 && cryptic>0) {
            double rate=s.params.at(Parameter::gamma);
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


void BuildSystem(SIRGSPN& bg, const std::vector<double>& pairwise,
      double cutoff, RandGen& rng) {
    using Edge=BuildGraph<SIRGSPN>::PlaceEdge;
    int64_t cnt=std::lround(std::sqrt(pairwise.size()));

    for (int64_t pcreate_idx=0; pcreate_idx<cnt; ++pcreate_idx) {
      int64_t house=0;
      bg.AddPlace( {pcreate_idx, house}, 0);
      int64_t cryptic=1;
      bg.AddPlace( {pcreate_idx, cryptic}, 0);
      int64_t notified=2;
      bg.AddPlace( {pcreate_idx, notified}, 0);
    }

    for (int64_t source_idx=0; source_idx<cnt; ++source_idx) {
        auto source=SIRPlace{source_idx, 0};

        bg.AddTransition({source_idx, source_idx, 5},
          {Edge{source, -1}, Edge{source, 1}},
          std::unique_ptr<SIRTransition>(new Birth()));

        bg.AddTransition({source_idx, source_idx, 6},
          {Edge{source, -1}},
          std::unique_ptr<SIRTransition>(new Death()));


        for (int64_t target_idx=0; target_idx<cnt; ++target_idx) {
            auto target=SIRPlace{target_idx, 0};
            if (pairwise[source_idx+cnt*target_idx]<cutoff) {
                bg.AddTransition({source_idx, target_idx, 1},
                    {Edge{source, -1}, Edge{target, 1}},
                    std::unique_ptr<SIRTransition>(new Move0()));
            } else {
                bg.AddTransition({source_idx, target_idx, 1},
                    {Edge{source, -1}, Edge{target, 1}},
                    std::unique_ptr<SIRTransition>(new Move1()));
            }
        }

        auto cryptic=SIRPlace{source_idx, 1};
        auto notified=SIRPlace{source_idx, 2};
        bg.AddTransition({source_idx, source_idx, 3},
          {Edge{source, -1}, Edge{cryptic, -1}, Edge{notified, 1}},
          std::unique_ptr<SIRTransition>(new Notify()));
    }
}



template<typename GSPN, typename SIRState>
struct SIROutput {
 public:
  SIROutput(const GSPN& gspn, int64_t individual_cnt,
    std::shared_ptr<TrajectoryObserver> observer)
  : gspn_(gspn), individual_cnt_(individual_cnt), observer_(observer),
    infected_cnt_(0) {}

  bool operator()(const SIRState& state) {
    auto transition=gspn_.VertexTransition(state.last_transition);
    switch (transition.kind) {
      case 1: // move
        {
          auto target_house=gspn_.PlaceVertex({transition.j, 0});
          if (Length<0>(state.marking, target_house)==1) {
            // infection event.
            observer_->event(state.CurrentTime(), 1, transition.j,
              transition.i);
            infected_cnt_+=1;
          } else {
            // ignore infection of already infected
          }
        }
        break;

      case 6:
        {
          auto perish_house=gspn_.PlaceVertex({transition.i, 0});
          if (Length<0>(state.marking, perish_house)==0) {
            // extinction event.
            observer_->event(state.CurrentTime(), 4, transition.i, 0);
            infected_cnt_-=1;
          } else {
            // ignore death of an individual if not last one.
          }
        }
        break;
                
      case 3:
        // notification event
        observer_->event(state.CurrentTime(), 3, transition.i, 0);
        break;

      default:
        assert(false);
        break;
    }
    return infected_cnt_<individual_cnt_;
  }

  void initial(const SIRState& state) {}

  void final(const SIRState& state) {}

 private:
  const GSPN& gspn_;
  int64_t individual_cnt_;
  int64_t infected_cnt_;
  std::shared_ptr<TrajectoryObserver> observer_;
};



int64_t SIR_run(std::map<std::string, boost::any> params,
    const std::vector<double>& pairwise_distance,
    std::shared_ptr<TrajectoryObserver> observer, RandGen& rng) {
  assert(params["individual_cnt"].type()==typeid(int64_t));
  using boost::any_cast;
  int64_t individual_cnt=boost::any_cast<int64_t>(params["individual_cnt"]);

  int64_t place_cnt=3*individual_cnt;
  int64_t transition_cnt=(individual_cnt*individual_cnt)*2 + individual_cnt;
  SIRGSPN gspn(place_cnt+transition_cnt);
  double cutoff=boost::any_cast<double>(params["cutoff"]);
  BuildSystem(gspn, pairwise_distance, cutoff, rng);

  using Mark=Marking<SIRGSPN::PlaceKey, Uncolored<AnonymousToken>>;
  using SIRState=GSPNState<Mark, SIRGSPN::TransitionKey,WithParams>;

  SIRState state;

  try {
    state.user.params[Parameter::birth]=any_cast<double>(params["birth"]);
    state.user.params[Parameter::death]=any_cast<double>(params["death"]);
    state.user.params[Parameter::carrying]=any_cast<double>(params["carrying"]);
    state.user.params[Parameter::move0]=any_cast<double>(params["move0"]);
    state.user.params[Parameter::move1]=any_cast<double>(params["move1"]);
    state.user.params[Parameter::gamma]=any_cast<double>(params["gamma"]);
  } catch (boost::bad_any_cast& e) {
    BOOST_LOG_TRIVIAL(error) << "Bad any_cast in bugs::SIR_Run";
  }

  int64_t initial_bug_count=boost::any_cast<int64_t>(params["initial_bug_cnt"]);

  std::uniform_int_distribution<int64_t> random_individual(0, individual_cnt-1);
  int64_t infected_start=random_individual(rng);
  for (int64_t init_idx=0; init_idx<initial_bug_count; ++init_idx) {
    Add<0>(state.marking, gspn.PlaceVertex({infected_start, 0}),
        AnonymousToken{});
  }

  //using Propagator=PropagateCompetingProcesses<int64_t,RandGen>;
  using Propagator=NonHomogeneousPoissonProcesses<int64_t,RandGen>;
  Propagator competing;
  using Dynamics=StochasticDynamics<SIRGSPN,SIRState,RandGen>;
  Dynamics dynamics(gspn, {&competing});

  SIROutput<SIRGSPN,SIRState> output_function(gspn, individual_cnt,
    observer);

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



} // namespace bugs
} // namespace hsb
