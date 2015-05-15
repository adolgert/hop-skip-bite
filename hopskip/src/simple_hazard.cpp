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
class ExcessInfect0 : public SIRTransition {
  public:
    ExcessInfect0() {}

    virtual std::pair<bool, std::unique_ptr<Dist>>
    Enabled(const UserState& s, const Local& lm,
        double te, double t0, RandGen& rng) override {
      BOOST_LOG_TRIVIAL(trace)<<"ExcessInfect0::Enabled()";
        int64_t I=lm.template Length<0>(0);
        int64_t S=lm.template Length<0>(1);
        if (S>0 && I>0) {
            return {true, std::unique_ptr<ExcessGrowth<RandGen>>(
              new ExcessGrowth<RandGen>(
                s.params.at(Parameter::N0),
                s.params.at(Parameter::carrying),
                s.params.at(Parameter::growthrate),
                s.params.at(Parameter::beta0),
                te))};
        } else {
            return {false, std::unique_ptr<Dist>(nullptr)};
        }
    }

    virtual void Fire(UserState& s, Local& lm, double t0,
        RandGen& rng) override {
      BOOST_LOG_TRIVIAL(trace)<<"ExcessInfect0::Fire()";
        lm.template Move<0,0>(1, 2, 1);
      BOOST_LOG_TRIVIAL(trace)<<"ExcessInfect0::~Fire()";
    }
};



// Infect a hop-distance away with one hazard.
class ExcessInfect1 : public SIRTransition {
  public:
    ExcessInfect1() {}

    virtual std::pair<bool, std::unique_ptr<Dist>>
    Enabled(const UserState& s, const Local& lm,
        double te, double t0, RandGen& rng) override {
      BOOST_LOG_TRIVIAL(trace)<<"ExcessInfect0::Enabled()";
        int64_t I=lm.template Length<0>(0);
        int64_t S=lm.template Length<0>(1);
        if (S>0 && I>0) {
            return {true, std::unique_ptr<ExcessGrowth<RandGen>>(
              new ExcessGrowth<RandGen>(
                s.params.at(Parameter::N0),
                s.params.at(Parameter::carrying),
                s.params.at(Parameter::growthrate),
                s.params.at(Parameter::beta1),
                te))};
        } else {
            return {false, std::unique_ptr<Dist>(nullptr)};
        }
    }

    virtual void Fire(UserState& s, Local& lm, double t0,
        RandGen& rng) override {
      BOOST_LOG_TRIVIAL(trace)<<"ExcessInfect0::Fire()";
        lm.template Move<0,0>(1, 2, 1);
      BOOST_LOG_TRIVIAL(trace)<<"ExcessInfect0::~Fire()";
    }
};

class Notify : public SIRTransition {
  public:
    Notify() {}

    virtual std::pair<bool, std::unique_ptr<Dist>>
    Enabled(const UserState& s, const Local& lm,
        double te, double t0, RandGen& rng) override {
        int64_t I=lm.template Length<0>(0);
        int64_t cryptic=lm.template Length<0>(1);
        if (I>0 && cryptic>0) {
            double rate=s.params.at(Parameter::gamma);
            return {true, std::unique_ptr<ExpDist>(new ExpDist(rate, te))};
        } else {
            return {false, std::unique_ptr<Dist>(nullptr)};
        }
    }

    virtual void Fire(UserState& s, Local& lm, double t0,
        RandGen& rng) override {
      BOOST_LOG_TRIVIAL(trace)<<"Notify::Fire()";
        lm.template Move<0,0>(1, 2, 1);
      BOOST_LOG_TRIVIAL(trace)<<"Notify::~Fire()";
    }
};




void BuildSystem(SIRGSPN& bg, const std::vector<double>& pairwise,
      double cutoff, RandGen& rng) {
  using Edge=BuildGraph<SIRGSPN>::PlaceEdge;
  int64_t cnt=std::lround(std::sqrt(pairwise.size()));

  // disease states are 
  // 0 susceptible, 1 infected, 2 undetected, 3 detected
  for (int64_t pcreate_idx=0; pcreate_idx<cnt; ++pcreate_idx) {
    for (int64_t dcreate_idx=0; dcreate_idx<4; ++dcreate_idx) {
      bg.AddPlace( {pcreate_idx, dcreate_idx}, 0);
    }
  }

  for (int64_t source_idx=0; source_idx<cnt; ++source_idx) {
    auto source=SIRPlace{source_idx, 1}; // infected
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
      {Edge{source, -1}, Edge{SIRPlace{source_idx, 2}, -1},
       Edge{SIRPlace{source_idx, 3}, 1}},
      std::unique_ptr<SIRTransition>(new Notify()));
  }
}



void BuildGrowthSystem(SIRGSPN& bg, const std::vector<double>& pairwise,
    double cutoff, RandGen& rng) {
  using Edge=BuildGraph<SIRGSPN>::PlaceEdge;
  int64_t cnt=std::lround(std::sqrt(pairwise.size()));

  for (int64_t pcreate_idx=0; pcreate_idx<cnt; ++pcreate_idx) {
    for (int64_t dcreate_idx=0; dcreate_idx<4; ++dcreate_idx) {
      bg.AddPlace( {pcreate_idx, dcreate_idx}, 0);
    }
  }

  for (int64_t source_idx=0; source_idx<cnt; ++source_idx) {
    auto source=SIRPlace{source_idx, 1}; // infected
    for (int64_t target_idx=0; target_idx<cnt; ++target_idx) {
      auto susceptible=SIRPlace{target_idx, 0};
      auto infected=SIRPlace{target_idx, 1};
      if (pairwise[source_idx+cnt*target_idx]<cutoff) {
        bg.AddTransition({source_idx, target_idx, 1},
          {Edge{source, -1}, Edge{susceptible, -1},
          Edge{infected, 1}},
          std::unique_ptr<SIRTransition>(new ExcessInfect0()));
      } else {
        bg.AddTransition({source_idx, target_idx, 1},
          {Edge{source, -1}, Edge{susceptible, -1},
          Edge{infected, 1}},
          std::unique_ptr<SIRTransition>(new ExcessInfect1()));
      }
    }

    bg.AddTransition({source_idx, source_idx, 3},
      {Edge{source, -1}, Edge{SIRPlace{source_idx, 2}, -1},
       Edge{SIRPlace{source_idx, 3}, 1}},
      std::unique_ptr<SIRTransition>(new Notify()));
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

  void final(const SIRState& state) {
    observer_->final();
  }

 private:
  const GSPN& gspn_;
  int64_t individual_cnt_;
  std::vector<double> times_;
  std::vector<std::array<int64_t,3>> sir_;
  std::shared_ptr<TrajectoryObserver> observer_;
};




SIRGSPN SimpleHazardGSPN(std::map<std::string, boost::any> params,
    const std::vector<double>& pairwise_distance, RandGen& rng) {
  BOOST_LOG_TRIVIAL(debug)<<"Entering SIR_run";
  assert(params["individual_cnt"].type()==typeid(int64_t));
  int64_t individual_cnt=boost::any_cast<int64_t>(params["individual_cnt"]);

  int64_t place_cnt=3*individual_cnt;
  int64_t transition_cnt=(individual_cnt*individual_cnt)*2 + individual_cnt;
  SIRGSPN gspn(place_cnt+transition_cnt);
  double cutoff=boost::any_cast<double>(params["cutoff"]);
  double growthrate=boost::any_cast<double>(params["growthrate"]);
  if (growthrate<=0) {
    BuildSystem(gspn, pairwise_distance, cutoff, rng);
  } else {
    BuildGrowthSystem(gspn, pairwise_distance, cutoff, rng);
  }
  return gspn;
}


int64_t SIR_run(std::map<std::string, boost::any> params,
    SIRGSPN& gspn, std::shared_ptr<TrajectoryObserver> observer,
    RandGen& rng) {
  BOOST_LOG_TRIVIAL(debug)<<"Entering SIR_run";
  assert(params["individual_cnt"].type()==typeid(int64_t));
  int64_t individual_cnt=boost::any_cast<int64_t>(params["individual_cnt"]);
  int64_t infected_start=boost::any_cast<int64_t>(params["initial"]);

  using Mark=Marking<SIRGSPN::PlaceKey, Uncolored<AnonymousToken>>;
  using SIRState=GSPNState<Mark, SIRGSPN::TransitionKey,WithParams>;

  static const std::map<Parameter,std::string> par_key= {
    {Parameter::N0, "N0"},
    {Parameter::beta0, "beta0"},
    {Parameter::beta1, "beta1"},
    {Parameter::beta2, "beta2"},
    {Parameter::gamma, "gamma"},
    {Parameter::growthrate, "growthrate"},
    {Parameter::carrying, "carrying"}
    };

  SIRState state;

  for (auto& kv : par_key) {
    try {
      state.user.params[kv.first]=boost::any_cast<double>(params[kv.second]);
      BOOST_LOG_TRIVIAL(debug)<<kv.second<<"="<<state.user.params[kv.first];
    } catch (boost::bad_any_cast) {
      BOOST_LOG_TRIVIAL(error)<<"Could not cast "<<kv.second;
    }
  }

  std::uniform_int_distribution<int64_t> random_individual(0, individual_cnt-1);
  //int64_t infected_start=random_individual(rng);
  for (int64_t init_idx=0; init_idx<individual_cnt; ++init_idx) {
    if (init_idx!=infected_start) {
      Add<0>(state.marking, gspn.PlaceVertex({init_idx, 0}), AnonymousToken{});
    } else {
      Add<0>(state.marking, gspn.PlaceVertex({init_idx, 1}), AnonymousToken{});
    }
    // Token in undetected place.
    Add<0>(state.marking, gspn.PlaceVertex({init_idx, 2}), AnonymousToken{});
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