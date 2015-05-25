#include <fstream>
#include <set>
#include "boost/any.hpp"
#include "Rcpp.h"
#include "simple_hazard.hpp"
#include "bugs.hpp"
#include "excess_growth.hpp"
#include "rng.hpp"
#include "smv.hpp"
#include "segment_intersect.hpp"


using namespace Rcpp;
using namespace hsb;


namespace hsb {
class SIRObserver : public hsb::TrajectoryObserver {
 public:
  virtual void event(double when, int what, int64_t who, int64_t who2) {
    when_.push_back(when);
    what_.push_back(what);
    who_.push_back(who);
    who2_.push_back(who2);
  }
  virtual void final() {}

  std::vector<double> when_;
  std::vector<int> what_;
  std::vector<int64_t> who_;
  std::vector<int64_t> who2_;
};



class CallbackEventObserver : public hsb::TrajectoryObserver {
  Rcpp::Function callback_;
  std::vector<double> when_;
  std::vector<int> what_;
  std::vector<int64_t> who_;
  std::vector<int64_t> who2_;
 public:
  CallbackEventObserver(Rcpp::Function callback) : callback_(callback) {}
  virtual ~CallbackEventObserver() {}
  virtual void event(double when, int what, int64_t who, int64_t who2) {
    when_.push_back(when);
    what_.push_back(what);
    who_.push_back(who);
    who2_.push_back(who2);
  }

  virtual void final() {
    NumericVector when(when_.begin(), when_.end());
    IntegerVector what(what_.begin(), what_.end());
    IntegerVector who(who_.begin(), who_.end());
    IntegerVector who2(who2_.begin(), who2_.end());

    auto df=DataFrame::create(Named("when")=when,
      Named("event")=what, Named("whom")=who, Named("who")=who2);
    BOOST_LOG_TRIVIAL(info)<<"Writing results to callback";
    callback_(df);
    when_.clear();
    what_.clear();
    who_.clear();
    who2_.clear();
  }
};
}

struct TKeyWriter {
  std::ofstream keyfile;
  TKeyWriter() : keyfile("keys.txt") {}
  void Write(int64_t id, hsb::simple_hazard::SIRTKey tid) {
      keyfile << id << "\t" << tid ;
  }
};



/*! How many times does the bisector of two units cross a street?
 *  
 *  sunits = matrix of (x,y) for unit locations.
 *  spoints = matrix of (x,y) of endpoints of streets
 *  sstreets = matrix of (p0, p1) index into endpoints defining a street.
 *
 *  \returns Number of crossings for each pair, ordered as a flattened
 *           two-dimensional matrix, suitable for adding to pairwise_distance. 
 */
// [[Rcpp::export]]
SEXP intersections(SEXP sunitsx, SEXP sunitsy, SEXP sendpointsx,
    SEXP sendpointsy, SEXP sstreetsp0, SEXP sstreetsp1) {
  afidd::LogInit("info");

  NumericVector unitsx(sunitsx);
  NumericVector unitsy(sunitsy);
  NumericVector endpointsx(sendpointsx);
  NumericVector endpointsy(sendpointsy);
  IntegerVector streetsp0(sstreetsp0);
  IntegerVector streetsp1(sstreetsp1);
  assert(unitsx.size()==unitsy.size());
  assert(endpointsx.size()==endpointsy.size());
  assert(streetsp0.size()==streetsp1.size());

  size_t unit_cnt=unitsx.size();
  size_t dist_cnt=unit_cnt*(unit_cnt-1)/2;
  size_t street_cnt=streetsp0.size();
  size_t street_pt_cnt=endpointsx.size();
  BOOST_LOG_TRIVIAL(trace) << "units "<<unit_cnt<<" dists "<<dist_cnt<<" streets "
    << street_cnt << " streetpt " << street_pt_cnt ;;
  std::vector<std::pair<double,double>> points(unit_cnt+street_pt_cnt);
  BOOST_LOG_TRIVIAL(trace) << "placing points: ";
  for (size_t uidx=0; uidx<unit_cnt; ++uidx) {
    BOOST_LOG_TRIVIAL(trace) << "("<<unitsx[uidx] << " " << unitsy[uidx] << ") ";
    points[uidx]=std::make_pair(unitsx[uidx], unitsy[uidx]);
  }
  BOOST_LOG_TRIVIAL(trace) ;
  BOOST_LOG_TRIVIAL(trace) << "street points ";
  for (size_t eidx=0; eidx<street_pt_cnt; ++eidx) {
    BOOST_LOG_TRIVIAL(trace) << "("<<endpointsx[eidx] << " " << endpointsy[eidx]<<") ";
    points[unit_cnt+eidx]=std::make_pair(endpointsx[eidx], endpointsy[eidx]);
  }
  // Make line segments
  std::vector<std::pair<size_t,size_t>> segments(dist_cnt+street_cnt);
  BOOST_LOG_TRIVIAL(trace) << "segments between pairs of units" ;
  // Add segments between all pairs of units.
  size_t seg_idx=0;
  for (size_t source_idx=0; source_idx<unit_cnt-1; ++source_idx) {
    for (size_t target_idx=source_idx+1; target_idx<unit_cnt; ++target_idx) {
      segments[seg_idx]=std::make_pair(source_idx, target_idx);
      ++seg_idx;
    }
  }
  // Add segments for streets.
  for (size_t street_idx=0; street_idx<street_cnt; ++street_idx) {
    // -1 to move to 0-based indexing
    segments[seg_idx]=std::make_pair(streetsp0[street_idx]-1+unit_cnt,
      streetsp1[street_idx]-1+unit_cnt);
    ++seg_idx;
  }
  if (seg_idx!=dist_cnt+street_cnt) {
    BOOST_LOG_TRIVIAL(trace) << "You cannot count!" ;
    assert(seg_idx==dist_cnt+street_cnt);
  }

  BOOST_LOG_TRIVIAL(debug) << "Calling segment_intersections";
  std::vector<std::pair<double,double>> intpoints;
    std::multimap<size_t,size_t> intverts;
    std::tie(intpoints, intverts)=segment_intersections_sweep(points, segments);

  BOOST_LOG_TRIVIAL(debug) << "Intersections to pass back ";
  for (auto mmi : intverts) {
    BOOST_LOG_TRIVIAL(trace) << '(' << mmi.first << ',' << mmi.second << ") ";
  }
  BOOST_LOG_TRIVIAL(trace) << "Points going back ";
  for (auto xy : intpoints) {
    BOOST_LOG_TRIVIAL(trace) << '(' << xy.first << ',' << xy.second << ") ";
  }

  IntegerVector crossings(unit_cnt*unit_cnt);

  std::multimap<size_t,size_t>::const_iterator vert_cursor=intverts.cbegin();
  while (vert_cursor!=intverts.cend()) {
    size_t point_idx=vert_cursor->first;
    BOOST_LOG_TRIVIAL(trace) << "looking at "<<point_idx;
    std::set<size_t> street;
    std::set<size_t> arc;
    while (vert_cursor!=intverts.cend() && vert_cursor->first==point_idx) {
      size_t segment_idx=vert_cursor->second;
      if (segment_idx<dist_cnt) {
        BOOST_LOG_TRIVIAL(trace) << "segment "<<segment_idx <<" in arc";
        arc.insert(segment_idx);
      } else {
        BOOST_LOG_TRIVIAL(trace) << "segment "<<segment_idx <<" in street";
        street.insert(segment_idx);
      }
      ++vert_cursor;
    }
    if (arc.size()>0 && street.size()>0) {
      BOOST_LOG_TRIVIAL(trace) << "Adding crossings ";
      for (size_t crossed : arc) {
        const auto& arcref=segments[crossed];
        // Add to i,j and j,i.
        auto idxl=arcref.first*unit_cnt + arcref.second;
        auto idxr=arcref.second*unit_cnt + arcref.first;
        BOOST_LOG_TRIVIAL(trace) << "("<<idxl << " " << idxr << ") " ;
        crossings[idxl]+=1;
        crossings[idxr]+=1;
      }
    }
  }

  return crossings;
}


void ConvertIntParams(Rcpp::List& p, std::map<std::string,boost::any>& params,
    const std::vector<std::string>& int_params) {
  for (auto iparam : int_params) {
    try {
      params[iparam]=int64_t{as<int>(p[iparam])};
    } catch (...) {
      BOOST_LOG_TRIVIAL(error)<<"Could not decode "<<iparam;
      std::stringstream pmsg;
      pmsg << "Could not decode " << iparam;
      throw std::runtime_error(pmsg.str());
    }
  }
}

void ConvertFloatParams(Rcpp::List& p,
    std::map<std::string,boost::any>& params,
    const std::vector<std::string>& float_params) {
  for (auto fparam : float_params) {
    try {
      params[fparam]=double{as<double>(p[fparam])};
    } catch (...) {
      BOOST_LOG_TRIVIAL(error)<<"Could not decode "<<fparam;
      std::stringstream pmsg;
      pmsg << "Could not decode " << fparam;
      throw std::runtime_error(pmsg.str());
    }
  }
}

// [[Rcpp::export]]
SEXP simple_hazard(SEXP pairwise_distanceS, SEXP street_matrixS,
    SEXP parametersS, Rcpp::Function callback) {
  afidd::LogInit("info");

  NumericVector pairwise(pairwise_distanceS);
  std::vector<double> distance(pairwise.begin(), pairwise.end());
  IntegerVector street_matrix(street_matrixS);
  std::vector<int> streets(street_matrix.begin(), street_matrix.end());

  BOOST_LOG_TRIVIAL(debug)<<"Unpacking parameters";

  static const std::vector<std::string> int_params={
    "individual_cnt", "initial", "seed"
  };
  static const std::vector<std::string> float_params={
    "N0", "beta0", "beta1", "beta2", "gamma", "streetfactor",
    "cutoff", "growthrate", "carrying"
  };
  List parameters(parametersS);
  std::map<std::string, boost::any> params;
  ConvertIntParams(parameters, params, int_params);
  ConvertFloatParams(parameters, params, float_params);

  int64_t rand_seed=int64_t{as<int>(parameters["seed"])};

  auto observer=std::make_shared<CallbackEventObserver>(callback);

  BOOST_LOG_TRIVIAL(debug)<<"Call SimpleHazardGSPN";
  RandGen rng(rand_seed);
  auto gspn=hsb::simple_hazard::SimpleHazardGSPN(params, distance,
    streets, rng);

  bool write_keys=false;
  if (write_keys) {
    TKeyWriter w;
    gspn.WriteKeys(w);
  }

  int run_cnt=1;
  if (parameters.containsElementNamed("runs")) {
    run_cnt=as<int>(parameters["runs"]);
    BOOST_LOG_TRIVIAL(trace)<<"has run cnt";
  } else {
    BOOST_LOG_TRIVIAL(trace)<<"doesn't have cnt";
  }

  BOOST_LOG_TRIVIAL(debug)<<"Start SIR_run loop";
  for (int run_idx=0; run_idx<run_cnt; ++run_idx) {
    hsb::simple_hazard::SIR_run(params, gspn, observer, rng);
  }
  return Rcpp::wrap(run_cnt);
}



// [[Rcpp::export]]
DataFrame bugs(SEXP pairwise_distanceS, SEXP street_matrixS,
    SEXP parametersS, Rcpp::Function callback) {
  afidd::LogInit("info");

  NumericVector pairwise(pairwise_distanceS);
  std::vector<double> distance(pairwise.begin(), pairwise.end());
  IntegerVector street_matrix(street_matrixS);
  std::vector<int> streets(street_matrix.begin(), street_matrix.end());

  BOOST_LOG_TRIVIAL(debug)<<"Unpacking parameters";

  static const std::vector<std::string> int_params={
    "individual_cnt", "initial_bug_cnt", "seed", "initial"
  };
  static const std::vector<std::string> float_params={
    "birth", "death", "carrying", "alpha1", "alpha2", "streetfactor",
    "cutoff", "move0", "move1", "gamma"
  };

  List parameters(parametersS);
  std::map<std::string, boost::any> params;
  ConvertIntParams(parameters, params, int_params);
  ConvertFloatParams(parameters, params, float_params);
  BOOST_LOG_TRIVIAL(debug)<<"Unpacked parameters";

  auto observer=std::make_shared<CallbackEventObserver>(callback);

  int64_t rand_seed=int64_t{as<int>(parameters["seed"])};
  RandGen rng(rand_seed);
  auto gspn=hsb::bugs::BugsGSPN(params, distance, streets, rng);

  int run_cnt=1;
  double max_time=-1;
  if (parameters.containsElementNamed("runs")) {
    run_cnt=as<int>(parameters["runs"]);
    BOOST_LOG_TRIVIAL(trace)<<"has run cnt";
  } else {
    BOOST_LOG_TRIVIAL(trace)<<"doesn't have cnt";
  }
  if (parameters.containsElementNamed("max_time")) {
    max_time=as<double>(parameters["max_time"]);
    BOOST_LOG_TRIVIAL(info)<<"has max time " << max_time;
  } else {
    BOOST_LOG_TRIVIAL(trace)<<"doesn't have max time";
  }
  params["max_time"]=max_time;

  for (int run_idx=0; run_idx<run_cnt; ++run_idx) {
    hsb::bugs::SIR_run(params, gspn, observer, rng);
  }
  return Rcpp::wrap(run_cnt);
}

// [[Rcpp::export]]
SEXP TestExcessGrowthDistribution() {
  return Rcpp::wrap(hsb::TestExcessGrowthDistribution());
}


template<typename B>
void fillit(B& callback, int idx) {
  std::vector<int> events { 1, 3, 5, 7, 9};
  Rcpp::IntegerVector sevents(events.size());
  for (size_t i=0; i<events.size(); ++i) {
    sevents[i]=events[i];
  }
  sevents[1]=idx;
  callback(sevents);
}

// [[Rcpp::export]]
SEXP TestCallback(Function callback) {
  fillit(callback, 4);
  fillit(callback, 7);
  return Rcpp::wrap(0);
}


// Hilbert code straight from Wikipedia.
// Hope it's right.
//rotate/flip a quadrant appropriately
void rot(int n, int *x, int *y, int rx, int ry) {
    if (ry == 0) {
        if (rx == 1) {
            *x = n-1 - *x;
            *y = n-1 - *y;
        }
 
        //Swap x and y
        int t  = *x;
        *x = *y;
        *y = t;
    }
}

//convert (x,y) to d
void xy2d (int const* const px, int const* const py, int *const pd, int s, int n) {
  for (int i=0; i<s; ++i) {
    int x=px[i];
    int y=py[i];
    int rx, ry, s, d=0;
    for (s=n/2; s>0; s/=2) {
        rx = (x & s) > 0;
        ry = (y & s) > 0;
        d += s * s * ((3 * rx) ^ ry);
        rot(s, &x, &y, rx, ry);
    }
    pd[i]=d;
  }
}
 
//convert d to (x,y)
void d2xy(int n, int d, int *x, int *y) {
    int rx, ry, s, t=d;
    *x = *y = 0;
    for (s=1; s<n; s*=2) {
        rx = 1 & (t/2);
        ry = 1 & (t ^ rx);
        rot(s, x, y, rx, ry);
        *x += s * rx;
        *y += s * ry;
        t /= 4;
    }
}

// [[Rcpp::export]]
SEXP hilbertXY2D(SEXP sx, SEXP sy, SEXP sn) {
  NumericVector x(sx);
  NumericVector y(sy);
  IntegerVector n(sn);
  std::vector<int> ix(x.size());
  std::vector<int> iy(x.size());
  for (int i=0; i<x.size(); ++i) {
    ix[i]=std::round(x[i]);
    iy[i]=std::round(y[i]);
  }
  std::vector<int> id(x.size());

  xy2d(ix.data(), iy.data(), id.data(), x.size(), n[0]);
  IntegerVector d(id.begin(), id.end());
  return d;
}
