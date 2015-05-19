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
    BOOST_LOG_TRIVIAL(debug)<<"Writing results to callback";
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
      keyfile << id << "\t" << tid << std::endl;
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
  afidd::LogInit("debug");

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
  std::cout << "units "<<unit_cnt<<" dists "<<dist_cnt<<" streets "
    << street_cnt << " streetpt " << street_pt_cnt << std::endl;;
  std::vector<std::pair<double,double>> points(unit_cnt+street_pt_cnt);
  std::cout << "placing points: ";
  for (size_t uidx=0; uidx<unit_cnt; ++uidx) {
    std::cout << "("<<unitsx[uidx] << " " << unitsy[uidx] << ") ";
    points[uidx]=std::make_pair(unitsx[uidx], unitsy[uidx]);
  }
  std::cout << std::endl;
  std::cout << "street points ";
  for (size_t eidx=0; eidx<street_pt_cnt; ++eidx) {
    std::cout << "("<<endpointsx[eidx] << " " << endpointsy[eidx]<<") ";
    points[unit_cnt+eidx]=std::make_pair(endpointsx[eidx], endpointsy[eidx]);
  }
  std::cout << std::endl;
  // Make line segments
  std::vector<std::pair<size_t,size_t>> segments(dist_cnt+street_cnt);
  std::cout << "segments between pairs of units" << std::endl;
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
    std::cout << "You cannot count!" << std::endl;
    assert(seg_idx==dist_cnt+street_cnt);
  }

  std::cout << "Calling segment_intersections"<<std::endl;
  std::vector<std::pair<double,double>> intpoints;
    std::multimap<size_t,size_t> intverts;
    std::tie(intpoints, intverts)=segment_intersections(points, segments);

  std::cout << "Intersections to pass back ";
  for (auto mmi : intverts) {
    std::cout << '(' << mmi.first << ',' << mmi.second << ") ";
  }
  std::cout << std::endl;
  std::cout << "Points going back ";
  for (auto xy : intpoints) {
    std::cout << '(' << xy.first << ',' << xy.second << ") ";
  }
  std::cout << std::endl;

  IntegerVector crossings(unit_cnt*unit_cnt);

  std::multimap<size_t,size_t>::const_iterator vert_cursor=intverts.cbegin();
  while (vert_cursor!=intverts.cend()) {
    size_t point_idx=vert_cursor->first;
    std::cout << "looking at "<<point_idx<<std::endl;
    std::set<size_t> street;
    std::set<size_t> arc;
    while (vert_cursor!=intverts.cend() && vert_cursor->first==point_idx) {
      size_t segment_idx=vert_cursor->second;
      if (segment_idx<dist_cnt) {
        std::cout << "segment "<<segment_idx <<" in arc"<<std::endl;
        arc.insert(segment_idx);
      } else {
        std::cout << "segment "<<segment_idx <<" in street"<<std::endl;
        street.insert(segment_idx);
      }
      ++vert_cursor;
    }
    if (arc.size()>0 && street.size()>0) {
      std::cout << "Adding crossings ";
      for (size_t crossed : arc) {
        const auto& arcref=segments[crossed];
        // Add to i,j and j,i.
        auto idxl=arcref.first*unit_cnt + arcref.second;
        auto idxr=arcref.second*unit_cnt + arcref.first;
        std::cout << "("<<idxl << " " << idxr << ") " << std::endl;
        crossings[idxl]+=1;
        crossings[idxr]+=1;
      }
      std::cout << std::endl;
    }
  }

  return crossings;
}


// [[Rcpp::export]]
SEXP simple_hazard(SEXP pairwise_distanceS, SEXP parametersS,
    Rcpp::Function callback) {
  afidd::LogInit("info");

  NumericVector pairwise(pairwise_distanceS);
  std::vector<double> distance(pairwise.begin(), pairwise.end());

  List parameters(parametersS);
  std::map<std::string, boost::any> params;
  params["individual_cnt"]=int64_t{as<int>(parameters["individual_cnt"])};
  params["initial"]=int64_t{as<int>(parameters["initial"])};
  params["N0"]=double{as<double>(parameters["N0"])};
  params["beta0"]=double{as<double>(parameters["beta0"])};
  params["beta1"]=double{as<double>(parameters["beta1"])};
  params["beta2"]=double{as<double>(parameters["beta2"])};
  params["gamma"]=double{as<double>(parameters["gamma"])};
  params["cutoff"]=double{as<double>(parameters["cutoff"])};
  params["growthrate"]=double{as<double>(parameters["growthrate"])};
  params["carrying"]=double{as<double>(parameters["carrying"])};
  int64_t rand_seed=int64_t{as<int>(parameters["seed"])};
  params["seed"]=rand_seed;

  auto observer=std::make_shared<CallbackEventObserver>(callback);

  RandGen rng(rand_seed);
  auto gspn=hsb::simple_hazard::SimpleHazardGSPN(params, distance, rng);

  bool write_keys=false;
  if (write_keys) {
    TKeyWriter w;
    gspn.WriteKeys(w);
  }

  int run_cnt=1;
  if (parameters.containsElementNamed("runs")) {
    run_cnt=as<int>(parameters["runs"]);
    BOOST_LOG_TRIVIAL(debug)<<"has run cnt";
  } else {
    BOOST_LOG_TRIVIAL(debug)<<"doesn't have cnt";
  }

  for (int run_idx=0; run_idx<run_cnt; ++run_idx) {
    hsb::simple_hazard::SIR_run(params, gspn, observer, rng);
  }
  return Rcpp::wrap(0);
}



// [[Rcpp::export]]
DataFrame bugs(SEXP pairwise_distanceS, SEXP parametersS) {
  afidd::LogInit("info");

  NumericVector pairwise(pairwise_distanceS);
  std::vector<double> distance(pairwise.begin(), pairwise.end());

  List parameters(parametersS);
  std::map<std::string, boost::any> params;
  for (auto iname : std::vector<std::string>{"individual_cnt", "seed",
      "initial_bug_cnt"}) {
    try {
      params[iname]=int64_t{as<int>(parameters[iname])};
    } catch (std::exception& e) {
      BOOST_LOG_TRIVIAL(error) << "Could not cast "<<iname<<
        " to an integer";
      throw;
    }
  }
  for (auto dname : std::vector<std::string>{"birth", "death", "carrying",
      "move0", "move1", "gamma", "cutoff"}) {
    try {
      params[dname]=double{as<double>(parameters[dname])};
    } catch (std::exception& e) {
      BOOST_LOG_TRIVIAL(error) << "Could not cast "<<dname<<" to a double";
      throw;
    }
  }

  auto observer=std::make_shared<SIRObserver>();

  int64_t rand_seed=boost::any_cast<int64_t>(params["seed"]);
  RandGen rng(rand_seed);
  hsb::bugs::SIR_run(params, distance, observer, rng);

  NumericVector when(observer->when_.begin(), observer->when_.end());
  IntegerVector what(observer->what_.begin(), observer->what_.end());
  IntegerVector who(observer->who_.begin(), observer->who_.end());
  IntegerVector who2(observer->who2_.begin(), observer->who2_.end());

  return DataFrame::create(Named("times")=when,
      Named("event")=what, Named("who")=who, Named("actor")=who2);
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
