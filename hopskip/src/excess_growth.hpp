#ifndef _EXCESS_GROWTH_H_
#define _EXCESS_GROWTH_H_ 1

#include <cmath>
#include <algorithm>
#include <memory>
#include "gsl/gsl_errno.h"
#include "gsl/gsl_math.h"
#include "gsl/gsl_roots.h"
#include "gspn_random.hpp"
#include "smv.hpp"

namespace hsb {

/*! Parameters for excess growth.
 */
struct ExcessGrowthParams {
  double xa;
};

double ExcessGrowthFunction(double y, void* params);


/*! Excess growth beyond logistic.
 *  The equation for logistic growth is
 *  \frac{ N_0 K e^{rt} }{ K + N_0 (e^{rt} - 1) }
 *  The excess is (N/K) * N. It assumes total growth rate
 *  is exponential, and logistic growth accounts for those who stay.
 *  This distribution is described in detail in accompanying
 *  documentation.
 *
 *  The scale parameter is a scaling of the hazard to account
 *  for loss of bugs.
 *
 *  There is a python file, growth.py, that demonstrates what
 *  this class does. It's allso in hopskip.tex.
 */
template<typename RNG>
class ExcessGrowth : public afidd::smv::TransitionDistribution<RNG> {
  std::unique_ptr<gsl_function> excess_function_;
  std::unique_ptr<ExcessGrowthParams> params_;
  const gsl_root_fsolver_type* solver_type_{gsl_root_fsolver_brent};
  std::unique_ptr<gsl_root_fsolver> solver_;
  double te_;
  double KP_; // K/N0
  double K_;  // K*scale
  double r_;
 public:
  ExcessGrowth(double N0, double K, double r, double scale, double te)
    : params_(new ExcessGrowthParams), te_(te),
      KP_(K/N0), K_(K*scale), r_(r) {
    assert(KP_>1); // Carrying capacity greater than initial seed.
    excess_function_.reset(new gsl_function());
    solver_.reset(gsl_root_fsolver_alloc(solver_type_));
    excess_function_->function=&ExcessGrowthFunction;
    excess_function_->params=params_.get();
  }
  virtual ~ExcessGrowth() {}

  /*! It calls the other version of sampling
   *  with a costly, but correct, transform.
   */
  virtual double Sample(double current_time, RNG& rng) const {
    return ImplicitHazardIntegral(-std::log(afidd::smv::uniform(rng)),
      current_time);
  }

  virtual double EnablingTime() const { return te_; }
  virtual bool BoundedHazard() const { return true; }

  virtual double HazardIntegral(double t0, double t1) const {
    assert(t0>=te_);
    assert(t1>=t0);
    double y0=y_of_t(t0-te_);
    double y1=y_of_t(t1-te_);
    return K_*( std::log((y1+1)/(y0+1)) - (y1-y0)/((y0+1)*(y1+1) ));
  }

  virtual double ImplicitHazardIntegral(double xa, double t0) const {
    if (t0<te_) {
      t0=te_;
    }
    if (xa/K_<1e-10) {
      // Avoid Brent when resolution is low.
      return this->SmallInverseHazardIntegral(xa, t0);
    }

    double absresolution=1e-9;
    double relresolution=1e-8;
    int iter_max=1000;
    assert(xa>0);
    double y0=y_of_t(t0-te_);
    assert(std::abs(te_+t_of_y(y0) - t0) < 1e-6);
    params_->xa=xa/K_ + std::log(y0+1) - y0/(y0+1);
    double low_bound=y0;
    // The 1.001 gives the high bound some space for numerical roundoff.
    double high_bound=1.001*(std::exp(params_->xa + 1) - 1);
    // auto f=[](double y, double xa)->double {
    //   ExcessGrowthParams p;
    //   p.xa=xa;
    //   return ExcessGrowthFunction(y, &p);
    // };
    // double lb=f(low_bound, params_->xa);
    // double hb=f(high_bound, params_->xa);
    // if ((lb<0 && hb<0) || (lb>0 && hb>0)) {
    //   BOOST_LOG_TRIVIAL(error)<<"roots on same side "<<lb<<" "<<hb;
    //   BOOST_LOG_TRIVIAL(error)<<"xa "<<xa<<" t0 "<< t0 << " y0 " << y0;
    //   BOOST_LOG_TRIVIAL(error)<<"te "<<te_<<" KP "<< KP_ << " K " << K_
    //     << " r " << r_;
    // }

    BOOST_LOG_TRIVIAL(debug)<<"xa "<<xa<<" t0 "<< t0 << " y0 " << y0;
    BOOST_LOG_TRIVIAL(debug)<<"te "<<te_<<" KP "<< KP_ << " K " << K_
      << " r " << r_;
    gsl_root_fsolver_set(solver_.get(), excess_function_.get(),
      low_bound, high_bound);

    int test_status=GSL_CONTINUE;
    int iter=0;
    double y1=0;
    while (test_status==GSL_CONTINUE && iter<iter_max) {
      int status=gsl_root_fsolver_iterate(solver_.get());
      switch (status) {
        case 0:
          break;
        case GSL_EBADFUNC:
          BOOST_LOG_TRIVIAL(error) << "Root solver bad function";
          break;
        case GSL_EZERODIV:
          BOOST_LOG_TRIVIAL(error) << "Root solver divide by zero";
          break;
        default:
          BOOST_LOG_TRIVIAL(error) << "Root solver error " << status;
        break;
      }
      assert(status==0);
      y1=gsl_root_fsolver_root(solver_.get());
      double t_low=gsl_root_fsolver_x_lower(solver_.get());
      double t_high=gsl_root_fsolver_x_upper(solver_.get());
      test_status=gsl_root_test_interval(t_low, t_high, absresolution,
        relresolution);
      ++iter;
    }
    if (test_status!=GSL_SUCCESS) {
      BOOST_LOG_TRIVIAL(error) << "test status " << test_status;
    }
    if (iter==iter_max) {
      BOOST_LOG_TRIVIAL(error)<< "Reached max iteration.";
    }
    double t_low=gsl_root_fsolver_x_lower(solver_.get());
    double t_high=gsl_root_fsolver_x_upper(solver_.get());
    BOOST_LOG_TRIVIAL(debug)<<"xa "<<xa<<" t0 "<< t0 << " y0 " << y0
      << " y1 " << y1;
    BOOST_LOG_TRIVIAL(debug)<<"iter "<<iter
    << " low "<<low_bound<<" high "<<high_bound;
    BOOST_LOG_TRIVIAL(debug)<< "tlow "<<t_low<<" thigh "<<t_high;
    return te_+this->t_of_y(y1);
  }


  double SmallInverseHazardIntegral(double xa, double t0) const {
    double y0=this->y_of_t(t0-te_);
    double dy=xa*std::pow(y0+1, 2)/(K_*y0);
    return te_+this->t_of_y(y0+dy);
  }

 private:
  inline double y_of_t(double t) const { return std::exp(r_*t)/(KP_-1); }
  inline double t_of_y(double y) const { return std::log(y*(KP_-1))/r_; }
};

int TestExcessGrowthDistribution();

} // namespace hsb

//_EXCESS_GROWTH_H_
#endif