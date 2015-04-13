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
 */
template<typename RNG>
class ExcessGrowth : public afidd::smv::TransitionDistribution<RNG> {
  std::unique_ptr<gsl_function> excess_function_;
  std::unique_ptr<ExcessGrowthParams> params_;
  const gsl_root_fsolver_type* solver_type_;
  std::unique_ptr<gsl_root_fsolver> solver_;
  double te_;
  double N0_;
  double K_;
  double r_;
 public:
  ExcessGrowth(double N0, double K, double r, double te)
    : params_(new ExcessGrowthParams), te_(te),
      N0_(N0), K_(K), r_(r),
      solver_type(gsl_root_fsolver_brent),
      excess_function_(new gsl_function()),
      solver_(gsl_root_fsolver_alloc(solver_type_)) {
    excess_function_->function=&ExcessGrowthHazard;
    excess_function_->params=params_.get();
  }
  virtual ~ExcessGrowth() {}

  /*! Don't intend to use this for production, so it calls
   *  the other version with a costly, but correct, transform.
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
    return (K_/r_)*( log((y1+1)/(y0+1)) - (y1-y0)/((y0+1)*(y1+1) )
  }

  virtual double ImplicitHazardIntegral(double xa, double t0) const {
    double resolution=1e9;
    double bounds_err=1e-13;
    int iter_max=1000;
    double y0=y_of_t(t0-te_);
    params_->xa=xa*r_/K_ + std::log(y0+1) - y0/(y0+1);
    double low_bound=y0;
    double high_bound=std::exp(params_->xa + 1) - 1;
  }

 private:
  double y_of_t(double t) { return (N0_/K_)*(std::exp(r*t)-1); }
  double t_of_y(double y) { return std::log(1+y*K_/N0_)/r_; }
};





//_EXCESS_GROWTH_H_
#endif