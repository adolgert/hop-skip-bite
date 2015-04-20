#include "excess_growth.hpp"
#include "rng.hpp"

namespace hsb {


double ExcessGrowthFunction(double y, void* params) {
  ExcessGrowthParams* egparams=static_cast<ExcessGrowthParams*>(params);
  return egparams->xa - std::log(y+1) - y/(y+1);
}


int TestExcessGrowthDistribution() {
  auto rng=RandGen(33333);
  int error_cnt=0;
  double eps=1e-7;
  for (int tidx=0; tidx<3; ++tidx) {
    double te=0.7*tidx;
    for (int nidx=0; nidx<3; ++nidx) {
      double N0=nidx;
      for (int kidx=0; kidx<3; ++kidx) {
        double K=kidx*100;
        for (int ridx=0; ridx<3; ++ridx) {
          double r=std::pow(10, -1-ridx);
          ExcessGrowth<RandGen> dist(N0, K, r, te);

          for (int t0idx=0; t0idx<5; ++t0idx) {
            double t0=te+t0idx*0.1;
            for (int t1idx=0; t1idx<5; ++t1idx) {
              double t1=t0+0.1*t1idx;

              double xa=dist.HazardIntegral(t0, t1);
              double t1p=dist.ImplicitHazardIntegral(xa, t0);
              if (std::abs(t1p-t1)>eps) {
                ++error_cnt;
                BOOST_LOG_TRIVIAL(error) << "mismatch " <<
                  N0 << " " << K << " " << r << " " << te << " " <<
                  t0 << " " << t1 << " " << xa << " " << t1p << std::endl;
              }
            }
          }
        }
      }
    }
  }
  return error_cnt;
}

} // namespace hsb