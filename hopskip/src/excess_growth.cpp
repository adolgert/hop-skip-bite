#include "excess_growth.hpp"


double ExcessGrowthFunction(double y, void* params) {
  ExcessGrowthParams* egparams=static_cast<ExcessGrowthParams*>(params);
  return egparams->xa - log(y+1) - y/(y+1);
}
