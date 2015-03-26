#ifndef _TRAJECTORY_OBSERVER_H_
#define _TRAJECTORY_OBSERVER_H_ 1

namespace hsb {

class TrajectoryObserver {
 public:
  virtual void event(double when, int what, int64_t who, int64_t who2)=0;
};



} // namespace hsb

#endif // _TRAJECTORY_OBSERVER_H_
