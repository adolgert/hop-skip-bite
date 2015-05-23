#include <iostream>
#include <fstream>

#include <set>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <cassert>
#include "boost/math/special_functions/fpclassify.hpp"
#include "segment_intersect.hpp"
#include "smv.hpp"
//#include "boost/math/constants"

struct Point {
  double x;
  double y;
  Point()=default;
  Point(double x, double y) : x(x), y(y) {}
  Point(const std::pair<double,double>& xy) : x(xy.first), y(xy.second) {}
  operator std::pair<double,double>() { return std::make_pair(x, y); }
  bool operator<(const Point& b) const {
    if (y!=b.y) {
      return y>b.y;
    }
    return x<b.x;
  }
  bool operator>(const Point& b) const {
    if (y!=b.y) {
      return y<b.y;
    }
    return x>b.x;
  }
};

/*! Counterclockwise.
 *  ab is (a,b) cd is (c,d)  ef is (e,f)=
 *  (d-b)/(c-a)  < (h-f)/(g-e) so (d-b)(g-e) < (h-f)(c-a)
 */
bool CCW(const Point& ab, const Point& cd, const Point& ef) {
  return (ef.y - ab.y)*(cd.x-ab.x) > (cd.y-ab.y)*(ef.x-ab.x);
}

bool Intersect(const Point& a, const Point& b, const Point& c, const Point& d) {
  if (CCW(a, c, d) == CCW(b, c, d)) {
    return false;
  } else if (CCW(a, b, c) == CCW(a, b, d)) {
    return false;
  } // else
  return true;
}

// Intersection of line from p0 to p1 and line from p2 to p3.
// http://en.wikipedia.org/wiki/Line%E2%80%93line_intersection
Point IntersectionRaw(const Point& p0, const Point& p1, const Point& p2,
    const Point& p3) {
  return Point{
    ((p0.x*p1.y-p0.y*p1.x)*(p2.x-p3.x) - (p0.x-p1.x)*(p2.x*p3.y-p2.y*p3.x))/
    ((p0.x-p1.x)*(p2.y-p3.y) - (p0.y-p1.y)*(p2.x-p3.x)),
    ((p0.x*p1.y-p0.y*p1.x)*(p2.y-p3.y) - (p0.y-p1.y)*(p2.x*p3.y-p2.y*p3.x))/
    ((p0.x-p1.x)*(p2.y-p3.y) - (p0.y-p1.y)*(p2.x-p3.x))
  };
}

void cross_product(const double* u, const double* v, double* w) {
  w[0]=u[1]*v[2]-u[2]*v[1];
  w[1]=u[2]*v[0]-u[0]*v[2];
  w[2]=u[0]*v[1]-u[1]*v[0];
}

std::tuple<Point,bool>
Intersection(const Point& p0, const Point& p1, const Point& p2,
    const Point& p3) {
  double l0[3];
  double l1[3];
  double u[3];
  double v[3];
  u[0]=p0.x;
  u[1]=p0.y;
  u[2]=1;
  v[0]=p1.x;
  v[1]=p1.y;
  v[2]=1;
  cross_product(u, v, l0);
  u[0]=p2.x;
  u[1]=p2.y;
  u[2]=1;
  v[0]=p3.x;
  v[1]=p3.y;
  v[2]=1;
  cross_product(u, v, l1);
  double w[3];
  cross_product(l0, l1, w);
  if (w[2]!=0) {
    double x=w[0]/w[2];
    double y=w[1]/w[2];
    // Check for denormalization. Call that intersection-not-found.
    if (boost::math::isnormal<double>(x) && boost::math::isnormal<double>(y)) {
      return std::tuple<Point,bool>({x, y}, true);
    }
  }
  return std::tuple<Point,bool>({0, 0}, false);
}


using PointIdx=size_t;
using SegmentIdx=size_t;
enum class IntLoc { None, Up, Down, Middle };

struct Segment {
  PointIdx p0;
  PointIdx p1;
  Segment()=default;
  Segment(PointIdx a, PointIdx b) : p0(a), p1(b) {}
  Segment(const std::pair<PointIdx,PointIdx>& ab) : p0(ab.first), p1(ab.second) {}
  bool Contains(PointIdx v) const {
    return p0==v || p1==v;
  }
};

PointIdx UpperPoint(const Segment& s, const std::vector<Point>& points) {
  if (points.at(s.p0) < points.at(s.p1)) {
    return s.p0;
  }
  return s.p1;
}

PointIdx LowerPoint(const Segment& s, const std::vector<Point>& points) {
  if (points.at(s.p0) > points.at(s.p1)) {
    return s.p0;
  }
  return s.p1;
}

struct IntersectQueueSort {
  const std::vector<Point>& points;
  IntersectQueueSort(const std::vector<Point>& p) : points(p) {}
  bool operator()(PointIdx a, PointIdx b) const {
    if (a==b) return false;
    const Point& pa=points.at(a);
    const Point& pb=points.at(b);
    if (pa<pb) {
      return true;
    } else if (pb<pa) {
      return false;
    }
    return a<b;
  }
};


/*! Left-to-right ordering just below the sweep line.
 */
struct StatusLoc {
  double x;
  double radians; // first and second quadrants, 0 to pi, excluding 0.
  // Constructor for a general line.
  StatusLoc(const Point& p0, const Point& p1, const Point& p) {
    Point intpt;
    bool bIntersects;
    std::tie(intpt, bIntersects)=Intersection(p0, p1,
      Point(0, p.y), Point(1, p.y));
    if (bIntersects) {
      x=intpt.x;
      double dy=p0.y-p1.y;
      double dx=p0.x-p1.x;
      // atan2 returns radians between -pi and pi.
      radians=std::atan2(dy, dx);
      // Or equal to 0 so that horizontal is the last in the list.
      if (radians <= 0) {
        radians+=M_PI; // boost::math::constants::pi<double>();
      }
    } else {
      if (p0.x<p1.x) {
        x=p0.x;
      } else {
        x=p1.x;
      }
      radians=M_PI;
    }
  }
  // Constructor for a line this point contains.
  // This enables an exact comparison of floating point where it counts.
  StatusLoc(const Point& p0, const Point& p1, const Point& p, int which) {
    x=p.x;
    double dy=p0.y-p1.y;
    double dx=p0.x-p1.x;
    // atan2 returns radians between -pi and pi.
    radians=std::atan2(dy, dx);
    bool p0lower=(which==0) && (p0.y<p1.y);
    bool p1lower=(which==1) && (p1.y<p0.y);
    // Lower endpoints come before upper endpoints.
    if (p0lower || p1lower) {
      if (radians<0 && radians>-M_PI) {
        radians=-M_PI-radians;
      } else if (radians>0) {
        radians=-radians;
      } else {
        ; // leave -M_PI alone
      }
    } else {
      if (radians<=0) {
        radians+=M_PI;
      }
    }
  }
  bool operator<(const StatusLoc& slb) const {
    if (x!=slb.x) {
      return x < slb.x;
    } else {
      return radians < slb.radians;
    }
  }
  bool operator>(const StatusLoc& slb) const {
    if (x!=slb.x) {
      return x > slb.x;
    } else {
      return radians > slb.radians;
    }
  }
  friend inline
  std::ostream& operator<<(std::ostream& os, const StatusLoc& sl) {
    return os << '(' << sl.x << ", " << sl.radians
      << ')';
  }
};


struct QEntry {
  SegmentIdx seg;
  IntLoc where;
  friend bool operator<(const QEntry& a, const QEntry& b) {
    return a.seg<b.seg;
  }
};

struct TEntry {
  SegmentIdx seg;
  PointIdx pt;
  IntLoc where;
};

/*! This asks if one point is higher or lefter just below the sweep line.
 *  Instead of introducing a small value just below the current sweep line,
 *  we look at the location of the point at the sweep line and then look
 *  at the angle of the segment. If it has lower slope, it will be more
 *  to the left just below the sweep line. A zero slope is defined as
 *  coming last in the list.
 *
 *  This sort keeps lines in order as the sweep line lowers. A mistake
 *  here for stability of the sort would be a real problem for the
 *  underlying algorithm. We remove segments that would switch places,
 *  change the current_point in this class, and then reinsert those
 *  same segments.
 *
 *  Horizontal lines have to have a stable sort, too.
 */
struct StatusSort {
  const std::vector<Point>& points;
  const std::vector<Segment>& lines;
  const PointIdx& current_point;
  static const std::map<IntLoc,int> locval;
  StatusSort(const std::vector<Point>& points,
      const std::vector<Segment>& lines, const PointIdx& p) :
    points(points), lines(lines), current_point(p) {}

  //! Decide which entry is greater or lesser.
  /*! Two with the same upper point sort according to which
   *  would have smaller x just below the current y.
   *  Two with same lower point sort according to which
   *  has smaller x just above current y. Horizontal lines
   *  sort according to their current intersection and are always
   *  after others at the same point.
   *
   *  Two horizontal lines intersect only at start or end points.
   */
  bool operator()(const TEntry& a, const TEntry& b) const {
    double ax;
    double atheta;
    double bx;
    double btheta;
    if (a.pt==current_point) {
      ax=points[current_point].x;
    } else {
      ax=x_intercept_of(a);
    }
    if (b.pt==current_point) {
      bx=points[current_point].x;
    } else {
      bx=x_intercept_of(b);
    }
    if (ax!=bx) {
      return ax<bx;
    }
    atheta=angle_of(a);
    btheta=angle_of(b);
    if (std::abs(atheta)>M_PI || std::abs(btheta)>M_PI) {
      throw std::runtime_error("Angle outside of +=pi.");
    }
    if (atheta!=btheta) {
      return atheta<btheta;
    }
    return locval.at(a.where)<locval.at(b.where);
  }

  double angle_of(const TEntry& t) const {
    const Point& p0=points[lines[t.seg].p0];
    const Point& p1=points[lines[t.seg].p1];
    double dy=p0.y-p1.y;
    double dx=p0.x-p1.x;
    // atan2 returns values from -pi to pi.
    double angle=std::atan2(dy, dx);
    if (t.where==IntLoc::Down) {
      // Lower points should return from -pi to 0.
      // -pi -> 0, -pi/2->-pi/2, 0->0
      if (angle<0 && angle>-M_PI) {
        angle=-M_PI-angle;
      } else if (angle>0) {
        // pi/2 -> -pi/2, pi->-pi
        angle=-angle;
      } else {
        ; // leave -M_PI alone. Horizontal end comes first.
      }
    } else {
      // Middle points should return from 0 to pi.
      if (angle<0) {
        angle+=M_PI;
      } else {
        ; // angle=angle
      }
    }
    return angle;
  }

  double x_intercept_of(const TEntry& t) const {
    const Point& p0=points[lines[t.seg].p0];
    const Point& p1=points[lines[t.seg].p1];
    const Point& p=points[current_point];
    bool bIntersects;
    Point intpt;
    std::tie(intpt, bIntersects)=Intersection(p0, p1,
      Point(0, p.y), Point(1, p.y));
    if (!bIntersects) {
      if (p0.x<p1.x) {
        return p0.x;
      } else {
        return p1.x;
      }
    }
    return intpt.x;
  }
};

const std::map<IntLoc,int> StatusSort::locval={{IntLoc::Up, 1},
  {IntLoc::Middle, 2}, {IntLoc::Down, 3}, {IntLoc::None, 4}};

void CheckStatusSort() {
  std::vector<Point> points;
  std::vector<Segment> lines;
  // the cases to test are:
  // not-horizontal, horizontal.
  // upper, middle, lower points.
  for (int i=0; i<4; ++i) {
    for (int j=0; j<4; ++j) {
      double di=i;
      double dj=j;
      points.push_back({di, dj});
    }
  }
  PointIdx curpt=0*4+0;
  StatusSort sorter(points, lines, curpt);
  lines.push_back({0*4+1, 1*4+0});
  lines.push_back({1*4+1, 1*4+0});
  if(sorter({0, curpt, IntLoc::Down}, {1, curpt, IntLoc::Down})!=true) {
    BOOST_LOG_TRIVIAL(error)<<"Lines 0 and 1 wrong order.";
  }
  lines.push_back({0*4+0, 1*4+1});
  lines.push_back({2*4+0, 1*4+1});
  curpt=0*4+1;
  if (sorter({2, curpt, IntLoc::Up}, {3, curpt, IntLoc::Up})!=true) {
    BOOST_LOG_TRIVIAL(error)<<"Lines 2-3 wrong order.";
  }
  lines.push_back({1*4+1, 2*4+1});
  if (sorter({2, curpt, IntLoc::Up}, {4, curpt, IntLoc::Up})!=true) {
    BOOST_LOG_TRIVIAL(error)<<"Lines 2-4 wrong order.";
  }
}


void FindNewEvent(SegmentIdx sl_idx, SegmentIdx sr_idx, PointIdx p,
    const std::vector<Segment>& lines, std::vector<Point>& points,
    std::multimap<PointIdx,QEntry,IntersectQueueSort>& Q,
    std::multimap<PointIdx,SegmentIdx>& intersections,
    std::vector<std::pair<double,double>>& intersection_points) {
  const Segment& sl=lines[sl_idx];
  const Segment& sr=lines[sr_idx];
  // If they share an endpoint, that endpoint is in Q already.
  std::set<SegmentIdx> same_points;
  same_points.insert(sl.p0);
  same_points.insert(sl.p1);
  same_points.insert(sr.p0);
  same_points.insert(sr.p1);
  if (same_points.size()<4) {
    std::cout << "FindNewEvent: Shared points" << std::endl;
    return;
  }
  if (Intersect(points[sl.p0], points[sl.p1], points[sr.p0], points[sr.p1])) {
    Point intersection;
    bool int_found;
    std::tie(intersection, int_found)=Intersection(
      points[sl.p0], points[sl.p1], points[sr.p0], points[sr.p1]);
    // If intersection is below or to the right of p.
    if (int_found && intersection > points.at(p)) {
      // If the intersection is not yet present in event queue.
      // Stay away from comparing doubles. Just look for any vertex which
      // contains an intersection of both segments.

      points.push_back(intersection);
      size_t added_idx=points.size()-1;
      auto sliter=Q.emplace(added_idx, QEntry{sl_idx, IntLoc::Middle});
      auto sriter=Q.emplace(added_idx, QEntry{sr_idx, IntLoc::Middle});

      bool bFound=false;
      auto slsearch=sliter;
      bool searching_left=true;
      while (searching_left && !bFound &&
          std::abs(points[slsearch->first].y-intersection.y)<0.001) {
        PointIdx qpt=slsearch->first;
        bool saw_l=false;
        bool saw_r=false;
        while (searching_left && slsearch->first==qpt) {
          if (slsearch->second.seg==sl_idx) {
            saw_l=true;
          }
          if (slsearch->second.seg==sr_idx) {
            saw_r=true;
          }
          if (saw_l && saw_r) {
            searching_left=false;
          }
          if (slsearch==Q.begin()) {
            searching_left=false;
          } else {
            --slsearch;
          }
        }
        if (saw_l && saw_r) {
          bFound=true;
        }
      }
      auto srsearch=sriter;
      while (srsearch!=Q.end() && !bFound &&
          std::abs(points[srsearch->first].y-intersection.y)<0.001) {
        PointIdx qpt=srsearch->first;
        bool saw_l=false;
        bool saw_r=false;
        while (srsearch!=Q.end() && srsearch->first==qpt) {
          if (srsearch->second.seg==sl_idx) {
            saw_l=true;
          }
          if (srsearch->second.seg==sr_idx) {
            saw_r=true;
          }
          ++srsearch;
        }
        if (saw_l && saw_r) {
          bFound=true;
        }
      }
      if (bFound) {
        BOOST_LOG_TRIVIAL(debug)<<"Point already in Q";
        Q.erase(sliter);
        Q.erase(sriter);
      }
    } else {
      std::cout << "Segment intersection above or left of p" << std::endl;
    }
  } else {
    std::cout << "Segments cannot intersect by CCW. "
      <<sl_idx <<" " << sr_idx << std::endl;
  }
}


void graphical(const std::vector<Point>& points,
    const std::vector<Segment>& lines,
    const std::map<StatusLoc,std::pair<PointIdx,SegmentIdx>>& T,
    PointIdx p, std::set<SegmentIdx>& L, std::set<SegmentIdx>& U,
    std::set<SegmentIdx>& C,
    const std::vector<std::pair<double,double>>& intpts) {
  std::ofstream slines("lines.txt", std::fstream::out | std::fstream::trunc);
  for (auto l : lines) {
    slines << l.p0 << " " << l.p1 << std::endl;
  }
  std::ofstream spoints("points.txt", std::fstream::out | std::fstream::trunc);
  for (auto pts : points) {
    spoints << pts.x << " " << pts.y << std::endl;
  }
  std::ofstream sipoints("intpoints.txt", std::fstream::out | std::fstream::trunc);
  for (auto ipts : intpts) {
    sipoints << ipts.first << " " << ipts.second << std::endl;
  }
  std::ofstream scurpoint("current.txt", std::fstream::out | std::fstream::app);
  scurpoint << p << std::endl;

  std::ofstream ststatus("tstatus.txt", std::fstream::out | std::fstream::app);
  ststatus << T.size();
  for (auto t : T) {
    ststatus << " " << t.first.x;
  }
  ststatus << std::endl;

  std::ofstream sl("l.txt", std::fstream::out | std::fstream::app);
  sl << L.size();
  for (auto l : L) {
    sl << " " << l;
  }
  sl << std::endl;
  std::ofstream su("u.txt", std::fstream::out | std::fstream::app);
  su << U.size();
  for (auto u : U) {
    su << " " << u;
  }
  su << std::endl;
  std::ofstream sc("c.txt", std::fstream::out | std::fstream::app);
  sc << C.size();
  for (auto c : C) {
    sc << " " << c;
  }
  sc << std::endl;
}

bool bSegDebug=true;

void CheckStatus(const std::set<TEntry,StatusSort>& T, const StatusSort& ss) {
  // Same segment shouldn't be in twice.
  std::set<SegmentIdx> segs;
  std::stringstream segstr;
  for (const auto& t : T) {
    segstr << t.seg;
    if (t.where==IntLoc::Up) segstr << "u ";
    else if (t.where==IntLoc::Middle) segstr << "c ";
    else if (t.where==IntLoc::Down) segstr << "l ";
    else segstr << "x ";
    auto res=segs.insert(t.seg);
    if (!res.second) {
      BOOST_LOG_TRIVIAL(error) << "Found segment twice " << t.seg;
    } else {
    }
  }
  BOOST_LOG_TRIVIAL(debug)<< "T "<< segstr.str();
  // Check that they are ordered according to status sort.
  auto statiter=T.begin();
  auto nextiter=statiter;
  bool ordered=true;
  if (nextiter!=T.end()) {
    ++nextiter;
    while (nextiter!=T.end()) {
      if (ss(*statiter, *nextiter)==false) {
        BOOST_LOG_TRIVIAL(error)<<"Segments out of order "
        << statiter->seg << "-" << nextiter->seg;
        ordered=false;
      } else {
        // BOOST_LOG_TRIVIAL(debug)<<"Segments in order "
        // << statiter->seg << "-" << nextiter->seg;
      }
      ++statiter;
      ++nextiter;
    }
  }
  if (!ordered) {
    throw std::runtime_error("Segments out of order");
  }
}

/*! Find intersections among line segments.
 *  Reading "Computational Geometry" by Berg et al., Chapter 2.
 */
std::tuple<std::vector<std::pair<double,double>>,
  std::multimap<PointIdx,SegmentIdx>>
segment_intersections_sweep(const std::vector<std::pair<double,double>>& inpoints,
    const std::vector<std::pair<size_t,size_t>>& inlines) {
  // Make a copy so that we can add points for intersections to the list.
  std::vector<Point> points(inpoints.begin(), inpoints.end());
  std::vector<Segment> lines(inlines.begin(), inlines.end());
  std::cout << "Copied inputs " << inpoints.size() << ' ' << inlines.size()
    << std::endl;

  std::cout << "Seen points ";
  size_t sh_pts=0;
  for (auto pp : points) {
    std::cout << sh_pts << "-(" << pp.x << ',' << pp.y << ") ";
    ++sh_pts;
  }
  std::cout << std::endl;
  std::cout << "Seen lines ";
  size_t sh_lines=0;
  for (auto lls : lines) {
    std::cout << sh_lines<<"-(" << lls.p0 << ',' << lls.p1 << ") ";
    ++sh_lines;
  }
  std::cout << std::endl;

  // Event queue starts with all segment endpoints but will
  // also include future intersections. It is sorted so that
  // largest y come first, and, for same y, smallest x come first.
  std::multimap<PointIdx,QEntry,IntersectQueueSort> Q(
      (IntersectQueueSort(points))); // Event queue.
  // Status of the algorithm.
  using StatusType=std::set<TEntry,StatusSort>;
  PointIdx sweep_p;
  StatusSort T_sort(points, lines, sweep_p);
  StatusType T(T_sort);
  // Outputs return to same datatypes passed in.
  std::multimap<PointIdx,SegmentIdx> intersections;
  std::vector<std::pair<double,double>> intersection_points;

  // Insert all segment endpoints into event queue.
  std::cout << "Insert all segment endpoints into queue" << std::endl;
  for (size_t ins_idx=0; ins_idx<lines.size(); ++ins_idx) {
    const auto& line=lines[ins_idx];
    if (points[line.p0] < points[line.p1]) {
      Q.emplace(line.p0, QEntry{ins_idx, IntLoc::Up});
      Q.emplace(line.p1, QEntry{ins_idx, IntLoc::Down});
    } else {
      Q.emplace(line.p0, QEntry{ins_idx, IntLoc::Down});
      Q.emplace(line.p1, QEntry{ins_idx, IntLoc::Up});
    }
  }
  std::cout <<"Starting Q size " << Q.size() << std::endl;

  auto next_event=Q.begin();
  while (next_event!=Q.end()) {
    if (bSegDebug) {
      CheckStatusSort();
      // Check invariant on Q: same pair insects only once in Q.
      // All intersections don't repeat a segment.
      std::set<std::pair<SegmentIdx,SegmentIdx>> segsect;
      auto qiter=Q.begin();
      while (qiter!=Q.end()) {
        auto qp=qiter->first;
        std::set<SegmentIdx> segs;
        while (qp==qiter->first) {
          auto seg=qiter->second.seg;
          if (segs.find(seg)!=segs.end()) {
            BOOST_LOG_TRIVIAL(error)<<"Same segment listed twice for one point"
              << " point " << qp << " seg " << seg;
          }
          segs.insert(seg);
          ++qiter;
        }
        auto segstart=segs.begin();
        while (segstart!=segs.end()) {
          auto segfinish=segstart;
          ++segfinish;
          while (segfinish!=segs.end()) {
            SegmentIdx l=std::min(*segstart, *segfinish);
            SegmentIdx m=std::max(*segstart, *segfinish);
            auto segpair=std::make_pair(l,m);
            auto res=segsect.insert(segpair);
            if (res.second==false) {
              BOOST_LOG_TRIVIAL(error) << "Same intersection twice for segs "
              << l << " and " << m;
            }
            ++segfinish;
          }
          ++segstart;
        }
      }
    }

    PointIdx p=next_event->first;
    using SegSet=std::set<SegmentIdx>;
    SegSet L;
    SegSet C;
    SegSet U;
    std::set<QEntry> UC;

    auto r=Q.equal_range(p);
    std::cout << "Line segments associated with " << p << ": ";
    for (auto ri=r.first; ri!=r.second; ++ri) {
      const QEntry& entry=ri->second;
      std::cout << entry.seg << ' ';
      if (entry.where==IntLoc::Up) {
        U.insert(entry.seg);
        UC.insert(entry);
      } else if (entry.where==IntLoc::Down) {
        L.insert(entry.seg);
      } else if (entry.where==IntLoc::Middle) {
        C.insert(entry.seg);
        UC.insert(entry);
      } else {
        assert(false);
        BOOST_LOG_TRIVIAL(error)<<"Segment lacks where classification.";
      }
    }
    std::cout << std::endl;
    Q.erase(r.first, r.second);
    std::cout << "next event p=" << p << " segs=";
    std::cout << std::endl << "queue length " << Q.size() << std::endl;

    std::set<SegmentIdx> LUC(L.begin(), L.end());
    LUC.insert(U.begin(), U.end());
    LUC.insert(C.begin(), C.end());
    std::set<SegmentIdx> LC(L.begin(), L.end());
    LC.insert(C.begin(), C.end());

    // 3-4. If L, U, C contain more than one segment. Report as intersection.
    // Insert them in a set in case L and U point to the same segment.
    if (LUC.size()>1) {
      PointIdx shared_vertex=intersection_points.size();
      intersection_points.push_back(std::make_pair(points[p].x, points[p].y));
      std::cout << "LUC intersection " << shared_vertex << " (" << points[p].x
        << ", " << points[p].y << ") " << LUC.size() << std::endl;
      for (auto ls : LUC) {
        intersections.emplace(shared_vertex, ls);
      }
    } // else no intersection to add.


    // Any midpoint intersection will change its ordering in T
    // the moment we set the StatusSort to the new current p.
    // So we have to remove those and then set the new current point.
    // We insert the segment, look nearby for it, then delete both.
    // The trick to using the sorted set is to insert a dummy point
    // and then search nearby.
    if (LC.size()>0) {
      BOOST_LOG_TRIVIAL(debug)<<"LC.size() "<<LC.size();
      TEntry tdel{*LC.begin(), sweep_p, IntLoc::None};
      auto tdel_ins=T.insert(tdel);
      auto tdel_save=tdel_ins.first;
      BOOST_LOG_TRIVIAL(debug)<<"LC.size() inserted "<<tdel_ins.second;
      auto low_iter=tdel_save;
      auto high_iter=tdel_save;
      ++high_iter;
      int search_cnt=0;
      std::vector<StatusType::iterator> to_erase;
      bool bLeftRun=true;
      bool bRightRun=true;
      while (to_erase.size()<LC.size() && (bLeftRun || bRightRun)) {
        BOOST_LOG_TRIVIAL(debug)<<"LC.size() Low iter "<<low_iter->seg;
        if (low_iter!=T.begin()) {
          --low_iter;
          auto low_found=LC.find(low_iter->seg);
          if (low_found!=LC.end()) {
            to_erase.push_back(low_iter);
          }
        } else {
          bLeftRun=false;
        }
        if (high_iter!=T.end()) {
          BOOST_LOG_TRIVIAL(debug)<<"LC.size() High iter "<<high_iter->seg;
          auto high_found=LC.find(high_iter->seg);
          if (high_found!=LC.end()) {
            to_erase.push_back(high_iter);
          }
          ++high_iter;
        } else {
          bRightRun=false;
        }
        ++search_cnt;
      }
      std::stringstream lcstr;
      for (auto lprint : LC) {
        lcstr << lprint << ' ';
      }
      if (to_erase.size()<LC.size()) {
        BOOST_LOG_TRIVIAL(error)<<"Didn't find all segs in LC: " << lcstr.str();
        CheckStatus(T, T_sort);
        throw std::runtime_error("Could not find segments to erase");
      }
      BOOST_LOG_TRIVIAL(debug)<<"LC.size() erasing "<<lcstr.str();
      for (auto erase_it : to_erase) {
        T.erase(erase_it);
      }
      BOOST_LOG_TRIVIAL(debug)<<"LC.size() erasing last";
      T.erase(tdel_save);
      BOOST_LOG_TRIVIAL(debug)<<"LC.size() search took steps "<<search_cnt
        <<" for "<<LC.size()<<" entries.";
    }

    // Now we can set the current point for the sort with impunity.
    sweep_p=p;

    //graphical(points, lines, T, p, L, U, C, intersection_points);

    // 7. if U+C==0
    std::cout << "UC size " << UC.size() << std::endl;
    if (UC.size()==0) {
      if (bSegDebug) CheckStatus(T, T_sort);
      std::cout << "size(UC)=0" << std::endl;
      // 8. let sl and sr be the left and right neighbors of p.
      auto Qlp=L.begin();
      auto probe_ins=T.insert(TEntry{*Qlp, p, IntLoc::None});
      auto probe_iter=probe_ins.first;
      if (probe_iter!=T.begin()) {
        auto sl=probe_iter;
        --sl;
        auto sr=probe_iter;
        ++sr;
        if (sr!=T.end()) {
          std::cout << "8. find event " << sl->seg << "-"
            << sr->seg << std::endl;
          FindNewEvent(sl->seg, sr->seg, p,
            lines, points, Q, intersections, intersection_points);
        } else {
          std::cout << "8. " << p << " has no right neighbor. no search" << std::endl;
        }
      } else {
        std::cout <<"8. " << p << " has no left neighbor. no search" << std::endl;
      }
      T.erase(probe_iter);
    } else {
      std::cout << "size(UC) " << UC.size() << std::endl;
      StatusType::iterator least_u;
      StatusType::iterator greatest_u;
      bool first=true;
      for (const QEntry& qe_ins : UC) {
        auto insiter=T.insert({qe_ins.seg, p, qe_ins.where});
        BOOST_LOG_TRIVIAL(debug)<<"Insert seg "<<qe_ins.seg;
        if (first) {
          least_u=insiter.first;
          greatest_u=insiter.first;
          first=false;
        } else {
          if (T_sort(*insiter.first, *least_u)) {
            least_u=insiter.first;
          }
          if (T_sort(*greatest_u, *insiter.first)) {
            greatest_u=insiter.first;
          }
        }
      }
      if (bSegDebug) CheckStatus(T, T_sort);
      std::cout << "\tLeast u "<<least_u->seg<<'-'<<least_u->pt
        <<" greatest u " <<greatest_u->seg <<'-'
        <<greatest_u->pt << std::endl;
      if (least_u!=T.begin()) {
        auto sl=least_u;
        --sl;
        std::cout << "8. Find event " << sl->seg << "-" <<
          least_u->seg << std::endl;
        FindNewEvent(sl->seg, least_u->seg, p, lines,
          points, Q, intersections, intersection_points);
      } else {
        std::cout << "8. Least_u is at beginning. No search." << std::endl;
      }
      auto sr=greatest_u;
      ++sr;
      if (sr!=T.end()) {
        FindNewEvent(greatest_u->seg, sr->seg, p, lines,
          points, Q, intersections, intersection_points);
      } else {
        std::cout << "8. Greatest_u is at end. No search." << std::endl;
      }
    }

    next_event=Q.begin();
  }
  return std::make_tuple(intersection_points, intersections);
}


std::tuple<std::vector<std::pair<double,double>>,
  std::multimap<PointIdx,SegmentIdx>>
segment_intersections(const std::vector<std::pair<double,double>>& inpoints,
    const std::vector<std::pair<size_t,size_t>>& inlines) {
  std::multimap<PointIdx,SegmentIdx> intersections;
  std::vector<std::pair<double,double>> intersection_points;
  // Make a copy so that we can add points for intersections to the list.
  for (size_t s0idx=0; s0idx<inlines.size()-1; ++s0idx) {
    size_t p0_idx=inlines[s0idx].first;
    size_t p1_idx=inlines[s0idx].second;
    Point p0(inpoints[p0_idx]);
    Point p1(inpoints[p1_idx]);
    for (size_t s1idx=s0idx+1; s1idx<inlines.size(); ++s1idx) {
      size_t p2_idx=inlines[s1idx].first;
      size_t p3_idx=inlines[s1idx].second;
      Point p2(inpoints[p2_idx]);
      Point p3(inpoints[p3_idx]);
      if (p0_idx!=p2_idx && p0_idx!=p3_idx && p1_idx!=p2_idx && p1_idx!=p3_idx){
        if (Intersect(p0, p1, p2, p3)) {
          Point intersection_pt;
          bool found;
          std::tie(intersection_pt, found)=Intersection(p0, p1, p2, p3);
          if (found) {
            intersection_points.push_back(std::make_pair(
              intersection_pt.x, intersection_pt.y));
            intersections.emplace(intersection_points.size()-1, s0idx);
            intersections.emplace(intersection_points.size()-1, s1idx);
          } else {
            std::cout << "Couldn't find intersection" << std::endl;
          }
        }
      }
    }
  }
  return std::make_tuple(intersection_points, intersections);
}
