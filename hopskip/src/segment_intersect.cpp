#include <iostream>
#include <fstream>

#include <set>
#include <algorithm>
#include <cmath>
#include <cassert>
#include "boost/math/special_functions/fpclassify.hpp"
#include "segment_intersect.hpp"
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
  if (w!=0) {
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
    return points.at(a) < points.at(b);
  }
};


struct StatusLoc {
  double x;
  double y;
  double radians; // first and second quadrants, 0 to pi, excluding 0.
  bool operator<(const StatusLoc& slb) const {
    if (y!=slb.y) {
      return y > slb.y;
    } else {
      if (x!=slb.x) {
        return x < slb.x;
      } else {
        return radians < slb.radians;
      }
    }
  }
  bool operator>(const StatusLoc& slb) const {
    if (y!=slb.y) {
      return y < slb.y;
    } else {
      if (x!=slb.x) {
        return x > slb.x;
      } else {
        return radians > slb.radians;
      }
    }
  }
  friend inline
  std::ostream& operator<<(std::ostream& os, const StatusLoc& sl) {
    return os << '(' << sl.x << ", " << sl.y << ", " << sl.radians
      << ')';
  }
};


/*! This asks if one point is higher or lefter just below the sweep line.
 *  Instead of introducing a small value just below the current sweep line,
 *  we look at the location of the point at the sweep line and then look
 *  at the angle of the segment. If it has lower slope, it will be more
 *  to the left just below the sweep line. A zero slope is defined as
 *  coming last in the list.
 */
struct StatusSort {
  const std::vector<Point>& points;
  const std::vector<Segment>& lines;
  StatusSort(const std::vector<Point>& points,
      const std::vector<Segment>& lines) :
    points(points), lines(lines) {}

  bool operator()(const StatusLoc& a, const StatusLoc& b) {
    bool higher_lefter= a.y>b.y;
    if (a.y==b.y) {
      higher_lefter=a.radians<b.radians;
    }
    return higher_lefter;
  }
};

void FindNewEvent(SegmentIdx sl_idx, SegmentIdx sr_idx, PointIdx p,
    const std::vector<Segment>& lines, std::vector<Point>& points,
    std::multimap<PointIdx,SegmentIdx,IntersectQueueSort> Q,
    std::multimap<PointIdx,SegmentIdx>& intersections,
    std::vector<std::pair<double,double>>& intersection_points) {
  const Segment& sl=lines[sl_idx];
  const Segment& sr=lines[sr_idx];
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
      bool already_entered=false;
      auto q_cursor=Q.cbegin();
      while (q_cursor!=Q.cend()) {
        auto found_one=std::find_if(q_cursor, Q.cend(),
          [sl_idx, sr_idx](const std::pair<PointIdx,SegmentIdx>& entry)->bool {
            return entry.second==sl_idx || entry.second==sr_idx;
          });
        if (found_one!=Q.cend()) {
          PointIdx v=found_one->first;
          ++found_one;
          while (found_one!=Q.cend() && found_one->first==v) {
            if (found_one->second==sl_idx || found_one->second==sr_idx) {
              already_entered=true;
            }
            ++found_one;
          }
        }
        q_cursor=found_one;
      }
      if (!already_entered) {
        std::cout << "Adding intersection" << std::endl;
        points.push_back(intersection);
        size_t added_idx=points.size()-1;
        Q.emplace(added_idx, sl_idx);
        Q.emplace(added_idx, sr_idx);
        intersection_points.emplace_back(intersection.x, intersection.y);
        size_t iadded_idx=intersection_points.size()-1;
        intersections.emplace(iadded_idx, sl_idx);
        intersections.emplace(iadded_idx, sr_idx);
      } else {
        std::cout << "Intersection already added" << std::endl;
      }
    } else {
      std::cout << "Segment intersection above or left of p" << std::endl;
    }
  } else {
    std::cout << "Segments cannot intersect by CCW." << std::endl;
  }
}

/*! Find intersections among line segments.
 *  Reading "Computational Geometry" by Berg et al., Chapter 2.
 */
std::tuple<std::vector<std::pair<double,double>>,
  std::multimap<PointIdx,SegmentIdx>>
segment_intersections(const std::vector<std::pair<double,double>>& inpoints,
    const std::vector<std::pair<size_t,size_t>>& inlines) {
  // Make a copy so that we can add points for intersections to the list.
  std::vector<Point> points(inpoints.begin(), inpoints.end());
  std::vector<Segment> lines(inlines.begin(), inlines.end());
  std::cout << "Copied inputs " << inpoints.size() << ' ' << inlines.size()
    << std::endl;

  // Event queue starts with all segment endpoints but will
  // also include future intersections. It is sorted so that
  // largest y come first, and, for same y, smallest x come first.
  std::multimap<PointIdx,SegmentIdx,IntersectQueueSort> Q(
      (IntersectQueueSort(points))); // Event queue.
  // Status of the algorithm.
  using StatusType=std::map<StatusLoc,SegmentIdx>;
  StatusType T;
  // Outputs return to same datatypes passed in.
  std::multimap<PointIdx,SegmentIdx> intersections;
  std::vector<std::pair<double,double>> intersection_points;

  // Insert all segment endpoints into event queue.
  std::cout << "Insert all segment endpoints into queue" << std::endl;
  for (size_t ins_idx=0; ins_idx<lines.size(); ++ins_idx) {
    const auto& line=lines[ins_idx];
    Q.emplace(line.p0, ins_idx);
    Q.emplace(line.p1, ins_idx);
  }
  std::cout <<"Starting Q size " << Q.size() << std::endl;

  auto next_event=Q.begin();
  while (next_event!=Q.end()) {
    auto p=next_event->first;
    std::set<SegmentIdx> p_segments;
    auto r=Q.equal_range(p);
    std::cout << "Line segments associated with " << p << std::endl;
    for (auto ri=r.first; ri!=r.second; ++ri) {
      std::cout << '\t' << ri->second << std::endl;
      p_segments.insert(ri->second);
    }
    Q.erase(r.first, r.second);
    std::cout << "next event p=" << p << " segs=";
    for (auto ps : p_segments) {
      std::cout << ps << ", ";
    }
    std::cout << "queue length " << Q.size() << std::endl;

    // 1. Let U(p) be all segments whose upper point is p.
    std::cout << "U(p) have p=" << p <<" as upper point" << std::endl;
    std::set<SegmentIdx> U;
    for (SegmentIdx si : p_segments) {
      const auto& segment=lines[si];
      std::cout << "In U? ("<<segment.p0<<", "<<segment.p1<<")"<<std::endl;
      if (UpperPoint(segment, points) == p) {
        std::cout << "Insert " << si << " into U(p)" << std::endl;
        U.insert(si);
      }
    }

    // 2. Find all line segments in T that contain p.
    // L are segments whose lower point is p. 
    // Store the iterator, too, for step 5.
    using LowerSet=std::set<SegmentIdx>;
    LowerSet L;
    StatusType::iterator Titer=std::find_if(T.begin(), T.end(),
      [p, &lines](const std::pair<StatusLoc,SegmentIdx>& tentry)->bool {
        const auto& l=lines.at(tentry.second);
        return l.Contains(p);
      });
    // These segments are guaranteed to be adjacent in T.
    while (Titer!=T.end() && lines.at(Titer->second).Contains(p)) {
      if (LowerPoint(lines.at(Titer->second), points)==p) {
        std::cout << "Insert " << Titer->second << " into L(p)" << std::endl;
        L.emplace(Titer->second);
      }
      ++Titer;
    }
    // C(p) are segments that contain p in their interior. No.

    // 3-4. If L, U, C contain more than one segment. Report as intersection.
    // Insert them in a set in case L and U point to the same segment.
    std::set<SegmentIdx> LU(U.begin(), U.end());
    LU.insert(L.begin(), L.end());
    if (LU.size()>1) {
      PointIdx shared_vertex=intersection_points.size();
      intersection_points.push_back(std::make_pair(points[p].x, points[p].y));
      std::cout << "intersection " << shared_vertex << "(" << points[p].x
        << ", " << points[p].y << ")" << std::endl;
      for (auto ls : LU) {
        intersections.emplace(shared_vertex, ls);
      }
    } // else no intersection to add.

    // 5. Delete L and C from T.
    for (auto del_l : L) {
      auto found_segment=std::find_if(T.begin(), T.end(),
        [del_l](const StatusType::value_type& v)->bool {
          return v.second==del_l;
        });
      if (found_segment!=T.end()) {
        std::cout << "Delete L from T "<< del_l << std::endl;
        T.erase(found_segment);
      }
    }

    // 6. 7. Insert U and C into T.
    // Store the leftmost and rightmost of U for use in 11.
    StatusLoc least_loc;
    StatusType::iterator least_u;
    StatusLoc greatest_loc;
    StatusType::iterator greatest_u;
    size_t t_before_cnt=T.size();
    bool first=true;
    for (auto segins : U) {
      const auto& pa=points[lines[segins].p0];
      const auto& pb=points[lines[segins].p1];
      double dy=pa.y-pb.y;
      double dx=pa.x-pb.x;
      // atan2 returns radians between -pi and pi.
      double angle=std::atan2(dy, dx);
      // Or equal to 0 so that horizontal is the last in the list.
      if (angle <= 0) {
        angle+=M_PI; // boost::math::constants::pi<double>();
      }
      StatusLoc loc{points[p].x, points[p].y, angle};
      std::cout << "Putting U point "<<segins<<" into T: ("
        << loc.x << ", " << loc.y << ", " << loc.radians << ")-"
        << segins << std::endl;
      bool success;
      if (first) {
        least_loc=loc;
        greatest_loc=loc;
        std::tie(least_u, success)=T.emplace(loc, segins);
        if (!success) {
          std::cout << "Point "<<p<<" is already in T" << std::endl;
          assert(success);
        }
        greatest_u=least_u;
        first=false;
      } else {
        StatusType::iterator ins_iter;
        std::tie(ins_iter, success)=T.emplace(loc, segins);
        if (!success) {
          std::cout << "Point "<<p<<" is already in T" << std::endl;
          assert(success);
        }
        if (loc<least_loc) {
          least_loc=loc;
          least_u=ins_iter;
        } else if (loc>greatest_loc) {
          greatest_loc=loc;
          greatest_u=ins_iter;
        }
      }
    }
    assert(t_before_cnt+U.size()==T.size());
    std::cout << "T ";
    for (auto tentry : T) {
      std::cout << tentry.first << '-' << tentry.second << ", ";
    }
    std::cout << std::endl;
    std::cout << "\tLeast u "<<least_u->first<<'-'<<least_u->second
      <<" greatest u " <<greatest_u->first <<'-'<<greatest_u->second<<std::endl;

    // 7. if U+C==0
    std::cout << "U size " << U.size() << std::endl;
    if (U.size()==0) {
      // 8. let sl and sr be the left and right neighbors of p.
      auto p_begin=std::find_if(T.begin(), T.end(),
        [p, &lines](const std::pair<StatusLoc,SegmentIdx>& tentry)->bool {
          const auto& l=lines.at(tentry.second);
          return l.Contains(p);
        });
      if (p_begin!=T.begin()) {
        auto sl=p_begin;
        --sl;
        while (p_begin!=T.end() && lines.at(p_begin->second).Contains(p)) {
          ++p_begin;
        }
        if (p_begin!=T.end()) {
          std::cout << "8. find event " << sl->second << "-" << p_begin->second
            << std::endl;
          FindNewEvent(sl->second, p_begin->second, p, lines, points, Q,
            intersections, intersection_points);
        } else {
          std::cout << "8. " << p << " has no right neighbor. no search" << std::endl;
        }
      } else {
        std::cout <<"8. " << p << " has no left neighbor. no search" << std::endl;
      }
    } else {
      if (least_u!=T.begin()) {
        auto sl=least_u;
        --sl;
        std::cout << "8. Find event " << sl->second << "-" <<
          least_u->second << std::endl;
        FindNewEvent(sl->second, least_u->second, p, lines, points, Q,
          intersections, intersection_points);
      } else {
        std::cout << "8. Least_u is at beginning. No search." << std::endl;
      }
      auto sr=greatest_u;
      ++sr;
      if (sr!=T.end()) {
        FindNewEvent(greatest_u->second, sr->second, p, lines, points, Q,
          intersections, intersection_points);
      } else {
        std::cout << "8. Greatest_u is at end. No search." << std::endl;
      }
    }

    next_event=Q.begin();
  }
  return std::make_tuple(intersection_points, intersections);
}

