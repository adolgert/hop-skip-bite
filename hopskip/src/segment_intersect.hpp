#ifndef _SEGMENT_INTERSECT_H_
#define _SEGMENT_INTERSECT_H_ 1

#include <map>
#include <utility>
#include <vector>

/*! Find intersections among line segments.
 *  Reading "Computational Geometry" by Berg et al., Chapter 2.
 */
std::tuple<std::vector<std::pair<double,double>>,
  std::multimap<size_t,size_t>>
segment_intersections(const std::vector<std::pair<double,double>>& inpoints,
    const std::vector<std::pair<size_t,size_t>>& inlines);

std::tuple<std::vector<std::pair<double,double>>,
  std::multimap<size_t,size_t>>
segment_intersections_sweep(
    const std::vector<std::pair<double,double>>& inpoints,
    const std::vector<std::pair<size_t,size_t>>& inlines);

std::vector<std::pair<size_t,size_t>>
segment_intersections_ab(const std::vector<std::pair<double,double>>& inpoints,
    const std::vector<std::pair<size_t,size_t>>& inlinesa,
    const std::vector<std::pair<size_t,size_t>>& inlinesb);
// _SEGMENT_INTERSECT_H_
#endif
