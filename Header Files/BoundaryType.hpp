//BoundaryType.hpp
#ifndef BoundaryType_HPP
#define BoundaryType_HPP

#include "../Header Files/CGAL_CUSTOM_CONSTRAINED_DELAUNAY_TRIANGULATION_2.h"
#include "../Header Files/triangulation_utils.hpp"
#include <CGAL/Polygon_2.h>
#include <vector>
#include <unordered_map>
#include <unordered_set>


bool convex_no_constrains(Polygon_2 boundary,std::vector<std::pair<int,int>> additional_constrains);

bool non_convex_parrallel(Polygon_2 boundary);
bool convex_cycle_constrains(const std::vector<int>& region_boundary, const std::vector<std::pair<int, int>>& additional_constraints,std::vector<int>& cycle);
#endif
