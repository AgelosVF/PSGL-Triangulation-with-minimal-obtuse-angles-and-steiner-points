//BoundaryType.hpp
#ifndef BoundaryType_HPP
#define BoundaryType_HPP

#include "../Header Files/CGAL_CUSTOM_CONSTRAINED_DELAUNAY_TRIANGULATION_2.h"
#include "../Header Files/triangulation_utils.hpp"
#include <CGAL/Polygon_2.h>
#include <vector>
#include <unordered_map>
#include <unordered_set>


int boundary_type(Polygon_2 boundary,const std::vector<int>& region_boundary, const std::vector<std::pair<int,int>>& additional_constrains, std::vector<int>& closed_p);
#endif
