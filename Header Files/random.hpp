#ifndef STEINERPOINTS_HPP
#define STEINERPOINTS_HPP

#include "triangulation_utils.hpp"

//inserts a random point around the centroid using gaussian distribution
int random_and_flips(Custom_CDT& cdt, Face_handle face, Polygon_2 region);

//inserts a point on the longest edge using the projection of a random point around the centroid
bool random_point_on_edge(Custom_CDT& cdt,Face_handle face) ;

#endif // RANDOM_HPP
