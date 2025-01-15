
#ifndef STEINERPOINTS_HPP
#define STEINERPOINTS_HPP
#include "CGAL_CUSTOM_CONSTRAINED_DELAUNAY_TRIANGULATION_2.h"
#include "triangulation_utils.hpp"
bool point_inside_triangle(Point p1, Point p2, Point p3,Point query);

//Returns the point located at the midpoint of the longest edge of a given face
Point steiner_longest_edge(Custom_CDT& cdel_tri, Face_handle face);
//testORass=1 add steiner on the midpoint of the longest edge if it reduces the obtuse count. testORadd!=1 return the number of obtuse triangles if the steiner would be added
int test_add_steiner_longest_edge(Custom_CDT& ccdt,Face_handle& face, const Polygon_2 & region_boundary, boost::associative_property_map<std::unordered_map<Face_handle, bool>> &in_domain,int testORadd);

//Returns the point located at the centroid of the face
Point steiner_centroid(Custom_CDT& cdel_tri, Face_handle face);
//testORass=1 add steiner on the centroid if it reduces the obtuse count. testORadd!=1 return the number of obtuse triangles if the steiner would be added
int test_add_steiner_centroid(Custom_CDT& ccdt,Face_handle& face, const Polygon_2 & region_boundary, boost::associative_property_map<std::unordered_map<Face_handle, bool>> &in_domain,int testORadd);

//Returns the point that is located at the projection of the obtuse angle on the oposite edge
Point steiner_projection(Custom_CDT& cdt, Face_handle F);
//testORass=1 add steiner on the projection of the obtuse angle if it reduces the obtuse count. testORadd!=1 return the number of obtuse triangles if the steiner would be added
int test_add_steiner_projection(Custom_CDT& ccdt,Face_handle& face, const Polygon_2 & region_boundary, boost::associative_property_map<std::unordered_map<Face_handle, bool>> &in_domain,int testORadd);

//Returns the point that is located at the circumcenter of the face
Point steiner_circumcenter(Custom_CDT& cdel_tri, Face_handle face);

//Returns the obtuse count after adding a a steiner and merging neighbor faces if possible.Simulation
//if testORadd !=1 calls test_merge_steiner else it performs the merge.
int test_add_steiner_circumcenter(Custom_CDT& ccdt,Face_handle& face, const Polygon_2 & region_boundary, boost::associative_property_map<std::unordered_map<Face_handle, bool>> &in_domain,int testORadd);

int simulate_merge_steiner(Custom_CDT& ccdt,Face_handle& ob_face,const Polygon_2& region_boundary,Point& neighbor_point,bool& found);
int test_add_steiner_merge(Custom_CDT& ccdt,Face_handle& face,Polygon_2& region_polygon,boost::associative_property_map<std::unordered_map<Face_handle, bool>> &in_domain, bool& found,Point& neighbor,int testORadd);
#endif // TRIANGULATION_UTILS_HPP
