
// triangulation_utils.hpp
#ifndef TRIANGULATION_UTILS_HPP
#define TRIANGULATION_UTILS_HPP

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <boost/property_map/property_map.hpp>
#include <string>
#include <unordered_map>
#include "CGAL_CUSTOM_CONSTRAINED_DELAUNAY_TRIANGULATION_2.h"
#include <CGAL/Polygon_2.h>

// CGAL Typedefs
typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Triangulation_vertex_base_2<Kernel> V_base;
typedef CGAL::Constrained_triangulation_face_base_2<Kernel> F_base;
typedef CGAL::Triangulation_data_structure_2<V_base, F_base> Tri_ds;
typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel, Tri_ds> CDT;
typedef CGAL::Constrained_triangulation_2<Kernel,Tri_ds> CT;
typedef CDT::Point Point;
typedef CDT::Vertex_handle Vertex_handle;
typedef CDT::Face_handle Face_handle;
typedef Custom_Constrained_Delaunay_triangulation_2<Kernel, Tri_ds> Custom_CDT;
typedef CGAL::Polygon_2<Kernel> Polygon_2;

//---------------------------------------------------------------------------------------------------------//
//In triangulation_utils.cpp
//returns true if face contains no constrain edge false if it does at least one
bool face_has_constrained_edge(const Custom_CDT& cdt, Face_handle face);
// fills the first cdt with the points,region boundary and additional_constrains. The region boundary is treated as a constrain
void fill_initial_cdt(Custom_CDT& cdel_tri, const std::vector<int>& points_x, const std::vector<int>& points_y, const std::vector<int>& region_boundary, const std::vector<std::pair<int, int>>& additional_constrains);
//return the number of obtuse triangles in the cdt that are inside the domain(region boundary)
int count_obtuse_faces(const Custom_CDT& cdel_tri, const boost::associative_property_map<std::unordered_map<Face_handle, bool>>& in_domain);
//takes a triangle returns true if it is obtuse false if it is not
bool is_obtuse_triangle(const CDT::Face_handle& f);
//takes three points and checks if the triangle they would create is obtuse or not
bool is_obtuse_simulated(const Point& p1, const Point& p2, const Point& p3);
// checks if the lines formed by p1 p2 and p3 p4 intersect
bool segments_intersect(const Point& p1, const Point& p2, const Point& p3, const Point& p4);
//checks if flip is valid by checking if it intersects with any other line
bool is_flip_valid(const CDT& cdel_tri, Face_handle f, Face_handle neighbor, int edge_index);
//simmulates fliping edge_index and if it reduces the obtuse triangles in the cdt it performs it
bool try_edge_flip(Custom_CDT& cdel_tri, Face_handle f, int edge_index, boost::associative_property_map<std::unordered_map<Face_handle, bool>>& in_domain);
//goes over all the obtuse triangles in the domain of the cdt and calls try_edge_flip for every of the non constrained edges of them to reduce the number
void reduce_obtuse_by_flips(Custom_CDT& cdel_tri,  boost::associative_property_map<std::unordered_map<Face_handle, bool>>& in_domain);

void print_point(const Point& point);
//find the face that is made from the three points in the ccdt
Face_handle find_face(Custom_CDT& ccdt,Point p1, Point p2, Point p3);
//find face that contains edge with vh1 vh2 save face at face and the index of the edge at index
bool find_face(Custom_CDT& ccdt,Vertex_handle vh1, Vertex_handle vh2,int& index, Face_handle& face);
//find the face handle of the face in a diffrent ccdt
Face_handle find_face(Custom_CDT& n_ccdt, Face_handle o_face);
//returns the vertex of the point
Vertex_handle find_vertex_by_point(const Custom_CDT& ccdt, const Point& point);
//returns yes if any of the vertices of the face is a part of a constrain edge
bool is_face_constrained(const Custom_CDT& ccdt, Face_handle face);
//---------------------------------------------------------------------------------------------------------//
//In SteinerPoints.cpp
bool point_inside_triangle(Point p1, Point p2, Point p3,Point query);

int find_longest_edge_index(Face_handle face);
//Returns the point located at the midpoint of the longest edge of a given face
Point steiner_longest_edge(Custom_CDT& cdel_tri, Face_handle face);
//testORass=1 add steiner on the midpoint of the longest edge if it reduces the obtuse count. testORadd!=1 return the number of obtuse triangles if the steiner would be added
int test_add_steiner_longest_edge(Custom_CDT& ccdt,Face_handle face, const Polygon_2 & region_boundary, boost::associative_property_map<std::unordered_map<Face_handle, bool>> &in_domain,int testORadd);

//Returns the point located at the centroid of the face
Point steiner_centroid(Custom_CDT& cdel_tri, Face_handle face);
//testORass=1 add steiner on the centroid if it reduces the obtuse count. testORadd!=1 return the number of obtuse triangles if the steiner would be added
int test_add_steiner_centroid(Custom_CDT& ccdt,Face_handle face, const Polygon_2 & region_boundary, boost::associative_property_map<std::unordered_map<Face_handle, bool>> &in_domain,int testORadd);

//Returns the point that is located at the projection of the obtuse angle on the oposite edge
Point steiner_projection(Custom_CDT& cdt, Face_handle F);
//testORass=1 add steiner on the projection of the obtuse angle if it reduces the obtuse count. testORadd!=1 return the number of obtuse triangles if the steiner would be added
int test_add_steiner_projection(Custom_CDT& ccdt,Face_handle face, const Polygon_2 & region_boundary, boost::associative_property_map<std::unordered_map<Face_handle, bool>> &in_domain,int testORadd);

//Returns the point that is located at the circumcenter of the face
Point steiner_circumcenter(Custom_CDT& cdel_tri, Face_handle face);

//Returns the obtuse count after adding a a steiner and merging neighbor faces if possible.Simulation
//if testORadd !=1 calls test_merge_steiner else it performs the merge.
int test_add_steiner_circumcenter(Custom_CDT& ccdt,Face_handle face, const Polygon_2 & region_boundary, boost::associative_property_map<std::unordered_map<Face_handle, bool>> &in_domain,int testORadd);

int simulate_merge_steiner(Custom_CDT& ccdt,Face_handle ob_face,const Polygon_2& region_boundary,int& best_merge);
int test_add_steiner_merge(Custom_CDT& ccdt,Face_handle face,Polygon_2& region_polygon,boost::associative_property_map<std::unordered_map<Face_handle, bool>> &in_domain, int& indx,int testORadd);


bool is_polygon_convex(const Point& v1, const Point& v2, const Point& v3, const Point& v4) ;
//------------------------------------------------------------------------------------------------------------//
//In output.cpp
//creates the output json file
void generate_output_json(const Custom_CDT& cdel_tri, const boost::optional<std::string>& inst_iud, const std::vector<Point>& initial_points,Polygon_2 region_polygon,std::string write_file);
//-----------------------------------------------------------------------------------------------------------//
//In local search
int local_search(Custom_CDT& ccdt, Polygon_2 region_polygon, int loops, boost::associative_property_map<std::unordered_map<Face_handle, bool>> &in_domain);
//In simmulate annealing
//
void simulated_annealing(Custom_CDT& ccdt, Polygon_2 region_polygon,boost::associative_property_map<std::unordered_map<Face_handle, bool>> &in_domain, int steiner_count, double a, double b, int L);
#endif // TRIANGULATION_UTILS_HPP
