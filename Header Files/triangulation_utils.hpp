
// triangulation_utils.hpp
#ifndef TRIANGULATION_UTILS_HPP
#define TRIANGULATION_UTILS_HPP

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include "boost_utils.hpp"
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
int count_obtuse_faces(const Custom_CDT& cdel_tri, const Polygon_2& region_polygon);
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
//used for debugging
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

bool is_in_region_polygon(Face_handle f, const Polygon_2& region_polygon);

bool face_still_exists(Custom_CDT cdt,Face_handle face);

bool is_polygon_convex(const Point& v1, const Point& v2, const Point& v3, const Point& v4) ;

bool point_inside_triangle(Point p1, Point p2, Point p3,Point query);

double calculate_convergence_rate(int steiner_points, int prev_obtuse, int current_obtuse);
//In local search
double compute_distance(const Point& p1, const Point p2) ;
void print_polygon_vertices(const Polygon_2& polygon);
//------------------------------------------------------------------------------------------------------------//
//In output.cpp
//creates the output json file
void generate_output_json(
    const Custom_CDT& cdel_tri, 
    const boost::optional<std::string>& inst_iud, 
    const std::vector<Point>& initial_points, 
    Polygon_2 region_polygon, 
    std::string write_file, 
    const std::string& method, 
    const boost::property_tree::ptree& parameters, 
    int obtuse_count,
    double final_rate

);
//-----------------------------------------------------------------------------------------------------------//
int local_search(Custom_CDT& ccdt, Polygon_2 region_polygon, int loops, boost::associative_property_map<std::unordered_map<Face_handle, bool>> &in_domain,double& c_rate);
int simulated_annealing(Custom_CDT& ccdt, Polygon_2 region_polygon,boost::associative_property_map<std::unordered_map<Face_handle, bool>> &in_domain, int steiner_count, double a, double b, int L,double& c_rate,int methods);
int ant_colony(Custom_CDT& cdt,Polygon_2 boundary,double alpha ,double beta, double xi, double ps, double lambda, int kappa, int L, int steiner_points, double& c_rate);
//In Previous_Project
int previous_triangulation(Custom_CDT& cdel_tri,const Polygon_2& region_polygon);
int find_obtuse_vertex_index(const Face_handle& Face);
#endif // TRIANGULATION_UTILS_HPP
