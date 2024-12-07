#include <CGAL/Qt/Basic_viewer_qt.h>
#include <boost/property_map/property_map.hpp>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <iterator>
#include <ostream>
#include <fstream>
#include <vector>
#include <map>
#include <unordered_map>
//custom
#include "../Header Files/boost_utils.hpp"
#include "../Header Files/triangulation_utils.hpp"
//boost libraries
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include "../Header Files/CGAL_CUSTOM_CONSTRAINED_DELAUNAY_TRIANGULATION_2.h"
//cgal libraries
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Triangulation_face_base_2.h>
#include <CGAL/draw_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_2.h>
#include <CGAL/mark_domain_in_triangulation.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>
#include <CGAL/Triangulation_data_structure_2.h>

#include <CGAL/polygon_function_objects.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/draw_polygon_2.h>
#include <vector>
#include <string>
#include <algorithm>


#include <random>

double getRandomDouble() {
// Create a random number generator
std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> dis(0.0, 1.0);

// Generate and return a random double in the range [0, 1)
return dis(gen);
}
void simulated_annealing(Custom_CDT& ccdt, Polygon_2 region_polygon,boost::associative_property_map<std::unordered_map<Face_handle, bool>> &in_domain, int steiner_count, double a, double b, int L){

}
