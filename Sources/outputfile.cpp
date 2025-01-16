//custom
#include "../Header Files/boost_utils.hpp"
#include "../Header Files/triangulation_utils.hpp"
#include "../Header Files//CGAL_CUSTOM_CONSTRAINED_DELAUNAY_TRIANGULATION_2.h"

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include<CGAL/Kernel/global_functions_2.h>
#include <CGAL/enum.h>
#include <CGAL/draw_triangulation_2.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <CGAL/Polygon_2.h>
#include <CGAL/polygon_function_objects.h>
#include <iostream>
#include <string>
#include <boost/algorithm/string/replace.hpp>

template <typename T>
std::string to_string(const T& number) {
	std::ostringstream oss;
	oss << number;
	return oss.str();
}

std::string rational_to_string(const Kernel::FT& coord) {
    const auto exact_coord = CGAL::exact(coord);
    std::ostringstream oss;
    oss << exact_coord.get_num() << "/" << exact_coord.get_den();
    return oss.str();
}
// Function to generate the output JSON file

// Function to generate the output JSON file
void generate_output_json(const Custom_CDT& cdel_tri, const boost::optional<std::string>& inst_iud, const std::vector<Point>& initial_points, Polygon_2 region_polygon, std::string write_file, const std::string& method, const boost::property_tree::ptree& parameters,int obtuse_count, double final_rate) {

    // Map to store indices for all points (initial + Steiner)
    std::map<Point, int> point_index_map;
    int index = 0;

    for (size_t i = 0; i < initial_points.size(); ++i) {
	point_index_map[initial_points[i]] = index;
	index++;
    }

    // Collect Steiner points
    std::vector<Point> steiner_points;
    for (auto vit = cdel_tri.finite_vertices_begin(); vit != cdel_tri.finite_vertices_end(); ++vit) {
	Point p = vit->point();
	if (point_index_map.find(p) == point_index_map.end()) {
	    point_index_map[p] = index++;
	    steiner_points.push_back(p);
	    }
    }


    // Collect edges
    std::vector<std::pair<int, int>> edges;
    for (auto eit = cdel_tri.finite_edges_begin(); eit != cdel_tri.finite_edges_end(); ++eit) {
	auto segment = cdel_tri.segment(*eit);
	int source_index = point_index_map[segment.source()];
	int target_index = point_index_map[segment.target()];

	auto midpoint = CGAL::midpoint(segment.source(), segment.target());
	if ((region_polygon.bounded_side(midpoint) == CGAL::ON_BOUNDED_SIDE) || (region_polygon.bounded_side(midpoint) == CGAL::ON_BOUNDARY)) {
	    edges.emplace_back(source_index, target_index);
	}
	else{
	}
    }

    // Create the JSON structure
    boost::property_tree::ptree root;
    root.put("content_type", "CG_SHOP_2025_Solution");
    root.put("instance_uid", inst_iud.get());
    root.put("method", method); // Add the method

    // Add parameters as a JSON object
    root.add_child("parameters", parameters);

    root.put("obtuse_count", obtuse_count); // Add the obtuse triangle count
    root.put("random", false);  // Add random

    boost::property_tree::ptree steiner_points_x;
    boost::property_tree::ptree steiner_points_y;

    for (const auto& steiner_point : steiner_points) {
	boost::property_tree::ptree x_node;
	x_node.put("", rational_to_string(steiner_point.x()));
	steiner_points_x.push_back(std::make_pair("", x_node));

	boost::property_tree::ptree y_node;
	y_node.put("", rational_to_string(steiner_point.y()));
	steiner_points_y.push_back(std::make_pair("", y_node));
    }
    root.add_child("steiner_points_x", steiner_points_x);
    root.add_child("steiner_points_y", steiner_points_y);

    boost::property_tree::ptree edges_node;
    for (const auto& edge : edges) {
	boost::property_tree::ptree edge_node;
	edge_node.push_back(std::make_pair("", boost::property_tree::ptree(std::to_string(edge.first))));
	edge_node.push_back(std::make_pair("", boost::property_tree::ptree(std::to_string(edge.second))));
	edges_node.push_back(std::make_pair("", edge_node));
    }
    root.add_child("edges", edges_node);

    // Write the JSON structure to a string
    std::ostringstream oss;
    boost::property_tree::write_json(oss, root);
    std::string json_string = oss.str();

    // Replace escaped backslashes with regular slashes
    boost::replace_all(json_string, "\\/", "/");

    // Write the modified JSON string to the output file
    std::ofstream output_file(write_file);
    output_file << json_string;
    output_file.close();
}



