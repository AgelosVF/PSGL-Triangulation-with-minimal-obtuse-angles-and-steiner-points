//custom
#include "../Header Files/boost_utils.hpp"
#include "../Header Files/triangulation_utils.hpp"
#include "../Header Files//CGAL_CUSTOM_CONSTRAINED_DELAUNAY_TRIANGULATION_2.h"

#include<CGAL/Kernel/global_functions_2.h>
#include <CGAL/enum.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <CGAL/Polygon_2.h>
#include <CGAL/polygon_function_objects.h>


template <typename T>
std::string to_string(const T& number) {
	std::ostringstream oss;
	oss << number;
	return oss.str();
}

// Function to generate the output JSON file
void generate_output_json(const Custom_CDT& cdel_tri, const boost::optional<std::string>& inst_iud, const std::vector<Point>& initial_points,Polygon_2 region_polygon) {
	// Collect Steiner points
	std::vector<Point> steiner_points;
	for (auto vit = cdel_tri.finite_vertices_begin(); vit != cdel_tri.finite_vertices_end(); ++vit) {
		Point p = vit->point();
		//if we cant find the point in the initial points it means its a steiener that was added after
		if (std::find(initial_points.begin(), initial_points.end(), p) == initial_points.end()) {
			steiner_points.push_back(p);
		}
	}
	// Collect edges
	std::vector<std::pair<int, int>> edges;
	std::map<Point, int> point_index_map;
	int index = 0;
	//go through all the vertices and map the indexes
	for (auto vit = cdel_tri.finite_vertices_begin(); vit != cdel_tri.finite_vertices_end(); ++vit) {
		point_index_map[vit->point()] = index++;
	}

	//go throught all the edges and map theyr source and target
	for (auto eit = cdel_tri.finite_edges_begin(); eit != cdel_tri.finite_edges_end(); ++eit) {
		auto segment = cdel_tri.segment(*eit);
		int source_index = point_index_map[segment.source()];
		int target_index = point_index_map[segment.target()];

		//calculate the midpoint of the edge and if it is inside the region boundary add the edge.
		auto midpoint= CGAL::midpoint(segment.source(), segment.target());
		if( (region_polygon.bounded_side(midpoint)== CGAL::ON_BOUNDED_SIDE) || (region_polygon.bounded_side(midpoint) == CGAL::ON_BOUNDARY) )
			edges.emplace_back(source_index, target_index);
	}

	//Create the JSON structure (boost tree)
	boost::property_tree::ptree root;
	//put name and instance iud
	root.put("content_type", "CG_SHOP_2025_Solution");
	root.put("instance_uid", inst_iud.get());

	boost::property_tree::ptree steiner_points_x;
	boost::property_tree::ptree steiner_points_y;

	//put the steiner point
	for (const auto& p : steiner_points) {
		boost::property_tree::ptree x_node;
		x_node.put("", to_string(p.x()));
		steiner_points_x.push_back(std::make_pair("", x_node));

		boost::property_tree::ptree y_node;
		y_node.put("", to_string(p.y()));
		steiner_points_y.push_back(std::make_pair("", y_node));
	}
	root.add_child("steiner_points_x", steiner_points_x);
	root.add_child("steiner_points_y", steiner_points_y);

	//make the edge tree
	boost::property_tree::ptree edges_node;
	for (const auto& edge : edges) {

		boost::property_tree::ptree edge_node;
		edge_node.push_back(std::make_pair("", boost::property_tree::ptree(std::to_string(edge.first))));
		edge_node.push_back(std::make_pair("", boost::property_tree::ptree(std::to_string(edge.second))));
		edges_node.push_back(std::make_pair("", edge_node));
	}
	//add it to the root
	root.add_child("edges", edges_node);

	std::string filename= "Results/"+inst_iud.get()+"_output.json";
	// Write the JSON to a file
	std::ofstream output_file(filename);
	boost::property_tree::write_json(output_file, root);
	output_file.close();
}


void generate_input_json(const Custom_CDT& cdel_tri, const boost::optional<std::string>& inst_iud,Polygon_2 region_polygon) {
	// Extract points and map them to indices
	std::vector<std::string> points_x, points_y;
	std::map<Point, int> point_index_map;
	int index = 0;


	//get the point coordinates
	for (auto vit = cdel_tri.finite_vertices_begin(); vit != cdel_tri.finite_vertices_end(); ++vit) {
		const auto& point = vit->point();
		point_index_map[point] = index++;
		points_x.push_back(to_string(CGAL::to_double(point.x())));
		points_y.push_back(to_string(CGAL::to_double(point.y())));
	}
//-------------------------------------------------    
//-----------------------------------------
	// Collect edges as additional constraints
	std::vector<std::pair<int, int>> additional_constraints;
	for (auto eit = cdel_tri.finite_edges_begin(); eit != cdel_tri.finite_edges_end(); ++eit) {
		auto segment = cdel_tri.segment(*eit);
		int source_index = point_index_map[segment.source()];
		int target_index = point_index_map[segment.target()];

		// Calculate the midpoint of the edge
		auto midpoint = CGAL::midpoint(segment.source(), segment.target());
		//make if the midpoint is outside of the region boundary we ignore the edge
		if (region_polygon.bounded_side(midpoint)== CGAL::ON_BOUNDED_SIDE)
			additional_constraints.emplace_back(source_index, target_index);
	}

	// Create the JSON structure (boost property tree)
	boost::property_tree::ptree root;

	// Add basic metadata
	root.put("instance_uid", inst_iud.get());
	root.put("num_points", points_x.size());

	// Add point coordinates
	boost::property_tree::ptree points_x_node, points_y_node;
	for (const auto& x : points_x) {
		boost::property_tree::ptree x_node;
		x_node.put("", x);
		points_x_node.push_back(std::make_pair("", x_node));
	}
	for (const auto& y : points_y) {
		boost::property_tree::ptree y_node;
		y_node.put("", y);
		points_y_node.push_back(std::make_pair("", y_node));
	}
	root.add_child("points_x", points_x_node);
	root.add_child("points_y", points_y_node);

	// Add additional constraints
	root.put("num_constraints", additional_constraints.size());
	boost::property_tree::ptree constraints_node;
	for (const auto& constraint : additional_constraints) {
		boost::property_tree::ptree constraint_node;
		constraint_node.push_back(std::make_pair("", boost::property_tree::ptree(std::to_string(constraint.first))));
		constraint_node.push_back(std::make_pair("", boost::property_tree::ptree(std::to_string(constraint.second))));
		constraints_node.push_back(std::make_pair("", constraint_node));
	}
	root.add_child("additional_constraints", constraints_node);

	// Write JSON to a file
	std::string filename = "Results/" + inst_iud.get() + "_input.json";
	std::ofstream output_file(filename);
	boost::property_tree::write_json(output_file, root);
	output_file.close();
}
