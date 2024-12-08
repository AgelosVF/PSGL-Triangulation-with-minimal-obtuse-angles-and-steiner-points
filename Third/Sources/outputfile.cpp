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
#include <string>


template <typename T>
std::string to_string(const T& number) {
	std::ostringstream oss;
	oss << number;
	return oss.str();
}

// Function to generate the output JSON file
void generate_output_json(const Custom_CDT& cdel_tri, const boost::optional<std::string>& inst_iud, const std::vector<Point>& initial_points, Polygon_2 region_polygon, std::string write_file) {
    // Helper lambda to convert a CGAL::FT to a rational string
    auto rational_to_string = [](const Kernel::FT& coord) {
        const auto exact_coord = CGAL::exact(coord); // Convert to exact representation
        std::ostringstream oss;
        oss << exact_coord.get_num() << "/" << exact_coord.get_den(); // Output as numerator/denominator
        return oss.str();
    };

    // Map to store indices for all points (initial + Steiner)
    std::map<Point, int> point_index_map;
    int index = 0;

    // Populate the map with initial points, assigning their original indices
    for (size_t i = 0; i < initial_points.size(); ++i) {
        point_index_map[initial_points[i]] = static_cast<int>(i); // Map point to its original index
    }

    // Collect Steiner points
    std::vector<Point> steiner_points;
    for (auto vit = cdel_tri.finite_vertices_begin(); vit != cdel_tri.finite_vertices_end(); ++vit) {
        Point p = vit->point();
        // Add point to map if it's not an initial point
        if (point_index_map.find(p) == point_index_map.end()) {
            point_index_map[p] = index++; // Assign new index for Steiner points
            steiner_points.push_back(p); // Add to Steiner points vector
        }
    }

    // Collect edges
    std::vector<std::pair<int, int>> edges;

    // Process all edges of the triangulation
    for (auto eit = cdel_tri.finite_edges_begin(); eit != cdel_tri.finite_edges_end(); ++eit) {
        auto segment = cdel_tri.segment(*eit); // Get edge segment and use the indexes of the map
        int source_index = point_index_map[segment.source()]; 
        int target_index = point_index_map[segment.target()]; 

        // Check if the midpoint of the edge lies inside the region boundary
        auto midpoint = CGAL::midpoint(segment.source(), segment.target());
        if ((region_polygon.bounded_side(midpoint) == CGAL::ON_BOUNDED_SIDE) || 
            (region_polygon.bounded_side(midpoint) == CGAL::ON_BOUNDARY)) {
            edges.emplace_back(source_index, target_index); // Add edge if inside the region
        }
    }

    // Create the JSON structure
    boost::property_tree::ptree root;
    root.put("content_type", "CG_SHOP_2025_Solution"); // Add content type
    root.put("instance_uid", inst_iud.get()); // Add instance UID

    boost::property_tree::ptree steiner_points_x; // JSON array for Steiner x-coordinates
    boost::property_tree::ptree steiner_points_y; // JSON array for Steiner y-coordinates

    // Add Steiner points in rational format
    for (const auto& steiner_point : steiner_points) {
        boost::property_tree::ptree x_node;
        x_node.put("", rational_to_string(steiner_point.x())); // Add x-coordinate in x/y form
        steiner_points_x.push_back(std::make_pair("", x_node));

        boost::property_tree::ptree y_node;
        y_node.put("", rational_to_string(steiner_point.y())); // Add y-coordinate in x/y form
        steiner_points_y.push_back(std::make_pair("", y_node));
    }
    root.add_child("steiner_points_x", steiner_points_x); // Attach x-coordinates to root
    root.add_child("steiner_points_y", steiner_points_y); // Attach y-coordinates to root

    // Add edges to the JSON
    boost::property_tree::ptree edges_node; // JSON array for edges
    for (const auto& edge : edges) {
        boost::property_tree::ptree edge_node;
        edge_node.push_back(std::make_pair("", boost::property_tree::ptree(std::to_string(edge.first)))); // Source index
        edge_node.push_back(std::make_pair("", boost::property_tree::ptree(std::to_string(edge.second)))); // Target index
        edges_node.push_back(std::make_pair("", edge_node));
    }
    root.add_child("edges", edges_node); // Attach edges to root

    // Write the JSON structure to the output file
    std::ofstream output_file(write_file);
    boost::property_tree::write_json(output_file, root);
    output_file.close();
}
