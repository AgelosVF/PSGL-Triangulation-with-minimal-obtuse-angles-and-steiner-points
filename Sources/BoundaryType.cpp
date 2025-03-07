#include "../Header Files/CGAL_CUSTOM_CONSTRAINED_DELAUNAY_TRIANGULATION_2.h"
#include "../Header Files/triangulation_utils.hpp"
#include "../Header Files/BoundaryType.hpp"
#include <CGAL/Polygon_2.h>
#include <algorithm>
#include <ostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <iostream>


bool is_convex(Polygon_2 boundary){
	if(boundary.is_convex())
		return true;
	else
		return false;
}

bool no_constrains(std::vector<std::pair<int,int>> additional_constrains){
	if(additional_constrains.empty())
		return true;
	else
		return false;
}

bool non_convex_parrallel(Polygon_2 boundary){

	if(is_convex(boundary))
		return false;
	for(auto  e  : boundary.edges()){
		
		Point p= e.start();
		Point p2=e.end();
		if(p.x()!=p2.x() && p.y() != p2.y()){
			return false;
		}
	}
	return true;
}


// Function to detect cycles in a graph using DFS and check if the target edge is part of the cycle
bool hasCycleDFS(int node, int parent, std::unordered_map<int, std::vector<int>>& graph, std::unordered_set<int>& visited, std::vector<int>& path, std::vector<int>& cycle, const std::pair<int, int>& target_edge) {

	visited.insert(node);
	path.push_back(node);
	for(unsigned int i=0;i<path.size();i++){
	}
	//loop throught the neighbors
	for (int neighbor : graph[node]) {
		if (neighbor == parent) {
			continue; // Skip the parent
		 }
		if (visited.count(neighbor)) {
			// Cycle detected - check if it includes the target edge
			auto it = std::find(path.begin(), path.end(), neighbor);
			if (it != path.end()) {
				cycle.clear();
				cycle.insert(cycle.end(), it, path.end());
				cycle.push_back(neighbor); // Close the cycle

				// Check if the target edge is part of the cycle
				for (size_t i = 0; i < cycle.size() - 1; ++i) {
					if ((cycle[i] == target_edge.first && cycle[i + 1] == target_edge.second) ||
					(cycle[i] == target_edge.second && cycle[i + 1] == target_edge.first)) {
					return true;
					}
				}
			}
		}
		else if (hasCycleDFS(neighbor, node, graph, visited, path, cycle, target_edge)) {
			return true;
		}
	}
	path.pop_back();
	return false;
}

bool convex_cycle_constrains(const std::vector<int>& region_boundary, const std::vector<std::pair<int, int>>& additional_constraints,std::vector<int>& cycle) {
	// Build the graph from the region boundary
	std::unordered_map<int, std::vector<int>> graph;
	for (size_t i = 0; i < region_boundary.size(); ++i) {
		int u = region_boundary[i];
		int v = region_boundary[(i + 1) % region_boundary.size()]; // Wrap around to form a closed polygon
		graph[u].push_back(v);
		graph[v].push_back(u);
	}

	// Add additional constraints to the graph and check for cycles
	for (const auto& constraint : additional_constraints) {
		int u = constraint.first;
		int v = constraint.second;

		//  Add the edge to the graph
		graph[u].push_back(v);
		graph[v].push_back(u);

		// Check for cycles with the edge
		std::unordered_set<int> visited;
		std::vector<int> path; // To track the current DFS path
		cycle.clear(); // Clear the cycle vector before starting DFS
		if (hasCycleDFS(u, -1, graph, visited, path, cycle, constraint)) {
			return true; // Cycle detected
		}

	// Remove the edge to restore the original graph
	//graph[u].pop_back();
	//graph[v].pop_back();
	}

	return false; // No cycles detected
}

int boundary_type(Polygon_2 boundary,const std::vector<int>& region_boundary, const std::vector<std::pair<int,int>>& additional_constrains, std::vector<int>& closed_p){

	if(is_convex(boundary)){
		if(no_constrains(additional_constrains)){
			std::cout<<"Is convex with no additional constrains\n";
			return 1;
		}
		else if (convex_cycle_constrains(region_boundary, additional_constrains, closed_p)){
			std::cout<<"Is convex with constrains \" ";
			for(unsigned int i=0;i<closed_p.size();i++){
				std::cout<<closed_p[i]<<"->";
			}
			std::cout<<" \" making a closed polygon\n";
			return 3;
		}
		else{
			std::cout<<"Is convex with open constrains\n";
			closed_p.clear(); // Clear the cycle vector before starting DFS
			return 2;
		}
	}
	else if (non_convex_parrallel(boundary)&&no_constrains(additional_constrains)) {
		std::cout<<"Is non convex with parrallel boundary edges to xx' or yy' and no additional constrains\n";
		return 4;
	
	}
	std::cout<<"Doesnt belong to any of the 4 predetermined categories\n";
	return 5;
}

