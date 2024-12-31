#include "../Header Files/CGAL_CUSTOM_CONSTRAINED_DELAUNAY_TRIANGULATION_2.h"
#include "../Header Files/triangulation_utils.hpp"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Kernel/global_functions_3.h>
#include <CGAL/enum.h>
#include <CGAL/squared_distance_2.h>
#include <algorithm>
#include <cmath> 
#include <iostream>
#include <random>

void ant_colony(Custom_CDT& cdt, Polygon_2 boundary, double alpha, double beta, double xi, double ps, double lambda, int kappa, int L, int steiner_points) {
    // Initialize variables
    double global_best_energy = std::numeric_limits<double>::max();
    Custom_CDT global_best_solution;

    // Initialize pheromone levels on edges
    std::map<CDT::Edge, double> pheromone_levels;
    for (auto edge_it = cdt.edges_begin(); edge_it != cdt.edges_end(); ++edge_it) {
        pheromone_levels[*edge_it] = 1.0; // Initial pheromone level
    }

    // Main loop for ant colony optimization
    for (int iteration = 0; iteration < kappa; ++iteration) {
        std::vector<Custom_CDT> ant_solutions;
        std::vector<double> ant_energies;

        // Each ant builds a solution
        for (int ant = 0; ant < L; ++ant) {
            Custom_CDT ant_cdt = cdt; // Start with the initial triangulation

            // Ant's construction phase: Add Steiner points
            for (int i = 0; i < steiner_points; ++i) {
                Face_handle selected_face = select_face_based_on_probability(ant_cdt, alpha, beta, pheromone_levels);
                if (!selected_face) {
                    std::cerr << "No valid face selected during construction phase.\n";
                    break;
                }

                Point new_steiner_point = compute_steiner_point(ant_cdt, selected_face, boundary);
                ant_cdt.insert(new_steiner_point);
            }

            // Calculate energy of the ant's solution
            double ant_energy = calculate_energy(ant_cdt);
            ant_solutions.push_back(ant_cdt);
            ant_energies.push_back(ant_energy);

            // Update global best solution if necessary
            if (ant_energy < global_best_energy) {
                global_best_energy = ant_energy;
                global_best_solution = ant_cdt;
            }
        }

        // Update pheromones
        for (auto& edge_pheromone : pheromone_levels) {
            edge_pheromone.second *= (1 - xi); // Evaporation
        }

        for (size_t ant = 0; ant < ant_solutions.size(); ++ant) {
            const Custom_CDT& ant_solution = ant_solutions[ant];
            double ant_energy = ant_energies[ant];

            for (auto edge_it = ant_solution.edges_begin(); edge_it != ant_solution.edges_end(); ++edge_it) {
                pheromone_levels[*edge_it] += ps / ant_energy; // Pheromone deposit proportional to solution quality
            }
        }

        // Optional: Print progress
        std::cout << "Iteration " << iteration + 1 << "/" << kappa << ", Best Energy: " << global_best_energy << "\n";
    }

    // Replace original CDT with the best solution found
    cdt = global_best_solution;
}

// Helper function to select a face based on probabilities
Face_handle select_face_based_on_probability(Custom_CDT& cdt, double alpha, double beta, const std::map<CDT::Edge, double>& pheromone_levels) {
    std::vector<std::pair<Face_handle, double>> face_probabilities;
    double total_probability = 0.0;

    for (auto face_it = cdt.all_faces_begin(); face_it != cdt.all_faces_end(); ++face_it) {
        if (cdt.is_infinite(face_it)) continue;

        double heuristic_value = compute_heuristic_value(cdt, face_it);
        double pheromone_value = compute_face_pheromone(face_it, pheromone_levels);
        double probability = pow(pheromone_value, alpha) * pow(heuristic_value, beta);

        face_probabilities.emplace_back(face_it, probability);
        total_probability += probability;
    }

    if (total_probability == 0.0) return Face_handle(); // No valid faces

    double random_value = ((double)rand() / RAND_MAX) * total_probability;
    for (const auto& [face, probability] : face_probabilities) {
        random_value -= probability;
        if (random_value <= 0.0) {
            return face;
        }
    }

    return face_probabilities.back().first; // Fallback
}

// Helper function to compute the heuristic value for a face
double compute_heuristic_value(Custom_CDT& cdt, Face_handle face) {
    // Example: Use the circumradius of the face as the heuristic
    return cdt.triangle(face).area();
}

// Helper function to compute pheromone for a face
double compute_face_pheromone(Face_handle face, const std::map<CDT::Edge, double>& pheromone_levels) {
    double total_pheromone = 0.0;
    for (int i = 0; i < 3; ++i) {
        CDT::Edge edge = CDT::Edge(face, i);
        auto it = pheromone_levels.find(edge);
        if (it != pheromone_levels.end()) {
            total_pheromone += it->second;
        }
    }
    return total_pheromone;
}

// Helper function to compute Steiner point within a face
Point compute_steiner_point(Custom_CDT& cdt, Face_handle face, const Polygon_2& boundary) {
    // Example: Use the circumcenter of the face
    Point circumcenter = CGAL::circumcenter(cdt.triangle(face));
    if (is_in_region_polygon(circumcenter, boundary)) {
        return circumcenter;
    }

    // Fallback: Use centroid if circumcenter is outside the boundary
    return CGAL::centroid(cdt.triangle(face));
}
