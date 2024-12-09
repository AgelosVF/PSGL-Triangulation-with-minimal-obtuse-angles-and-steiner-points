#include "../Header Files/CGAL_CUSTOM_CONSTRAINED_DELAUNAY_TRIANGULATION_2.h"
#include "../Header Files/triangulation_utils.hpp"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Kernel/global_functions_3.h>
#include <CGAL/squared_distance_2.h>
#include <algorithm>
#include <cmath> // for std::sqrt
#include <iostream>
#include <random>

// Define the CGAL Kernel

// Function to compute the distance between two points
double compute_distance(const Point& p1, const Point p2) {
    // Compute the squared distance
    Kernel::FT sq_dist = CGAL::squared_distance(p1, p2);

    // Return the square root of the squared distance for the actual distance
    return std::sqrt(CGAL::to_double(sq_dist));
}

double calculate_r(const Face_handle& Face){

	Point p0=Face->vertex(0)->point();
	Point p1=Face->vertex(1)->point();
	Point p2=Face->vertex(2)->point();

	Point circumcenter=CGAL::circumcenter(p0,p1,p2);

	double R=compute_distance(p0,circumcenter);

	int l_side_i=find_longest_edge_index(Face);
	Point obtuse_vert=Face->vertex(l_side_i)->point();
	Kernel::Line_2 longest_edge(Face->vertex((l_side_i+1)%3)->point(),Face->vertex((l_side_i+2)%3)->point());
	Point projection=longest_edge.projection(obtuse_vert);
	double hfle=compute_distance(projection, obtuse_vert);
	double r= R/hfle;
	return r;
}

double h_projection(double r){
	return std::max(0.0,((r-1)/r));
}

double h_circumcenter(double r){
	return (r/(r+2));
}

double h_longest_edge(double r){
	return std::max(0.0,((3-2*r)/3));
}
double h_merge(const Custom_CDT& ccdt,const Face_handle& Face,Polygon_2& region_polygon){

	double h=0.0;
	
	for(int i=0;i<3;i++){
		Face_handle neighbor=Face->neighbor(i);
		//skip  constrained ,non obtuse, outside of region_boundary neighbors
		if( ccdt.is_infinite(neighbor) ){
			continue;
		}
		
		if( !(is_obtuse_triangle(neighbor)) ){
			continue;
		}
		if( !(is_in_region_polygon(neighbor,region_polygon)) ){
			continue;
		}
		if(is_face_constrained(ccdt,neighbor)){
			continue;
		}
		h=1.0;

	}
	return h;

}

//make sure there are obtuse faces before calling
Face_handle select_random_obtuse_face(const Custom_CDT& cdt, Polygon_2 boundary){

	std::vector<Face_handle> obtuse_faces;

	for(auto face= cdt.finite_faces_begin(); face!= cdt.finite_faces_end(); face++){
		if(is_obtuse_triangle(face) && is_in_region_polygon(face,boundary) ){
			obtuse_faces.push_back(face);
		}
	}
	if(obtuse_faces.empty()){
		std::cerr<<"Error tried to select obtuse face when they are none\n";
		exit(-1);
	}

	std::srand(static_cast<unsigned>(std::time(nullptr)));
	int random_index= (std::rand() %obtuse_faces.size());

	Face_handle o_face=obtuse_faces[random_index];

	if( o_face == nullptr){
		std::cerr<<"Error: couldnt locate original ccdt face in the copy cccdt"<<std::endl;
	}
	return o_face;



}



Face_handle impruveTriangulation(Custom_CDT& cdt,Polygon_2 boundary,double alpha,double beta,int steiner_count,double xi,double ps, const std::vector<double>& pheromones,double& ant_energy,int& selected_method,double best_E){

	//choose random face
	Face_handle face=select_random_obtuse_face(cdt, boundary);
	
	//choose random method
	
	//find heuristics
	double r=calculate_r(face);
	double h_mrg=h_merge(cdt,face,boundary);
	double h_circ=h_circumcenter(r);
	double h_le=h_longest_edge(r);
	double h_proj=h_projection(r);

	std::vector<double> heuristics= {h_proj,h_le,h_circ,h_mrg};
	std::vector<double> probabilities(pheromones.size());


	//chose random with probability p(k)= (t^xi * h^psi)/ Σ_{i} t^x_{i}h^ps_{i} 
	double denominator = 0.0;
	for (size_t i = 0; i < pheromones.size(); i++) {
		probabilities[i] = std::pow(pheromones[i], xi) * std::pow(heuristics[i], ps);
		denominator += probabilities[i];
	}

	//normalize the probabilities
	for (double& p :probabilities){
		p/=denominator;
	}

	// Perform weighted random selection
	std::random_device rd;
	std::mt19937 gen(rd());
	std::discrete_distribution<> dist(probabilities.begin(), probabilities.end());

	selected_method = dist(gen); // Randomly choose a method index based on probabilities
	
	int obtuse_faces=99999;
	switch (selected_method){
		case 0:
			std::cout<<"Selected projection\n";
			obtuse_faces=0;
			break;
		case 1:
			std::cout<<"Selected longest edge\n";
			obtuse_faces=1;
			break;
		case 2:
			std::cout<<"Selected circumcenter\n";
			obtuse_faces=2;
			break;
		case 3:
			std::cout<<"Selected merge\n";
			obtuse_faces=3;
			break;
		default:
			std::cerr<<"Error: chose invalid method\n";
			exit(-1);
	}

	if(obtuse_faces==0){
		ant_energy=0;
	}
	else{
		ant_energy=alpha*obtuse_faces +beta*(steiner_count+1);
	}

	return face;


}
double update_pheromon(double t, double lambda,double alpha, int obtuse_count, double beta, int steiner_count){
	double denom=1+alpha*obtuse_count + beta* steiner_count;
	return ( (1-lambda) * t + (1/denom) );

}

bool face_still_exists(Custom_CDT cdt,Face_handle face){
	
	Point p0 = face->vertex(0)->point();
	Point p1 = face->vertex(1)->point();
	Point p2 = face->vertex(2)->point();
	Point centroid= CGAL::centroid(p0,p1,p2);
	Face_handle n_face = cdt.locate(centroid);
	if( n_face == nullptr){
		return false;
	}
	Point np0=n_face->vertex(0)->point();
	Point np1=n_face->vertex(1)->point();
	Point np2=n_face->vertex(2)->point();

	std::set<Point> og_points={p0,p1,p2};
	std::set<Point> new_points={np0,np1,np2};

	return og_points==new_points;
}

void ant_colony(Custom_CDT& cdt,Polygon_2 boundary,double alpha ,double beta, double xi, double ps, double lambda, int kappa, int L, int steiner_points){


	std ::vector<double> pheromones(5,1.0);
	int current_steiner=steiner_points;
	double best_Energy=alpha*count_obtuse_faces(cdt,boundary) + beta*steiner_points;
	std::vector<double> ant_energy(kappa,best_Energy);
	std::vector<int> ant_method(kappa,-1);

	Custom_CDT best_cdt(cdt);
	for(int t=0;t<L;t++){

		//vector that will have each face each ant worked on
		std::vector<Face_handle> ant_faces(kappa);
		for(int ant=0;ant<kappa;ant++){
		
			//each ant should have its own cdt
			Face_handle edited_face = impruveTriangulation( cdt, boundary, alpha, beta, current_steiner, xi, ps, pheromones, ant_energy[ant], ant_method[ant],best_Energy);
			//then we count the obtuse faces of the "improved" cdt
			ant_faces[ant] = edited_face;
		}

		for(int ant=0;ant<kappa;ant++){
			if(ant_energy[ant]==0){
				//apply the method i to the face i and break since we reduced obtuse to 0
			}

		}
		//create index energy pair vector
		std::vector<std::pair<double, int>> energy_index_pairs;
		for(int i=0;i<kappa;i++){
			energy_index_pairs.emplace_back(ant_energy[i],i);
		}

		//sort them from smallest to largest
		std::sort(energy_index_pairs.begin(), energy_index_pairs.end());
		//choose the ants that reduced the energy of the best cdt
		//first apply the one that reduced the energy the most
		//then by always choosing the one that reduced the most 
		//to filter conflicts check if its the same face or 1 is the edited neighbor of the other
		for(int ant=0;ant<kappa;ant++){
			double current_ant_energy=energy_index_pairs[ant].first;
			int current_ant_index=energy_index_pairs[ant].second;
			if(current_ant_energy<best_Energy){
				Face_handle c_ant_face=ant_faces[current_ant_index];
				int c_meth=ant_method[current_ant_index];
				if(face_still_exists(best_cdt, c_ant_face)){
					//apply the steiner
					//update the pheramon[method]
				}

			}
			else{
				//the ant increased the energy so we stop looking
				break;
			}

			
		}



	}
}

