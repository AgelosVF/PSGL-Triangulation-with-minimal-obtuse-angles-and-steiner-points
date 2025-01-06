
#include "../Header Files/CGAL_CUSTOM_CONSTRAINED_DELAUNAY_TRIANGULATION_2.h"
#include "../Header Files/triangulation_utils.hpp"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Kernel/global_functions_3.h>
#include <CGAL/enum.h>
#include <CGAL/squared_distance_2.h>
#include <algorithm>
#include <cmath> 
#include <iostream>
#include <ostream>
#include <random>


//-------------------------------Helping functions-----------------------------------------------------//

bool face_still_exists(Custom_CDT cdt, Face_handle face){

	Point p0 = face->vertex(0)->point();
	Point p1 = face->vertex(1)->point();
	Point p2 = face->vertex(2)->point();
	Point centroid = CGAL::centroid(p0,p1,p2);
	Face_handle n_face = cdt.locate(centroid);

	if( n_face == nullptr){
		std::cout<<"Got Null ptr when trying to locate face: ("<<p0.x()<<","<<p0.y()<<")  "<<"("<<p1.x()<<","<<p1.y()<<")  "<<"("<<p2.x()<<","<<p2.y()<<")\n";
		return false;
	}


	Point np0=n_face->vertex(0)->point();
	Point np1=n_face->vertex(1)->point();
	Point np2=n_face->vertex(2)->point();

	std::set<Point> og_points={p0,p1,p2};
	std::set<Point> new_points={np0,np1,np2};


	if(!(og_points==new_points)){
		
		std::cout<<"NEW PTR ("<<np0.x()<<","<<np0.y()<<")  "<<"("<<np1.x()<<","<<np1.y()<<")  "<<"("<<np2.x()<<","<<np2.y()<<")\n";
		std::cout<<"NULL PTR ("<<p0.x()<<","<<p0.y()<<")  "<<"("<<p1.x()<<","<<p1.y()<<")  "<<"("<<p2.x()<<","<<p2.y()<<")\n";

	}

	return og_points==new_points;
}

// Function to compute the distance between two points
double compute_distance(const Point& p1, const Point p2) {
    // Compute the squared distance
    Kernel::FT sq_dist = CGAL::squared_distance(p1, p2);

    // Return the square root of the squared distance for the actual distance
    return std::sqrt(CGAL::to_double(sq_dist));
}

//---------------------------------------------------------------------------------------------------------//
//--Heuristic Functions---//
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
//-----------------------------------------------------------------------------------------------------//

//-----------------------------Steiner Functions------------------------------------------------------//

Face_handle ant_longest_edge(Custom_CDT& cdt, Face_handle face, Polygon_2 region){

	
	Point p0 = face->vertex(0)->point();
	Point p1 = face->vertex(1)->point();
	Point p2 = face->vertex(2)->point();
	Face_handle neighbor;
	// Calculate the squared lengths of the triangle's sides

	Kernel::FT side0 = CGAL::squared_distance(p1, p2);
	Kernel::FT side1 = CGAL::squared_distance(p2, p0);
	Kernel::FT side2 = CGAL::squared_distance(p0, p1);
	// Determine the longest edge and calculate its midpoint
	Point midpoint;
	if( (side0>side2) && (side0>side1) ){
		midpoint = CGAL::midpoint(p1, p2);
		neighbor=face->neighbor(0);
	}
	else if ( (side1>side2) && (side1>side0)){
		midpoint= CGAL::midpoint(p2,p0);
		neighbor=face->neighbor(1);
	}
	else if ( (side2>side1) && (side2>side0)){
		midpoint=CGAL::midpoint(p0,p1);
		neighbor=face->neighbor(2);
	}
	else{
		std::cerr<<"Could determine longest edge between: "<<side1<<" "<<side2<<" "<<side0<<std::endl;
	}
	// Insert the Steiner point without flipping
	if(cdt.is_infinite(neighbor)||!(is_in_region_polygon(neighbor, region))){
		//we dont care about the neighbor since its outside the region boundary
		Face_handle ignore=face;
		cdt.insert_no_flip(midpoint,face);
		return ignore;
	}
	cdt.insert_no_flip(midpoint,face);
	return neighbor;

}

Face_handle ant_projection(Custom_CDT& cdt, Face_handle F, Polygon_2 region){
	//get the vertices of the face
	Point v0=F->vertex(0)->point();
	Point v1=F->vertex(1)->point();
	Point v2=F->vertex(2)->point();

	//calculate the length of the edges
	Kernel::FT length0 = CGAL::squared_distance(v1, v2);
	Kernel::FT length1 = CGAL::squared_distance(v2, v0);
	Kernel::FT length2 = CGAL::squared_distance(v0, v1);

	//FInd the longest edge and the oposite vertex
	Face_handle neighbor;
	Point le_start,le_end,obtuse_vertex;
	if( (length0 > length2) && (length0>length1)){
		le_start=v1;
		le_end=v2;
		obtuse_vertex=v0;
		neighbor=F->neighbor(0);
	}
	else if ((length1>length2)&& (length1>length0) ) {
		le_start=v2;
		le_end=v0;
		obtuse_vertex=v1;
		neighbor=F->neighbor(1);
	}
	else if ((length2>length1) &&(length2>length0) ) {
		le_start=v1;
		le_end=v0;
		obtuse_vertex=v2;
		neighbor=F->neighbor(2);
	
	}
	Point p0 = neighbor->vertex(0)->point();
	Point p1 = neighbor->vertex(1)->point();
	Point p2 = neighbor->vertex(2)->point();

	//calculate the projection of the obtuse vertex on the longest edge
	Kernel::Line_2 longest_edge(le_start,le_end);
	Point projection= longest_edge.projection(obtuse_vertex);
	if(!(is_in_region_polygon(neighbor, region))||cdt.is_infinite(neighbor)){
		Face_handle ignore=F;
		cdt.insert_no_flip(projection,F);
	//	std::cout<<"Will ignore\n";
		return ignore;
	}
	cdt.insert_no_flip(projection,F);
	return neighbor;
}


Face_handle ant_circumcenter(Custom_CDT& cdt, Face_handle face, const Polygon_2& boundary){

	// Get the vertices of the face
	int i=find_longest_edge_index(face);
	Point obP0=face->vertex(i)->point();
	Point sP1=face->vertex((i+1)%3)->point();
	Point sP2=face->vertex((i+2)%3)->point();

	// Calculate the circumcenter of the triangle formed by p0, p1, and p2
	Point circumcenter = CGAL::circumcenter(obP0, sP1, sP2);
	if(boundary.bounded_side(circumcenter)==CGAL::ON_UNBOUNDED_SIDE){
		Point centroid=steiner_centroid(cdt,face);
		cdt.insert_no_flip(centroid,face);
		return face;
	}
	if(is_face_constrained(cdt, face)){
		Point centroid=steiner_centroid(cdt,face);
		cdt.insert_no_flip(centroid,face);
		return face;
	}
	

	Face_handle neighbor=face->neighbor(i);
	Point nonShared=neighbor->vertex((neighbor->index(face)))->point();
	Point n0=neighbor->vertex(0)->point();
	Point n1=neighbor->vertex(1)->point();
	Point n2=neighbor->vertex(2)->point();

	//point not in neighbor
	if(!point_inside_triangle(n0,n1,n2,circumcenter)){
		Point centroid=steiner_centroid(cdt,face);
		cdt.insert_no_flip(centroid,face);
		return face;
	}
	cdt.remove_no_flip( (cdt.insert_no_flip(obP0)) );
	cdt.remove_no_flip( (cdt.insert_no_flip(sP1)) );
	cdt.remove_no_flip( (cdt.insert_no_flip(sP2)) );

	cdt.insert_points_and_constraint_no_flip(obP0, circumcenter);
	cdt.insert_no_flip(sP1);
	cdt.insert_no_flip(sP2);

	Vertex_handle vh0=cdt.insert_no_flip(obP0);
	Vertex_handle vhC=cdt.insert_no_flip(circumcenter);
	Face_handle face_p0c;
	int edge_indx;
	if(!(find_face(cdt,vh0,vhC,edge_indx, face_p0c))){
		std::cout<<"Couldnt locate vh1 vh3 face\n";
		exit(-1);
	}
	cdt.remove_constraint_no_flip(face_p0c,edge_indx);
	return neighbor;

};
Face_handle ant_merge(Custom_CDT& cdt,Face_handle face,Polygon_2& region_polygon){
	int test_obtuse_count;
	int index=-1;
	simulate_merge_steiner(cdt, face, region_polygon, index);
	if(index==-1){
		Point centroid=steiner_centroid(cdt,face);
		cdt.insert_no_flip(centroid,face);
		return face;
	}
	auto neighbor=face->neighbor(index);
	//skip  constrained ,non obtuse, outside of region_boundary neighbors
	if( cdt.is_infinite(neighbor) ){
		Point centroid=steiner_centroid(cdt,face);
		cdt.insert_no_flip(centroid,face);
		return face;
	}
	
	if( !(is_obtuse_triangle(neighbor)) ){
		Point centroid=steiner_centroid(cdt,face);
		cdt.insert_no_flip(centroid,face);
		return face;
	}
	if( !(is_in_region_polygon(neighbor, region_polygon)) ){
		Point centroid=steiner_centroid(cdt,face);
		cdt.insert_no_flip(centroid,face);
		return face;
	}
	if(is_face_constrained(cdt,neighbor)){
		Point centroid=steiner_centroid(cdt,face);
		cdt.insert_no_flip(centroid,face);
		return face;
	}

	Point v1 = face->vertex((index + 1) % 3)->point();
	Point v2 = face->vertex((index + 2) % 3)->point();
	Point v3 = neighbor->vertex(neighbor->index(face))->point();
	Point v4 = face->vertex(index)->point();

	if( !( is_polygon_convex(v1, v2, v3,v4) ) ){
		Point centroid=steiner_centroid(cdt,face);
		cdt.insert_no_flip(centroid,face);
		return face;
	}
	Point centroid=CGAL::centroid(v1,v2,v3,v4);

	//remove the points
	//v1
	cdt.remove_no_flip(face->vertex((index+1) % 3));
	//v2
	cdt.remove_no_flip((cdt.insert_no_flip(v2)));
	//v3
	cdt.remove_no_flip((cdt.insert_no_flip(v3)));
	//v4
	cdt.remove_no_flip((cdt.insert_no_flip(v4)));

	//add the centroid insert_no_flip
	auto vh_cent=cdt.insert_no_flip(centroid);
	
	//add back v1 v2 v3 v4 insert_no_flip
	//make (v1,v3) (v3,v2) (v2,v4) (v4,v1) constrains
	auto vh1=cdt.insert_no_flip(v1);
	auto vh3=cdt.insert_no_flip(v3);
	cdt.insert_constraint_no_flip(vh1,vh3);
	auto vh2=cdt.insert_no_flip(v2);
	auto vh4=cdt.insert_no_flip(v4);
	cdt.insert_constraint_no_flip(vh2,vh4);
	cdt.insert_constraint_no_flip(vh4,vh1);
	cdt.insert_constraint_no_flip(vh3, vh2);
	
	
	//remove the constrains remove_constrain_no_flip
	int edge_indx;
	Face_handle face13;
	if(!(find_face(cdt,vh1,vh3,edge_indx, face13))){
		std::cout<<"Couldnt locate vh1 vh3 face\n";
		exit(-1);
	}
	cdt.remove_constraint_no_flip(face13,edge_indx);
	Face_handle face24;
	if(!find_face(cdt,vh2,vh4,edge_indx, face24)){
		std::cout<<"Couldnt locate vh2 vh4 face\n";
		exit(-1);
	}
	cdt.remove_constraint_no_flip(face24,edge_indx);
	Face_handle face41;
	if(!find_face(cdt,vh4,vh1,edge_indx, face41)){
		std::cout<<"Couldnt locate vh1 vh4 face\n";
		exit(-1);
	}
	cdt.remove_constraint_no_flip(face41,edge_indx);
	Face_handle face32;
	if(!find_face(cdt,vh3,vh2,edge_indx, face32)){
		std::cout<<"Couldnt locate vh2 vh3 face\n";
		exit(-1);
	}
	cdt.remove_constraint_no_flip(face32,edge_indx);

return neighbor;
}



//----------------------------------------------------------------------------------------------------//

std::random_device rd; // Seed for random number generator
std::mt19937 gen(rd()); // Mersenne Twister engine
//select random face
Face_handle select_random_obtuse_face(const Custom_CDT& cdt, Polygon_2 boundary){
	
	std::vector<Face_handle> obtuse_faces;

	for(auto face=cdt.finite_faces_begin(); face!=cdt.finite_faces_end();face++){
		if(is_obtuse_triangle(face) && is_in_region_polygon(face, boundary) ){
		     obtuse_faces.push_back(face);
		}
	}
	if(obtuse_faces.empty()){
		std::cerr<<"Error tried to select obtuse face when there are none\n";
		exit(-1);
	}
	std::uniform_int_distribution<> distr(0,obtuse_faces.size()-1);
	int random_index=distr(gen);
	Face_handle o_face = obtuse_faces[random_index];
	if( o_face == nullptr){
		std::cerr<<"Error: couldnt locate original ccdt face in the copy cccdt"<<std::endl;
	}
	return o_face;
}

//-----------------------------------------------------------------------------------------------------//
std::pair <Face_handle, Face_handle> ImproveTriangulation(Custom_CDT& cdt, Polygon_2 boundary, double alpha, double beta,double xi,double psi, int steiner_count, double& ant_energy , const std::vector<double>& pheromones ,int& selected_method ){
	
	Face_handle face=select_random_obtuse_face(cdt,boundary);
	
	//calculate heuristics
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
		probabilities[i] = std::pow(pheromones[i], xi) * std::pow(heuristics[i], psi);
		denominator += probabilities[i];
	}
	//normalize the probabilities
	for (double& p :probabilities){
		p/=denominator;
	}

	std::discrete_distribution<> distr(probabilities.begin(),probabilities.end());
	selected_method=distr(gen);

	Custom_CDT cdt_copy(cdt);
	Face_handle face_copy=find_face(cdt_copy,face);
	Face_handle neighbor;

	switch (selected_method){
		case 0: {
			//make copy, add the steiner, save the edited face and neighbor . What if neigbor doesn't exist or is outside. Do i care?
			neighbor=ant_projection(cdt_copy, face_copy,boundary);
			break;
		}
		case 1:{
			//make copy, add the steiner, save the edited face and neighbor . What if neigbor doesn't exist or is outside. Do i care?
			neighbor=ant_longest_edge(cdt_copy, face_copy,boundary);
			break;
		}
		case 2:{
			ant_circumcenter(cdt_copy, face_copy,boundary);
			break;
		}
		case 3:{
			neighbor=ant_merge(cdt_copy, face_copy, boundary);
			break;
		}
		default:
			std::cerr<<"Error: chose invalid method\n";
			exit(-1);
	}

	int obtuse_faces=count_obtuse_faces(cdt_copy,boundary);

	if(obtuse_faces==0){
		ant_energy=0;
	}
	else{
		ant_energy=alpha*obtuse_faces +beta*(steiner_count+1);
	}
	return {face,neighbor};

}


double decrease_pheromon(double t, double lambda){
	return ((1-lambda)*t);
}
double increase_pheromon(double t,double lambda, double ant_energy){
	double denom=1+ant_energy;
	return (1-lambda)*t + 1/denom;
}



void ant_colony(Custom_CDT& cdt,Polygon_2 boundary,double alpha ,double beta, double xi, double psi, double lambda, int kappa, int L, int steiner_points){

	int obt_count=count_obtuse_faces(cdt,boundary);
	double best_Energy=alpha*obt_count+beta*steiner_points;

	std::vector<double> pheromones(4,0.5);
	std::vector<double> ant_energy(kappa,best_Energy);
	std::vector<int> ant_method(kappa,-1);

	Face_handle c_ant_face,c_ant_neighbor;
	double current_ant_energy;
	int current_ant_index;
	int cycle_steiners_added;
	for(int t=0;t<L;t++){
		cycle_steiners_added=0;
		//vector for each ant with the face it worked on and neighbor it changed (potentialy)
		std::vector<std::pair<Face_handle, Face_handle>> ant_faces(kappa);
		//each ant chooses a random method, applies it to its own cdt and returns the results
		for(int ant=0;ant<kappa;ant++){
			ant_faces[ant]=ImproveTriangulation(cdt, boundary, alpha, beta, xi, psi, steiner_points, ant_energy[ant], pheromones, ant_method[ant]);
		}
		//vector with [energy,ant_index] that will sort by energy
		std::vector<std::pair<double,int>> energy_index_pairs;
		for(int i=0;i<kappa;i++){
			energy_index_pairs.emplace_back(ant_energy[i],i);
		}
		std::sort(energy_index_pairs.begin(),energy_index_pairs.end());
		//for each ant
		for(int ant=0;ant<kappa;ant++){
			current_ant_energy=energy_index_pairs[ant].first;
			current_ant_index=energy_index_pairs[ant].second;
			//if the ant reduces the energy of the previous cycle
			if(current_ant_energy<best_Energy){
				//if both of the faces that will change still exist
				if(face_still_exists(cdt,c_ant_face)){
					if(face_still_exists(cdt, c_ant_neighbor)){
						Face_handle apply_face=find_face(cdt,c_ant_face);
						switch (ant_method[current_ant_index]){
							case 0: {
								//std::cout<<"Selected projection edge\n";
								ant_projection(cdt, apply_face,boundary);
								break;
							}
							case 1:{
								std::cout<<"Selected longest edge\n";
								ant_longest_edge(cdt,apply_face,boundary);
								break;
							}
							case 2:{
								std::cout<<"Selected circumcenter\n";
								ant_circumcenter(cdt,apply_face,boundary);
								break;
							}
							case 3:{
								std::cout<<"Selected merge\n";
								ant_merge(cdt,apply_face, boundary);
								break;
							}
							default:
								std::cerr<<"Error: chose invalid method\n";
								exit(-1);
							}
						cycle_steiners_added++;
						//reinforce pheromon
						pheromones[ant_method[current_ant_index]]=increase_pheromon(t, lambda, current_ant_energy);
						
					}
					else{
						std::cout<<"Face exists but neighbor doesnt exist "<<ant_method[current_ant_index]<<std::endl;
						pheromones[ant_method[current_ant_index]]=decrease_pheromon( t, lambda);
					}
				}
				else{
					std::cout<<"Face  exist "<<ant_method[current_ant_index]<<std::endl;
					pheromones[ant_method[current_ant_index]]=decrease_pheromon( t, lambda);
				}
			}
			else{
				//didnt decrease energy
				pheromones[ant_method[current_ant_index]]=decrease_pheromon( t, lambda);
			}
		}
		//for each ant
		//if changes were made
		if(cycle_steiners_added){
			obt_count=count_obtuse_faces(cdt,boundary);
			steiner_points+=cycle_steiners_added;
			best_Energy=obt_count*alpha + steiner_points*beta;
		}
	}
}

