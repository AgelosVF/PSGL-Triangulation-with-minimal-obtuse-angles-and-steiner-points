#include "../Header Files/CGAL_CUSTOM_CONSTRAINED_DELAUNAY_TRIANGULATION_2.h"
#include "../Header Files/triangulation_utils.hpp"
#include "../Header Files/SteinerPoints.hpp"
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Kernel/global_functions_3.h>
#include <CGAL/enum.h>
#include <CGAL/squared_distance_2.h>
#include <algorithm>
#include <cmath> 
#include <iostream>
#include <ostream>
#include <random>
#include <strings.h>
#include <thread>


//-------------------------------Helping functions-----------------------------------------------------//

bool same_face(Face_handle face1, Face_handle face2){
	if( face1==nullptr){
		std::cerr<<"Got null face for face 1\n";
		return false;
	}
	else if( face2==nullptr){
		std::cerr<<"Got null face for face 2\n";
		return false;
	}
	Point p0=face1->vertex(0)->point();
	Point p1=face1->vertex(1)->point();
	Point p2=face1->vertex(2)->point();

	Point np0=face2->vertex(0)->point();
	Point np1=face2->vertex(1)->point();
	Point np2=face2->vertex(2)->point();

	std::set<Point> og_points={p0,p1,p2};
	std::set<Point> new_points={np0,np1,np2};
	
	return og_points==new_points;
}
// Helper function to mark duplicate faces
void resolve_conflicts(std::vector<std::pair<double, int>>& energy_index_pairs, std::vector<std::pair<Face_handle, Face_handle>>& ant_faces) {
	for(size_t i = 0; i < energy_index_pairs.size(); i++) {
		// Skip already marked faces
		if(energy_index_pairs[i].second == -1) continue;

		int current_ant_index = energy_index_pairs[i].second;
		Face_handle current_face = ant_faces[current_ant_index].first;
		Face_handle current_neighbor = ant_faces[current_ant_index].second;

		// Compare with all subsequent faces
		for(size_t j = i + 1; j < energy_index_pairs.size(); j++) {
			if(energy_index_pairs[j].second == -1) continue;

			int compare_ant_index = energy_index_pairs[j].second;
			Face_handle compare_face = ant_faces[compare_ant_index].first;
			Face_handle compare_neighbor = ant_faces[compare_ant_index].second;


			// Check if either the face or neighbor matches
			if(same_face(current_face, compare_face)){
				//std::cout<<"Same face\n";
				energy_index_pairs[j].second = -1;
			}
			else if(compare_neighbor != nullptr){
				//std::cout<<"compare neighbor exists\n";

				if(same_face(current_face, compare_neighbor)){
					//std::cout<<"Same face and neighbor\n";
					energy_index_pairs[j].second = -1;
					}

			}
			else if(current_neighbor!=nullptr){
				//std::cout<<"Current neighbor exists\n";
				if(same_face(current_neighbor, compare_face)){
					energy_index_pairs[j].second=-1;
					//std::cout<<"Same neighbor and face\n";
				}
				if(compare_neighbor!=nullptr){
					if( same_face(current_neighbor,compare_neighbor)){
						//std::cout<<"Same neighbor\n";
						energy_index_pairs[j].second = -1;
					}
				}
			}
		}
	}
}



//---------------------------------------------------------------------------------------------------------//
//--Heuristic Functions---//
double calculate_r(const Face_handle& Face){

	Point p0=Face->vertex(0)->point();
	Point p1=Face->vertex(1)->point();
	Point p2=Face->vertex(2)->point();

	Point circumcenter=CGAL::circumcenter(p0,p1,p2);

	double R=compute_distance(p0,circumcenter);

	int l_side_i=find_obtuse_vertex_index(Face);
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


Face_handle ant_circumcenter(Custom_CDT& cdt, Face_handle& face, const Polygon_2& boundary){

	int obtuse;

	Point circ=steiner_circumcenter(cdt,face);
	if( boundary.bounded_side(circ)!=CGAL::ON_BOUNDED_SIDE && boundary.bounded_side(circ)!=CGAL::ON_BOUNDARY){
		Point centroid=steiner_centroid(cdt,face);
		cdt.insert_no_flip(centroid,face);
		return face;
	}


	//get the vertices of the face
	Point p0=face->vertex(0)->point();
	Point p1=face->vertex(1)->point();
	Point p2=face->vertex(2)->point();

	//calculate the length of the edges
	//
	Kernel::FT length1 = CGAL::squared_distance(p0, p1);
	Kernel::FT length2 = CGAL::squared_distance(p1, p2);
	Kernel::FT length3 = CGAL::squared_distance(p2, p0);

	//FInd the longest edge and the oposite vertex
	Point le_start,le_end,obtuse_point;
	if( (length1 > length2) && (length1>length3)){
		le_start=p0;
		le_end=p1;
		obtuse_point=p2;
	}
	else if ((length2>length1)&& (length2>length3) ) {
		le_start=p1;
		le_end=p2;
		obtuse_point=p0;
	}
	else if ((length3>length1) &&(length3>length2) ) {
		le_start=p2;
		le_end=p0;
		obtuse_point=p1;
	
	}


	Vertex_handle obtuse_vertex=cdt.insert_no_flip(obtuse_point,face);
	Face_handle neighbor=face->neighbor(face->index(obtuse_vertex));
	

	Point n1=neighbor->vertex(0)->point();
	Point n2=neighbor->vertex(1)->point();
	Point n3=neighbor->vertex(2)->point();
	if(!point_inside_triangle(n1, n2, n3, circ)){
		Point centroid=steiner_centroid(cdt,face);
		cdt.insert_no_flip(centroid,face);
		return face;
	}
	
	//neighbor i has the same index as the opposite point
	Point nonShared=neighbor->vertex(neighbor->index(face))->point();
	//check if the common edge is constrained
	if(face->is_constrained(face->index(neighbor))){
		Point centroid=steiner_centroid(cdt,face);
		cdt.insert_no_flip(centroid,face);
		return face;
	}
	Vertex_handle circ_vh=cdt.insert_no_flip(circ,neighbor);
	//make the obtuse point circ a constain to remove the edge in between
	obtuse_vertex=cdt.insert_no_flip(obtuse_point);
	cdt.insert_constraint_no_flip(circ_vh, obtuse_vertex);
	//after remove the constrain without flipping
	obtuse_vertex=cdt.insert_no_flip(obtuse_point);
	circ_vh=cdt.insert_no_flip(circ,neighbor);
	Face_handle temp_face;
	int edge_index=-1;
	cdt.is_edge(obtuse_vertex, circ_vh, temp_face, edge_index);
	cdt.remove_constraint_no_flip(temp_face, edge_index);
	return neighbor;
};

Face_handle ant_merge(Custom_CDT& cdt,Face_handle face,Polygon_2& region_polygon){
	int test_obtuse_count;
	bool found=false;
	Point neighbor_point;
	simulate_merge_steiner(cdt, face, region_polygon,neighbor_point,found );
	if(found==false){
		Point centroid=steiner_centroid(cdt,face);
		cdt.insert_no_flip(centroid,face);
		return face;
	}
	int index=-1;
	for(int i=0;i<3;i++){
		Point test=face->vertex(i)->point();
		if(test==neighbor_point){
			index=i;
			found=true;
			break;
		}
	}
	if(index==-1){
		std::cout<<"Ant merge didnt find index\n";
		exit(-1);
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
	Face_handle neighbor_copy;

	switch (selected_method){
		case 0: {
			//make copy, add the steiner, save the edited face and neighbor . What if neigbor doesn't exist or is outside. Do i care?
			neighbor_copy=ant_projection(cdt_copy, face_copy,boundary);
			break;
		}
		case 1:{
			//make copy, add the steiner, save the edited face and neighbor_copy . What if neigbor doesn't exist or is outside. Do i care?
			neighbor_copy=ant_longest_edge(cdt_copy, face_copy,boundary);
			break;
		}
		case 2:{
			neighbor_copy=ant_circumcenter(cdt_copy, face_copy,boundary);
			break;
		}
		case 3:{
			neighbor_copy=ant_merge(cdt_copy, face_copy, boundary);
			break;
		}
		default:
			std::cerr<<"Error: chose invalid method\n";
			exit(-1);
	}

	int obtuse_faces=count_obtuse_faces(cdt_copy,boundary);
	Face_handle neighbor=find_face(cdt,neighbor_copy);

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


int ant_colony(Custom_CDT& cdt, Polygon_2 boundary, double alpha, double beta, double xi, double psi, double lambda, int kappa, int L, int steiner_points,double& c_rate) {
	int obt_count = count_obtuse_faces(cdt, boundary);
	int prev_count=obt_count;
	if(obt_count==0){
		return steiner_points;
	}
	double best_Energy = alpha * obt_count + beta * steiner_points;

	std::vector<double> pheromones(4, 0.5);
	std::vector<double> ant_energy(kappa, best_Energy);
	std::vector<int> ant_method(kappa, -1);

	Face_handle c_ant_face, c_ant_neighbor;
	double current_ant_energy;
	int current_ant_index;
	int cycle_steiners_added;
	int merge_used = 0, proj_used = 0, longest_used = 0, circum_used = 0;

	// Function to be executed by each thread
	auto ant_task = [&](int ant_id, std::vector<std::pair<Face_handle, Face_handle>>& ant_faces) {
		// These variables are unique per ant/thread, so no mutex needed
		ant_faces[ant_id] = ImproveTriangulation(cdt, boundary, alpha, beta, xi, psi, steiner_points, ant_energy[ant_id], pheromones, ant_method[ant_id]);
	};

	for (int t = 0; t < L; t++) {
		cycle_steiners_added = 0;
		std::vector<std::pair<Face_handle, Face_handle>> ant_faces(kappa);
		std::vector<std::thread> threads;

		// Launch threads for each ant
		for (int ant = 0; ant < kappa; ant++) {
			threads.emplace_back(ant_task, ant, std::ref(ant_faces));
		}

		// Wait for all threads to complete
		for (auto& thread : threads) {
			thread.join();
		}

		// Sort results by energy
		std::vector<std::pair<double, int>> energy_index_pairs;
		for (int i = 0; i < kappa; i++) {
			energy_index_pairs.emplace_back(ant_energy[i], i);
		}
		std::sort(energy_index_pairs.begin(), energy_index_pairs.end());

		// Resolve conflicts
		resolve_conflicts(energy_index_pairs, ant_faces);

		// Apply changes from successful ants
		for (int ant = 0; ant < kappa; ant++) {
			if (energy_index_pairs[ant].second == -1) continue;

			current_ant_energy = energy_index_pairs[ant].first;
			current_ant_index = energy_index_pairs[ant].second;

			if (current_ant_energy < best_Energy) {
				c_ant_face = ant_faces[current_ant_index].first;
				c_ant_neighbor = ant_faces[current_ant_index].second;
				Face_handle apply_face = find_face(cdt, c_ant_face);

				switch (ant_method[current_ant_index]) {
				case 0: {
					ant_projection(cdt, apply_face, boundary);
					proj_used++;
					obt_count = count_obtuse_faces(cdt, boundary);
					if(steiner_points + cycle_steiners_added > 1) {
						c_rate += calculate_convergence_rate(steiner_points + cycle_steiners_added, prev_count, obt_count);
					}
					prev_count = obt_count;
					cycle_steiners_added++;
					break;
				}
				case 1: {
					ant_longest_edge(cdt, apply_face, boundary);
					longest_used++;
					obt_count = count_obtuse_faces(cdt, boundary);
					if(steiner_points + cycle_steiners_added > 1) {
						c_rate += calculate_convergence_rate(steiner_points + cycle_steiners_added, prev_count, obt_count);
					}
					prev_count = obt_count;
					cycle_steiners_added++;
					break;
				}
				case 2: {
					ant_circumcenter(cdt, apply_face, boundary);
					circum_used++;
					obt_count = count_obtuse_faces(cdt, boundary);
					if(steiner_points + cycle_steiners_added > 1) {
						c_rate += calculate_convergence_rate(steiner_points + cycle_steiners_added, prev_count, obt_count);
					}
					prev_count = obt_count;
					cycle_steiners_added++;
					break;
				}
				case 3: {
					ant_merge(cdt, apply_face, boundary);
					merge_used++;
					obt_count = count_obtuse_faces(cdt, boundary);
					if(steiner_points + cycle_steiners_added > 1) {
						c_rate += calculate_convergence_rate(steiner_points + cycle_steiners_added, prev_count, obt_count);
					}
					prev_count = obt_count;
					cycle_steiners_added++;
					break;
				}
				default:
					std::cerr << "Error: chose invalid method\n";
					exit(-1);
				}

				if (current_ant_energy == 0) {
					obt_count = count_obtuse_faces(cdt, boundary);
					if (obt_count != 0) {
						std::cerr << "Ant energy is 0 but there are still obtuse faces in the cdt\n";
					}
					steiner_points += cycle_steiners_added;
					best_Energy = obt_count * alpha + steiner_points * beta;
					std::cout << "Used:\n\tMerges:" << merge_used 
					<< "\n\tCircumcenter:" << circum_used 
					<< "\n\tLongest edge:" << longest_used 
					<< "\n\tProjections:" << proj_used << std::endl;
					return steiner_points;
				}

				pheromones[ant_method[current_ant_index]] = 
				increase_pheromon(t, lambda, current_ant_energy);
			}
			else {
				pheromones[ant_method[current_ant_index]] =decrease_pheromon(t, lambda);
			}
		}

		if (cycle_steiners_added) {

			obt_count = count_obtuse_faces(cdt, boundary);
			steiner_points += cycle_steiners_added;
			best_Energy = obt_count * alpha + steiner_points * beta;
		}
	}

	std::cout<< "Used:\n\tMerges:" << merge_used 
		<< "\n\tCircumcenter:" << circum_used 
		<< "\n\tLongest edge:" << longest_used 
		<< "\n\tProjections:" << proj_used << std::endl;
	return steiner_points;
}

