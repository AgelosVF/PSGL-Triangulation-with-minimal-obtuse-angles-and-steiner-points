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
	std::random_device rd; // Seed for random number generator
	std::mt19937 gen(rd()); // Mersenne Twister engine
	std::uniform_int_distribution<> dis(0, obtuse_faces.size() - 1);
	int random_index = dis(gen);
	Face_handle o_face = obtuse_faces[random_index];
	if( o_face == nullptr){
		std::cerr<<"Error: couldnt locate original ccdt face in the copy cccdt"<<std::endl;
	}
	return o_face;
}

Face_handle ant_longest_edge(Custom_CDT& ccdt, Face_handle face,Polygon_2 region) {
    // Get the vertices of the face
	// Get the three vertices of the triangle
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
	if(ccdt.is_infinite(neighbor)||!(is_in_region_polygon(neighbor, region))){
		Face_handle ignore=face;
		ccdt.insert_no_flip(midpoint,face);
		return ignore;
	}
	ccdt.insert_no_flip(midpoint,face);
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
	//std::cout<<"PROJECTION PTR ("<<p0.x()<<","<<p0.y()<<")  "<<"("<<p1.x()<<","<<p1.y()<<")  "<<"("<<p2.x()<<","<<p2.y()<<")\n";
	//std::cout<<"Is infinite: "<<cdt.is_infinite(neighbor)<<std::endl;

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
	Point p0 = face->vertex(0)->point();
	Point p1 = face->vertex(1)->point();
	Point p2 = face->vertex(2)->point();

	// Calculate the circumcenter of the triangle formed by p0, p1, and p2
	Point circumcenter = CGAL::circumcenter(p0, p1, p2);
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
	
	int i=find_longest_edge_index(face);
	Point obP0=face->vertex(i)->point();
	Point sP1=face->vertex((i+1)%3)->point();
	Point sP2=face->vertex((i+2)%3)->point();

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

static int reduced=0;
//if the obtuse count didnt get reduced
double update_pheromon_0(double t, double lambda){
	if(reduced<5)
		std::cout<<"Reduced pheromon from "<<t<<" to "<< (1-lambda)*t<<std::endl;
	reduced++;
	return ( (1-lambda) * t );

}
//if the obtuse count got reduced
double update_pheromon_1(double t, double lambda,double alpha, int obtuse_count, double beta, int steiner_count){
	double denom=1+alpha*obtuse_count + beta* steiner_count;
	std::cout<<"Reinforce pheromon from "<<t<<" to "<<(1-lambda)*t + (1/denom)<<std::endl;
	return ( (1-lambda) * t + (1/denom) );

}



std::pair<Face_handle, Face_handle> impruveTriangulation(Custom_CDT& cdt,Polygon_2 boundary,double alpha,double beta,int steiner_count,double xi,double ps, const std::vector<double>& pheromones,double& ant_energy,int& selected_method,double best_E){

	//choose random face
	Face_handle face=select_random_obtuse_face(cdt, boundary);
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
	selected_method=0;
	
	Custom_CDT cdt_copy(cdt);
	Face_handle face_copy=find_face(cdt_copy,face);
	Face_handle neighbor;
	int obtuse_faces=99999;
	switch (selected_method){
		case 0: {
			//make copy, add the steiner, save the edited face and neighbor . What if neigbor doesn't exist or is outside. Do i care?
			neighbor=ant_projection(cdt_copy, face_copy,boundary);
			obtuse_faces=count_obtuse_faces(cdt_copy,boundary);
			break;
		}
		case 1:{
			//make copy, add the steiner, save the edited face and neighbor . What if neigbor doesn't exist or is outside. Do i care?
			neighbor=ant_longest_edge(cdt_copy, face_copy,boundary);
			obtuse_faces=count_obtuse_faces(cdt_copy,boundary);
			break;
		}
		case 2:{
			ant_circumcenter(cdt_copy, face_copy,boundary);
			obtuse_faces=count_obtuse_faces(cdt_copy,boundary);
			break;
		}
		case 3:{
			neighbor=ant_merge(cdt_copy, face_copy, boundary);
			obtuse_faces=count_obtuse_faces(cdt_copy,boundary);
			break;
		}
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
	return {face,neighbor};
}

bool face_still_exists(Custom_CDT cdt,Face_handle face){
	
	Point p0 = face->vertex(0)->point();
	Point p1 = face->vertex(1)->point();
	Point p2 = face->vertex(2)->point();
	Point centroid= CGAL::centroid(p0,p1,p2);
	Face_handle n_face = cdt.locate(centroid);
	if( n_face == nullptr){
		std::cout<<"NULL PTR ("<<p0.x()<<","<<p0.y()<<")  "<<"("<<p1.x()<<","<<p1.y()<<")  "<<"("<<p2.x()<<","<<p2.y()<<")\n";
		return false;
	}
	Point np0=n_face->vertex(0)->point();
	Point np1=n_face->vertex(1)->point();
	Point np2=n_face->vertex(2)->point();

	std::set<Point> og_points={p0,p1,p2};
	std::set<Point> new_points={np0,np1,np2};
	std::cout<<"\tNEW PTR ("<<np0.x()<<","<<np0.y()<<")  "<<"("<<np1.x()<<","<<np1.y()<<")  "<<"("<<np2.x()<<","<<np2.y()<<")\n";
	std::cout<<"\tNULL PTR ("<<p0.x()<<","<<p0.y()<<")  "<<"("<<p1.x()<<","<<p1.y()<<")  "<<"("<<p2.x()<<","<<p2.y()<<")\n";
	if(!(og_points==new_points)){
		
		std::cout<<"NEW PTR ("<<np0.x()<<","<<np0.y()<<")  "<<"("<<np1.x()<<","<<np1.y()<<")  "<<"("<<np2.x()<<","<<np2.y()<<")\n";
		std::cout<<"NULL PTR ("<<p0.x()<<","<<p0.y()<<")  "<<"("<<p1.x()<<","<<p1.y()<<")  "<<"("<<p2.x()<<","<<p2.y()<<")\n";

	}

	return og_points==new_points;
}

void ant_colony(Custom_CDT& cdt,Polygon_2 boundary,double alpha ,double beta, double xi, double ps, double lambda, int kappa, int L, int steiner_points){

	std ::vector<double> pheromones(4,0.5);
	int current_steiner=steiner_points;
	int temp_steiner=steiner_points;
	int obt_count=count_obtuse_faces(cdt,boundary);
	int temp_obt_count=obt_count;
	double best_Energy=alpha*obt_count + beta*steiner_points;
	std::vector<double> ant_energy(kappa,best_Energy);
	std::vector<int> ant_method(kappa,-1);
	int merge=0,circ=0,proj=0,longE=0;
	int Tmerge=0,Tcirc=0,Tproj=0,TlongE=0;
	Face_handle c_ant_face,c_ant_neigbor;
	for(int t=0;t<L;t++){
		//std::cout<<"Starting energy: "<<best_Energy<<std::endl;
		//vector that will have each face each ant worked on
		std::vector<std::pair<Face_handle, Face_handle>> ant_faces(kappa);
		for(int ant=0;ant<kappa;ant++){
			//each ant should have its own cdt
			ant_faces[ant]= impruveTriangulation( cdt, boundary, alpha, beta, current_steiner, xi, ps, pheromones, ant_energy[ant], ant_method[ant],best_Energy);
		}
		//extra can skip for now
		for(int ant=0;ant<kappa;ant++){
			if(ant_energy[ant]==0){
				std::cout<<"Ant reduced energy to 0\n";
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
		for(int ant=0;ant<kappa;ant++){
			double current_ant_energy=energy_index_pairs[ant].first;
			int current_ant_index=energy_index_pairs[ant].second;
			//std::cout<<"Current ant: "<<current_ant_index<<" : "<<current_ant_energy<<std::endl;
			c_ant_face=ant_faces[current_ant_index].first;
			c_ant_neigbor=ant_faces[current_ant_index].second;
			int c_meth=ant_method[current_ant_index];
			switch(c_meth){
				case 0:{
					Tproj++;
					break;
				}
				case 1:{
					TlongE++;
					break;
				}
				case 2:{
					Tcirc++;
					break;
				}
				case 3:{
					Tmerge++;
					break;
				}
				default:
					std::cerr<<"Error: chose invalid method\n";

			}
			if(current_ant_energy<best_Energy){
				std::cout<<"Better energy:" <<current_ant_energy<<std::endl;
				if(face_still_exists(cdt, c_ant_face)){
					if(face_still_exists(cdt, c_ant_neigbor)){
						Face_handle apply_face=find_face(cdt,c_ant_face);
						switch (ant_method[current_ant_index]){
							case 0: {
								//std::cout<<"Selected projection edge\n";
								ant_projection(cdt, apply_face,boundary);
								proj++;
								break;
							}
							case 1:{
								std::cout<<"Selected longest edge\n";
								ant_longest_edge(cdt,apply_face,boundary);
								longE++;
								break;
							}
							case 2:{
								std::cout<<"Selected circumcenter\n";
								ant_circumcenter(cdt,apply_face,boundary);
								circ++;
								break;
							}
							case 3:{
								std::cout<<"Selected merge\n";
								ant_merge(cdt,apply_face, boundary);
								merge++;
								break;
							}
							default:
								std::cerr<<"Error: chose invalid method\n";
								exit(-1);
						}
						current_steiner++;
						temp_obt_count=count_obtuse_faces(cdt,boundary);
						if(temp_obt_count<obt_count){
							obt_count=temp_obt_count;
							std::cout<<"Current ant index:"<<current_ant_index<<std::endl;
							pheromones[ant_method[current_ant_index]]=update_pheromon_1( pheromones[ant_method[current_ant_index]], lambda, alpha, obt_count, beta, current_steiner);
						}
						else{
							//we changed the obt_count to temp obt count before we should
							std::cout<<"NEVER GONNA GIVE YOU UP\n";
							obt_count=temp_obt_count;
							pheromones[ant_method[current_ant_index]]=update_pheromon_0( pheromones[ant_method[current_ant_index]],  lambda);
						}
					}
					else{
						std::cout<<"Neighbor doesnt exist "<<ant_method[current_ant_index]<<std::endl;
						pheromones[ant_method[current_ant_index]]=update_pheromon_0( pheromones[ant_method[current_ant_index]],  lambda);

					}
				}
				else{
					std::cout<<"Face doesnt exist\n";
				}
			}
			else{
				//std::cout<<"Didnt improve energy\n";
				pheromones[ant_method[current_ant_index]]=update_pheromon_0( pheromones[ant_method[current_ant_index]],  lambda);
			}
		}
		obt_count=count_obtuse_faces(cdt,boundary);
		best_Energy=obt_count*alpha+current_steiner*beta;
		for(int i=0;i<kappa;i++){
			std::cout<<"\t"<<energy_index_pairs[i].first<<" "<<energy_index_pairs[i].second<<std::endl;
		}
		std::cout<<std::endl;
	
	}
	std::cout<<"Pheromones: \n";
	for(int i=0;i<4;i++){
		std::cout<<"\t"<<i<<": "<<pheromones[i]<<std::endl;
	}
	std::cout<<"Points used: "<<"\n\t projections:"<<proj<<"\n\tLongest edge:"<<longE<<"\n\tCircumcenter:"<<circ<<"\n\tMerge:"<<merge<<std::endl;
	std::cout<<"Points tested: "<<"\n\t projections:"<<Tproj<<"\n\tLongest edge:"<<TlongE<<"\n\tCircumcenter:"<<Tcirc<<"\n\tMerge:"<<Tmerge<<std::endl;
	std::cout<<"Final energy: "<<best_Energy<<std::endl;
}

