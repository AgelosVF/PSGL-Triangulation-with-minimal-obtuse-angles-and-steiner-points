#include <CGAL/Qt/Basic_viewer_qt.h>
#include <CGAL/Uncertain.h>
#include <boost/property_map/property_map.hpp>
#include <cmath>
#include <cstddef>
#include <cstdlib>
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
#include "../Header Files/SteinerPoints.hpp"
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

//-----------------------------------------------------------------------------------------------------------------------------------------//
int annealing_test_add_steiner_centroid(Custom_CDT& ccdt,Face_handle face, const Polygon_2 & region_boundary,int testORadd){
	
	int og_obtuse_count=count_obtuse_faces(ccdt,region_boundary);
	//create a copy of the ccdt and add the steiner there to check if the obtuse triangles are reduces
	Custom_CDT test_ccdt(ccdt);
	Face_handle test_face=find_face(test_ccdt, face);
	Point test_steiner=steiner_centroid(test_ccdt, test_face);
	
	test_ccdt.insert_no_flip(test_steiner);

	int test_obutse_count=count_obtuse_faces(test_ccdt, region_boundary);

	if(testORadd==1){
		Point steiner=steiner_centroid(ccdt, face);
		ccdt.insert_no_flip(steiner);
		
	}

	return test_obutse_count;
}

int annealing_simulate_merge_steiner(Custom_CDT& ccdt,Face_handle ob_face,const Polygon_2& region_boundary,Point& neighbor_point, bool& found){
	
	found=false;
	//check if the vertices of the face are constrained
	if(is_face_constrained(ccdt,ob_face) ){
		int test_obtuse_count=annealing_test_add_steiner_centroid(ccdt, ob_face, region_boundary, 0);
		return test_obtuse_count;
	}
	//create a copy and locate the face in it
	Custom_CDT ccdt_copy(ccdt);

	Face_handle face_copy;
	Face_handle neighbor;
	Point temp_point;
	int ob_count=annealing_test_add_steiner_centroid(ccdt, ob_face, region_boundary, 0);

	for(int i=0;i<3;i++){
		face_copy=find_face(ccdt_copy,ob_face);
		neighbor=face_copy->neighbor(i);
	temp_point=face_copy->vertex(i)->point();
		//skip  constrained ,non obtuse, outside of region_boundary neighbors
		if( ccdt_copy.is_infinite(neighbor) ){
			continue;
		}
		
		if( !(is_obtuse_triangle(neighbor)) ){
			continue;
		}
		if( !(is_in_region_polygon(face_copy, region_boundary)) ){
			continue;
		}
		if(is_face_constrained(ccdt_copy,neighbor)){
			continue;
		}

		Point v1 = face_copy->vertex((i + 1) % 3)->point();
		Point v2 = face_copy->vertex((i + 2) % 3)->point();
		Point v3 = neighbor->vertex(neighbor->index(face_copy))->point();
		Point v4 = face_copy->vertex(i)->point();

		if( !( is_polygon_convex(v1, v2, v3,v4) ) ){
			continue;
		}
		Point centroid=CGAL::centroid(v1,v2,v3,v4);

		//remove the points
		//v1
		ccdt_copy.remove_no_flip(face_copy->vertex((i+1) % 3));
		//v2
		ccdt_copy.remove_no_flip((ccdt_copy.insert_no_flip(v2)));
		//v3
		ccdt_copy.remove_no_flip((ccdt_copy.insert_no_flip(v3)));
		//v4
		ccdt_copy.remove_no_flip((ccdt_copy.insert_no_flip(v4)));

		//add the centroid insert_no_flip
		auto vh_cent=ccdt_copy.insert_no_flip(centroid);
		
		//add back v1 v2 v3 v4 insert_no_flip
		//make (v1,v3) (v3,v2) (v2,v4) (v4,v1) constrains
		auto vh1=ccdt_copy.insert_no_flip(v1);
		auto vh3=ccdt_copy.insert_no_flip(v3);
		ccdt_copy.insert_constraint_no_flip(vh1,vh3);
		auto vh2=ccdt_copy.insert_no_flip(v2);
		auto vh4=ccdt_copy.insert_no_flip(v4);
		ccdt_copy.insert_constraint_no_flip(vh2,vh4);
		ccdt_copy.insert_constraint_no_flip(vh4,vh1);
		ccdt_copy.insert_constraint_no_flip(vh3, vh2);
		
		
		//remove the constrains remove_constrain_no_flip
		int edge_indx;
		Face_handle face13;
		if(!(find_face(ccdt_copy,vh1,vh3,edge_indx, face13))){
			std::cout<<"Couldnt locate vh1 vh3 face\n";
			exit(-1);
		}
		ccdt_copy.remove_constraint_no_flip(face13,edge_indx);
		Face_handle face24;
		if(!find_face(ccdt_copy,vh2,vh4,edge_indx, face24)){
			std::cout<<"Couldnt locate vh2 vh4 face\n";
			exit(-1);
		}
		ccdt_copy.remove_constraint_no_flip(face24,edge_indx);
		Face_handle face41;
		if(!find_face(ccdt_copy,vh4,vh1,edge_indx, face41)){
			std::cout<<"Couldnt locate vh1 vh4 face\n";
			exit(-1);
		}
		ccdt_copy.remove_constraint_no_flip(face41,edge_indx);
		Face_handle face32;
		if(!find_face(ccdt_copy,vh3,vh2,edge_indx, face32)){
			std::cout<<"Couldnt locate vh2 vh3 face\n";
			exit(-1);
		}
		ccdt_copy.remove_constraint_no_flip(face32,edge_indx);

		//count obtuses and reset by removing the steiner, forcing the edge that was removed
		int temp_count=count_obtuse_faces(ccdt_copy, region_boundary);
		if(temp_count<ob_count){
			ob_count=temp_count;
			found=true;
			neighbor_point=temp_point;
		}
		ccdt_copy.remove_no_flip(vh_cent);


		vh1=ccdt_copy.insert_no_flip(v1);
		vh2=ccdt_copy.insert_no_flip(v2);
		ccdt_copy.insert_constraint_no_flip(vh1, vh2);
		Face_handle face12;
		if(!find_face(ccdt_copy,vh1,vh2,edge_indx, face12)){
			std::cout<<"Couldnt locate vh2 vh2 face\n";
			exit(-1);
		}
		ccdt_copy.remove_constraint_no_flip(face12,edge_indx);
		
	}
	return ob_count;
}
int annealing_test_add_steiner_merge(Custom_CDT& ccdt,Face_handle face,Polygon_2& region_polygon, bool&found, Point& neighbor_point,int testORadd){
	int test_obtuse_count;
	if(testORadd!=1){
		found=false;
		test_obtuse_count=annealing_simulate_merge_steiner(ccdt,face,region_polygon,neighbor_point,found);
		return test_obtuse_count;
	}
	

	//skip  constrained ,non obtuse, outside of region_boundary neighbors
	if(found==false){
		test_obtuse_count=annealing_test_add_steiner_centroid(ccdt, face, region_polygon, 1);
		return test_obtuse_count;

	}

	int indx;
	found=false;
	for(int i=0;i<3;i++){
		Point test=face->vertex(i)->point();
		if(test==neighbor_point){
			indx=i;
			found=true;
			break;
		}
	}
	if(!found){
		std::cout<<"didnt find neighbor\n";
	}
	auto face_2=face->neighbor(indx);
	if( ccdt.is_infinite(face_2) ){
		test_obtuse_count=annealing_test_add_steiner_centroid(ccdt, face, region_polygon, 1);
		return test_obtuse_count;
	}
	
	if( !(is_obtuse_triangle(face_2)) ){
		test_obtuse_count=annealing_test_add_steiner_centroid(ccdt, face, region_polygon, 1);
		return test_obtuse_count;
	}
	if( !(is_in_region_polygon(face_2, region_polygon)) ){
		test_obtuse_count=annealing_test_add_steiner_centroid(ccdt, face, region_polygon, 1);
		return test_obtuse_count;
	}
	if(is_face_constrained(ccdt,face_2)){
		test_obtuse_count=annealing_test_add_steiner_centroid(ccdt, face, region_polygon, 1);
		return test_obtuse_count;
	}

	Point v1 = face->vertex((indx + 1) % 3)->point();
	Point v2 = face->vertex((indx + 2) % 3)->point();
	Point v3 = face_2->vertex(face_2->index(face))->point();
	Point v4 = face->vertex(indx)->point();

	if( !( is_polygon_convex(v1, v2, v3,v4) ) ){
		test_obtuse_count=annealing_test_add_steiner_centroid(ccdt, face, region_polygon, 1);
		return test_obtuse_count;
	}
	Point centroid=CGAL::centroid(v1,v2,v3,v4);

	//remove the points
	//v1
	ccdt.remove_no_flip(face->vertex((indx+1) % 3));
	//v2
	ccdt.remove_no_flip((ccdt.insert_no_flip(v2)));
	//v3
	ccdt.remove_no_flip((ccdt.insert_no_flip(v3)));
	//v4
	ccdt.remove_no_flip((ccdt.insert_no_flip(v4)));

	//add the centroid insert_no_flip
	auto vh_cent=ccdt.insert_no_flip(centroid);
	
	//add back v1 v2 v3 v4 insert_no_flip
	//make (v1,v3) (v3,v2) (v2,v4) (v4,v1) constrains
	auto vh1=ccdt.insert_no_flip(v1);
	auto vh3=ccdt.insert_no_flip(v3);
	ccdt.insert_constraint_no_flip(vh1,vh3);
	auto vh2=ccdt.insert_no_flip(v2);
	auto vh4=ccdt.insert_no_flip(v4);
	ccdt.insert_constraint_no_flip(vh2,vh4);
	ccdt.insert_constraint_no_flip(vh4,vh1);
	ccdt.insert_constraint_no_flip(vh3, vh2);
	
	
	//remove the constrains remove_constrain_no_flip
	int edge_indx;
	Face_handle face13;
	if(!(find_face(ccdt,vh1,vh3,edge_indx, face13))){
		std::cout<<"Couldnt locate vh1 vh3 face\n";
		exit(-1);
	}
	ccdt.remove_constraint_no_flip(face13,edge_indx);
	Face_handle face24;
	if(!find_face(ccdt,vh2,vh4,edge_indx, face24)){
		std::cout<<"Couldnt locate vh2 vh4 face\n";
		exit(-1);
	}
	ccdt.remove_constraint_no_flip(face24,edge_indx);
	Face_handle face41;
	if(!find_face(ccdt,vh4,vh1,edge_indx, face41)){
		std::cout<<"Couldnt locate vh1 vh4 face\n";
		exit(-1);
	}
	ccdt.remove_constraint_no_flip(face41,edge_indx);
	Face_handle face32;
	if(!find_face(ccdt,vh3,vh2,edge_indx, face32)){
		std::cout<<"Couldnt locate vh2 vh3 face\n";
		exit(-1);
	}
	ccdt.remove_constraint_no_flip(face32,edge_indx);

	test_obtuse_count=count_obtuse_faces(ccdt, region_polygon);


return test_obtuse_count;
}

int annealing_test_add_steiner_projection(Custom_CDT& ccdt,Face_handle face, const Polygon_2 & region_boundary, int testORadd){
	
	int og_obtuse_count=count_obtuse_faces(ccdt,region_boundary);
	//create a copy of the ccdt and add the steiner there to check if the obtuse triangles are reduces
	Custom_CDT test_ccdt(ccdt);
	Face_handle test_face=find_face(test_ccdt, face);
	Point test_steiner=steiner_projection(test_ccdt, test_face);
	
	test_ccdt.insert_no_flip(test_steiner);


	int test_obutse_count=count_obtuse_faces(test_ccdt, region_boundary);

	if(testORadd==1){
			Point steiner=steiner_projection(ccdt, face);
			
			ccdt.insert_no_flip(steiner,face);
	}

	return test_obutse_count;
}


int annealing_test_add_steiner_longest_edge(Custom_CDT& ccdt,Face_handle face, const Polygon_2 & region_boundary,int testORadd){
	
	int og_obtuse_count=count_obtuse_faces(ccdt,region_boundary);
	//create a copy of the ccdt and add the steiner there to check if the obtuse triangles are reduces
	Custom_CDT test_ccdt(ccdt);
	Face_handle test_face=find_face(test_ccdt, face);
	Point test_steiner=steiner_longest_edge(test_ccdt, test_face);
	
	test_ccdt.insert_no_flip(test_steiner);

	int test_obutse_count=count_obtuse_faces(test_ccdt, region_boundary);


	if(testORadd==1){
		Point steiner=steiner_longest_edge(ccdt, face);
		ccdt.insert_no_flip(steiner,face);
	}

	return test_obutse_count;
}
//------------------------------
int annealing_simulate_circumcenter(const Custom_CDT& cdt,Face_handle face,const Polygon_2& region_boundary){

	Custom_CDT cdt_copy(cdt);
	Face_handle face_copy=find_face(cdt_copy,face);

	std::unordered_map<Face_handle, bool> in_domain_map;
	boost::associative_property_map<std::unordered_map<Face_handle, bool>> in_domain(in_domain_map);
	CGAL::mark_domain_in_triangulation(cdt_copy, in_domain);

	Point circ=steiner_circumcenter(cdt_copy,face_copy);
	if( region_boundary.bounded_side(circ)!=CGAL::ON_BOUNDED_SIDE && region_boundary.bounded_side(circ)!=CGAL::ON_BOUNDARY){
		int test_obtuse_count=annealing_test_add_steiner_centroid(cdt_copy, face_copy, region_boundary, 0);
		return test_obtuse_count;
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


	Vertex_handle obtuse_vertex=cdt_copy.insert_no_flip(obtuse_point,face_copy);
	Face_handle neighbor=face_copy->neighbor(face_copy->index(obtuse_vertex));
	

	Point n1=neighbor->vertex(0)->point();
	Point n2=neighbor->vertex(1)->point();
	Point n3=neighbor->vertex(2)->point();
	if(!point_inside_triangle(n1, n2, n3, circ)){
		int test_obtuse_count=annealing_test_add_steiner_centroid(cdt_copy, face_copy, region_boundary, 0);
		return test_obtuse_count;
	}
	
	//neighbor i has the same index as the opposite point
	Point nonShared=neighbor->vertex(neighbor->index(face_copy))->point();
	//check if the common edge is constrained
	if(face_copy->is_constrained(face_copy->index(neighbor))){
		int test_obtuse_count=annealing_test_add_steiner_centroid(cdt_copy, face_copy, region_boundary, 0);
		return test_obtuse_count;
	}

	
	Vertex_handle circ_vh=cdt_copy.insert_no_flip(circ,neighbor);

	//make the obtutse point circ a constain to remove the edge in between
	obtuse_vertex=cdt_copy.insert_no_flip(obtuse_point);
	cdt_copy.insert_constraint_no_flip(circ_vh, obtuse_vertex);
	CGAL::mark_domain_in_triangulation(cdt_copy, in_domain);

	//after remove the constrain without flipping
	obtuse_vertex=cdt_copy.insert_no_flip(obtuse_point);
	circ_vh=cdt_copy.insert_no_flip(circ,neighbor);
	Face_handle temp_face;
	int edge_index=-1;
	cdt_copy.is_edge(obtuse_vertex, circ_vh, temp_face, edge_index);
	cdt_copy.remove_constraint_no_flip(temp_face, edge_index);

	CGAL::mark_domain_in_triangulation(cdt_copy, in_domain);

	int new_count=count_obtuse_faces(cdt_copy,region_boundary);
	return new_count;
}

//add centroid of polygon latter
int annealing_test_add_steiner_circumcenter(Custom_CDT& ccdt,Face_handle face, const Polygon_2 & region_boundary,int testORadd){

	int obtuse;
	if(testORadd!=1){
		obtuse=annealing_simulate_circumcenter(ccdt, face, region_boundary);
		return obtuse;
	}

	Point circ=steiner_circumcenter(ccdt,face);
	if( region_boundary.bounded_side(circ)!=CGAL::ON_BOUNDED_SIDE && region_boundary.bounded_side(circ)!=CGAL::ON_BOUNDARY){
		int test_obtuse_count=annealing_test_add_steiner_centroid(ccdt, face, region_boundary, 1);
		return test_obtuse_count;
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


	Vertex_handle obtuse_vertex=ccdt.insert_no_flip(obtuse_point,face);
	Face_handle neighbor=face->neighbor(face->index(obtuse_vertex));
	

	Point n1=neighbor->vertex(0)->point();
	Point n2=neighbor->vertex(1)->point();
	Point n3=neighbor->vertex(2)->point();
	if(!point_inside_triangle(n1, n2, n3, circ)){
		int test_obtuse_count=annealing_test_add_steiner_centroid(ccdt, face, region_boundary, 1);
		return test_obtuse_count;
	}
	
	//neighbor i has the same index as the opposite point
	Point nonShared=neighbor->vertex(neighbor->index(face))->point();
	//check if the common edge is constrained
	if(face->is_constrained(face->index(neighbor))){
		int test_obtuse_count=annealing_test_add_steiner_centroid(ccdt, face, region_boundary, 1);
		return test_obtuse_count;
	}

	
	Vertex_handle circ_vh=ccdt.insert_no_flip(circ,neighbor);

	//make the obtutse point circ a constain to remove the edge in between
	obtuse_vertex=ccdt.insert_no_flip(obtuse_point);
	ccdt.insert_constraint_no_flip(circ_vh, obtuse_vertex);

	//after remove the constrain without flipping
	obtuse_vertex=ccdt.insert_no_flip(obtuse_point);
	circ_vh=ccdt.insert_no_flip(circ,neighbor);
	Face_handle temp_face;
	int edge_index=-1;
	ccdt.is_edge(obtuse_vertex, circ_vh, temp_face, edge_index);
	ccdt.remove_constraint_no_flip(temp_face, edge_index);


	int new_count=count_obtuse_faces(ccdt,region_boundary);
	return new_count;

}
//------------------------------




int simulated_annealing(Custom_CDT& ccdt, Polygon_2 region_polygon,boost::associative_property_map<std::unordered_map<Face_handle, bool>> &in_domain, int steiner_count, double a, double b, int L,double& c_rate,int methods){
	int current_steiner=steiner_count;
	int current_obtuse=count_obtuse_faces(ccdt,in_domain);
	int prev_obtuse=current_obtuse;
	if(current_obtuse==0){
		return steiner_count;
	}
	int best_steiner=current_steiner;

	int test_obtuse=current_obtuse;

	double cooling_rate=1.0/L;
	double Temperature=1.0;
	double Best_E=current_steiner * b + current_obtuse * a;
	double current_E=Best_E;
	double test_E=current_E;
	double delta_E;
	//create random number between 0 and 4
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> func_selector(0,methods); 
	int choise= func_selector(gen);
	bool merge_index=false;
	Point merge_point;

	int centroid_added=0,circ_added=0,projection_added=0,midpoint_added=0,merged=0;
	Custom_CDT current_cdt(ccdt);
	int if_loop=0;
	for(int i=0;i<L;i++){
		if(Temperature<=0){
			break;;
		}
		for(Face_handle face:current_cdt.finite_face_handles()){
			if(is_in_region_polygon(face,region_polygon) && is_obtuse_triangle(face)){

				choise=func_selector(gen);
				switch (choise) {
					case 0:
						test_obtuse=annealing_test_add_steiner_projection(current_cdt, face,region_polygon, 0);
						break;
					case 1:
						test_obtuse=annealing_test_add_steiner_circumcenter(current_cdt, face,region_polygon, 0);
						break;

					case 2:
						test_obtuse=annealing_test_add_steiner_longest_edge(current_cdt, face,region_polygon, 0);
						break;
					case 3:
						test_obtuse=annealing_test_add_steiner_merge(current_cdt,face,region_polygon,merge_index,merge_point,0);
						break;
					case 4:
					
						test_obtuse=annealing_test_add_steiner_centroid(current_cdt, face,region_polygon, 0);
						break;
					default:
						std::cerr<<"Error anealing 1 switch invalid random choise\n";
						exit(-1);

				
				}

				test_E=test_obtuse*a + (current_steiner+1)*b;
				delta_E=test_E-Best_E;
				double random=getRandomDouble();
				if( (delta_E<0) || (random<=exp(-1*delta_E/Temperature)) ){
					switch (choise) {
						case 0:
							current_obtuse=annealing_test_add_steiner_projection(current_cdt, face,region_polygon, 1);
							projection_added++;
							break;
						case 1:
							current_obtuse=annealing_test_add_steiner_circumcenter(current_cdt, face,region_polygon, 1);
							circ_added++;
							break;
						case 2:
							current_obtuse=annealing_test_add_steiner_longest_edge(current_cdt, face,region_polygon, 1);
							midpoint_added++;
							break;
						case 3:
							test_obtuse=annealing_test_add_steiner_merge(current_cdt,face,region_polygon,merge_index,merge_point,1);
							merged++;
							break;
						case 4:
							current_obtuse=annealing_test_add_steiner_centroid(current_cdt, face,region_polygon, 1);
							centroid_added++;
							break;
						
						default:
							std::cerr<<"Error anealing 1 switch invalid random choise\n";
							exit(-1);

					
					}
					current_steiner++;
					current_E=test_E;

					if(current_E<Best_E){
						ccdt.clear();
						ccdt=current_cdt;
						Best_E=current_E;
						current_obtuse=count_obtuse_faces(ccdt,region_polygon);
						if(best_steiner > 1) {  // Only calculate rate after first steiner point
							c_rate += calculate_convergence_rate(best_steiner, prev_obtuse, current_obtuse);
						}
						prev_obtuse = current_obtuse;  // Update previous obtuse count
						best_steiner = current_steiner;
					}
					break;
				}
			}
		}
		Temperature=Temperature*(1-cooling_rate);
		if(prev_obtuse==0){
			//reached 0 obtuse no need to keep looking
			break;
		}
	}

	current_cdt.clear();
	std::cout<<"\nAfter SA we got:\n\tSteiners:"<<best_steiner<<"\n\t Obtuse count:"<<count_obtuse_faces(ccdt,region_polygon)
		<<"\nUsed:\n\tCentroid:"<<centroid_added<<"\n\tProjections:"<<projection_added<<"\n\tCircumcenter:"<<circ_added
		<<"\n\tMerged faces:"<<merged<<"\n\tLongest edge: "<<midpoint_added<<std::endl;

	return best_steiner;
}
