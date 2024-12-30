#include "../Header Files/CGAL_CUSTOM_CONSTRAINED_DELAUNAY_TRIANGULATION_2.h"
#include "../Header Files/triangulation_utils.hpp"
#include <CGAL/Kernel/global_functions_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/draw_polygon_2.h>
#include <CGAL/Qt/Basic_viewer_qt.h>
#include <CGAL/enum.h>
#include <CGAL/squared_distance_2.h> // For calculating squared distances
#include <boost/property_map/property_map.hpp>
#include <iostream>
#include <iterator>
#include <unistd.h>
#include <unordered_map>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/mark_domain_in_triangulation.h>
#include <CGAL/draw_constrained_triangulation_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>


// steiner.cpp
// sometimes it is better to insert in projections/longest edge
//------------------------------------------------------------------------------------

void print_polygon_vertices(const Polygon_2& polygon) {
	if (polygon.is_empty()) {
		std::cout << "Polygon is empty." << std::endl;
		return;
	}

	std::cout << "Polygon vertices:" << std::endl;
	for (const auto& vertex : polygon.vertices()) {
		std::cout << "(" << vertex.x() << ", " << vertex.y() << ")" << std::endl;
	}
}
bool point_inside_triangle(Point p1, Point p2, Point p3,Point query){

	Polygon_2 triangle;
	triangle.push_back(p1);

	triangle.push_back(p2);
	triangle.push_back(p3);
	if(triangle.bounded_side(query)==CGAL::ON_BOUNDED_SIDE || triangle.bounded_side(query)==CGAL::ON_BOUNDARY){
		return true;
	}
	return false;

}
// Function to check if an edge flip is valid by testing if the quadrilateral is convex
bool is_polygon_convex(const Point& v1, const Point& v2, const Point& v3, const Point& v4) {
	Polygon_2 pol;
	pol.push_back(v1);
	pol.push_back(v3);
	pol.push_back(v2);
	pol.push_back(v4);

	
	
	// Check if the quadrilateral is convex
	return pol.is_convex();
}

//---------------------------------------------------------------------------------------

// Adds a Steiner point at the midpoint of the longest edge of a given face
Point steiner_longest_edge(Custom_CDT& cdel_tri, Face_handle face) {
    // Get the vertices of the face
	// Get the three vertices of the triangle
	Point p1 = face->vertex(0)->point();
	Point p2 = face->vertex(1)->point();
	Point p3 = face->vertex(2)->point();

	point_inside_triangle( p1,  p2,  p3,p1);
	// Calculate the squared lengths of the triangle's sides

	Kernel::FT side1 = CGAL::squared_distance(p1, p2);
	Kernel::FT side2 = CGAL::squared_distance(p2, p3);
	Kernel::FT side3 = CGAL::squared_distance(p3, p1);
	// Determine the longest edge and calculate its midpoint
	Point midpoint;
	if( (side1>side2) && (side1>side3) )
		midpoint = CGAL::midpoint(p1, p2);
	else if ( (side2>side1) && (side2>side3))
		midpoint= CGAL::midpoint(p2,p3);
	else if ( (side3>side1) && (side3>side2))
		midpoint=CGAL::midpoint(p3,p1);
	else
		std::cerr<<"Could determine longest edge between: "<<side1<<" "<<side2<<" "<<side3<<std::endl;

	// Insert the Steiner point without flipping
	return midpoint;
}

//adds a steiner point in a copy of ccdt and counts the new obtuse triangle count. If testOradd==1 if the count is reduced it adds the point to the ccdt and returns 1 else it returns 0. If testORadd==0 it returns the number of obtuse triangles that inserting the steiner would produce.
int test_add_steiner_longest_edge(Custom_CDT& ccdt,Face_handle face, const Polygon_2 & region_boundary, boost::associative_property_map<std::unordered_map<Face_handle, bool>> &in_domain,int testORadd){
	
	int og_obtuse_count=count_obtuse_faces(ccdt,in_domain);
	//create a copy of the ccdt and add the steiner there to check if the obtuse triangles are reduces
	Custom_CDT test_ccdt(ccdt);
	Face_handle test_face=find_face(test_ccdt, face);
	Point test_steiner=steiner_longest_edge(test_ccdt, test_face);
	
	test_ccdt.insert_no_flip(test_steiner);

	std::unordered_map<Face_handle, bool> in_domain_map_copy;
	boost::associative_property_map<std::unordered_map<Face_handle, bool>> in_domain_copy(in_domain_map_copy);
	CGAL::mark_domain_in_triangulation(test_ccdt, in_domain_copy);
	int test_obutse_count=count_obtuse_faces(test_ccdt, in_domain_copy);


	if(testORadd==1){
		if(test_obutse_count<og_obtuse_count){
			Point steiner=steiner_longest_edge(ccdt, face);
			ccdt.insert_no_flip(steiner,face);
			return 1;
		}
		else{
			return -1;
		}
	}

	return test_obutse_count;
}

//------------------------------------------------------------------------------------------------//

Point steiner_centroid(Custom_CDT& cdel_tri, Face_handle face){

	// Get the vertices of the face
	Point p0 = face->vertex(0)->point();
	Point p1 = face->vertex(1)->point();
	Point p2 = face->vertex(2)->point();

	Point centroid=CGAL::centroid(p0,p1,p2);
	return centroid;
}


//adds a steiner point in a copy of ccdt and counts the new obtuse triangle count. If testOradd==1 if the count is reduced it adds the point to the ccdt and returns 1 else it returns 0. If testORadd==0 it returns the number of obtuse triangles that inserting the steiner would produce.
int test_add_steiner_centroid(Custom_CDT& ccdt,Face_handle face, const Polygon_2 & region_boundary, boost::associative_property_map<std::unordered_map<Face_handle, bool>> &in_domain,int testORadd){
	
	int og_obtuse_count=count_obtuse_faces(ccdt,in_domain);
	//create a copy of the ccdt and add the steiner there to check if the obtuse triangles are reduces
	Custom_CDT test_ccdt(ccdt);
	Face_handle test_face=find_face(test_ccdt, face);
	Point test_steiner=steiner_centroid(test_ccdt, test_face);
	
	test_ccdt.insert_no_flip(test_steiner);

	std::unordered_map<Face_handle, bool> in_domain_map_copy;
	boost::associative_property_map<std::unordered_map<Face_handle, bool>> in_domain_copy(in_domain_map_copy);
	CGAL::mark_domain_in_triangulation(test_ccdt, in_domain_copy);
	//CGAL::draw(test_ccdt,in_domain_copy);
	int test_obutse_count=count_obtuse_faces(test_ccdt, in_domain_copy);

	if(testORadd==1){
		if(test_obutse_count<og_obtuse_count){
			Point steiner=steiner_centroid(ccdt, face);
			ccdt.insert_no_flip(steiner,face);
			return 1;
		}
		else{
			return -1;
		}
	}

	return test_obutse_count;
}


//--------------------------------------------//
Point steiner_projection(Custom_CDT& cdt, Face_handle F){
	
	//get the vertices of the face
	Point v1=F->vertex(1)->point();
	Point v2=F->vertex(2)->point();
	Point v0=F->vertex(0)->point();

	//calculate the length of the edges
	//
	Kernel::FT length1 = CGAL::squared_distance(v0, v1);
	Kernel::FT length2 = CGAL::squared_distance(v1, v2);
	Kernel::FT length3 = CGAL::squared_distance(v2, v0);

	//FInd the longest edge and the oposite vertex
	Point le_start,le_end,obtuse_vertex;
	if( (length1 > length2) && (length1>length3)){
		le_start=v0;
		le_end=v1;
		obtuse_vertex=v2;
	}
	else if ((length2>length1)&& (length2>length3) ) {
		le_start=v1;
		le_end=v2;
		obtuse_vertex=v0;
	}
	else if ((length3>length1) &&(length3>length2) ) {
		le_start=v2;
		le_end=v0;
		obtuse_vertex=v1;
	
	}

	//calculate the projection of the obtuse vertex on the longest edge
	Kernel::Line_2 longest_edge(le_start,le_end);
	Point projection= longest_edge.projection(obtuse_vertex);
	return projection;
}

int test_add_steiner_projection(Custom_CDT& ccdt,Face_handle face, const Polygon_2 & region_boundary, boost::associative_property_map<std::unordered_map<Face_handle, bool>> &in_domain,int testORadd){
	
	int og_obtuse_count=count_obtuse_faces(ccdt,in_domain);
	//create a copy of the ccdt and add the steiner there to check if the obtuse triangles are reduces
	Custom_CDT test_ccdt(ccdt);
	Face_handle test_face=find_face(test_ccdt, face);
	Point test_steiner=steiner_projection(test_ccdt, test_face);
	
	test_ccdt.insert_no_flip(test_steiner);

	std::unordered_map<Face_handle, bool> in_domain_map_copy;
	boost::associative_property_map<std::unordered_map<Face_handle, bool>> in_domain_copy(in_domain_map_copy);
	CGAL::mark_domain_in_triangulation(test_ccdt, in_domain_copy);
	int test_obutse_count=count_obtuse_faces(test_ccdt, in_domain_copy);

	if(testORadd==1){
		if(test_obutse_count<og_obtuse_count){
			Point steiner=steiner_projection(ccdt, face);
			
			ccdt.insert_no_flip(steiner,face);
			return 1;
		}
		else{
			return -1;
		}
	}

	return test_obutse_count;
}

//---------------------------------------------------------------------------------------------------------//
Point steiner_circumcenter(Custom_CDT& cdel_tri, Face_handle face){

	// Get the vertices of the face
	Point p0 = face->vertex(0)->point();
	Point p1 = face->vertex(1)->point();
	Point p2 = face->vertex(2)->point();

	// Calculate the circumcenter of the triangle formed by p0, p1, and p2
	Point circumcenter = CGAL::circumcenter(p0, p1, p2);

	return circumcenter;
};

int simulate_merge_steiner(Custom_CDT& ccdt,Face_handle ob_face,const Polygon_2& region_boundary,int& best_merge){
	
	//check if the vertices of the face are constrained
	if(is_face_constrained(ccdt,ob_face) ){
		return 21474836;
	}
	//create a copy and locate the face in it
	Custom_CDT ccdt_copy(ccdt);

	Face_handle face_copy;
	Face_handle neighbor;
	std::unordered_map<Face_handle, bool> in_domain_map_copy;
	boost::associative_property_map<std::unordered_map<Face_handle, bool>> in_domain_copy(in_domain_map_copy);
	CGAL::mark_domain_in_triangulation(ccdt_copy, in_domain_copy);
	int new_ob_count=21474836;
	best_merge=-1;
	for(int i=0;i<3;i++){
		face_copy=find_face(ccdt_copy,ob_face);
		neighbor=face_copy->neighbor(i);

		//skip  constrained ,non obtuse, outside of region_boundary neighbors
		if( ccdt_copy.is_infinite(neighbor) ){
			continue;
		}
		
		if( !(is_obtuse_triangle(neighbor)) ){
			continue;
		}
		if( !(get(in_domain_copy,neighbor)) ){
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
		CGAL::mark_domain_in_triangulation(ccdt_copy, in_domain_copy);
		int temp_count=count_obtuse_faces(ccdt_copy, in_domain_copy);
		if(temp_count<new_ob_count){
			new_ob_count=temp_count;
			best_merge=i;
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
		//remove constraint vh1 vh2
		

		CGAL::mark_domain_in_triangulation(ccdt_copy, in_domain_copy);
	}
	return new_ob_count;
}

int test_add_steiner_merge(Custom_CDT& ccdt,Face_handle face,Polygon_2& region_polygon,boost::associative_property_map<std::unordered_map<Face_handle, bool>> &in_domain, int& indx,int testORadd){
	int obtuse;
	if(testORadd!=1){
		obtuse=simulate_merge_steiner(ccdt,face,region_polygon,indx);
		return obtuse;
	}
	
	auto face_2=face->neighbor(indx);

	//skip  constrained ,non obtuse, outside of region_boundary neighbors
	if( ccdt.is_infinite(face_2) ){
		std::cerr<<"Tried to merge infinite face\n";
		exit(-1);
	}
	
	if( !(is_obtuse_triangle(face_2)) ){
		std::cerr<<"Tried to merge with non obtuse face\n";
		exit(-1);
	}
	if( !(get(in_domain,face_2)) ){
		std::cerr<<"Tried to merge with outside of the domain\n";
		exit(-1);
	}
	if(is_face_constrained(ccdt,face_2)){
		std::cerr<<"Tried to merge with constrained face\n";
		exit(-1);
	}

	Point v1 = face->vertex((indx + 1) % 3)->point();
	Point v2 = face->vertex((indx + 2) % 3)->point();
	Point v3 = face_2->vertex(face_2->index(face))->point();
	Point v4 = face->vertex(indx)->point();

	if( !( is_polygon_convex(v1, v2, v3,v4) ) ){
		std::cerr<<"Tried to merge non convex polygon\n";
		exit(-1);
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

	CGAL::mark_domain_in_triangulation(ccdt, in_domain);

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
	
	CGAL::mark_domain_in_triangulation(ccdt, in_domain);
	
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

	CGAL::mark_domain_in_triangulation(ccdt, in_domain);
	obtuse=count_obtuse_faces(ccdt, in_domain);


return obtuse;
}



int find_longest_edge_index(Face_handle face){
	// Get the three vertices of the triangle
	Point p1 = face->vertex(0)->point();
	Point p2 = face->vertex(1)->point();
	Point p3 = face->vertex(2)->point();
	// Calculate the squared lengths of the triangle's sides

	Kernel::FT side1 = CGAL::squared_distance(p1, p2);
	Kernel::FT side2 = CGAL::squared_distance(p2, p3);
	Kernel::FT side3 = CGAL::squared_distance(p3, p1);
	// Determine the longest edge and calculate its midpoint
	Point midpoint;
	if( (side1>side2) && (side1>side3) )
		return 0;
	else if ( (side2>side1) && (side2>side3))
		return 1;
	else if ( (side3>side1) && (side3>side2))
		return 2;
	else
		std::cerr<<"Could determine longest edge between: "<<side1<<" "<<side2<<" "<<side3<<std::endl;

	return -1;

}
//add centroid of polygon latter
int test_add_steiner_circumcenter(Custom_CDT& ccdt,Face_handle face, const Polygon_2 & region_boundary, boost::associative_property_map<std::unordered_map<Face_handle, bool>> &in_domain,int testORadd){

	Point circ=steiner_circumcenter(ccdt,face);
	CGAL::mark_domain_in_triangulation(ccdt, in_domain);
	//check if the point is inside the boundary
	if( region_boundary.bounded_side(circ)!=CGAL::ON_BOUNDED_SIDE && region_boundary.bounded_side(circ)!=CGAL::ON_BOUNDARY){
		return 214748364;
	}
	if(is_face_constrained(ccdt,face)){
		
		return 214748364;
	}

	//oposite neighbor i shares edge i with face and the oposite vertex is i
	int i=find_longest_edge_index(face);


	Point obp0=face->vertex(i)->point();
	Point cp1=face->vertex((i+1)%3)->point();
	Point cp2=face->vertex((i+2)%3)->point();
	
	Face_handle neighbor=face->neighbor(i);
	Point n1=neighbor->vertex(0)->point();
	Point n2=neighbor->vertex(1)->point();
	Point n3=neighbor->vertex(2)->point();

	Point nonShared=neighbor->vertex(neighbor->index(face))->point();
	//check if the circ is inside the neighbor
	if(!point_inside_triangle(n1, n2, n3, circ)){
		//std::cout<<"\tCircumcenter outside neighbor\n";
		return 214748364;
	}

	//remove then add steiner then make constrain then return point then unconstrain
	
	ccdt.remove_no_flip( (ccdt.insert_no_flip(obp0)) );
	ccdt.remove_no_flip( (ccdt.insert_no_flip(cp1)) );
	ccdt.remove_no_flip( (ccdt.insert_no_flip(cp2)) );

	std::cout<<"\tRemoved triangle\n";
	CGAL::mark_domain_in_triangulation(ccdt, in_domain);
	ccdt.insert_points_and_constraint_no_flip(obp0, circ);
	ccdt.insert_no_flip(cp1);
	ccdt.insert_no_flip(cp2);

	//get the vertex handles of circ and the old obtuse
	Vertex_handle vh0=ccdt.insert_no_flip(obp0);
	Vertex_handle vhC=ccdt.insert_no_flip(circ);
	
	Face_handle face_p0c;
	int edge_indx;
	if(!(find_face(ccdt,vh0,vhC,edge_indx, face_p0c))){
		std::cout<<"Couldnt locate vh1 vh3 face\n";
		exit(-1);
	}
	ccdt.remove_constraint_no_flip(face_p0c,edge_indx);

	if(testORadd==1){
		CGAL::mark_domain_in_triangulation(ccdt, in_domain);
		return count_obtuse_faces(ccdt,in_domain);
	}
	else{
		
		CGAL::mark_domain_in_triangulation(ccdt, in_domain);
		//we need to reset the cdt to its previous state remove the steiner force the old edges as constrains
		ccdt.remove_no_flip(vhC);
		Vertex_handle vh1=ccdt.insert_no_flip(cp1);
		Vertex_handle vhn=ccdt.insert_no_flip(nonShared);
		ccdt.insert_constraint_no_flip(vh1 , vhn );
		Vertex_handle vh2=ccdt.insert_no_flip(cp2);
		ccdt.insert_constraint_no_flip( vh2, vhn);
		Vertex_handle vh0=ccdt.insert_no_flip(obp0);

		ccdt.insert_constraint_no_flip(vh1, vh0);
		ccdt.insert_constraint_no_flip(vh2, vh0);


		Face_handle face01;
		if(!find_face(ccdt,vh0,vh1,edge_indx, face01)){
			std::cout<<"Couldnt locate vh1 vh0 face\n";
			exit(-1);
		}
		ccdt.remove_constraint_no_flip(face01,edge_indx);
		
		Face_handle face02;
		if(!find_face(ccdt,vh0,vh2,edge_indx, face02)){
			std::cout<<"Couldnt locate vh0 vh2 face\n";
			exit(-1);
		}
		ccdt.remove_constraint_no_flip(face02,edge_indx);

		Face_handle face2n;
		if(!find_face(ccdt,vhn,vh2,edge_indx, face2n)){
			std::cout<<"Couldnt locate vhn vh2 face\n";
			exit(-1);
		}
		ccdt.remove_constraint_no_flip(face2n,edge_indx);


		Face_handle face1n;
		if(!find_face(ccdt,vhn,vh1,edge_indx, face1n)){
			std::cout<<"Couldnt locate vhn vh1 face\n";
			exit(-1);
		}
		ccdt.remove_constraint_no_flip(face1n,edge_indx);



		std::cout<<"\n Circumcenter:Reset triangulation\n---------------------------------\n";
		CGAL::mark_domain_in_triangulation(ccdt, in_domain);
		return count_obtuse_faces(ccdt,in_domain);
	}
}
