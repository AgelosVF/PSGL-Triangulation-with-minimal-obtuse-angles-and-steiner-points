#include <CGAL/Distance_2/Point_2_Segment_2.h>
#include <CGAL/Distance_2/Segment_2_Ray_2.h>
#include <CGAL/Distance_2/Segment_2_Segment_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include "../Header Files/SteinerPoints.hpp"
#include "../Header Files/triangulation_utils.hpp"
#include "../Header Files/random.hpp"
#include <CGAL/Line_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Qt/Basic_viewer_qt.h>
#include <random>
#include "../Header Files/CGAL_CUSTOM_CONSTRAINED_DELAUNAY_TRIANGULATION_2.h"
#include <CGAL/mark_domain_in_triangulation.h>
#include <unordered_map>
#include "CGAL/draw_triangulation_2.h"

bool flip_valid(const Point& v1, const Point& v2, const Point& v3, const Point& v4) {
	Polygon_2 pol;
	pol.push_back(v1);
	pol.push_back(v3);
	pol.push_back(v2);
	pol.push_back(v4);

	
	//CGAL::draw(pol);
	// Check if the quadrilateral is convex
	return pol.is_convex();
}

bool try_flip(Custom_CDT& cdt, Vertex_handle other, Face_handle face, Polygon_2 region){
    
    int edge_index=face->index(other);
    Face_handle neighbor=face->neighbor(edge_index);

    
    
    if (!is_in_region_polygon(neighbor, region) || cdt.is_constrained({face, edge_index})){
        return false;
    }
    Point v1 = face->vertex((edge_index + 1) % 3)->point();
    Point v2 = face->vertex((edge_index + 2) % 3)->point();
    Point v3 = neighbor->vertex(neighbor->index(face))->point();
    Point v4 = face->vertex(edge_index)->point();

  
    if (!flip_valid(v1,v2,v3,v4)) {
        return false;
    }

    cdt.flip(face, edge_index);                        
    return true;
}

Point random_point_in_face_gaussian(Face_handle face) {
    // Calculate centroid of the face
    Point p1 = face->vertex(0)->point();
    Point p2 = face->vertex(1)->point();
    Point p3 = face->vertex(2)->point();
    // Calculate the centroid of the face
    Point centroid=CGAL::centroid(p1,p2,p3);
    double cent_x=(CGAL::to_double(centroid.x()));
    double cent_y=(CGAL::to_double(centroid.y()));

    // Create random number generator
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::normal_distribution<> dist_x(cent_x, 0.1);
    std::normal_distribution<> dist_y(cent_y, 0.1);

    // Generate random point using Gaussian distribution
    double x = dist_x(gen);
    double y = dist_y(gen);
    Point random(x,y);
    // Check if the point is inside the face
    if (!point_inside_triangle(p1, p2, p3, centroid)) {
        //if the point is not inside the face, try again
        return random_point_in_face_gaussian(face);
    } 

    return Point(x, y);
}

bool random_point_on_edge(Custom_CDT& cdt,Face_handle face) {
    
    Face_handle c_face=find_face(cdt,face);
    if(face_still_exists(cdt,c_face)){
	Point random=random_point_in_face_gaussian(face);
	Point start,end;
	Point p1 = face->vertex(0)->point();
	Point p2 = face->vertex(1)->point();
	Point p3 = face->vertex(2)->point();
	// Calculate the squared lengths of the triangle's sides

	Kernel::FT side1 = CGAL::squared_distance(p1, p2);
	Kernel::FT side2 = CGAL::squared_distance(p2, p3);
	Kernel::FT side3 = CGAL::squared_distance(p3, p1);
	// Determine the longest edge and calculate its midpoint
	if( (side1>side2) && (side1>side3) ){
	    start=p1;
	    end=p2;

	}
	else if ( (side2>side1) && (side2>side3)){
	    start=p2;
	    end=p3;
	}
	else if ( (side3>side1) && (side3>side2)){
	    start=p3;
	    end=p1;

	}
	Kernel::Line_2 longest_edge(start,end);
	
	Point midpoint=CGAL::midpoint(start,end);

	Point random_on_edge=longest_edge.projection(random);
	cdt.insert_no_flip(random_on_edge);
	return true;
    }
    return false;
}

bool random_point_with_merge(Custom_CDT& cdt, Face_handle face){

    Face_handle c_face=find_face(cdt,face);
    
    if(face_still_exists(cdt,c_face)){
    
	//get the vertices of the face
	Point p0=c_face->vertex(0)->point();
	Point p1=c_face->vertex(1)->point();
	Point p2=c_face->vertex(2)->point();
	

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

	Vertex_handle obtuse_vertex=cdt.insert_no_flip(obtuse_point,c_face);
	Face_handle neighbor=c_face->neighbor(c_face->index(obtuse_vertex));

	Point nonShared=neighbor->vertex(neighbor->index(c_face))->point();
	int nonShared_idx = neighbor->index(c_face);
	Point shared1 = neighbor->vertex(cdt.cw(nonShared_idx))->point();
	Point shared2 = neighbor->vertex(cdt.ccw(nonShared_idx))->point();
	std::cout<<nonShared<<" "<<shared2<<" "<<shared1<<std::endl;

	if(c_face->is_constrained(c_face->index(neighbor)))
	    return false;
	//after finding the neighbor we try to find a neighbor of his so we can merge 3 faces
	length1=CGAL::squared_distance(nonShared, shared1);
	length2=CGAL::squared_distance(nonShared, shared2);
	Point neighbor_2,le_start_2,le_end_2;

	if(length1>=length2){
	    le_start_2=nonShared;
	    le_end_2=shared1;
	    neighbor_2=shared2;
	    
	    Vertex_handle neighbors_neighbor_vertex=cdt.insert_no_flip(neighbor_2);
	    Face_handle neighbors_neighbor=neighbor->neighbor(neighbor->index(neighbors_neighbor_vertex));
	    int nonShared_idx2=neighbors_neighbor->index(neighbor);
	    Vertex_handle NonShared2=neighbor->vertex(nonShared_idx2);

	    Polygon_2 merge;
	    merge.push_back(obtuse_point);



	}
	


	return true;
    }
    return false;
}


int random_and_flips(Custom_CDT& cdt, Face_handle face, Polygon_2 region){
    int og_count=count_obtuse_faces(cdt,region);
    Point random;
    int flag=0;
    int steiners_added=0;
    for(int i=0;i<400;i++){
	Custom_CDT cdt_copy(cdt);

	random=random_point_in_face_gaussian(face);
	Face_handle face_copy=find_face(cdt_copy,face);
	Vertex_handle v_copy=cdt_copy.insert_no_flip(random,face_copy);
	std::unordered_map<Face_handle, bool> in_domain_map_copy;
	boost::associative_property_map<std::unordered_map<Face_handle, bool>> in_domain_copy(in_domain_map_copy);
	CGAL::mark_domain_in_triangulation(cdt_copy, in_domain_copy);
	reduce_obtuse_by_flips(cdt_copy,in_domain_copy);
	int new_count=count_obtuse_faces(cdt_copy,in_domain_copy);
	if(new_count<=og_count){
	    flag=1;
	    break;
	}
    }
    if(flag==1){
	Vertex_handle v=cdt.insert_no_flip(random,face);
	std::unordered_map<Face_handle, bool> in_domain_map;
	boost::associative_property_map<std::unordered_map<Face_handle, bool>> in_domain(in_domain_map);
	CGAL::mark_domain_in_triangulation(cdt, in_domain);
	reduce_obtuse_by_flips(cdt,in_domain);
	steiners_added++;
    }
    return steiners_added;
}
int reduce_random_local(Custom_CDT& cdt,int starting_obt_count,Polygon_2 region_polygon,bool& randomed,double& c_rate){

    int c_obtuse_count=starting_obt_count;
    int steiner_count=0;
    int current_steiners=0;
    Custom_CDT c_cdt(cdt);
    std::unordered_map<Face_handle, bool> in_domain_map;
    boost::associative_property_map<std::unordered_map<Face_handle, bool>> in_domain(in_domain_map);
    CGAL::mark_domain_in_triangulation(c_cdt, in_domain);
    

    for(int i=0;i<10;i++){
	//Add random steiners to 1/4 of the obtuse faces to move around the obtuse angles
	std::vector<Face_handle> obtuse_faces;
	for(Face_handle f : c_cdt.finite_face_handles()) {
		if(get(in_domain,f) && is_obtuse_triangle(f)) {
			obtuse_faces.push_back(f);
		}
	}
	size_t quarter_size = (obtuse_faces.size() + 3) / 4;  // This rounds up division
	for(size_t i = 0; i < quarter_size; i++) {
		if(random_point_on_edge(c_cdt, obtuse_faces[i]))
			current_steiners++;
	}

	CGAL::mark_domain_in_triangulation(c_cdt, in_domain);
	reduce_obtuse_by_flips(c_cdt, in_domain);
	current_steiners+=local_search(c_cdt, region_polygon, 1000, in_domain, c_rate);

	CGAL::mark_domain_in_triangulation(c_cdt, in_domain);
	reduce_obtuse_by_flips(c_cdt, in_domain);
	c_obtuse_count=count_obtuse_faces(c_cdt,in_domain);
	if(c_obtuse_count<starting_obt_count){
		randomed=true;
		cdt.clear();
		cdt=c_cdt;
		steiner_count+=current_steiners;
		current_steiners=0;
		starting_obt_count=c_obtuse_count;
	}
	if(c_obtuse_count==0){
	    break;
	}
    }
    return steiner_count;
}

