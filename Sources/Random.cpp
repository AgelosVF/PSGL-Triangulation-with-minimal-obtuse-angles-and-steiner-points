#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include "../Header Files/triangulation_utils.hpp"
#include "../Header Files/random.hpp"
#include <CGAL/Line_2.h>
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

