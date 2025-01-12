//custom libraries
#include "../Header Files/triangulation_utils.hpp"
#include "../Header Files/CGAL_CUSTOM_CONSTRAINED_DELAUNAY_TRIANGULATION_2.h"
//libraries
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Qt/Basic_viewer_qt.h>
#include <CGAL/Triangulation_face_base_2.h>
#include <CGAL/draw_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_2.h>
#include <CGAL/mark_domain_in_triangulation.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>
#include <CGAL/Triangulation_data_structure_2.h>
//other libs
#include <cmath>
#include <iostream>
#include <cstddef>
#include <iostream>
#include <ostream>
#include <vector>
#include <unordered_map>
#include <CGAL/polygon_function_objects.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/draw_polygon_2.h>

int find_obtuse_vertex_index(const Face_handle& Face){

	Point p0=Face->vertex(0)->point();
	Point p1=Face->vertex(1)->point();
	Point p2=Face->vertex(2)->point();

	Kernel::FT length1 = CGAL::squared_distance(p0, p1);
	Kernel::FT length2 = CGAL::squared_distance(p1, p2);
	Kernel::FT length3 = CGAL::squared_distance(p2, p0);
	//FInd the longest edge and the oposite vertex
	if( (length1 > length2) && (length1>length3)){
		return 2;

	}
	else if ((length2>length1)&& (length2>length3) ) {
		return 0;
	}
	else if ((length3>length1) &&(length3>length2) ) {
		return 1;
	}
	else
		std::cerr<<"Got isosceles face\n";
		exit(-1);

}
double calculate_convergence_rate(int steiner_points, int prev_obtuse, int current_obtuse){
//we calculate p(n-1)	
	if(steiner_points==0 ||prev_obtuse==0){
		std::cerr<<"Called convergence rate with "<<steiner_points<<" steiner points, "<<prev_obtuse<<" previous obtuse count and "<<current_obtuse<<" current obtuse"<<"\n";
		exit(-1);
	}
	if(current_obtuse==0){
		return 0;
	}
	double rate = log((double)current_obtuse/prev_obtuse) / log((double)steiner_points/(steiner_points-1));
	std::cout<<"Steiner points: "<<steiner_points<<"\tPrev_count: "<<prev_obtuse<<"\tCurrent count: "<<current_obtuse<< "\nRate: "<<rate<<std::endl;
	return  rate;
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

// Function to check if a face is inside or on the boundary of the polygon
bool is_in_region_polygon(Face_handle f, const Polygon_2& region_polygon) {
	// Compute the centroid of the face
	Point centroid = CGAL::centroid(
	f->vertex(0)->point(), 
	f->vertex(1)->point(), 
	f->vertex(2)->point()
	);

	// Check if the centroid is inside or on the boundary of the polygon
	return (region_polygon.bounded_side(centroid) == CGAL::ON_BOUNDED_SIDE || region_polygon.bounded_side(centroid) == CGAL::ON_BOUNDARY);
}

// Function to count obtuse triangles in the triangulation within the polygon
int count_obtuse_faces(const Custom_CDT& cdel_tri, const Polygon_2& region_polygon) {
	int obtuse_count = 0;

	// Iterate over all finite faces in the triangulation
	for (Face_handle f : cdel_tri.finite_face_handles()) {
		// Check if the face is inside the polygon and is obtuse
		if (is_in_region_polygon(f, region_polygon) && is_obtuse_triangle(f)) {
		    obtuse_count++;
		}
	}

return obtuse_count;
}

// Function to compute the distance between two points
double compute_distance(const Point& p1, const Point p2) {
    // Compute the squared distance
    Kernel::FT sq_dist = CGAL::squared_distance(p1, p2);

    // Return the square root of the squared distance for the actual distance
    return std::sqrt(CGAL::to_double(sq_dist));
}


//get the face handle of the same face in a diffrent ccdt
Face_handle find_face(Custom_CDT& n_ccdt, Face_handle o_face){

	
	Point p0 = o_face->vertex(0)->point();
	Point p1 = o_face->vertex(1)->point();
	Point p2 = o_face->vertex(2)->point();
	Point centroid= CGAL::centroid(p0,p1,p2);
	Face_handle n_face = n_ccdt.locate(centroid);
	if( n_face == nullptr){
		std::cerr<<"Error: couldnt locate original ccdt face in the copy cccdt"<<std::endl;
	}
	return n_face;
}
//get the face that contains the three points
Face_handle find_face(Custom_CDT& ccdt,Point p1, Point p2, Point p3){


	Point centroid= CGAL::centroid(p1,p2,p3);
	Face_handle face = ccdt.locate(centroid);
	if( face == nullptr){
		std::cerr<<"Error: couldnt locate original ccdt face in the copy cccdt"<<std::endl;
	}
	return face;


}
//finds the face that has both vh1 and vh2
bool find_face(Custom_CDT& ccdt,Vertex_handle vh1, Vertex_handle vh2,int& index, Face_handle& face){

	//go throught the faces that have vh1
	for(auto fc=vh1->incident_faces();fc != nullptr;++fc){
		//if one has alse vh2
		if(fc->has_vertex(vh2)){
			for(int i=0;i<3;i++){
				//find the index of the edge that has both vertices
				if( ((fc->vertex(i)==vh1) && (fc->vertex( (i+1)%3)==vh2)) || ((fc->vertex(i)==vh2) && (fc->vertex( (i+1)%3)==vh1)) ){
					face=fc;
					index=(i+2)%3;
					return true;
				}

			}
		}
	}
	return false;}


//helper function for debugging
void print_point(const Point& point) {
	std::cout << "Point coordinates: (" << point.x() << ", " << point.y() << ")" << std::endl;
}

//returns yes if face is constrained
bool is_face_constrained(const Custom_CDT& ccdt, Face_handle face) {
	for (int i = 0; i < 3; ++i) {
		Vertex_handle vh = face->vertex(i);
		if (ccdt.are_there_incident_constraints(vh)) {
			return true;
		}
	}
return false;
}

Vertex_handle find_vertex_by_point(const Custom_CDT& ccdt, const Point& point) {
	for (auto vit = ccdt.finite_vertices_begin(); vit != ccdt.finite_vertices_end(); ++vit) {
		if (vit->point() == point) {
			return vit;
		}
	}
return Vertex_handle(); // Return an invalid handle if not found
}

bool face_has_constrained_edge(const Custom_CDT& cdt, Face_handle face) {
    for (int i = 0; i < 3; ++i) {
        if (cdt.is_constrained(CDT::Edge(face, i))) {
            return true; // At least one edge is constrained
        }
    }
    return false; // No constrained edges
}

bool is_obtuse_simulated(const Point& p1, const Point& p2, const Point& p3) {
	// Calculate squared side lengths
	Kernel::FT side1 = CGAL::squared_distance(p1, p2);
	Kernel::FT side2 = CGAL::squared_distance(p2, p3);
	Kernel::FT side3 = CGAL::squared_distance(p3, p1);

	// Check obtuseness
	return (side1 + side2 < side3) || (side2 + side3 < side1) || (side3 + side1 < side2);
}

// Function to check if two segments intersect
bool segments_intersect(const Point& p1, const Point& p2, const Point& p3, const Point& p4) {
	CGAL::Segment_2<Kernel> seg1(p1, p2);
	CGAL::Segment_2<Kernel> seg2(p3, p4);
	return CGAL::do_intersect(seg1, seg2);
}


// Function to check if an edge flip is valid by testing if the quadrilateral is convex
bool is_flip_valid(const Point& v1, const Point& v2, const Point& v3, const Point& v4) {
	Polygon_2 pol;
	pol.push_back(v1);
	pol.push_back(v3);
	pol.push_back(v2);
	pol.push_back(v4);

	
	// Check if the quadrilateral is convex
	return pol.is_convex();
}

//takes the face edge we are trying to flip and in_domain map
bool try_edge_flip(Custom_CDT& cdel_tri, Face_handle f, int edge_index, boost::associative_property_map<std::unordered_map<Face_handle, bool>>& in_domain) {
	Face_handle neighbor = f->neighbor(edge_index);

	// Only flip if the neighbor face is also in-domain and the edge is not constrained
	if (!get(in_domain, neighbor) || cdel_tri.is_constrained({f, edge_index})) {
		return false;
	}

	// Get the vertices of the face and its neighbor before flipping
	Point v1 = f->vertex((edge_index + 1) % 3)->point();
	Point v2 = f->vertex((edge_index + 2) % 3)->point();
	Point v3 = neighbor->vertex(neighbor->index(f))->point();
	Point v4 = f->vertex(edge_index)->point();


	int obtuse_pre_flip=0;
	int obtuse_post_flip=0;
	//count obtuse pre flip
	if(is_obtuse_triangle(f))
		obtuse_pre_flip++;
	if(is_obtuse_triangle(neighbor))
		obtuse_pre_flip++;
	//count obtuse post flip by simulation the triangles
	if(is_obtuse_simulated(v1, v3, v4))
		obtuse_post_flip++;
	if(is_obtuse_simulated(v2, v3, v4))
		obtuse_post_flip++;

	// Apply the flip to the original triangulation only if it reduces obtuse triangles
	if (obtuse_pre_flip > obtuse_post_flip) {
		if (is_flip_valid(v1,v2,v3,v4) ){
			cdel_tri.flip(f, edge_index);                        
			
			// Reset and re-fill the in_domain map after each successful flip
			CGAL::mark_domain_in_triangulation(cdel_tri, in_domain);
			return true;
		}
	}

	return false;
}


void reduce_obtuse_by_flips(Custom_CDT& cdel_tri, boost::associative_property_map<std::unordered_map<Face_handle, bool>>& in_domain) {
	bool flip_made = true;

	while (flip_made) {
		flip_made = false;

		// Iterate over all faces in the triangulation
		for (Face_handle f : cdel_tri.finite_face_handles()) {
			if (get(in_domain, f) && is_obtuse_triangle(f)) {
				// Try flipping each edge of the obtuse triangle
				for (int i = 0; i < 3; ++i) {
					if (try_edge_flip(cdel_tri, f, i, in_domain)) {
						flip_made = true;
						break;  // Stop further flips on this face, as it has been modified
					}
				}
				if (flip_made) break;  // Recalibrate face handles after a successful flip
			}
		}
	}
}



void fill_initial_cdt(
    Custom_CDT& cdel_tri,
    const std::vector<int>& points_x,
    const std::vector<int>& points_y,
    const std::vector<int>& region_boundary,
    const std::vector<std::pair<int, int>>& additional_constrains)
{
	std::vector<CDT::Vertex_handle> vertices;
	for(int i=0; i< points_y.size();i++){
		Point p(points_x[i],points_y[i]);
		vertices.push_back(cdel_tri.insert(p));
	}

	//add the region boundary as constrain
	for(size_t i=0;i<region_boundary.size();i++){
		int start=region_boundary[i];
		int end=region_boundary[(i+1)%region_boundary.size()]; //mod it so that the last point gets connected to the first
		if( ( start >= vertices.size() ) || ( end >= vertices.size() ) ){
			std::cerr<<"Error: invalid region boundary: ("<<start<<", "<<end<<std::endl;
			exit(-1);
		}
		else if ( ( start<0 ) || ( end<0 ) ) {
			std::cerr<<"Error: invalid region boundary: ("<<start<<", "<<end<<std::endl;
			exit(-1);
		}
		else
			cdel_tri.insert_constraint(vertices[start],vertices[end]);
	}

	
	//add the additional constrains
	for(size_t i=0; i<additional_constrains.size();i++){
		int start=additional_constrains[i].first;
		int end=additional_constrains[i].second;
		if( ( start>=0 ) && ( start < vertices.size() ) && ( end >=0 ) && (end < vertices.size() ) )
			cdel_tri.insert_constraint(vertices[start],vertices[end]);
		else{
			std::cerr<<"Error: invalid additional constrain: ("<<start<<", "<<end<<")"<<std::endl;
			exit(-1);
		}
	}
}


int count_obtuse_faces(const Custom_CDT& cdel_tri, const boost::associative_property_map<std::unordered_map<Face_handle, bool>>& in_domain) {
	int obtuse_count = 0;

	// Iterate over all faces in the triangulation
	for (Face_handle f : cdel_tri.finite_face_handles()) {
		// Check if the face is inside the domain
		if (get(in_domain, f)) {
			if(is_obtuse_triangle(f))
				obtuse_count++;
		}
	}

	return obtuse_count;
}


void print_polygon_vertices(const Polygon_2& polygon) {
	if (polygon.is_empty()) {
		std::cout << "Polygon is empty." << std::endl;
		return;
	}

	std::cout << "Polygon vertices:" << std::endl;
	for (const auto& vertex : polygon.vertices()) {
		std::cout << "(" << vertex.x() << ", " << vertex.y() << ")" << std::endl;
	}
	draw(polygon);
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

bool is_obtuse_triangle(const typename CDT::Face_handle& f) {

	// Get the three vertices of the triangle
	Point v1 = f->vertex(0)->point();
	Point v2 = f->vertex(1)->point();
	Point v3 = f->vertex(2)->point();

	// Calculate the squared lengths of the triangle's sides

	Kernel::FT side1 = CGAL::squared_distance(v1, v2);
	Kernel::FT side2 = CGAL::squared_distance(v2, v3);
	Kernel::FT side3 = CGAL::squared_distance(v3, v1);

	// Check if it's an obtuse triangle using the inequality x^2 + y^2 < z^2
	return (side1 + side2 < side3) || (side2 + side3 < side1) || (side3 + side1 < side2);
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
