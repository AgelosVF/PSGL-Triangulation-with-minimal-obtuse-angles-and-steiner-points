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
#include <iostream>
#include <cstddef>
#include <iostream>
#include <ostream>
#include <vector>
#include <unordered_map>
#include <CGAL/polygon_function_objects.h>
#include <CGAL/Polygon_2.h>
//#include <CGAL/draw_polygon_2.h>

Point old_add_steiner_point_longest_edge(Custom_CDT& cdel_tri, Face_handle face) {
    // Get the vertices of the face
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
//-----------------------------------------------------------------------------------------
int old_reduce_obtuse_midpoint_steiner(Custom_CDT &cdel_tri, const Polygon_2 &region_boundary, boost::associative_property_map<std::unordered_map<Face_handle, bool>> &in_domain){
	bool steiner_added=true;
	int steiner_count=0;

	while (steiner_added){
		steiner_added=false;
		int og_obtuse_count=count_obtuse_faces(cdel_tri,in_domain);

		for(Face_handle f:cdel_tri.finite_face_handles()){

			//if the face is obtuse and in the region boundary
			if(get(in_domain,f) && is_obtuse_triangle(f) ){
				//get the point that is in the midpoint of the longest edge
				Custom_CDT local_cdt(cdel_tri);

				Point steiner=old_add_steiner_point_longest_edge(cdel_tri,f);
				// Retrieve the points of face `f` in `cdel_tri`
				Point p0 = f->vertex(0)->point();
				Point p1 = f->vertex(1)->point();
				Point p2 = f->vertex(2)->point();
				Point centroid= CGAL::centroid(p0,p1,p2);
				Face_handle local_f = local_cdt.locate(centroid);
				if( local_f != Face_handle() ){

					local_cdt.insert_no_flip(steiner,local_f);
					std::unordered_map<Face_handle, bool> in_domain_map_copy;
					boost::associative_property_map<std::unordered_map<Face_handle, bool>> in_domain_copy(in_domain_map_copy);
					CGAL::mark_domain_in_triangulation(local_cdt, in_domain_copy);

					int new_count=count_obtuse_faces(local_cdt, in_domain_copy);

					if(new_count<og_obtuse_count){
						
						cdel_tri.insert_no_flip(steiner,f);
						steiner_count++;
						steiner_added=true;
						CGAL::mark_domain_in_triangulation(cdel_tri,in_domain);
						break;
					}
				}

			}
		}
	}

	return steiner_count;
};
//------------------------------------------------------------------------------------------
Point old_add_steiner_point_circumcenter(Custom_CDT& cdel_tri, Face_handle face, const Polygon_2& region_boundary){

	// Get the vertices of the face
	Point p0 = face->vertex(0)->point();
	Point p1 = face->vertex(1)->point();
	Point p2 = face->vertex(2)->point();

	// Calculate the circumcenter of the triangle formed by p0, p1, and p2
	Point circumcenter = CGAL::circumcenter(p0, p1, p2);

	//check if the circumcenter is inside the region boundary polygon
	CGAL::Bounded_side check= CGAL::bounded_side_2(region_boundary.vertices_begin(),region_boundary.vertices_end(),circumcenter);
	if( (check==CGAL::ON_BOUNDARY) || (check==CGAL::ON_BOUNDED_SIDE) ){
		return circumcenter;
	}
	else{
		//if the circumcenter is outside the boundary add point in centroid
		Point centroid=CGAL::centroid(p0,p1,p2);
		return centroid;
	}
};
int old_reduce_obtuse_circumcenter_steiner(Custom_CDT &cdel_tri, const Polygon_2 &region_boundary, boost::associative_property_map<std::unordered_map<Face_handle, bool>> &in_domain){
	bool steiner_added=true;
	int steiner_count=0;

	while (steiner_added){
		steiner_added=false;
		int og_obtuse_count=count_obtuse_faces(cdel_tri,in_domain);

		for(Face_handle f:cdel_tri.finite_face_handles()){

			if(get(in_domain,f) && is_obtuse_triangle(f) ){
				Point steiner=old_add_steiner_point_circumcenter(cdel_tri, f,region_boundary);
				Custom_CDT cdel_copy(cdel_tri);
				Face_handle insert_face=cdel_copy.locate(steiner);//if point in face a face hander is returned

				if( insert_face != Face_handle() ){
					cdel_copy.insert_no_flip(steiner,insert_face);
					std::unordered_map<Face_handle, bool> in_domain_map_copy;
					boost::associative_property_map<std::unordered_map<Face_handle, bool>> in_domain_copy(in_domain_map_copy);
					CGAL::mark_domain_in_triangulation(cdel_copy, in_domain_copy);

					int new_obtuse_count=count_obtuse_faces(cdel_copy,in_domain_copy);

					if (new_obtuse_count<og_obtuse_count){
						cdel_tri.insert_no_flip(steiner,f);
						steiner_count++;
						steiner_added=true;
						CGAL::mark_domain_in_triangulation(cdel_tri, in_domain); // Reset and update the in_domain map
						break;
					}

				}
			}
		}


	}
	return steiner_count;
};
//--------------------------------------------//
Point old_add_steiner_point_projection(CDT& cdt, Face_handle F){
	
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
int old_reduce_obtuse_projection_steiner(Custom_CDT &cdel_tri, const Polygon_2 &region_boundary, boost::associative_property_map<std::unordered_map<Face_handle, bool>> &in_domain){
	bool steiner_added=true;
	int steiner_count=0;

	while (steiner_added){
		steiner_added=false;
		int og_obtuse_count=count_obtuse_faces(cdel_tri,in_domain);

		for(Face_handle f:cdel_tri.finite_face_handles()){

			//if the face is obtuse and in the region boundary
			if(get(in_domain,f) && is_obtuse_triangle(f) ){
				//get the point that is in the midpoint of the longest edge
				Custom_CDT local_cdt(cdel_tri);

				Point steiner=old_add_steiner_point_projection(cdel_tri,f);
				// Retrieve the points of face `f` in `cdel_tri`
				Point p0 = f->vertex(0)->point();
				Point p1 = f->vertex(1)->point();
				Point p2 = f->vertex(2)->point();
				Point centroid= CGAL::centroid(p0,p1,p2);
				Face_handle local_f = local_cdt.locate(centroid);
				if( local_f != Face_handle() ){

					local_cdt.insert_no_flip(steiner,local_f);
					std::unordered_map<Face_handle, bool> in_domain_map_copy;
					boost::associative_property_map<std::unordered_map<Face_handle, bool>> in_domain_copy(in_domain_map_copy);
					CGAL::mark_domain_in_triangulation(local_cdt, in_domain_copy);

					int new_count=count_obtuse_faces(local_cdt, in_domain_copy);

					
					if(new_count<og_obtuse_count){
						cdel_tri.insert(steiner);
						steiner_count++;
						steiner_added=true;
						CGAL::mark_domain_in_triangulation(cdel_tri, in_domain);
						break;
					}
				}
			}
		}
	}
	return steiner_count;
};
int previous_triangulation(Custom_CDT& cdel_tri,const Polygon_2& region_polygon){
	
	// Mark the domain inside the region boundary
	std::unordered_map<Face_handle, bool> in_domain_map;
	boost::associative_property_map<std::unordered_map<Face_handle, bool>> in_domain(in_domain_map);
	CGAL::mark_domain_in_triangulation(cdel_tri, in_domain);//marks faces connected with non constrained edges as inside of the domain based on the nesting level.
	
	unsigned int num_obtuse=count_obtuse_faces(cdel_tri,in_domain);
	unsigned int new_obtuse=count_obtuse_faces(cdel_tri, in_domain);
	int steiner_count=0;

	bool reduced=true;
	while(reduced){
		reduced=false;
		num_obtuse=new_obtuse;
		//try edge flips
		reduce_obtuse_by_flips(cdel_tri, in_domain);
		CGAL::mark_domain_in_triangulation(cdel_tri, in_domain);
		new_obtuse=count_obtuse_faces(cdel_tri,in_domain);
		if(new_obtuse<num_obtuse)
			reduced=true;
		//try projection steiner
		steiner_count+=old_reduce_obtuse_projection_steiner(cdel_tri, region_polygon, in_domain);
		CGAL::mark_domain_in_triangulation(cdel_tri, in_domain);
		new_obtuse=count_obtuse_faces(cdel_tri,in_domain);
		if(new_obtuse<num_obtuse)
			reduced=true;
		//try circumcenter steiner
		steiner_count+=old_reduce_obtuse_circumcenter_steiner(cdel_tri,region_polygon,in_domain);
		CGAL::mark_domain_in_triangulation(cdel_tri, in_domain);
		new_obtuse=count_obtuse_faces(cdel_tri,in_domain);
		if(new_obtuse<num_obtuse)
			reduced=true;
		//try midpoint of largest edge steiners
		steiner_count+=old_reduce_obtuse_midpoint_steiner(cdel_tri,region_polygon, in_domain);
		new_obtuse=count_obtuse_faces(cdel_tri,in_domain);
		if(new_obtuse<num_obtuse)
			reduced=true;
		CGAL::mark_domain_in_triangulation(cdel_tri, in_domain);
	}
	return steiner_count;
}


