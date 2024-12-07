#include <CGAL/Qt/Basic_viewer_qt.h>
#include <CGAL/Uncertain.h>
#include <boost/property_map/property_map.hpp>
#include <cmath>
#include <cstddef>
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
		ccdt.insert_no_flip(steiner,face);
	}

	return test_obutse_count;
}
void simulated_annealing(Custom_CDT& ccdt, Polygon_2 region_polygon,boost::associative_property_map<std::unordered_map<Face_handle, bool>> &in_domain, int steiner_count, double a, double b, int L){
	int current_steiner=steiner_count;
	int current_obtuse=count_obtuse_faces(ccdt,in_domain);
	int best_steiner=current_steiner;

	// not needed will be current +1 int test_steiner=current_steiner;
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
	std::uniform_int_distribution<> func_selector(0, 2); 
	int choise= func_selector(gen);
	std::cout<<func_selector(gen)<<" "<<func_selector(gen)<<" "<<func_selector(gen)<<" "<<func_selector(gen);

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
				std::cout<<"Choose:"<<choise<<std::endl;
				switch (choise) {
					case 0:
						test_obtuse=annealing_test_add_steiner_projection(current_cdt, face,region_polygon, 0);
						break;
					case 1:
						test_obtuse=annealing_test_add_steiner_longest_edge(current_cdt, face,region_polygon, 0);
						break;

					case 2:
						test_obtuse=annealing_test_add_steiner_centroid(current_cdt, face,region_polygon, 0);
						break;
					case 3:
						break;
					case 4:
						break;
					
					default:
						std::cerr<<"Error anealing 1 switch invalid random choise\n";
						exit(-1);

				
				}

				test_E=test_obtuse*a + (current_steiner+1)*b;
				delta_E=test_E-Best_E;

				if( (delta_E<0) || (getRandomDouble()<exp(-1*delta_E/Temperature)) ){
					switch (choise) {
						case 0:
							current_obtuse=annealing_test_add_steiner_projection(current_cdt, face,region_polygon, 1);
							projection_added++;
							break;
						case 1:
							current_obtuse=annealing_test_add_steiner_longest_edge(current_cdt, face,region_polygon, 1);
							midpoint_added++;
							break;
						case 2:
							current_obtuse=annealing_test_add_steiner_centroid(current_cdt, face,region_polygon, 1);
							centroid_added++;
							break;
						case 3:
							break;
						case 4:
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
						best_steiner=current_steiner;
						std::cout<<"Eddited best!!!!\n";
//						CGAL::draw(ccdt);
					}
					break;
				}
			}
		}
		Temperature=Temperature*(1-cooling_rate);
	}
	//std::cout<<"\n-----------------------------------------\n-----------------------------------------Final current\n";
	//CGAL::draw(current_cdt);
	current_cdt.clear();
	//CGAL::draw(ccdt);
	std::cout<<"\nAfter SA we got:\n\tSteiners:"<<best_steiner<<"\n\t Obtuse count:"<<count_obtuse_faces(ccdt,region_polygon)
		<<"\nUsed:\n\tCentroid:"<<centroid_added<<"\n\tProjections:"<<projection_added<<"\n\tCircumcenter:"<<circ_added<<std::endl;


}
