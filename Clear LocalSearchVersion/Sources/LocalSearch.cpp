#include <CGAL/Qt/Basic_viewer_qt.h>
#include <boost/property_map/property_map.hpp>
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

int choose_smallest(int a, int b, int c, int d, int e) {
    std::vector<int> values = {a, b, c, d, e};
    int min_value = *std::min_element(values.begin(), values.end());
    
    // Collect indices of all occurrences of the minimum value
    std::string result = "";
    for (size_t i = 0; i < values.size(); ++i) {
        if (values[i] == min_value) {
            return i+1; // Indices are 1-based
        }
    }
    
    // Convert the concatenated indices back to an integer
    return -1;
}

void local_search(Custom_CDT& ccdt, Polygon_2 region_polygon, int L, boost::associative_property_map<std::unordered_map<Face_handle, bool>> &in_domain){
	std::cout<<"Cool\n";
	bool steiner_added=true;
	//for some reason i need to re-mark_domain else it is blank
	CGAL::mark_domain_in_triangulation(ccdt, in_domain);
	int og_obtuse_count=count_obtuse_faces(ccdt,in_domain);
	int added_centroid=0,added_projections=0,added_longest_edge=0,added_merge=0,added_circ=0;
	int e=0;
	int best_i=-1;
	int best_steiner=-1;
	int merge_flag=-1;
	int obt_count=og_obtuse_count;
	int steiner_le_count=99999,steiner_ctr_count=99999,steiner_proj_count=99999,steiner_merge_count=99999,steiner_circ_count=999999;
	std::cout<<"Original obtuse triangles: "<<og_obtuse_count<<"\t"<<steiner_le_count<<std::endl;
	while((steiner_added==true) &&(L>0) ){
		steiner_added=false;
		L--;
		//std::cout<<"Loop:"<<L<<std::endl;
		for(Face_handle f:ccdt.finite_face_handles()){
			if(get(in_domain,f) && is_obtuse_triangle(f) )	{
				steiner_ctr_count=test_add_steiner_centroid(ccdt, f, region_polygon,in_domain, 0);
				steiner_le_count=test_add_steiner_longest_edge(ccdt, f, region_polygon, in_domain, 0);
				steiner_proj_count=test_add_steiner_projection(ccdt, f, region_polygon, in_domain, 0);
				steiner_merge_count=test_add_steiner_merge(ccdt,f,region_polygon,in_domain,merge_flag,0);
				steiner_circ_count=test_add_steiner_circumcenter(ccdt, f, region_polygon, in_domain, 0);

				//choose the best centroid , projection , longest edge , circumcenter , merge i case of ties the first in the previous order will be chosen.
				int best_i=choose_smallest(steiner_ctr_count,steiner_proj_count,steiner_le_count,steiner_circ_count,steiner_merge_count);
				//std::cout<<"Best choice="<<best_i<<std::endl;
				if( best_i==1){
					if(steiner_ctr_count<obt_count){
						e=test_add_steiner_centroid(ccdt,f,region_polygon,in_domain,1);
						steiner_added=true;
						if(e==-1){
							std::cerr<<"Tried to add point on centroid that wouldnt reduce the count in 1 L="<<L<<std::endl;
							steiner_added=false;
							break;
						}
						added_centroid++;
						break;
					}
				}
				if( best_i==2){
					if(steiner_proj_count<obt_count){
						e=test_add_steiner_projection(ccdt,f,region_polygon,in_domain,1);
						steiner_added=true;
						if(e==-1){
							std::cerr<<"Tried to add point on projection that wouldnt reduce the count in 3 L="<<L<<std::endl;
							steiner_added=false;
							break;
						}
						added_projections++;
						break;
					}
				}
				if( best_i==3){
					if(steiner_le_count<obt_count){
						e=test_add_steiner_longest_edge(ccdt,f,region_polygon,in_domain,1);
						steiner_added=true;
						if(e==-1){
							std::cerr<<"Tried to add point on longest edge that wouldnt reduce the count in 2 L="<<L<<std::endl;
							steiner_added=false;
							break;
						}
						added_longest_edge++;
						break;
					}
				}
				if( best_i==4){
					if(steiner_circ_count<obt_count){
						e=test_add_steiner_circumcenter(ccdt, f, region_polygon, in_domain, 1);
						steiner_added=true;
						if(e==214748364){
						std::cerr<<"Tried to add point on circumcenter that wasnt valid in 4 L="<<L<<std::endl;
						steiner_added=false;
						break;
						}
					added_circ++;
					break;
						
					}
				}
				if(best_i==5){
					if(steiner_merge_count<obt_count){
						e=test_add_steiner_merge(ccdt, f, region_polygon, in_domain,merge_flag, 1);
						steiner_added=true;
						if(e==214748364){
							std::cerr<<"Tried to add point on invalid merge of triangles in 5 L="<<L<<std::endl;
							steiner_added=false;
							break;
						}
					added_merge++;
					break;
					}
				}
			}
		}
		CGAL::mark_domain_in_triangulation(ccdt, in_domain); // Reset and update the in_domain map
		obt_count=count_obtuse_faces(ccdt, in_domain);
	}

	std::cout<<"After local search i added:\n\t"<<"Steiners centroid:"<<added_centroid<<"\n\tSteiners projection:"<<added_projections<<"\n\tSteiners on longest edge:"<<added_longest_edge<<"\n\tSteiners added on circumcenter:"<<added_circ<<"\n\tSteiners added by merging neighbors:"<<added_merge<<std::endl;
	CGAL::mark_domain_in_triangulation(ccdt, in_domain); // Reset and update the in_domain map
	obt_count=count_obtuse_faces(ccdt, in_domain);
	std::cout<<"Total number of obtuse triangles remaining:"<<obt_count<<std::endl;
}
