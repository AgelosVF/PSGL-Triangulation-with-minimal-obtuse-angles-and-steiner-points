#include "../Header Files/CGAL_CUSTOM_CONSTRAINED_DELAUNAY_TRIANGULATION_2.h"
#include "../Header Files/triangulation_utils.hpp"
#include <CGAL/Polygon_2.h>


bool is_convex(Polygon_2 boundary){
	if(boundary.is_convex())
		return true;
	else
		return false;
}

bool convex_no_constrains(Polygon_2 boundary,std::vector<std::pair<int,int>> additional_constrains){
	if(is_convex(boundary) && additional_constrains.empty())
		return true;
	else
		return false;
}

bool non_convex_parrallel(Polygon_2 boundary){

	for(auto  e  : boundary.edges()){
		
		Point p= e.start();
		Point p2=e.end();
		//std::cout<<"Start: "<<p<<" End: "<<p2<<std::endl;
		if(p.x()!=p2.x() && p.y() != p2.y()){
			std::cout<<"Start: "<<p<<" End: "<<p2<<std::endl;
			return false;
		}
	}
	return true;
}
