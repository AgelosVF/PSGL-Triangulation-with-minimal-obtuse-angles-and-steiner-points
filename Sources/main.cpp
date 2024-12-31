#include <CGAL/Qt/Basic_viewer_qt.h>
#include <boost/property_map/property_map.hpp>
#include <cstddef>
#include <cstring>
#include <iostream>
#include <ostream>
#include <fstream>
#include <vector>
#include <map>
#include <unordered_map>
//custom
#include "../Header Files/boost_utils.hpp"
#include "../Header Files/triangulation_utils.hpp"
#include "../Header Files/CGAL_CUSTOM_CONSTRAINED_DELAUNAY_TRIANGULATION_2.h"
//boost libraries
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
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



int main(int argc,char *argv[]){
	if(argc != 5){
		std::cout<<"Usage: "<<argv[0]<<" -i /path/to/input.json -o /path/to/output.json"<<std::endl;
		return 1;
	}
	std::string input_file,output_file;

	for(int i=1;i<argc;i+=2){
		if(strcmp(argv[i], "-i")==0){
			input_file=argv[i+1];
		}else if (strcmp(argv[i], "-o")==0){
			output_file=argv[i+1];
		}else {
			std::cerr<<"Error: Uknown flag"<<argv[i]<<std::endl;
			return 0;
		}

	}

	if(input_file.empty() || output_file.empty()){
		std::cerr<< "Error: Missing input or output file path."<<std::endl;
		return 0;
	}

	boost::property_tree::ptree pt;
	try{
		boost::property_tree::read_json(input_file,pt);
	} catch (const boost::property_tree::json_parser_error &e){
		std::cerr<<"Error: reading JSON file: "<<e.what()<<std::endl;
		exit(-1);
	}
	//extract the info from the file using the utils functions
	std::vector<int> points_x=extract_x_coordinates(pt);
	std::vector<int> points_y=extract_y_coordinates(pt);
	std::vector<int> region_boundary=extract_region_boundary(pt); 
	std::vector<std::pair<int,int>> additional_constrains= extract_additional_constrains(pt);
	boost::optional<std::string> inst_iud = pt.get_optional<std::string>("instance_uid");

	std::string method=extract_method(pt);
	boost::property_tree::ptree parameters=extract_parameters(pt);
	if(!inst_iud){
		std::cerr<<"Error: missing inst_iud from JSON file"<<std::endl;
		exit(-1);
	}
	//make sure that we have the same number of x and y coordinates
	if(points_x.size() != points_y.size()){
		std::cerr<<"Error: mismatch in the number of x and y coordinates."<<std::endl;
		exit(-1);
	}

	Custom_CDT cdel_tri;
	//use custom function to put the points and constrains in the trianglulation
	fill_initial_cdt(cdel_tri, points_x, points_y, region_boundary, additional_constrains);

	// Mark the domain inside the region boundary
	std::unordered_map<Face_handle, bool> in_domain_map;
	boost::associative_property_map<std::unordered_map<Face_handle, bool>> in_domain(in_domain_map);
	CGAL::mark_domain_in_triangulation(cdel_tri, in_domain);
	//create the polygon of the region boundry
	Polygon_2 region_polygon;
	for (int i : region_boundary) {
		region_polygon.push_back(Point(points_x[i], points_y[i]));
	}

	int steiner_count=0;
/*
	if(!(extract_delaunay(pt))){
		std::cout<<"Starting by using applying the previous Project to the triangulation\n";
		steiner_count+=previous_triangulation(cdel_tri, region_polygon);
	}
	*/
	unsigned int obtuse_count=count_obtuse_faces(cdel_tri,in_domain);

	std::cout<<obtuse_count<<std::endl;
	CGAL::mark_domain_in_triangulation(cdel_tri, in_domain);	
	CGAL::draw(cdel_tri,in_domain);
	ant_colony(cdel_tri,region_polygon,4.0 ,0.01, 1.0, 3.0, 0.5, 3, 1000, 0);
	CGAL::mark_domain_in_triangulation(cdel_tri, in_domain);	
	CGAL::draw(cdel_tri,in_domain);
	obtuse_count=count_obtuse_faces(cdel_tri,in_domain);
	std::cout<<"Final obtuse count:"<<obtuse_count<<std::endl;

	/*
	if(method=="local"){
		int L;
		extract_local_parameters(parameters, L);
		steiner_count+=local_search(cdel_tri, region_polygon, 5000, in_domain);
	}
	else if(method=="SA"){
		double alpha,beta;
		int L;
		extract_sa_parameters(parameters,alpha, beta, L);
		simulated_annealing(cdel_tri,region_polygon,in_domain, steiner_count, alpha, beta, L);
	}
	else if (method=="ant") {
		double alpha,beta,xi,psi,lambda;
		int kappa,L;
		extract_ant_parameters(parameters,alpha, beta, xi, psi, lambda, kappa, L);
		std::cout<<"Ant colony still under constraction the code is in AntCollony.cpp but isnt linked to program.\n";
	
	}
	CGAL::mark_domain_in_triangulation(cdel_tri, in_domain);
	CGAL::draw(cdel_tri,in_domain);
	
	/*
	obtuse_count=count_obtuse_faces(cdel_tri, in_domain);
	ant_colony(cdel_tri,region_polygon,4.0 ,2.0, 1.0, 3.0, 0.5, 4, 30, obtuse_count);
	std::cout<<"Final obtuse count:"<<obtuse_count<<std::endl;
	
	*/
	// Generate the output JSON file
	std::vector<Point> initial_points;
	for (size_t i = 0; i < points_x.size(); ++i) {
		initial_points.emplace_back(points_x[i], points_y[i]);
	}
	generate_output_json(cdel_tri, inst_iud, initial_points, region_polygon, output_file,  method, parameters, obtuse_count);
	return 0;
}



