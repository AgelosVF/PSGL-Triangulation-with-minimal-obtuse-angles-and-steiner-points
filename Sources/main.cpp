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
#include "../Header Files/BoundaryType.hpp"
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
	if( (argc != 5) && (argc!=6 )){
		std::cout<<"Usage: "<<argv[0]<<" -i /path/to/input.json -o /path/to/output.json (optional -preselected_params) "<<argc<<std::endl;
		return 1;
	}
	std::string input_file,output_file;
	bool preselect=false;
	bool flag1=false,flag2=false;
	for(int i=1;i<argc;i++){
		if(strcmp(argv[i], "-i")==0){
			input_file=argv[i+1];
			i++;
			flag1=true;
		}else if (strcmp(argv[i], "-o")==0){
			output_file=argv[i+1];
			i++;
			flag2=true;
		}else if (strcmp(argv[i], "-preselected_params")==0){
			preselect=true;
		}
	}
	if( !flag1 || !flag2 || (preselect==false && argc!=5)){
		std::cerr<<"Error:invalid parameters. Run ./"<<argv[0]<<" if you want to see expected input.\n";
		if(!flag1)
			std::cerr<<"Missing input flag\n";
		if(!flag2)
			std::cerr<<"Missing output flag\n";
		return 0;
	}
	if(input_file.empty()){
		std::cerr<< "Error: Missing input file path."<<std::endl;
		return 0;
	}
	if(output_file.empty()){
		std::cerr<< "Error: Missing output file path."<<std::endl;
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


	/*
	std::cout<<"Region boundary\n";
	for(unsigned int i=0;i<region_boundary.size();i++){
		std::cout<<region_boundary[i]<<std::endl;

	}
	std::cout<<"Additional Constrains\n";
	for(unsigned int i=0;i<additional_constrains.size();i++){
		std::cout<<additional_constrains[i].first<<", "<<additional_constrains[i].second<<std::endl;
	}
	*/
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

	Polygon_2 region_polygon;
	for (int i : region_boundary) {
		region_polygon.push_back(Point(points_x[i], points_y[i]));
	}
	CGAL::draw(region_polygon);

	std::vector<int> closed_p;
	Polygon_2 Pcycle;
	int type=boundary_type(region_polygon,region_boundary, additional_constrains, closed_p);
	if(closed_p.size()>2){
		for(unsigned int i=0;i<closed_p.size();i++){
			std::cout<<closed_p[i]<<"->";
			int x=closed_p[i];
			Pcycle.push_back(Point(points_x[x],points_y[x]));
		}
		CGAL::draw(Pcycle);
	}
	
	// Mark the domain inside the region boundary
	std::unordered_map<Face_handle, bool> in_domain_map;
	boost::associative_property_map<std::unordered_map<Face_handle, bool>> in_domain(in_domain_map);
	CGAL::mark_domain_in_triangulation(cdel_tri, in_domain);
	CGAL::draw(cdel_tri,in_domain);

	unsigned int obtuse_count=count_obtuse_faces(cdel_tri,in_domain);
	std::cout<<obtuse_count<<std::endl;
	//create the polygon of the region boundry
	int steiner_count=0;

	if(preselect==true){
		std::string method=extract_method(pt);
		boost::property_tree::ptree parameters=extract_parameters(pt);

		if(!(extract_delaunay(pt))){
			std::cout<<"Starting by using applying the previous Project to the triangulation\n";
			steiner_count+=previous_triangulation(cdel_tri, region_polygon);
		}
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
	}
/*
	//MOVED

	CGAL::mark_domain_in_triangulation(cdel_tri, in_domain);	
	//void ant_colony(Custom_CDT& cdt,Polygon_2 boundary,double alpha ,double beta, double xi, double ps, double lambda, int kappa, int L, int steiner_points
	ant_colony(cdel_tri,region_polygon,4.0 ,0.01, 1.0, 3.0, 0.5, 11, 4, 0);
	CGAL::mark_domain_in_triangulation(cdel_tri, in_domain);	
	CGAL::draw(cdel_tri,in_domain);
	obtuse_count=count_obtuse_faces(cdel_tri,in_domain);
	std::cout<<"Final obtuse count:"<<obtuse_count<<std::endl;


	CGAL::mark_domain_in_triangulation(cdel_tri, in_domain);
	CGAL::draw(cdel_tri,in_domain);
	
	/*
	obtuse_count=count_obtuse_faces(cdel_tri, in_domain);
	ant_colony(cdel_tri,region_polygon,4.0 ,2.0, 1.0, 3.0, 0.5, 4, 30, obtuse_count);
	std::cout<<"Final obtuse count:"<<obtuse_count<<std::endl;
	
	
	// Generate the output JSON file
	std::vector<Point> initial_points;
	for (size_t i = 0; i < points_x.size(); ++i) {
		initial_points.emplace_back(points_x[i], points_y[i]);
	}
	
	generate_output_json(cdel_tri, inst_iud, initial_points, region_polygon, output_file,  method, parameters, obtuse_count);
	*/
	return 0;
}



