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
#include "../Header Files/random.hpp"
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

struct InputParameters {
	std::string input_file;
	std::string output_file;
	bool preselect;
	boost::property_tree::ptree pt;
};

InputParameters read_inputs(int argc, char* argv[]);

int main(int argc, char* argv[]) {
	auto params = read_inputs(argc, argv);
	
	// Extract info from the file using the utils functions
	std::vector<int> points_x = extract_x_coordinates(params.pt);
	std::vector<int> points_y = extract_y_coordinates(params.pt);
	std::vector<int> region_boundary = extract_region_boundary(params.pt); 
	std::vector<std::pair<int,int>> additional_constrains = extract_additional_constrains(params.pt);
	boost::optional<std::string> inst_iud = params.pt.get_optional<std::string>("instance_uid");

	if (!inst_iud) {
		std::cerr << "Error: missing inst_iud from JSON file" << std::endl;
		exit(-1);
	}

	if (points_x.size() != points_y.size()) {
		std::cerr << "Error: mismatch in the number of x and y coordinates." << std::endl;
		exit(-1);
	}

	Custom_CDT ccdt;
	fill_initial_cdt(ccdt, points_x, points_y, region_boundary, additional_constrains);

	Polygon_2 region_polygon;
	for (int i : region_boundary) {
		region_polygon.push_back(Point(points_x[i], points_y[i]));
	}

	
	std::unordered_map<Face_handle, bool> in_domain_map;
	boost::associative_property_map<std::unordered_map<Face_handle, bool>> in_domain(in_domain_map);
	CGAL::mark_domain_in_triangulation(ccdt, in_domain);
	unsigned int obtuse_count = count_obtuse_faces(ccdt, in_domain);
	std::cout << "Original obtuse count: " << obtuse_count << std::endl;

	int steiner_count = 0;
	std::string method;
	boost::property_tree::ptree parameters;
	double c_rate=0;
	double final_rate=0;

	if (!params.preselect) {
		method = extract_method(params.pt);
		parameters = extract_parameters(params.pt);
		
		if (!extract_delaunay(params.pt)) {
			std::cout << "Starting by using applying the previous Project to the triangulation\n";
			steiner_count += previous_triangulation(ccdt, region_polygon);
		}

		if (method == "local") {
			int L;
			extract_local_parameters(parameters, L);
			steiner_count += local_search(ccdt, region_polygon, 5000, in_domain,c_rate);
		}
		else if (method == "SA") {
			double alpha, beta;
			int L;
			extract_sa_parameters(parameters, alpha, beta, L);
			simulated_annealing(ccdt, region_polygon, in_domain, steiner_count, alpha, beta, L,c_rate);
		}
		else if (method == "ant") {
			double alpha, beta, xi, psi, lambda;
			int kappa, L;
			extract_ant_parameters(parameters, alpha, beta, xi, psi, lambda, kappa, L);
			ant_colony(ccdt, region_polygon, alpha, beta, xi, psi, lambda, kappa, L, 0, c_rate);
		}
	}
	else{
		std::vector<int> closed_p;
		Polygon_2 Pcycle;
		int type = boundary_type(region_polygon, region_boundary, additional_constrains, closed_p);
	

		steiner_count+=simulated_annealing(ccdt, region_polygon, in_domain, steiner_count, 3.0, 0.2, 200,c_rate);

		CGAL::mark_domain_in_triangulation(ccdt, in_domain);
		obtuse_count=count_obtuse_faces(ccdt,in_domain);
		int saved=obtuse_count;
		bool randomed=false;
		//steiner_count+=reduce_random_local(ccdt,obtuse_count,region_polygon,randomed, c_rate);
		CGAL::mark_domain_in_triangulation(ccdt, in_domain);
		obtuse_count=count_obtuse_faces(ccdt,in_domain);
		if(randomed){
			std::cout<<"Random did reduce the obtuse triangle count from: "<<saved<<" to "<<obtuse_count<<std::endl; 
		}




		CGAL::mark_domain_in_triangulation(ccdt, in_domain);
		reduce_obtuse_by_flips(ccdt, in_domain);

		/*

		steiner_count+=local_search(ccdt, region_polygon, 1000, in_domain,c_rate);
		steiner_count=ant_colony(ccdt, region_polygon, 3.0, 0.2, 1.5, 0.8, 0.2, 3, 600, 0,c_rate);       
int ant_colony(Custom_CDT& cdt, Polygon_2 boundary, double alpha, double beta, double xi, double psi, double lambda, int kappa, int L, int steiner_points,double& c_rate) {
		steiner_count=ant_colony(ccdt, region_polygon, 4.0, 0.01, 1.0, 3.0, 0.2, 3, 1000, 0,c_rate);       
		bool randomed=false;
		steiner_count+=reduce_random_local(ccdt,obtuse_count,region_polygon,randomed, c_rate);
		CGAL::mark_domain_in_triangulation(ccdt, in_domain);
		obtuse_count = count_obtuse_faces(ccdt, in_domain);
	
	*/
	}

	CGAL::mark_domain_in_triangulation(ccdt, in_domain);
	reduce_obtuse_by_flips(ccdt, in_domain);
	obtuse_count = count_obtuse_faces(ccdt, in_domain);
		
	final_rate = (steiner_count > 0) ? abs(c_rate/steiner_count) : 0.0;

	std::cout<<"Convergence rate: "<<final_rate<<" Acumilated rate: "<<c_rate<<std::endl;
	std::cout << "Final obtuse count:" << obtuse_count <<"\nUsed "<<steiner_count<<" steiner points."<< std::endl;
	//---------------------------------------
	//CGAL::draw(ccdt,in_domain);

	//--------------------------------------
	std::vector<Point> initial_points;
				steiner_count+=1;
	for (size_t i = 0; i < points_x.size(); ++i) {
		initial_points.emplace_back(points_x[i], points_y[i]);
	}
	
	generate_output_json(ccdt, inst_iud, initial_points, region_polygon, params.output_file, method, parameters, obtuse_count,final_rate);
	return 0;
}



InputParameters read_inputs(int argc, char* argv[]) {
	if ((argc != 5) && (argc != 6)) {
		std::cout << "Usage: " << argv[0] << " -i /path/to/input.json -o /path/to/output.json (optional -preselected_params) " << argc << std::endl;
		exit(1);
	}

	InputParameters params;
	params.preselect = false;
	bool flag1 = false, flag2 = false;

	for (int i = 1; i < argc; i++) {
		if (strcmp(argv[i], "-i") == 0) {
			params.input_file = argv[i + 1];
			i++;
			flag1 = true;
		} else if (strcmp(argv[i], "-o") == 0) {
			params.output_file = argv[i + 1];
			i++;
			flag2 = true;
		} else if (strcmp(argv[i], "-preselected_params") == 0) {
			params.preselect = true;
		}
	}

	if (!flag1 || !flag2 || (params.preselect == false && argc != 5)) {
		std::cerr << "Error: invalid parameters. Run ./" << argv[0] << " if you want to see expected input.\n";
		if (!flag1) std::cerr << "Missing input flag\n";
		if (!flag2) std::cerr << "Missing output flag\n";
		exit(1);
	}

	if (params.input_file.empty() || params.output_file.empty()) {
		std::cerr << "Error: Missing file path(s)." << std::endl;
		exit(1);
	}

	try {
		boost::property_tree::read_json(params.input_file, params.pt);
	}
	catch (const boost::property_tree::json_parser_error& e) {
		std::cerr << "Error: reading JSON file: " << e.what() << std::endl;
		exit(1);
	}

	return params;
}
