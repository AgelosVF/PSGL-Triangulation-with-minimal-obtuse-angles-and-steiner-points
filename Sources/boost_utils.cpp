#include "../Header Files/boost_utils.hpp"
#include <iostream>
#include <vector>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>
#include <boost/optional/optional.hpp>
#include <boost/property_tree/ptree_fwd.hpp>

//Function that extracts points_x from the tree that has been created from the JSON file using boost
std::vector<int> extract_x_coordinates(const boost::property_tree::ptree &pt){
	std::vector<int> x_coordinates;
	boost::optional<const boost::property_tree::ptree&> points_x_opt = pt.get_child_optional("points_x");
	if (!points_x_opt){
		std::cerr<<"Error: Missing 'points_x' from JSON file"<<std::endl;
		exit(-1);
	}
	BOOST_FOREACH(const boost::property_tree::ptree::value_type &v, *points_x_opt)
		x_coordinates.push_back(v.second.get_value<int>());
	return x_coordinates;
}

std::string extract_method(const boost::property_tree::ptree &pt) {
    boost::optional<std::string> method_opt = pt.get_optional<std::string>("method");
    if (!method_opt) {
        std::cerr << "Error: Missing 'method' field from JSON file" << std::endl;
        exit(-1);
    }
    return method_opt.get();
}
bool extract_delaunay(const boost::property_tree::ptree &pt) {
    boost::optional<bool> delaunay_opt = pt.get_optional<bool>("delaunay");
    if (!delaunay_opt) {
        std::cerr << "Error: Missing 'delaunay' field from JSON file" << std::endl;
        exit(-1);
    }
    return delaunay_opt.get();
}

boost::property_tree::ptree extract_parameters(const boost::property_tree::ptree &pt) {
    boost::optional<const boost::property_tree::ptree&> params_opt = pt.get_child_optional("parameters");
    if (!params_opt) {
        std::cerr << "Error: Missing 'parameters' field from JSON file" << std::endl;
        exit(-1);
    }
    return params_opt.get();
}


void extract_sa_parameters(const boost::property_tree::ptree &parameters, double &alpha, double &beta, int &L) {
    try {
        alpha = parameters.get<double>("alpha");
        beta = parameters.get<double>("beta");
        L = parameters.get<int>("L");
    } catch (const boost::property_tree::ptree_error &e) {
        std::cerr << "Error extracting 'sa' parameters: " << e.what() << std::endl;
        exit(-1);
    }
}

void extract_local_parameters(const boost::property_tree::ptree &parameters, int &L) {
    try {
        L = parameters.get<int>("L");
    } catch (const boost::property_tree::ptree_error &e) {
        std::cerr << "Error extracting 'local' parameters: " << e.what() << std::endl;
        exit(-1);
    }
}
void extract_ant_parameters(const boost::property_tree::ptree &parameters, double &alpha, double &beta, double &xi, double &psi, double &lambda, int &kappa, int &L) {
    try {
        alpha = parameters.get<double>("alpha");
        beta = parameters.get<double>("beta");
        xi = parameters.get<double>("xi");
        psi = parameters.get<double>("psi");
        lambda = parameters.get<double>("lambda");
        kappa = parameters.get<int>("kappa");
        L = parameters.get<int>("L");
    } catch (const boost::property_tree::ptree_error &e) {
        std::cerr << "Error extracting 'ant' parameters: " << e.what() << std::endl;
        exit(-1);
    }
}
//Function that extracts points_y from the tree that has been created from the JSON file using boost
std::vector<int> extract_y_coordinates(const boost::property_tree::ptree &pt){
	std::vector<int> y_coordinates;
	boost::optional<const boost::property_tree::ptree&> points_y_opt = pt.get_child_optional("points_y");
	if (!points_y_opt){
		std::cerr<<"Error: Missing 'points_y' from JSON file"<<std::endl;
		exit(-1);
	}
	BOOST_FOREACH(const boost::property_tree::ptree::value_type &v, *points_y_opt)
		y_coordinates.push_back(v.second.get_value<int>());
	return y_coordinates;
}

//Function that extracts the region_boundary from the JSON file
std::vector<int> extract_region_boundary(const boost::property_tree::ptree &pt){
	std::vector<int> region_boundary;
	boost::optional<const boost::property_tree::ptree&> bound_opt=pt.get_child_optional("region_boundary");
	if(!bound_opt){
		std::cerr<<"Error Mising 'region boundary' from JSON file"<<std::endl;
		exit(-1);
	}
	BOOST_FOREACH(const boost::property_tree::ptree::value_type &v, *bound_opt){
		region_boundary.push_back(v.second.get_value<int>());
	}

	return region_boundary;// the region boundary is in the form of indices to the points that form the polygon (1,2,4,6). Still need to get the cordinates and make them to pairs.
}

std::vector<std::pair<int,int>> extract_additional_constrains(const boost::property_tree::ptree &pt){
	boost::optional<const boost::property_tree::ptree&> opt_constr=pt.get_child_optional("additional_constraints");
	if(!opt_constr){
		std::cerr<<"Error: missing \"additional_constraints\" field from JSON file"<<std::endl;
		exit(-1);
	}
	std::vector<std::pair<int,int>> add_constr;
	int first,second;
	BOOST_FOREACH(const boost::property_tree::ptree::value_type &v,*opt_constr){
		first=v.second.front().second.get_value<int>();
		second=v.second.back().second.get_value<int>();
		add_constr.emplace_back(first,second);
	}
	return add_constr;
}


