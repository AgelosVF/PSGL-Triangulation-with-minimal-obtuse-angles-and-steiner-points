#ifndef BOOST_UTILS_HPP
#define BOOST_UTILS_HPP

#include <boost/property_tree/ptree_fwd.hpp>
#include <vector>
#include <boost/property_tree/ptree.hpp>

std::vector<int> extract_x_coordinates(const boost::property_tree::ptree &pt);
std::vector<int> extract_y_coordinates(const boost::property_tree::ptree &pt);
std::vector<int> extract_region_boundary(const boost::property_tree::ptree &pt);
std::vector<std::pair<int,int>> extract_additional_constrains(const boost::property_tree::ptree &pt);

std::string extract_method(const boost::property_tree::ptree &pt);
bool extract_delaunay(const boost::property_tree::ptree &pt);
boost::property_tree::ptree extract_parameters(const boost::property_tree::ptree &pt);
void extract_sa_parameters(const boost::property_tree::ptree &parameters, double &alpha, double &beta, int &L);

void extract_local_parameters(const boost::property_tree::ptree &parameters, int &L);
void extract_ant_parameters(const boost::property_tree::ptree &parameters, double &alpha, double &beta, double &xi, double &psi, double &lambda, int &kappa, int &L);
#endif
