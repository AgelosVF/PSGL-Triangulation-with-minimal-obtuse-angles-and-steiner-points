#ifndef BOOST_UTILS_HPP
#define BOOST_UTILS_HPP

#include <boost/property_tree/ptree_fwd.hpp>
#include <vector>
#include <boost/property_tree/ptree.hpp>

std::vector<int> extract_x_coordinates(const boost::property_tree::ptree &pt);
std::vector<int> extract_y_coordinates(const boost::property_tree::ptree &pt);
std::vector<int> extract_region_boundary(const boost::property_tree::ptree &pt);
std::vector<std::pair<int,int>> extract_additional_constrains(const boost::property_tree::ptree &pt);

#endif
