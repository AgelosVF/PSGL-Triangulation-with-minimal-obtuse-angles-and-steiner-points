#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include "../Header Files/triangulation_utils.hpp"
#include <random>
#include "../Header Files/CGAL_CUSTOM_CONSTRAINED_DELAUNAY_TRIANGULATION_2.h"
#include "CGAL/draw_triangulation_2.h"

Point random_point_in_face_gaussian(Face_handle face) {
    // Calculate centroid of the face
    Point p1 = face->vertex(0)->point();
    Point p2 = face->vertex(1)->point();
    Point p3 = face->vertex(2)->point();
    // Calculate the centroid of the face
    Point centroid=CGAL::centroid(p1,p2,p3);
    double cent_x=(CGAL::to_double(centroid.x()));
    double cent_y=(CGAL::to_double(centroid.y()));

    // Create random number generator
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::normal_distribution<> dist_x(cent_x, 0.1);
    std::normal_distribution<> dist_y(cent_y, 0.1);

    // Generate random point using Gaussian distribution
    double x = dist_x(gen);
    double y = dist_y(gen);
    Point random(x,y);
    // Check if the point is inside the face
    if (!point_inside_triangle(p1, p2, p3, centroid)) {
        //if the point is not inside the face, try again
        return random_point_in_face_gaussian(face);
    } 

    return Point(x, y);
}

void random_and_flips(Custom_CDT& cdt, Face_handle face, boost::associative_property_map<std::unordered_map<Face_handle, bool>>& in_domain){
    Point random=random_point_in_face_gaussian(face);
    std::cout<<"Random point="<<random.x()<<" , "<<random.y()<<std::endl;
    Point p0=face->vertex(0)->point();
    Point p1=face->vertex(1)->point();
    Point p2=face->vertex(2)->point();
    std::cout<<"p0="<<p0.x()<<" , "<<p0.y()<<std::endl;
    std::cout<<"p1="<<p1.x()<<" , "<<p1.y()<<std::endl;
    std::cout<<"p2="<<p2.x()<<" , "<<p2.y()<<std::endl;
    Face_handle s_face=find_face(cdt,face);
    Vertex_handle v=cdt.insert_no_flip(random,face);
    reduce_obtuse_by_flips(cdt,in_domain); 

}


