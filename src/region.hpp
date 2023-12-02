#ifndef REGION_HPP
#define REGION_HPP

#include <array>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>
#include <unordered_map>
#include <map>
#include <set>
#include <chrono>
#include <iomanip>
#include <iterator>
#include <filesystem>
#include <list>
#include "triangulation.hpp"
#include "criterion.hpp"

typedef std::pair<double, double> point;

// Base class to define the regions for the refinement
class Region {
        
        public:

        Triangulation &tr;
        std::set<int>::iterator it;
        int number_vertices = 0;
        Criterion &criterion;

        public:

        Region(Triangulation &triangulation, Criterion &criterion) : tr(triangulation), it(triangulation.triangle_list.begin()),
            number_vertices(triangulation.n_vertices), criterion(criterion) {}


        // Returns a pair containing a triangle with its metric. If there is no triangle left in the region, returns {-1, 0}
        virtual std::pair<int, double> get_pair() {
            std::pair<bool, double> toBeRefined_metric_pair;
            std::pair<int, double> pair_to_return;
            while (it != tr.triangle_list.end()) {
                if (in_region(*it)) {
                    toBeRefined_metric_pair = criterion.has_to_be_refined(*it);
                    if (toBeRefined_metric_pair.first) {
                        pair_to_return = {*it, toBeRefined_metric_pair.second};
                        it++;
                        return pair_to_return;
                    }
                }
                it++;
            }
            return {-1, 0};
        }


        // Method to be overriden by the specific region. Determines if a triangle is inside the refinement region.
        virtual bool in_region(int triangle) {
            return true;
        }


        // Print the refinement region
        virtual void print(std::string filename) {
            std::ofstream out(filename);

            //  out<<"{ appearance  {+edge +face linewidth 2} LIST\n";
            out<<"OFF"<<std::endl;
            //num_vertices num_polygons 0
            out<<std::setprecision(15)<<tr.n_vertices<<" "<<tr.n_faces<<" 0"<<std::endl;
            //print nodes
            for(int i = 0; i < tr.n_vertices; i++)
                out<<tr.get_PointX(i)<<" "<<tr.get_PointY(i)<<" 0"<<std::endl; 
            //print triangles
            for(int i: tr.triangle_list){
                out<<3;
                int curr_edge = i;
                for (int j = 0; j < 3; j++) {
                    out<<" "<<tr.origin(curr_edge);
                    curr_edge = tr.next(curr_edge);
                }
                out<<std::endl;
            }
            out.close();
        }
    };


    class Circle: public Region {

        point center;
        double radius;

        public:

        Circle(point center, double radius, Triangulation &triangulation, Criterion &criterion) 
            : Region(triangulation, criterion), center(center), radius(radius) {number_vertices = 50;}

        bool in_region(int triangle) {
            point p;
            int current_edge = triangle;
            for (int i=0; i<3; i++) {
                int v = tr.origin(current_edge);
                p = {tr.get_PointX(v), tr.get_PointY(v)};

                // Checks if the vertex is inside the circle
                if (tr.distance(center, p) <= radius) {
                    return true;
                }
                current_edge = tr.next(current_edge);
            }
            return false;
        }

        // Print the refinement region
        void print(std::string filename) {
            std::ofstream out(filename);
            int clockwise = -1;
            
            out<<"OFF"<<std::endl;
            //num_vertices num_polygons 0
            out<<std::setprecision(15)<<number_vertices<<" "<<1<<" 0"<<std::endl;
            double angle = 0;
            for (int i = 0; i<number_vertices; i++) {
                out << std::to_string(radius*std::cos(angle) + center.first) << " "
                    << std::to_string(radius*std::sin(angle) + center.second)
                    << " 0" << std::endl;
                angle += clockwise * M_PI * 2/number_vertices;
            }

            out <<number_vertices;
            for (int i = 0; i<number_vertices; i++) {
                out << " " << i;
            }
            out<<std::endl;

            out.close();
        }
    };


    class Rectangle: public Region {

        point center;
        double min_x, max_x;
        double min_y, max_y;

        public:

        Rectangle(point center, double width, double height, Triangulation &triangulation, Criterion &criterion) 
            : Region(triangulation, criterion), center(center) {
                min_x = center.first - width/2;
                max_x = center.first + width/2;
                min_y = center.second - height/2;
                max_y = center.second + height/2;
                number_vertices = 4;
            }

        bool in_region(int triangle) {
            int current_edge = triangle;
            for (int i=0; i<3; i++) {
                int v = tr.origin(current_edge);
                double x = tr.get_PointX(v);
                double y = tr.get_PointY(v);

                // Checks if the vertex is inside the rectangle defined by the coordinates.
                if (min_x <= x && x <= max_x && min_y <= y && y <= max_y) {
                    return true;
                }
                current_edge = tr.next(current_edge);
            }
            return false;
        }

        // Print the rectangle
        void print(std::string filename) {
            std::ofstream out(filename);
            bool clockwise = false;

            //  out<<"{ appearance  {+edge +face linewidth 2} LIST\n";
            out<<"OFF"<<std::endl;
            //num_vertices num_polygons 0
            out<<std::setprecision(15)<<"4 1 0"<<std::endl;
            //print nodes
            out << std::to_string(min_x) << " " << std::to_string(min_y) << " 0" << std::endl;
            out << std::to_string(max_x) << " " << std::to_string(min_y) << " 0" << std::endl;
            out << std::to_string(max_x) << " " << std::to_string(max_y) << " 0" << std::endl;
            out << std::to_string(min_x) << " " << std::to_string(max_y) << " 0" << std::endl;
            if (clockwise) {
                out << "4 3 2 1 0" << std::endl;
            }
            else {
                out << "4 0 1 2 3" << std::endl;
            }
            out.close();
        }
    };


    class Polygon: public Region {

        std::vector<point> points;
        std::vector<int> edges;
        std::string filename;
        double min_x, max_x;
        double min_y, max_y;

        public:

        Polygon(std::string name, Triangulation &triangulation, Criterion &criterion) 
            : Region(triangulation, criterion), filename(name) {

                // Reads the polygon defined in the given file (with .off format).
                int number_vertices = 0;
                int number_faces;
                std::string line;
                std::ifstream offfile(name);
                double a1, a2, a3;
                std::string tmp;
                if (offfile.is_open())
                {
                    //Check first line is a OFF file
                    while (std::getline(offfile, line)){ //add check boundary vertices flag
                        std::istringstream(line) >> tmp;
                        //std::cout<<"tmp: "<<tmp<<std::endl;
                        if (tmp[0] != '#' && !isWhitespace(line))
                        {
                            if(tmp[0] == 'O' && tmp[1] == 'F' && tmp[2] == 'F') //Check if the format is OFF
                                break;
                            else{
                                throw std::runtime_error("The file " + name + " is not an OFF file");
                            }
                        }
                    }

                    //Read the number of vertices and faces
                    while (std::getline(offfile, line)){ //add check boundary vertices flag
                        std::istringstream(line) >> tmp;
                    // std::cout<<"tmp: "<<tmp<<std::endl;
                        if (tmp[0] != '#' && !isWhitespace(line))
                        { 
                            std::istringstream(line) >> number_vertices >> number_faces;

                            points.reserve(number_vertices);
                            if (number_faces != 1) {
                                throw std::runtime_error("Only one polygon can be specified.");
                            }
                            break;       
                        }
                    }

                    //Read vertices
                    int index = 0;
                    while (index < number_vertices && std::getline(offfile, line) )
                    {
                        std::istringstream(line) >> tmp;
                        // std::cout<<"tmp: "<<tmp<<std::endl;
                        if (tmp[0] != '#' && !isWhitespace(line))
                        {
                            std::istringstream(line) >> a1 >> a2 >> a3;
                            if (points.size() == 0) {
                                min_x = a1;
                                max_x = a1;
                                min_y = a2;
                                max_y = a2;
                            }
                            else {
                                if (a1 < min_x) {
                                    min_x = a1;
                                }
                                else if (a1 > max_x) {
                                    max_x = a1;
                                }
                                if (a2 < min_y) {
                                    min_y = a2;
                                }
                                else if (a2 > max_y) {
                                    max_y = a2;
                                }
                            }
                            points.push_back({a1, a2});
                            index++;
                        }
                    }
                    //Read faces
                    
                    int length, e1;
                    index = 0;
                    while (index < 1 && std::getline(offfile, line) )
                    {
                        std::istringstream(line) >> tmp;
                        // std::cout<<"tmp: "<<tmp<<std::endl;
                        if (tmp[0] != '#' && !isWhitespace(line))
                        {
                            std::istringstream list(line);
                            list >> length;
                            for (int i = 0; i < length; i++) {
                                list >> e1;
                                edges.push_back(e1);
                            }
                            index++;
                        }
                    }
                }
            }


        bool in_region(int triangle) {
            int current_edge = triangle;
            for (int i=0; i<3; i++) {
                int v = tr.origin(current_edge);
                double x = tr.get_PointX(v);
                double y = tr.get_PointY(v);
                // Determines if the the vertex is inside the region using the ray crossings algorithm
                if (min_x <= x && x <= max_x && min_y <= y && y <= max_y) {
                    int count = count_intersection(x, y);
                    if (count % 2 == 1) {
                        return true;
                    }
                }
                current_edge = tr.next(current_edge);
            }
            return false;
        }

        
        // Counts the intersections between the horizontal ray formed with the point given (x, y) and the polygon. 
        // This is the main part of the ray-crossings algorithm, to check if the point (x, y) lies inside the polygon.
        int count_intersection(double x, double y) {
            int count = 0;
            std::pair<bool, bool> intersects_edge = tr.intersects({x, y}, {max_x + 10, y}, points.at(edges.front()), points.at(edges.back()));
            if (intersects_edge.second) {
                return 1;
            }
            if (intersects_edge.first) {
                if (std::max(points.at(edges.front()).second, points.at(edges.back()).second) > y) {
                    count++;
                }
            }
            for (int index = 0; index < edges.size() - 1; index ++) {
                intersects_edge = tr.intersects({x, y}, {max_x + 10, y}, points.at(index), points.at(index+1));
                if (intersects_edge.second) {
                    return 1;
                }
                if (intersects_edge.first) {
                    if (std::max(points.at(index).second, points.at(index+1).second) > y) {
                        count++;
                    }
                }
            }
            return count;
        }


        void print(std::string new_filename) {
            std::filesystem::copy(filename, new_filename);
        }
    };

#endif