#include <algorithm>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <chrono>

#include "src/polylla.hpp"
#include "src/triangulation.hpp"
#include "src/parser.hpp"
#include "src/queue.hpp"
#include "src/criterion.hpp"
#include "src/refiner.hpp"
#include "src/region.hpp"


// Main program to parse and run the algorithms
int main(int argc, char **argv) {

    try {
        std::chrono::high_resolution_clock::time_point start_time = std::chrono::high_resolution_clock::now();;
        Parser parser(argc, argv);
        
        std::string output = std::any_cast<std::string>(parser.getOption("-output"));
        if (output.back() != '/') {
            output += "/";
        }

        std::vector<std::string_view> args = parser.getArgs();
        Triangulation *triangulation;
        // Get the number of files specified to instantiate the proper Triangulation object
        if (parser.getNumberFiles() == 1) {
            triangulation = new Triangulation(args.at(0).data());
        }
        else {
            triangulation = new Triangulation(args.at(0).data(), args.at(1).data(), args.at(2).data());
        }

        // Stores if the triangulation was modified in order to know if the new triangulation has to be written
        bool modified_triangulation = false;

        // If the triangulation has to be refined
        bool refine = false;

        // Checks if the given triangulation is Delaunay, if working with Delaunay triangulations
        bool is_delaunay = false;
        bool delaunay = !parser.hasOption("-not_delaunay");
        if (delaunay) {
            is_delaunay = triangulation->is_delaunay_triangulation();
            if (!is_delaunay) {
                std::cout << "The given triangulation is not Delaunay." << std::endl;
                modified_triangulation = true;
                triangulation->repairDelaunay();
            }
            else {
                std::cout << "The triangulation is Delaunay." << std::endl;
            }
        }

        // In case the user explicitly wants to check if the triangulation is Delaunay
        bool check_delaunay = parser.hasOption("-check_delaunay");
        if (parser.hasOption("-check_delaunay") && !delaunay) {
            is_delaunay = triangulation->is_delaunay_triangulation();
            if (is_delaunay) {
                std::cout << "The triangulation is Delaunay." << std::endl;
            }
            else {
                std::cout << "The triangulation is Not Delaunay." << std::endl;
            }
        }

        // The type of insertion to perform
        bool bisection = !parser.hasOption("-centroid");

        // The criterion to use
        Criterion *criterion;
        if (parser.hasOption("-factor")) {
            refine = true;
            modified_triangulation = true;
            criterion = new AverageAreaCriterion(*triangulation);
        }
        if (parser.hasOption("-max_area")) {
            refine = true;
            modified_triangulation = true;
            double max_area = std::any_cast<double>(parser.getOption("-max_area"));
            criterion = new MaxAreaCriterion(*triangulation, max_area);
        }

        // The region specified
        Region *region;
        if (parser.hasOption("-region")) {
            std::string regionType = std::any_cast<std::string>(parser.getOption("-region"));
            
            if (regionType.compare("-polygon") == 0) {
                std::string filename = std::any_cast<std::string>(parser.getOption(regionType));
                region = new Polygon(filename, *triangulation, *criterion);
            }
            else { //The region has to be a rectangle or a circle
                std::vector<double> values = std::any_cast<std::vector<double>>(parser.getOption(regionType));
                std::pair<double, double> center = {values.at(0), values.at(1)};

                if (regionType.compare("-rectangle") == 0) {
                    region = new Rectangle(center, values.at(2), values.at(3), *triangulation, *criterion);
                }
                else {
                    region = new Circle(center, values.at(2), *triangulation, *criterion);
                }
            }
            
            if (parser.hasOption("-print_region")) {
                region->print(output + "region" + regionType + ".off");
            }
        }
        else {
            region = new Region(*triangulation, *criterion);
        }

        std::chrono::high_resolution_clock::time_point start_without_queue;
        std::chrono::high_resolution_clock::time_point start_with_queue;
        std::chrono::high_resolution_clock::time_point stop;
        
        // If the triangulation has to be refined
        if (refine) {
            start_with_queue = std::chrono::high_resolution_clock::now();
            if (parser.hasOption("-factor")) {
                
                QueueAverage queue(*region, std::any_cast<double>(parser.getOption("-factor")));
                Refiner refiner(*triangulation, queue, delaunay, bisection);

                // Sets the directory to print every iteration of the algorithm, if specified
                if (parser.hasOption("-print_iter")) {
                    std::string directory = std::any_cast<std::string>(parser.getOption("-print_iter"));
                    if (directory.back() != '/') {
                        directory += "/";
                    }
                    refiner.setDirectory(directory);
                }
                start_without_queue = std::chrono::high_resolution_clock::now();
                refiner.refine_area();
                stop = std::chrono::high_resolution_clock::now();
            }
            else {
                QueueThreshold queue(*region);
                Refiner refiner(*triangulation, queue, delaunay, bisection);
                if (parser.hasOption("-print_iter")) {
                    std::string directory = std::any_cast<std::string>(parser.getOption("-print_iter"));
                    if (directory.back() != '/') {
                        directory += "/";
                    }
                    refiner.setDirectory(directory);
                }
                start_without_queue = std::chrono::high_resolution_clock::now();
                refiner.refine_area();
                stop = std::chrono::high_resolution_clock::now();
            }
        }
        
        // If the triangulation was modified, is written in the output directory
        if (modified_triangulation) {
            auto t = std::time(nullptr);
            auto tm = *std::localtime(&t);

            std::ostringstream oss;
            oss << std::put_time(&tm, "%d-%m-%Y %H-%M-%S");
            auto str = oss.str();

            auto duration_with_queue = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start_with_queue);
            auto duration_without_queue = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start_without_queue);

            if (!parser.hasOption("-dont_write")) {
                triangulation->print_OFF(output + "refined-triangulation" + str + ".off");
            }

            std::cout << "  Refinement time: " << duration_with_queue.count() << " milliseconds." << std::endl;
            //triangulation->print_OFF(output+"refined-triangulation" + str + ".off");
        }

        int number_triangles1 = triangulation->get_number_triangles();
        int number_v1 = triangulation->get_number_vertices();
        int number_edges = triangulation->get_number_edges();

        std::vector<double> stats = triangulation->get_metrics();
        auto total_stop = std::chrono::high_resolution_clock::now();
        auto total_time = std::chrono::duration_cast<std::chrono::milliseconds>(total_stop - start_time);
        std::cout << "  Number of triangles: " << number_triangles1 << std::endl;
        std::cout << "  Number of vertices: " << number_v1 << std::endl;
        std::cout << "  Number of edges: " << number_edges << std::endl;
        std::cout << "  Max area: " << stats.at(0) << std::endl;
        std::cout << "  Min area: " << stats.at(1) << std::endl;
        std::cout << "  Average area: " << stats.at(2) << std::endl;
        std::cout << "  Max angle (degrees): " << stats.at(3) * 180/M_PI << std::endl;
        std::cout << "  Min angle (degrees): " << stats.at(4) * 180/M_PI << std::endl;
        std::cout << "  Total time: " << total_time.count() << " milliseconds." << std::endl;
        std::cout << "\n";
        
        // Polylla algorithm
        if (parser.hasOption("-polylla")) {
            Polylla mesh(triangulation);
            
            std::string filename = "";
            if (parser.hasOption("-filename")) {
                filename = std::any_cast<std::string>(parser.getOption("-filename"));
                filename += "-";
            }

            mesh.print_OFF(output + filename + "polylla-mesh.off");

            mesh.find_non_convex();

            mesh.print_reflex(output + filename + "polylla-reflex.off");

            std::chrono::high_resolution_clock::time_point start;
            std::chrono::high_resolution_clock::time_point stop;
            start = std::chrono::high_resolution_clock::now();
            mesh.repair_reflex();
            stop = std::chrono::high_resolution_clock::now();

            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
            std::cout << "  Duration millisecond: " << duration.count() << std::endl;

            mesh.print_OFF(output + filename + "polylla-mesh-repaired.off");
        }
        return EXIT_SUCCESS;
        
    }
    catch (const std::exception &x) {
        std::cout << x.what() << std::endl;
        return EXIT_FAILURE;    
    }
}