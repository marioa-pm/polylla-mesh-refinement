/* Polygon mesh generator
//POSIBLE BUG: el algoritmo no viaja por todos los halfedges dentro de un poligono, 
    //por lo que pueden haber semillas que no se borren y tener poligonos repetidos de output
*/

#ifndef POLYLLA_HPP
#define POLYLLA_HPP
#define _USE_MATH_DEFINES
#include <algorithm>


#include <array>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include "triangulation.hpp"
#include <chrono>
#include <iomanip>

#define print_e(eddddge) eddddge<<" ( "<<mesh_input->origin(eddddge)<<" - "<<mesh_input->target(eddddge)<<") "


class Polylla
{
private:
    typedef std::vector<int> _polygon; 
    typedef std::vector<char> bit_vector; 

    struct isReflex {
        bool is_reflex;
        double angle;
    };

    Triangulation *mesh_input; // Halfedge triangulation
    Triangulation *mesh_output;
    std::vector<int> output_seeds; //Seeds of the polygon

    //std::vector<int> triangles; //True if the edge generated a triangle CHANGE!!!!

    bit_vector max_edges; //True if the edge i is a max edge
    bit_vector frontier_edges; //True if the edge i is a frontier edge
    std::vector<int> seed_edges; //Seed edges that generate polygon simple and non-simple
    std::vector<int> seed_to_repair_reflex;

    // Auxiliary array used during the barrier-edge elimination
    std::vector<int> triangle_list;
    bit_vector seed_bet_mark;

    //Statistics
    int m_polygons = 0; //Number of polygons
    int n_frontier_edges = 0; //Number of frontier edges
    int n_barrier_edge_tips = 0; //Number of barrier edge tips
    int n_polygons_to_repair = 0;
    int n_polygons_added_after_repair = 0;
    int n_polygons_to_repair_reflex = 0;
    int n_polygons_added_after_repair_reflex = 0;

    // Times
    double t_label_max_edges = 0;
    double t_label_frontier_edges = 0;
    double t_label_seed_edges = 0;
    double t_traversal_and_repair = 0;
    double t_traversal = 0;
    double t_repair = 0;
    
public:

    Polylla() {}; //Default constructor

    //Constructor with triangulation
    Polylla(Triangulation *input_mesh){
        this->mesh_input = input_mesh;
        mesh_output = new Triangulation(*mesh_input);
        construct_Polylla();
    }


    //Constructor from a OFF file
    Polylla(std::string off_file){

        this->mesh_input = new Triangulation(off_file);
        mesh_output = new Triangulation(*mesh_input);
        construct_Polylla();
    }

    //Constructor from a node_file, ele_file and neigh_file
    Polylla(std::string node_file, std::string ele_file, std::string neigh_file){
        this->mesh_input = new Triangulation(node_file, ele_file, neigh_file);
        //call copy constructor
        mesh_output = new Triangulation(*mesh_input);
        construct_Polylla();
    }


//Constructor random data construictor
    Polylla(int size){
        this->mesh_input = new Triangulation(size);
        mesh_output = new Triangulation(*mesh_input);
        construct_Polylla();
    }

    ~Polylla() {
        //triangles.clear(); 
        max_edges.clear(); 
        frontier_edges.clear();
        seed_edges.clear(); 
        seed_bet_mark.clear();
        triangle_list.clear();
        delete mesh_input;
        delete mesh_output;
    }

    void construct_Polylla(){

        max_edges = bit_vector(mesh_input->halfEdges(), false);
        frontier_edges = bit_vector(mesh_input->halfEdges(), false);
        //triangles = mesh_input->get_Triangles(); //Change by triangle list
        seed_bet_mark = bit_vector(this->mesh_input->halfEdges(), false);

        //terminal_edges = bit_vector(mesh_input->halfEdges(), false);
        //seed_edges = bit_vector(mesh_input->halfEdges(), false);
        
        std::cout<<"Creating Polylla..."<<std::endl;
        //Label max edges of each triangle
        int largo_triangle_list = mesh_input->triangle_list.size();
        auto t_start = std::chrono::high_resolution_clock::now();
        int numero_max_edges = 0;
        for(int i: mesh_input->triangle_list){
            numero_max_edges++;
            max_edges[label_max_edge(i)] = true;
        }

        auto t_end = std::chrono::high_resolution_clock::now();
        t_label_max_edges = std::chrono::duration<double, std::milli>(t_end-t_start).count();
        std::cout<<"Labered max edges in "<<t_label_max_edges<<" ms"<<std::endl;

        t_start = std::chrono::high_resolution_clock::now();
        //Label frontier edges
        int largo_halfedges = mesh_input->halfEdges();
        for (std::size_t e = 0; e < mesh_input->halfEdges(); e++){

            halfEdge curr = mesh_input->HalfEdges.at(e);

            int ori = mesh_input->origin(e);
            double orx = mesh_input->get_PointX(ori);
            double ory = mesh_input->get_PointY(ori);
            int dest = mesh_input->target(e);
            double destx = mesh_input->get_PointX(dest);
            double desty = mesh_input->get_PointY(dest);

            if(is_frontier_edge(e)){
                frontier_edges[e] = true;
                n_frontier_edges++;
            }
        }

        t_end = std::chrono::high_resolution_clock::now();
        t_label_frontier_edges = std::chrono::duration<double, std::milli>(t_end-t_start).count();
        std::cout<<"Labeled frontier edges in "<<t_label_frontier_edges<<" ms"<<std::endl;
        
        t_start = std::chrono::high_resolution_clock::now();
        //label seeds edges,
        for (std::size_t e = 0; e < mesh_input->halfEdges(); e++){

            if(mesh_input->is_interior_face(e) && is_seed_edge(e)) {
                seed_edges.push_back(e);
            }
        }

        int numero_seed = seed_edges.size();
        t_end = std::chrono::high_resolution_clock::now();
        t_label_seed_edges = std::chrono::duration<double, std::milli>(t_end-t_start).count();
        std::cout<<"Labeled seed edges in "<<t_label_seed_edges<<" ms"<<std::endl; 

        //Travel phase: Generate polygon mesh
        int polygon_seed;
        //Foreach seed edge generate polygon
        t_start = std::chrono::high_resolution_clock::now();
        int iter = 0;
        for(auto &e : seed_edges){
            int ori = mesh_input->origin(e);
            int dest = mesh_input->target(e);
            polygon_seed = travel_triangles(e);
            //output_seeds.push_back(polygon_seed);
            if(!has_BarrierEdgeTip(polygon_seed)){ //If the polygon is a simple polygon then is part of the mesh
                output_seeds.push_back(polygon_seed);
            }else{ //Else, the polygon is send to reparation phase
                auto t_start_repair = std::chrono::high_resolution_clock::now();
                barrieredge_tip_reparation(polygon_seed);
                auto t_end_repair = std::chrono::high_resolution_clock::now();
                t_repair += std::chrono::duration<double, std::milli>(t_end_repair-t_start_repair).count();
            }
            //print_OFF("output/final-4/iter/step" + std::to_string(iter) + ".off");
            iter++;        
        }    
        t_end = std::chrono::high_resolution_clock::now();
        t_traversal_and_repair = std::chrono::duration<double, std::milli>(t_end-t_start).count();
        t_traversal = t_traversal_and_repair - t_repair;
        
        this->m_polygons = output_seeds.size();

        std::cout<<"Mesh with "<<m_polygons<<" polygons "<<n_frontier_edges/2<<" edges and "<<n_barrier_edge_tips<<" barrier-edge tips."<<std::endl;
        //mesh_input->print_pg(std::to_string(mesh_input->vertices()) + ".pg");             
    }


    void print_stats(std::string filename){
        //Time
        std::cout<<"Time to generate Triangulation: "<<mesh_input->get_triangulation_generation_time()<<" ms"<<std::endl;
        std::cout<<"Time to label max edges "<<t_label_max_edges<<" ms"<<std::endl;
        std::cout<<"Time to label frontier edges "<<t_label_frontier_edges<<" ms"<<std::endl;
        std::cout<<"Time to label seed edges "<<t_label_seed_edges<<" ms"<<std::endl;
        std::cout<<"Time to label total"<<t_label_max_edges+t_label_frontier_edges+t_label_seed_edges<<" ms"<<std::endl;
        std::cout<<"Time to traversal and repair "<<t_traversal_and_repair<<" ms"<<std::endl;
        std::cout<<"Time to traversal "<<t_traversal<<" ms"<<std::endl;
        std::cout<<"Time to repair "<<t_repair<<" ms"<<std::endl;
        std::cout<<"Time to generate polygonal mesh "<<t_label_max_edges + t_label_frontier_edges + t_label_seed_edges + t_traversal_and_repair<<" ms"<<std::endl;

        //Memory
        long long m_max_edges =  sizeof(decltype(max_edges.back())) * max_edges.capacity();
        long long m_frontier_edge = sizeof(decltype(frontier_edges.back())) * frontier_edges.capacity();
        long long m_seed_edges = sizeof(decltype(seed_edges.back())) * seed_edges.capacity();
        long long m_seed_bet_mar = sizeof(decltype(seed_bet_mark.back())) * seed_bet_mark.capacity();
        long long m_triangle_list = sizeof(decltype(triangle_list.back())) * triangle_list.capacity();
        long long m_mesh_input = mesh_input->get_size_vertex_half_edge();
        long long m_mesh_output = mesh_output->get_size_vertex_half_edge();
        long long m_vertices_input = mesh_input->get_size_vertex_struct();
        long long m_vertices_output = mesh_output->get_size_vertex_struct();

        std::ofstream out(filename);
        std::cout<<"Printing JSON file as "<<filename<<std::endl;
        out<<"{"<<std::endl;
        out<<"\"n_polygons\": "<<m_polygons<<","<<std::endl;
        out<<"\"n_frontier_edges\": "<<n_frontier_edges/2<<","<<std::endl;
        out<<"\"n_barrier_edge_tips\": "<<n_barrier_edge_tips<<","<<std::endl;
        out<<"\"n_half_edges\": "<<mesh_input->halfEdges()<<","<<std::endl;
        out<<"\"n_faces\": "<<mesh_input->faces()<<","<<std::endl;
        out<<"\"n_vertices\": "<<mesh_input->vertices()<<","<<std::endl;
        out<<"\"n_polygons_to_repair\": "<<n_polygons_to_repair<<","<<std::endl;
        out<<"\"n_polygons_added_after_repair\": "<<n_polygons_added_after_repair<<","<<std::endl;
        out<<"\"time_triangulation_generation\": "<<mesh_input->get_triangulation_generation_time()<<","<<std::endl;
        out<<"\"time_to_label_max_edges\": "<<t_label_max_edges<<","<<std::endl;
        out<<"\"time_to_label_frontier_edges\": "<<t_label_frontier_edges<<","<<std::endl;
        out<<"\"time_to_label_seed_edges\": "<<t_label_seed_edges<<","<<std::endl;
        out<<"\"time_to_label_total\": "<<t_label_max_edges+t_label_frontier_edges+t_label_seed_edges<<","<<std::endl;
        out<<"\"time_to_traversal_and_repair\": "<<t_traversal_and_repair<<","<<std::endl;
        out<<"\"time_to_traversal\": "<<t_traversal<<","<<std::endl;
        out<<"\"time_to_repair\": "<<t_repair<<","<<std::endl;
        out<<"\"time_to_generate_polygonal_mesh\": "<<t_label_max_edges + t_label_frontier_edges + t_label_seed_edges + t_traversal_and_repair<<","<<std::endl;
        out<<"\t\"memory_max_edges\": "<<m_max_edges<<","<<std::endl;
        out<<"\t\"memory_frontier_edge\": "<<m_frontier_edge<<","<<std::endl;
        out<<"\t\"memory_seed_edges\": "<<m_seed_edges<<","<<std::endl;
        out<<"\t\"memory_seed_bet_mar\": "<<m_seed_bet_mar<<","<<std::endl;
        out<<"\t\"memory_triangle_list\": "<<m_triangle_list<<","<<std::endl;
        out<<"\t\"memory_mesh_input\": "<<m_mesh_input<<","<<std::endl;
        out<<"\t\"memory_mesh_output\": "<<m_mesh_output<<","<<std::endl;
        out<<"\t\"memory_vertices_input\": "<<m_vertices_input<<","<<std::endl;
        out<<"\t\"memory_vertices_output\": "<<m_vertices_output<<","<<std::endl;
        out<<"\t\"memory_total\": "<<m_max_edges + m_frontier_edge + m_seed_edges + m_seed_bet_mar + m_triangle_list + m_mesh_input + m_mesh_output + m_vertices_input + m_vertices_output<<std::endl;
        out<<"}"<<std::endl;
        out.close();
    }


    //Print ale file of the polylla mesh
    void print_ALE(std::string filename){
        std::ofstream out(filename);
        _polygon poly;
        out<<"# domain type\nCustom\n";
        out<<"# nodal coordinates: number of nodes followed by the coordinates \n";
        out<<mesh_input->vertices()<<std::endl;
        //print nodes
        for(std::size_t v = 0; v < mesh_input->vertices(); v++)
            out<<std::setprecision(15)<<mesh_input->get_PointX(v)<<" "<<mesh_input->get_PointY(v)<<std::endl; 
        out<<"# element connectivity: number of elements followed by the elements\n";
        out<<this->m_polygons<<std::endl;
        //print polygons
        int size_poly;
        int e_curr;
        for(auto &e_init : output_seeds){
            size_poly = 1;
            e_curr = mesh_output->next(e_init);
            while(e_init != e_curr){
                size_poly++;
                e_curr = mesh_output->next(e_curr);
            }
            out<<size_poly<<" ";            

            out<<mesh_output->origin(e_init)<<" ";
            e_curr = mesh_output->next(e_init);
            while(e_init != e_curr){
                out<<mesh_output->origin(e_curr)<<" ";
                e_curr = mesh_output->next(e_curr);
            }
            out<<std::endl; 
        }
        //Print borderedges
        out<<"# indices of nodes located on the Dirichlet boundary\n";
        ///Find borderedges
        int b_curr, b_init = 0;
        for(std::size_t i = mesh_input->halfEdges()-1; i != 0; i--){
            if(mesh_input->is_border_edge(i)){
                b_init = i;
                break;
            }
        }
        out<<mesh_input->origin(b_init)<<" ";
        b_curr = mesh_input->prev(b_init);
        while(b_init != b_curr){
            out<<mesh_input->origin(b_curr)<<" ";
            b_curr = mesh_input->prev(b_curr);
        }
        out<<std::endl;
        out<<"# indices of nodes located on the Neumann boundary\n0\n";
        out<<"# xmin, xmax, ymin, ymax of the bounding box\n";
        double xmax = mesh_input->get_PointX(0);
        double xmin = mesh_input->get_PointX(0);
        double ymax = mesh_input->get_PointY(0);
        double ymin = mesh_input->get_PointY(0);
        //Search min and max coordinates
        for(std::size_t v = 0; v < mesh_input->vertices(); v++){
            //search range x
            if(mesh_input->get_PointX(v) > xmax )
                xmax = mesh_input->get_PointX(v);
            if(mesh_input->get_PointX(v) < xmin )
                xmin = mesh_input->get_PointX(v);
            //search range y
            if(mesh_input->get_PointY(v) > ymax )
                ymax = mesh_input->get_PointY(v);
            if(mesh_input->get_PointY(v) < ymin )
                ymin = mesh_input->get_PointY(v);
        }
        out<<xmin<<" "<<xmax<<" "<<ymin<<" "<<ymax<<std::endl;
        out.close();
    }

    //Print off file of the polylla mesh
    void print_OFF(std::string filename){
        std::ofstream out(filename);

      //  out<<"{ appearance  {+edge +face linewidth 2} LIST\n";
        out<<"OFF"<<std::endl;
        //num_vertices num_polygons 0
        out<<std::setprecision(15)<<mesh_input->vertices()<<" "<<output_seeds.size()<<" 0"<<std::endl;
        //print nodes
        for(std::size_t v = 0; v < mesh_input->vertices(); v++)
            out<<mesh_input->get_PointX(v)<<" "<<mesh_input->get_PointY(v)<<" 0"<<std::endl; 
        //print polygons
        int size_poly;
        int e_curr;
        for(auto &e_init : output_seeds){
            size_poly = 1;
            e_curr = mesh_output->next(e_init);
            while(e_init != e_curr){
                size_poly++;
                e_curr = mesh_output->next(e_curr);
            }
            out<<size_poly<<" ";            

            out<<mesh_output->origin(e_init)<<" ";
            e_curr = mesh_output->next(e_init);
            while(e_init != e_curr){
                out<<mesh_output->origin(e_curr)<<" ";
                e_curr = mesh_output->next(e_curr);
            }
            out<<std::endl; 
        }
      //  out<<"}"<<std::endl;
        out.close();
    }

    // Print only non-convex polygons.
    void print_reflex(std::string filename){
        std::ofstream out(filename);
        std::stringstream strfile;

        //print nodes
        strfile<<std::setprecision(15);
        for(std::size_t v = 0; v < mesh_input->vertices(); v++)
            strfile<<mesh_input->get_PointX(v)<<" "<<mesh_input->get_PointY(v)<<" 0"<<std::endl; 
        //print polygons
        int n_reflex_polygons = 0;
        for(auto &e_init : seed_to_repair_reflex){
            int e = e_init;
            int size_poly = 0;
            do {
                size_poly++;
                e = mesh_output->next(e);
            } while (e != e_init);

            n_reflex_polygons++;
            strfile<<size_poly<<" ";
            do {
                strfile<<mesh_output->origin(e)<<" ";
                e = mesh_output->next(e);
            } while (e_init != e);
            strfile<<std::endl;
        }
      //  out<<"}"<<std::endl;

        out<<"OFF"<<std::endl;
        //num_vertices num_polygons 0
        out<<std::setprecision(15)<<mesh_input->vertices()<<" "<<n_reflex_polygons<<" 0"<<std::endl;
        out<<strfile.rdbuf();
        out.close();
    }

    // Indicates if the polygon is in counterclockwise order
    // Input: any edge of the polygon
    bool is_polygon_counterclockwise(int e_init) {
        int v_init = mesh_output->origin(e_init);
        int v_min = v_init;
        int e_min = e_init;

        for (int e = mesh_output->next(e_init); e != e_init; e = mesh_output->next(e)) {
            int v_curr = mesh_output->origin(e);
            double x_curr = mesh_output->get_PointX(v_curr);
            double y_curr = mesh_output->get_PointY(v_curr);
            double y_min = mesh_output->get_PointY(v_min);
            if (y_curr < y_min || (y_curr == y_min && mesh_output->get_PointX(v_curr) < mesh_output->get_PointX(v_min))) {
                v_min = v_curr;
                e_min = e;
            }
        }

        return mesh_output->is_counterclockwise({mesh_output->origin(mesh_output->prev(e_min)), v_min, mesh_output->target(e_min)});
    }

    // Determines all the non convex polygons and stores the seeds
    void find_non_convex() {
        for(auto &e_init : output_seeds) {

            // Determines the orientation of the polygon
            bool counterclockwise = is_polygon_counterclockwise(e_init);
            int e = e_init;

            // Loop through every vertex in the polygon checking if the angle is reflex
            do {
                if (is_reflex_angle(e, counterclockwise)) {
                    seed_to_repair_reflex.push_back(e_init);
                    n_polygons_to_repair_reflex++;
                    break;
                }
                e = mesh_output->next(e);
            } while (e != e_init);
        }
    }

    // Loop through every non convex polygon repairing the reflex angles
    void repair_reflex(){
        for(auto &e : seed_to_repair_reflex){
            reflex_reparation(e);
            this->m_polygons = output_seeds.size();
            //this->print_OFF("output/reflex-cambio/iter/repaired" + std::to_string(i) + ".off");
        }
        std::cout << "  " << n_polygons_to_repair_reflex << " concave polygons." << std::endl;
        std::cout << "  " << n_polygons_added_after_repair_reflex << " repaired polygons." << std::endl;
    }
    

private:

    //Return true if the edge is terminal-edge or terminal border edge, 
    //but only selects one halfedge as terminal-edge, the halfedge with lowest index is selected
    bool is_seed_edge(int e){
        int twin = mesh_input->twin(e);

        bool is_terminal_edge = (mesh_input->is_interior_face(twin) &&  (max_edges[e] && max_edges[twin]) );
        bool is_terminal_border_edge = (mesh_input->is_border_edge(twin) && max_edges[e]);

        if( (is_terminal_edge && e < twin ) || is_terminal_border_edge){
            return true;
        }

        return false;
    }

    int Equality(double a, double b, double epsilon) {
        return fabs(a - b) < epsilon;
    }
    
    int GreaterEqualthan(double a, double b, double epsilon){
            return Equality(a,b,epsilon) || a > b;
    }

    //Label max edges of all triangles in the triangulation
    //input: edge e indicent to a triangle t
    //output: position of edge e in max_edges[e] is labeled as true
    int label_max_edge(const int e)
    {
        //Calculates the size of each edge of a triangle 
        double dist0 = mesh_input->distance(e);
        double dist1 = mesh_input->distance(mesh_input->next(e));
        double dist2 = mesh_input->distance(mesh_input->prev(e));
        //Find the longest edge of the triangle
        if(std::max({dist0, dist1, dist2}) == dist0)
            return e;
        else if(std::max({dist0, dist1, dist2}) == dist1)
            return mesh_input->next(e);
        else
            return mesh_input->prev(e);
        return -1;
    }

 
    //Return true if the edge e is the lowest edge both triangles incident to e
    //in case of border edges, they are always labeled as frontier-edge
    bool is_frontier_edge(const int e)
    {
        int twin = mesh_input->twin(e);
        bool is_border_edge = mesh_input->is_border_edge(e) || mesh_input->is_border_edge(twin);
        bool is_not_max_edge = !(max_edges[e] || max_edges[twin]);
        if(is_border_edge || is_not_max_edge)
            return true;
        else
            return false;
    }

    //Travel in CCW order around the edges of vertex v from the edge e looking for the next frontier edge
    int search_frontier_edge(const int e)
    {
        int nxt = e;
        int orie = mesh_input->origin(e);
        int deste = mesh_input->target(e);
        while(!frontier_edges[nxt]){
            nxt = mesh_input->CW_edge_to_vertex(nxt);
            int ori = mesh_input->origin(nxt);
            int dest = mesh_input->target(nxt);
        }
        return nxt;
    }

    //return true if the polygon is not simple
    bool has_BarrierEdgeTip(int e_init){

        int e_curr = mesh_output->next(e_init);
        //travel inside frontier-edges of polygon
        while(e_curr != e_init){   
            //if the twin of the next halfedge is the current halfedge, then the polygon is not simple
            if( mesh_output->twin(mesh_output->next(e_curr)) == e_curr)
                return true;
            //travel to next half-edge
            e_curr = mesh_output->next(e_curr);
        }
        return false;
    }   

    //generate a polygon from a seed edge
    //input: Seed-edge
    //Output: seed frontier-edge of new popygon
    int travel_triangles(const int e)
    {   
        int ori = mesh_input->origin(e);
        int dest = mesh_input->target(e);
        //search next frontier-edge
        int e_init = search_frontier_edge(e);
        int ori_e_init = mesh_input->origin(e_init);
        int dest_e_init = mesh_input->target(e_init);

        int e_curr = mesh_input->next(e_init);
        int ori_e_curr = mesh_input->origin(e_curr);
        int dest_e_curr = mesh_input->target(e_curr);

        int e_fe = e_init;
        
        //travel inside frontier-edges of polygon
        do{   
            e_curr = search_frontier_edge(e_curr);
            //update next of previous frontier-edge
            mesh_output->set_next(e_fe, e_curr);  
            //update prev of current frontier-edge
            mesh_output->set_prev(e_curr, e_fe);

            //travel to next half-edge
            e_fe = e_curr;
            e_curr = mesh_input->next(e_curr);
        }while(e_fe != e_init);
        return e_init;
    }
    
    //Given a barrier-edge tip v, return the middle edge incident to v
    //The function first calculate the degree of v - 1 and then divide it by 2, after travel to until the middle-edge
    //input: vertex v
    //output: edge incident to v
    int calculate_middle_edge(const int v){
        int frontieredge_with_bet = this->search_frontier_edge(mesh_input->edge_of_vertex(v));
        int internal_edges =mesh_input->degree(v) - 1; //internal-edges incident to v
        int adv = (internal_edges%2 == 0) ? internal_edges/2 - 1 : internal_edges/2 ;
        int nxt = mesh_input->CW_edge_to_vertex(frontieredge_with_bet);
        //back to traversing the edges of v_bet until select the middle-edge
        while (adv != 0){
            nxt = mesh_input->CW_edge_to_vertex(nxt);
            adv--;
        }
        return nxt;
    }


    // Return the best edges to insert in the angle formed by origin(e)
    // ccw indicates if the polygon is in counterclockwise order.
    std::vector<int> get_best_edges_reflex(int e, bool ccw) {

        int e_first = mesh_output->twin(mesh_output->prev(e));
        int e_prev = e_first;
        int e_curr = e_first;
        double angle = 0;
        std::vector<int> best_edges;
        std::vector<int> tmp_edges;

        // Iterate over the edges finding the already inserted edges
        while (true) {

            // Depending on the orientation of the polygon, the edges have to be iterated in clockwise or counter-clockwise order
            e_curr = ccw ? mesh_input->CW_edge_to_vertex(e_prev, ccw) : mesh_input->CCW_edge_to_vertex(e_prev, ccw);

            // The angle subtended by the edges visited so far
            angle += mesh_input->angle(e_prev, e_curr);
            
            // If the edge has already been inserted
            if (frontier_edges[e_curr]) {
                if (angle > M_PI) {
                    tmp_edges = get_best_edges_between(e_first, e_curr, angle, ccw);
                    best_edges.insert(best_edges.end(), tmp_edges.begin(), tmp_edges.end());
                }
                e_first = e_curr;
                angle = 0;
            }

            e_prev = e_curr;

            if (e_curr == e) {
                return best_edges;
            }
        }
    }


    // Returns the best edges to insert between edges e1 and e2.
    // ccw indicates if the angle belongs to a counterclockwise polygon or not.
    //input: edge e1, edge e2, angle formed by e1 and e2, orientation of polygon ccw
    //output: vector with the best edges to insert
    std::vector<int> get_best_edges_between(int e1, int e2, double angle, bool ccw) {

        int best_edge_half = e1;
        int best_edge_first_third = e1;
        int best_edge_second_third = e1;
        double best_angle = 0;
        double total_angle = 0;
        int e_curr = e1;
        int next_edge = ccw ? mesh_input->CW_edge_to_vertex(e_curr, ccw) : mesh_input->CCW_edge_to_vertex(e_curr, ccw);

        // Iterates over the edges between e1 and e2
        do {

            double curr_angle = mesh_input->angle(e_curr, next_edge);
            double new_total_angle = total_angle + curr_angle;
            
            // If the edge repairs two reflex angles (its two ends), is the best edge to insert (following the
            // minimum number of edges added criterion).
            if (new_total_angle < M_PI && (angle - new_total_angle) < M_PI && fixes_angle(next_edge)) {
                return std::vector<int>({next_edge});
            }
            
            // The edge splits the angle evenly-ish
            if (std::abs(angle/2 - new_total_angle) < std::abs(angle/2 - total_angle)) {
                best_edge_half = next_edge;
                best_angle = new_total_angle;
            }

            // The edge splits the angle in three parts (or so)
            if (std::abs(angle/3 - new_total_angle) < std::abs(angle/3 - total_angle)) {
                best_edge_first_third = next_edge;
            }
            if (std::abs(2*angle/3 - new_total_angle) < std::abs(2*angle/3 - total_angle)) {
                best_edge_second_third = next_edge;
            }
            total_angle = new_total_angle;
            e_curr = next_edge;
            next_edge = ccw ? mesh_input->CW_edge_to_vertex(e_curr, ccw) : mesh_input->CCW_edge_to_vertex(e_curr, ccw);
        } while (next_edge != e2);

        // If only needs one edge to insert
        if (best_angle < M_PI && (angle - best_angle) < M_PI) {
            return std::vector({best_edge_half});
        }

        return std::vector({best_edge_first_third, best_edge_second_third});
    }

    // Determines if the target of edge e repairs a reflex angle
    bool fixes_angle(int e) {
        int e_ccw = mesh_input->twin(e);
        int e_cw = e_ccw;

        double angle_ccw = 0;
        double angle_cw = 0;

        // Go through the edges in counterclockwise order until the next frontier edge is found, computing the total angle
        do {
            angle_ccw += mesh_input->angle(e_ccw, mesh_input->CCW_edge_to_vertex(e_ccw));
            e_ccw = mesh_input->CCW_edge_to_vertex(e_ccw);
        } while (!frontier_edges.at(e_ccw));

        if (angle_ccw > M_PI) {
            return false;
        }

        // Same as before, but in clockwise order
        do {
            angle_cw += mesh_input->angle(e_cw, mesh_input->CW_edge_to_vertex(e_cw));
            e_cw = mesh_input->CW_edge_to_vertex(e_cw);
        } while (!frontier_edges.at(e_cw));

        if (angle_cw > M_PI) {
            return false;
        }

        return angle_ccw + angle_cw > M_PI;
    }


    //Given a seed edge e that generated polygon, split the polygon until remove al barrier-edge tips
    //input: seed edge e, polygon poly
    //output: polygon without barrier-edge tips
    void barrieredge_tip_reparation(const int e)
    {
        this->n_polygons_to_repair++;
        int x, y, i;
        int t1, t2;
        int middle_edge, v_bet;

        int e_init = e;
        int e_curr = mesh_output->next(e_init);
        //search by barrier-edge tips
        while(e_curr != e_init){   
            //if the twin of the next halfedge is the current halfedge, then the polygon is not simple
            if( mesh_output->twin(mesh_output->next(e_curr)) == e_curr){
                //std::cout<<"e_curr "<<e_curr<<" e_next "<<mesh_output->next(e_curr)<<" next del next "<<mesh_output->next(mesh_output->next(e_curr))<<" twin curr "<<mesh_output->twin(e_curr)<<" twin next "<<mesh_output->twin(mesh_output->next(e_curr))<<std::endl;

                n_barrier_edge_tips++;
                n_frontier_edges+=2;

                //select edge with bet
                v_bet = mesh_output->target(e_curr);
                middle_edge = calculate_middle_edge(v_bet);

                //middle edge that contains v_bet
                t1 = middle_edge;
                t2 = mesh_output->twin(middle_edge);
                
                //edges of middle-edge are labeled as frontier-edge
                this->frontier_edges[t1] = true;
                this->frontier_edges[t2] = true;

                //edges are use as seed edges and saves in a list
                triangle_list.push_back(t1);
                triangle_list.push_back(t2);

                seed_bet_mark[t1] = true;
                seed_bet_mark[t2] = true;
            }
                
            //travel to next half-edge
            e_curr = mesh_output->next(e_curr);
        }

        int t_curr;
        //generate polygons from seeds,
        //two seeds can generate the same polygon
        //so the bit_vector seed_bet_mark is used to label as false the edges that are already used
        int new_polygon_seed;
        while (!triangle_list.empty()){
            t_curr = triangle_list.back();
            triangle_list.pop_back();
            if(seed_bet_mark[t_curr]){
                this->n_polygons_added_after_repair++;
                seed_bet_mark[t_curr] = false;
                new_polygon_seed = generate_repaired_polygon(t_curr, seed_bet_mark);
                //Store the polygon in the as part of the mesh
                output_seeds.push_back(new_polygon_seed);
            }
        }

    }

    // Returns if the origin of edge e is a reflex angle
    // input: edge e
    bool is_reflex_angle(int e, bool counterclockwise_polygon) {
        int v1 = mesh_output->origin(mesh_output->prev(e));
        int v2 = mesh_output->origin(e);
        int v3 = mesh_output->target(e);
        int ori = mesh_output->orientation({
            mesh_output->get_PointX(v1), mesh_output->get_PointY(v1)},
            {mesh_output->get_PointX(v2), mesh_output->get_PointY(v2)},
            {mesh_output->get_PointX(v3), mesh_output->get_PointY(v3)});
        if (ori == 0) {
            return false;
        }
        bool is_counter = (ori == 1);
        return is_counter != counterclockwise_polygon;
    }

    // Adds edges until the polygon is convex.
    // very similar to the function barrieredge_tip_reparation
    void reflex_reparation(const int e) {
        bool counterclockwise = is_polygon_counterclockwise(e);
        int t1, t2;

        std::vector<int> best_edges;
        int n_convex_polygons_added = 0;
        int e_curr = e;
        int n_reflex_angles = 0;
        int n_edges_added = 0;
        do {
            double x = mesh_output->get_PointX(mesh_output->origin(e_curr));
            double y = mesh_output->get_PointY(mesh_output->origin(e_curr));
            double xn = mesh_output->get_PointX(mesh_output->origin(mesh_output->next(e_curr)));
            double yn = mesh_output->get_PointY(mesh_output->origin(mesh_output->next(e_curr)));
            bool is_reflex = is_reflex_angle(e_curr, counterclockwise);
            if(is_reflex){

                n_reflex_angles++;

                best_edges = get_best_edges_reflex(e_curr, counterclockwise);

                n_edges_added += best_edges.size();

                n_frontier_edges += best_edges.size();
                
                for (int edge: best_edges) {
                    int t1 = edge, t2 = mesh_input->twin(edge);
                    this->frontier_edges[t1] = true;
                    this->frontier_edges[t2] = true;
                    triangle_list.push_back(t1);
                    triangle_list.push_back(t2);
                    seed_bet_mark[t1] = true;
                    seed_bet_mark[t2] = true;
                }
            }
                
            //travel to next half-edge
            e_curr = mesh_output->next(e_curr);
        } while(e_curr != e);

        int t_curr;
        //generate polygons from seeds,
        //two seeds can generate the same polygon
        //so the bit_vector seed_bet_mark is used to label as false the edges that are already used
        int new_polygon_seed;
        while (!triangle_list.empty())
        {
            t_curr = triangle_list.back();
            triangle_list.pop_back();
            if(seed_bet_mark[t_curr]){
                seed_bet_mark[t_curr] = false;
                new_polygon_seed = generate_repaired_polygon(t_curr, seed_bet_mark);
                //Store the polygon in the as part of the mesh
                output_seeds.push_back(new_polygon_seed);
            }
        }
        n_polygons_added_after_repair_reflex += n_edges_added;

    }

    /*
    //Generate a polygon from a seed-edge and remove repeated seed from seed_list
    //POSIBLE BUG: el algoritmo no viaja por todos los halfedges dentro de un poligono, 
    //por lo que pueden haber semillas que no se borren y tener poligonos repetidos de output
    int generate_repaired_polygon(const int e, bit_vector &seed_list)
    {   
        int e_init = e;
        //search next frontier-edge
        while(!frontier_edges[e_init]){
            e_init = mesh_input->CW_edge_to_vertex(e_init);
            seed_list[e_init] = false; 
            //seed_list[mesh_input->twin(e_init)] = false;
        }        
        //first frontier-edge is store to calculate the prev of next frontier-edfge
        int e_prev = e_init; 
        int v_init = mesh_input->origin(e_init);

        int e_curr = mesh_input->next(e_init);
        int v_curr = mesh_input->origin(e_curr);
        seed_list[e_curr] = false;

        //travel inside frontier-edges of polygon
        while(e_curr != e_init && v_curr != v_init){   
            while(!frontier_edges[e_curr])
            {
                e_curr = mesh_input->CW_edge_to_vertex(e_curr);
                seed_list[e_curr] = false;
          //      seed_list[mesh_input->twin(e_curr)] = false;
            } 

            //update next of previous frontier-edge
            mesh_output->set_next(e_prev, e_curr);  
            //update prev of current frontier-edge
            mesh_output->set_prev(e_curr, e_prev);

            //travel to next half-edge
            e_prev = e_curr;        
            e_curr = mesh_input->next(e_curr);
            v_curr = mesh_input->origin(e_curr);
            seed_list[e_curr] = false;
            //seed_list[mesh_input->twin(e_curr)] = false;
        }
        mesh_output->set_next(e_prev, e_init);
        mesh_output->set_prev(e_init, e_prev);
        return e_init;
    }
*/

    //Generate a polygon from a seed-edge and remove repeated seed from seed_list
    //POSIBLE BUG: el algoritmo no viaja por todos los halfedges dentro de un poligono, 
    //por lo que pueden haber semillas que no se borren y tener poligonos repetidos de output
    int generate_repaired_polygon(const int e, bit_vector &seed_list)
    {   
        int e_init = e;

        //search next frontier-edge
        while(!frontier_edges[e_init]){
            e_init = mesh_input->CW_edge_to_vertex(e_init);
            seed_list[e_init] = false; 
            //seed_list[mesh_input->twin(e_init)] = false;
        }   
        int e_curr = mesh_input->next(e_init);    
        seed_list[e_curr] = false;
    
        int e_fe = e_init; 

        //travel inside frontier-edges of polygon
        do{   
            while(!frontier_edges[e_curr])
            {
                e_curr = mesh_input->CW_edge_to_vertex(e_curr);
                seed_list[e_curr] = false;
          //      seed_list[mesh_input->twin(e_curr)] = false;
            } 
            //update next of previous frontier-edge
            mesh_output->set_next(e_fe, e_curr);  
            //update prev of current frontier-edge
            mesh_output->set_prev(e_curr, e_fe);

            //travel to next half-edge
            e_fe = e_curr;
            e_curr = mesh_input->next(e_curr);
            seed_list[e_curr] = false;

        }while(e_fe != e_init);
        return e_init;
    }
    
};

#endif