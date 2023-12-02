#ifndef REFINER_HPP
#define REFINER_HPP

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
#include "queue.hpp"

class Refiner {

    public:

    Triangulation &tr;
    BaseQueue &queue;
    bool delaunay = true;
    bool bisection = true;
    std::string directory_refinement = "";
    bool print_iter = false;
    int number_lepp_bisection = 0;


    Refiner(Triangulation &triangulation, BaseQueue &queue, bool delaunay, bool bisection) : tr(triangulation), queue(queue), delaunay(delaunay), bisection(bisection) {}
    Refiner(Triangulation &triangulation, BaseQueue &queue) : tr(triangulation), queue(queue), delaunay(true), bisection(true) {}

    
    void refine_area() {
        while (!queue.is_refined()) {
            /*
            double current_factor = queue.get_average_area()/originalAverageArea;
            if (last_factor - current_factor > 0.00001) {
                std::cout << "factor: " << current_factor << std::endl;
                last_factor = current_factor;
            }
            */
            leppBisectionAuxRefinement(queue.get_first().second);
            
        }
    }


    // Set the directory where the iterations must be saved
    void setDirectory(std::string directory) {
        directory_refinement = directory;
        print_iter = true;
    }


    // Performs the lepp refinement algorithm over the triangle that has e as edge
    void leppBisectionAuxRefinement(int e) {
        std::vector<int> leppTriangleE;

        leppTriangleE = tr.lepp(e);
        
        // The face of the triangle to be refined (could be the same as e).
        int targetTriangle = leppTriangleE.front();

        // Get terminal-triangle
        int toRefine = leppTriangleE.back();

        // Vector to hold the triangles that should be updated in the queue
        std::vector<int> trianglesToUpdate;

        // List to hold the edges to legalize.
        std::list<int> edgesToFlip;

        while (true) {
            
            // Inserts the Terminal-triangle to be updated
            trianglesToUpdate.push_back(tr.face(toRefine));

            // If the edge is not border (meaning there are two terminal-triangles)
            if (!tr.is_border_edge(tr.twin(toRefine))) {
                trianglesToUpdate.push_back(tr.face(tr.twin(toRefine)));
            }

            // Get the position of the created triangles to be inserted in the queue
            std::vector<std::set<int>::iterator> triangle_inserted;

            // Determines the type of insertion to perform
            if (bisection || tr.is_border_edge(tr.twin(toRefine))) {
                triangle_inserted = insertVertexMiddleEdge(toRefine);
            }
            else {
                triangle_inserted = insertVertexCentroid(toRefine);
            }
            
            // Iterate over the modified triangles to update them in the queue
            for (int i = 0; i < trianglesToUpdate.size(); i++) {

                // Creates a vector of pairs
                std::pair<int, int> tup = {trianglesToUpdate.at(i), trianglesToUpdate.at(i)};
                std::vector<std::pair<int, int>> v = {tup};

                // Updates the pairs
                queue.update(v);
                queue.insert(*(triangle_inserted.at(i)));
            }
            
            trianglesToUpdate.clear();

            // If the legallization of edges destroys the triangle to refine, the algorithm stops
            bool stopBisection = false;

            if (delaunay) {
                // The edges that could no longer be legal due to the insertion
                edgesToFlip = {tr.twin(tr.next(toRefine)), tr.twin(tr.next(tr.twin(tr.prev(toRefine))))}; 
                if (!tr.is_border_edge(tr.twin(toRefine))) {
                    edgesToFlip.push_back(tr.twin(tr.prev(tr.twin(toRefine))));
                    edgesToFlip.push_back(tr.twin(tr.prev(tr.twin(tr.next(tr.twin(toRefine))))));
                }

                // The vertex inserted is the last in the list of vertices
                int inserted_vertex = tr.Vertices.size() - 1;
                stopBisection = legalize_edges(inserted_vertex, edgesToFlip, tr.face(targetTriangle));
                edgesToFlip.clear();
            }

            // Saves all iterations of the algorithm
            if (print_iter) {
                tr.print_OFF(directory_refinement + "step" + std::to_string(number_lepp_bisection) + ".off");
                number_lepp_bisection++;
            }

            // If the triangle to refine was refined or was destroyed by the legallization of edges.
            if ((toRefine == targetTriangle) || stopBisection) {
                break;
            }

            leppTriangleE = tr.lepp(targetTriangle);
            toRefine = leppTriangleE.back();
        }
    }


    bool legalize_edges(int inserted_vertex, std::list<int> edges, int target_triangle) {
        std::vector<std::pair<int, int>> triangles_to_update;
        std::list<int> current_edges;
        bool stop_lepp_refinement = false;

        // Loop through every edge that could be illegal
        for (int edge: edges) {

            // List to store the new vertices that need to be checked (the propagation part)
            current_edges = {edge};
            
            while (!current_edges.empty()) {
                int e = current_edges.back();
                current_edges.pop_back();
                if (tr.isIllegal(e)) {
                    
                    // If the target triangle to be refined is destroyed by the edge flips, the algorithm stops
                    if ((target_triangle == tr.face(edge)) || (target_triangle == tr.face(tr.twin(edge)))) {
                        stop_lepp_refinement = true;
                    }

                    // Triangles to update in the queue
                    triangles_to_update = edgeFlip(e);
                    queue.update(triangles_to_update);

                    // Depending on the inserted vertex, continues the propagation
                    if (tr.origin(e) == inserted_vertex) {
                        current_edges.push_back(tr.twin(tr.next(e)));
                        current_edges.push_back(tr.twin(tr.prev(tr.twin(e))));
                    }
                    else if (tr.target(e) == inserted_vertex) {
                        current_edges.push_back(tr.twin(tr.prev(e)));
                        current_edges.push_back(tr.twin(tr.next(tr.twin(e))));
                    }
                    else {
                        throw std::runtime_error("Ninguno de los extremos es el vértice insertado (por que?)");
                    }
                }
            }
        }
        return stop_lepp_refinement;
    }


    // Flips the edge e. e must be a diagonal of a cuadrilateral (cannot be a border edge)
    std::vector<std::pair<int, int>> edgeFlip(int e) {
        
        std::vector<std::pair<int, int>> trianglesModified;

        // Changes the incident halfedges that may no longer be valid
        tr.Vertices.at(tr.origin(e)).incident_halfedge = tr.next(tr.twin(e));
        tr.Vertices.at(tr.target(e)).incident_halfedge = tr.next(e);

        // The new halfedges
        halfEdge e1 = tr.HalfEdges.at(e);
        halfEdge e2 = tr.HalfEdges.at(tr.twin(e));

        // Changing the adjacencies
        e1.next = tr.next(tr.next(e));
        e1.prev = tr.next(tr.twin(e));
        e1.origin = tr.target(tr.next(tr.twin(e)));

        e2.next = tr.prev(tr.twin(e));
        e2.prev = tr.next(e);
        e2.origin = tr.target(tr.next(e));

        tr.HalfEdges.at(tr.prev(e)).next = tr.next(tr.twin(e));
        tr.HalfEdges.at(tr.prev(e)).prev = e;

        tr.HalfEdges.at(tr.next(e)).next = tr.twin(e);
        tr.HalfEdges.at(tr.next(e)).prev = tr.prev(tr.twin(e));

        // The case where one of the faces was determined by the edge to be flipped
        if (e2.face == tr.next(tr.twin(e))) {

            // Changes the face of the halfedges
            tr.HalfEdges.at(tr.prev(tr.twin(e))).face = tr.twin(e);
            tr.HalfEdges.at(tr.next(e)).face = tr.twin(e);
            std::set<int>::iterator position;

            // Finds the face in the list of faces
            position = tr.triangle_list.find(e2.face);

            // If not found (this should never happen)
            if (position == tr.triangle_list.end()) {
                std::cout << "The face was not found (This should never happen)." << std::endl;
            }
            else {
                // Delete the old face and insert the new one
                tr.triangle_list.erase(position);
                tr.triangle_list.insert(tr.twin(e));
            }
            trianglesModified.push_back({e2.face, tr.twin(e)});
            e2.face = tr.twin(e);
        }
        // No face is modified
        else {
            trianglesModified.push_back({e2.face, e2.face});
            tr.HalfEdges.at(tr.next(e)).face = e2.face;
        }

        tr.HalfEdges.at(tr.next(tr.twin(e))).next = e;
        tr.HalfEdges.at(tr.next(tr.twin(e))).prev = tr.prev(e);

        // Same as before, but with the other face
        if (e1.face == tr.next(e)) {
            tr.HalfEdges.at(tr.prev(e)).face = e;
            tr.HalfEdges.at(tr.next(tr.twin(e))).face = e;
            std::set<int>::iterator position;
            position = tr.triangle_list.find(e1.face);
            if (position == tr.triangle_list.end()) {
                std::cout << "The face was not found (why?)." << std::endl;
            }
            else {
                tr.triangle_list.erase(position);
                tr.triangle_list.insert(e);
            }
            trianglesModified.push_back({e1.face, e});
            e1.face = e;
        }
        else {
            trianglesModified.push_back({e1.face, e1.face});
            tr.HalfEdges.at(tr.next(tr.twin(e))).face = e1.face;
        }

        tr.HalfEdges.at(tr.prev(tr.twin(e))).next = tr.next(e);
        tr.HalfEdges.at(tr.prev(tr.twin(e))).prev = tr.twin(e);

        tr.HalfEdges.at(tr.twin(e)) = e2;
        tr.HalfEdges.at(e) = e1;

        return trianglesModified;
    }


    /*
    Insert a vertex in the middle of edge e.
    */
    std::vector<std::set<int>::iterator> insertVertexMiddleEdge(int e) {
        
        std::vector<std::set<int>::iterator> inserted;

        // The case where e is border edge is handled later, so we start with the twin
        if (tr.is_border_edge(e)) {
            e = tr.twin(e);
        }

        // Creates the new vertex to be inserted
        double x_mid = (tr.get_PointX(tr.origin(e)) + tr.get_PointX(tr.target(e))) / 2;
        double y_mid = (tr.get_PointY(tr.origin(e)) + tr.get_PointY(tr.target(e))) / 2;
        vertex v = {x_mid, y_mid, tr.is_border_edge(e) || tr.is_border_edge(tr.twin(e)), e};

        // Adds the new vertex to the list
        tr.Vertices.push_back(v);
        tr.n_vertices++;

        // Creates the new halfedges for the triangle to be split (most of these are taken as copies of the old
        // halfedges)
        halfEdge heMedianPrev;
        halfEdge heMedianNext = tr.HalfEdges.at(tr.prev(e));
        halfEdge heNext = tr.HalfEdges.at(e), hePrev = tr.HalfEdges.at(e), heOtherPrev = tr.HalfEdges.at(tr.prev(e));

        // Insert the median of one triangle, and stores the position of the new created face in the variable 'inserted'
        inserted.push_back(insertMedian(e, heNext, hePrev, heMedianNext, heMedianPrev, heOtherPrev));
        
        // If is not a border edge
        if (!tr.is_border_edge(tr.twin(e))) {
            
            // Same as before, inserting the median for the other triangle
            halfEdge heMedianPrev2;
            halfEdge heMedianNext2 = tr.HalfEdges.at(tr.prev(tr.twin(e)));
            halfEdge heNext2 = tr.HalfEdges.at(tr.twin(e)), hePrev2 = tr.HalfEdges.at(tr.twin(e)), heOtherPrev2 = tr.HalfEdges.at(tr.prev(tr.twin(e)));

            inserted.push_back(insertMedian(tr.twin(e), heNext2, hePrev2, heMedianNext2, heMedianPrev2, heOtherPrev2));

            // Updates the adjacencies with the new edges
            heNext.twin = tr.n_halfedges - 3;
            tr.HalfEdges.at(tr.n_halfedges - 6).twin = tr.twin(e);
            heNext2.twin = tr.n_halfedges - 6;
            tr.HalfEdges.at(tr.n_halfedges - 3).twin = e;

            // Replaces the old edges by the updated ones
            tr.HalfEdges.at(tr.prev(e)) = heMedianNext;
            tr.HalfEdges.at(tr.prev(tr.twin(e))) = heMedianNext2;
            tr.HalfEdges.at(tr.twin(e)) = heNext2;
            tr.HalfEdges.at(e) = heNext;

        }
        // Is a border edge
        else {

            // Creates just one new edge (split the border edge in two)
            halfEdge heBorderPrev = tr.HalfEdges.at(tr.twin(e));

            // Assign the adjacencies of the new edge
            heBorderPrev.next = tr.twin(e);
            heBorderPrev.twin = e;
            tr.HalfEdges.push_back(heBorderPrev);
            tr.n_halfedges++;

            heNext.twin = tr.n_halfedges - 1;
            tr.HalfEdges.at(tr.n_halfedges - 4).twin = tr.twin(e);
            
            // Changes the adjacencies of the splitted edge
            tr.HalfEdges.at(tr.twin(e)).prev = tr.n_halfedges - 1;
            tr.HalfEdges.at(tr.twin(e)).origin = tr.n_vertices - 1;
            tr.HalfEdges.at(tr.twin(e)).twin = tr.n_halfedges - 4;

            // Replaces the old edges by the updated ones
            tr.HalfEdges.at(tr.prev(e)) = heMedianNext;
            tr.HalfEdges.at(e) = heNext;
            tr.n_border_edges++;
        }
        return inserted;
    }

    /*
    Insert a vertex in the centroid of the quadrilateral that has e as diagonal.
    */
    std::vector<std::set<int>::iterator> insertVertexCentroid(int e) {
        
        std::vector<std::set<int>::iterator> inserted;

        // The case where e is border edge is handled later, so we start with the twin
        if (tr.is_border_edge(e)) {
            e = tr.twin(e);
        }

        // Creates the new vertex to be inserted
        double x_mid = (tr.get_PointX(tr.origin(e)) + tr.get_PointX(tr.target(e))) / 2;
        double y_mid = (tr.get_PointY(tr.origin(e)) + tr.get_PointY(tr.target(e))) / 2;
        point centroid = tr.get_centroid(e);
        vertex v = {centroid.first, centroid.second, tr.is_border_edge(e) || tr.is_border_edge(tr.twin(e)), e};

        // Adds the new vertex to the list
        tr.Vertices.push_back(v);
        tr.n_vertices++;

        // Creates the new halfedges for the triangle to be split (most of these are taken as copies of the old
        // halfedges)
        halfEdge heMedianPrev;
        halfEdge heMedianNext = tr.HalfEdges.at(tr.prev(e));
        halfEdge heNext = tr.HalfEdges.at(e), hePrev = tr.HalfEdges.at(e), heOtherPrev = tr.HalfEdges.at(tr.prev(e));

        // Insert the median of one triangle, and stores the position of the new created face in the variable 'inserted'
        inserted.push_back(insertMedian(e, heNext, hePrev, heMedianNext, heMedianPrev, heOtherPrev));
        
        // If is not a border edge
        if (!tr.is_border_edge(tr.twin(e))) {
            
            // Same as before, inserting the median for the other triangle
            halfEdge heMedianPrev2;
            halfEdge heMedianNext2 = tr.HalfEdges.at(tr.prev(tr.twin(e)));
            halfEdge heNext2 = tr.HalfEdges.at(tr.twin(e)), hePrev2 = tr.HalfEdges.at(tr.twin(e)), heOtherPrev2 = tr.HalfEdges.at(tr.prev(tr.twin(e)));

            inserted.push_back(insertMedian(tr.twin(e), heNext2, hePrev2, heMedianNext2, heMedianPrev2, heOtherPrev2));

            // Updates the adjacencies with the new edges
            heNext.twin = tr.n_halfedges - 3;
            tr.HalfEdges.at(tr.n_halfedges - 6).twin = tr.twin(e);
            heNext2.twin = tr.n_halfedges - 6;
            tr.HalfEdges.at(tr.n_halfedges - 3).twin = e;

            // Replaces the old edges by the updated ones
            tr.HalfEdges.at(tr.prev(e)) = heMedianNext;
            tr.HalfEdges.at(tr.prev(tr.twin(e))) = heMedianNext2;
            tr.HalfEdges.at(tr.twin(e)) = heNext2;
            tr.HalfEdges.at(e) = heNext;

        }
        // Is a border edge
        else {

            // Creates just one new edge (split the border edge in two)
            halfEdge heBorderPrev = tr.HalfEdges.at(tr.twin(e));

            // Assign the adjacencies of the new edge
            heBorderPrev.next = tr.twin(e);
            heBorderPrev.twin = e;
            tr.HalfEdges.push_back(heBorderPrev);
            tr.n_halfedges++;

            // Changes the incident halfedge that may no longer be valid
            tr.Vertices.at(tr.origin(tr.twin(e))).incident_halfedge = tr.n_halfedges - 1;

            // Changes the adjacencies of the splitted edge
            heNext.twin = tr.n_halfedges - 1;
            tr.HalfEdges.at(tr.n_halfedges - 4).twin = tr.twin(e);
            tr.HalfEdges.at(tr.twin(e)).prev = tr.n_halfedges - 1;
            tr.HalfEdges.at(tr.twin(e)).origin = tr.n_vertices - 1;
            tr.HalfEdges.at(tr.twin(e)).twin = tr.n_halfedges - 4;

            // Replaces the old edges by the updated ones
            tr.HalfEdges.at(tr.prev(e)) = heMedianNext;
            tr.HalfEdges.at(e) = heNext;

            tr.n_border_edges++;
        }
        return inserted;
    }

    // Inserts a new edge splitting the edge e in two edges (therefore, splitting the triangle that has e as an edge)
    // Returns the position of the new created triangle
    // heOtherPrev and heOtherNext must be copies of the originals edges, this is, copies of prev(e) and next(e), respectively.
    // The old triangle becomes the triangle that has as vertex the target of the original edge e to be divided (IMPORTANT).
    std::set<int>::iterator insertMedian(int e, halfEdge &heNext, halfEdge &hePrev, halfEdge &heMedianNext, halfEdge &heMedianPrev, halfEdge &heOtherPrev) {
            
            // Updates the number of faces and halfedges
            tr.n_halfedges += 3;
            tr.n_faces++;

            // Stores the position of the new inserted triangle
            std::pair<std::set<int>::iterator, bool> inserted = tr.triangle_list.insert(tr.n_halfedges - 1);

            // The indexes of the new vertex and halfedges (to be inserted in the structure later)
            int newVertexIndex = tr.n_vertices - 1;
            int hePrevIndex = tr.n_halfedges - 3;
            int heMedianPrevIndex = tr.n_halfedges - 2;
            int heOtherPrevIndex = tr.n_halfedges - 1;

            // Updates the incident halfedge of the vertex that may no longer be valid
            tr.Vertices.at(tr.origin(e)).incident_halfedge = hePrevIndex;
            
            // Assign the adjacencies of the new halfedges
            hePrev.origin = tr.origin(e);
            hePrev.next = heMedianPrevIndex;
            hePrev.prev = heOtherPrevIndex;
            hePrev.face = heOtherPrevIndex;;

            heMedianPrev.origin = newVertexIndex;
            heMedianPrev.twin = tr.prev(e);
            heMedianPrev.next = heOtherPrevIndex;
            heMedianPrev.prev = hePrevIndex;
            heMedianPrev.is_border = false;
            heMedianPrev.face = heOtherPrevIndex;

            heOtherPrev.next = hePrevIndex;
            heOtherPrev.prev = heMedianPrevIndex;
            heOtherPrev.face = heOtherPrevIndex;

            heNext.origin = newVertexIndex;
            heNext.next = tr.next(e);
            heNext.prev = tr.prev(e);

            heMedianNext.origin = tr.target(tr.next(e));
            heMedianNext.twin = heMedianPrevIndex;
            heMedianNext.next = e;
            heMedianNext.prev = tr.next(e);
            heMedianNext.is_border = false;

            tr.HalfEdges.at(heOtherPrev.twin).twin = heOtherPrevIndex;

            // Inserts the new haldfedges in the vectors
            tr.HalfEdges.push_back(hePrev);
            tr.HalfEdges.push_back(heMedianPrev);
            tr.HalfEdges.push_back(heOtherPrev);

            return inserted.first;
    }

};

#endif