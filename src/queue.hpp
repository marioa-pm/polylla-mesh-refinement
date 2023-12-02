#ifndef QUEUE_HPP
#define QUEUE_HPP

#include <vector>
#include <map>
#include <any>
#include "triangulation.hpp"
#include "region.hpp"

// Represents a queue to store triangles with their respective metric (the metric is the criterion to order it)
class BaseQueue {
        
        protected:

        typedef std::multimap<double, int> metricMultimap;
        typedef std::unordered_map<int, std::multimap<double, int>::iterator> triangleMap;
        typedef std::pair<double, int> metricTrianglePair;
        
        // Dictionary with the metric as key and the triangle as value (multimap is an ordered container,
        // so is ordered by the metric)
        metricMultimap metrics;

        // Dictionary with the triangles as keys and its respective position in "metrics" as values
        triangleMap triangles;

        int size = 0;
        Region &region;

        public:

        BaseQueue(Region &region) : region(region) {}
        
        ~BaseQueue() {
            metrics.clear();
            triangles.clear();
        }
        
        // Returns the triangle with higher priority (if the queue is empty returns {-1, -1})
        metricTrianglePair get_first() {
            if (!metrics.empty()) {
                // Depending on the criterion, the triangle with the higher priority could be
                // the one with the higher or lower metric 
                if (region.criterion.max) {
                    metricMultimap::reverse_iterator it = metrics.rbegin();
                    metricTrianglePair pair = *it;
                    return pair;
                }
                else {
                    metricMultimap::iterator it = metrics.begin();
                    metricTrianglePair pair = *it;
                    return pair;
                }
            }
            return {-1, -1};
        }

        // Returns a pair with the mertric and its triangle
        metricTrianglePair get(int triangle) {
            metricMultimap::iterator it_metrics = (*triangles.find(triangle)).second;
            return *it_metrics;
        }

        metricMultimap iterator() {
            return metrics;
        }

        // Inserts a triangle and its already calculated metric in the queue
        virtual void insert(int triangle, double metric) {
            metricMultimap::iterator inserted;

            // Checks if the metric doesn't satisfie the criterion, otherwise is not inserted
            if (region.criterion.has_to_be_refined(triangle, metric)) {
                inserted = metrics.insert(metricTrianglePair(metric, triangle));
                triangles.insert(std::pair<int, metricMultimap::iterator>(triangle, inserted));
                size++;
            }
        }

        // Inserts a triangle in the queue
        virtual void insert(int triangle) {
            metricMultimap::iterator inserted;
            // Checks if the metric doesn't satisfie the criterion, otherwise is not inserted
            std::pair<bool, double> criteria_data = region.criterion.has_to_be_refined(triangle);
            if (criteria_data.first && region.in_region(triangle)) {
                inserted = metrics.insert(metricTrianglePair(criteria_data.second, triangle));
                triangles.insert(std::pair<int, metricMultimap::iterator>(triangle, inserted));
                size++;
            }
        }

        // Removes a triangle.
        // Returns the removed pair, if the triangle is not in the queue, returns {0, -1}
        virtual metricTrianglePair remove(int triangle) {
            triangleMap::iterator it_triangles = triangles.find(triangle);
            if (it_triangles == triangles.end()) {
                return {0, -1};
            }
            metricMultimap::iterator it_metrics = (*it_triangles).second;
            metricTrianglePair removed = *it_metrics;
            triangles.erase(it_triangles);
            metrics.erase(it_metrics);
            size--;
            return removed;
        }

        // Updates the triangles in the queue.
        // triangle_pairs: vector of pairs, with each pair having as its first element the old triangle
        // to be updated and as second element the new triangle (the metric of a triangle can be updated
        // by just providing a pair with the same triangle)
        void update(std::vector<std::pair<int, int>> triangle_pairs) {
            for (std::pair<int, int> triangle_pair: triangle_pairs) {
                int old_triangle = triangle_pair.first;
                int new_triangle = triangle_pair.second;

                metricTrianglePair removed = remove(old_triangle);
                insert(new_triangle);
            }
        }

        // Checks if the queue is empty (which means there is no triangle left to refine)
        virtual bool is_refined() {
            return triangles.size() == 0;
        }
    };

    
    // Represents a queue whose criterion is the average value of the metric
    class QueueAverage: public BaseQueue {
        
        protected:

        double original_average_area = 0;
        double current_average_area = 0;
        double factor = 0;

        public:

        // In this case, all the triangles of the region are inserted immediately, as
        // all of them could potentially be refined (this can be done in a better way)
        QueueAverage(Region &region, double factor) : BaseQueue(region), factor(factor) {
            std::pair<int, double> triangle_area = region.get_pair();

            // While the region still contains a triangle
            while (triangle_area.first != -1) {
                insert(triangle_area.first, triangle_area.second);
                original_average_area += triangle_area.second;
                triangle_area = region.get_pair();
            }
            if (triangles.size() == 0) {
                throw std::runtime_error("The region defined doesn't contain any triangle.");
            }
            original_average_area /= triangles.size();
            current_average_area = original_average_area;
        }

        // Inserts a triangle and its already calculated metric in the queue
        void insert(int triangle, double area) {
            metricMultimap::iterator inserted;
            inserted = metrics.insert(metricTrianglePair(area, triangle));
            triangles.insert(std::pair<int, metricMultimap::iterator>(triangle, inserted));
            current_average_area *= (triangles.size() - 1) / (double) triangles.size();
            size++;
        }

        // Inserts a triangle in the queue
        void insert(int triangle) {
            metricMultimap::iterator inserted;
            std::pair<bool, double> criteria_data = region.criterion.has_to_be_refined(triangle);
            if (criteria_data.first && region.in_region(triangle)) {
                inserted = metrics.insert(metricTrianglePair(criteria_data.second, triangle));
                triangles.insert(std::pair<int, metricMultimap::iterator>(triangle, inserted));
                current_average_area *= (triangles.size() - 1) / (double) triangles.size();
                size++;
            }
        }

        // Removes a triangle.
        // Returns the removed pair, if the triangle is not in the queue, returns {0, -1}
        metricTrianglePair remove(int triangle) {
            triangleMap::iterator it_triangles = triangles.find(triangle);
            if (it_triangles == triangles.end()) {
                return {0, -1};
            }
            metricMultimap::iterator it_metrics = (*it_triangles).second;
            metricTrianglePair removed = *it_metrics;
            triangles.erase(it_triangles);
            metrics.erase(it_metrics);
            current_average_area *= (triangles.size() + 1) / (double) triangles.size();
            size--;
            return removed;
        }

        // Indicates if the average value satisfies the factor imposed, in which case the refinement ends
        bool is_refined() {
            return current_average_area / original_average_area <= factor;
        }


        double get_average_area() {
            return current_average_area;
        }
    };


    // Represents a queue whose criterion is a threshold value for a given metric for every triangle
    class QueueThreshold: public BaseQueue {
        
        public:

        QueueThreshold(Region &region) : BaseQueue(region) {
            std::pair<int, double> triangle_metric = region.get_pair();
            while (triangle_metric.first != -1) {
                insert(triangle_metric.first, triangle_metric.second);
                triangle_metric = region.get_pair();
            }
            if (triangles.size() == 0) {
                throw std::runtime_error("The region defined doesn't contain any triangle.");
            }
        }
    };

#endif