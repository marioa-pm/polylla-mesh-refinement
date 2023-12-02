#ifndef CRITERION_HPP
#define CRITERION_HPP

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

/*
    Represents the criteria to determine if a triangle is bad, according to the
    given metric and the given threshold.
    max: indicates if the threshold is a maximum value for the metric or not (or is the minimum value).
    */
    class Criterion {

        public:

        Triangulation &triangulation;
        bool max;
        double threshold;

        public:

        Criterion(Triangulation &triangulation, bool max, double threshold) : triangulation(triangulation), max(max), threshold(threshold) {}
        Criterion(Triangulation &triangulation) : triangulation(triangulation), max(true), threshold(0) {}

        virtual double get_metric(int triangle) {
            return 0;
        }

        virtual std::pair<bool, double> has_to_be_refined(int triangle) {
            return {true, get_metric(triangle)};
        }

        virtual bool has_to_be_refined(int triangle, double metric) {
            return true;
        }

    };


    // Represents the refinement criterion of maximum area.
    class MaxAreaCriterion: public Criterion {

        public:

        MaxAreaCriterion(Triangulation &triangulation, double threshold) : Criterion(triangulation, true, threshold) {}

        double get_metric(int triangle) {
            return triangulation.get_area(triangle);
        }
        

        // Returns a pair indicating if the triangle has to be refined and its area
        std::pair<bool, double> has_to_be_refined(int triangle) {
            double metric = get_metric(triangle);
            if (max) {
                return {metric > threshold, metric};
            }
            return {metric < threshold, metric};
        }


        // Returns if the triangle has to be refined, according to the given metric.
        bool has_to_be_refined(int triangle, double metric) {
            if (max) {
                return metric > threshold;
            }
            return metric < threshold;
        }
    };


    // Represents the refinement criterion of average area
    class AverageAreaCriterion: public Criterion {

        public:

        AverageAreaCriterion(Triangulation &triangulation) : Criterion(triangulation) {}

        double get_metric(int triangle) {
            return triangulation.get_area(triangle);
        }
    };

#endif