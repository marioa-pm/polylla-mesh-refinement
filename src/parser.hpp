#ifndef PARSER_HPP
#define PARSER_HPP

#include <vector>
#include <map>
#include <any>


class Parser {
    
    std::vector<std::string_view> args;
    std::map<std::string, std::any> optArgs;

    public:

    Parser(int argc, char **argv) {

        if (argc == 1) {
            std::stringstream msg;
            msg << "Usage: " << argv[0] << " <off_file .off> <output name>" << std::endl;
            msg << "Usage: " << argv[0] << " <node_file .node> <ele_file .ele> <neigh_file .neigh> <output name>";
            throw std::runtime_error(msg.str());
        }

        std::vector<std::string_view> argsc(argv + 1, argv + argc);
        auto it = argsc.begin();
        auto end = argsc.end();
        while (it < end) {
            if ((*it)[0] != '-') {
                args.push_back(*it);
            }
            else {
                if (args.size() != 2 && args.size() != 4) {
                    std::stringstream msg;
                    msg << "Incorrect number of arguments." << std::endl;
                    msg << "Usage: " << argv[0] << " <off file .off> <output name>" << std::endl;
                    msg << "Usage: " << argv[0] << " <node_file .node> <ele_file .ele> <neigh_file .neigh> <output name>";
                    throw std::runtime_error(msg.str());
                }
                break;
            }
            it++;
        }

        optArgs.insert({"-output", std::string(args.back())});
        args.pop_back();
        
        while (it < end) {
            if ((*it) == "-factor") {
                if (optArgs.count("-max_area") > 0) {
                    throw std::runtime_error("-max_area option already specified.");
                }
                if (optArgs.count("-factor") > 0) {
                    throw std::runtime_error("-factor option already specified.");
                }
                if (it + 1 >= end) {
                    throw std::runtime_error("Missing value for -factor option (float value).");
                }
                std::string_view v = *(++it);
                double factor = parseToDouble(v);
                if (factor > 1 || factor <= 0) {
                    throw std::runtime_error("Factor value must be greater than zero and less than one.");
                }
                optArgs.insert({"-factor", factor});
                ++it;
            }
            else if ((*it) == "-max_area") {
                if (optArgs.count("-factor") > 0) {
                    throw std::runtime_error("-factor option already specified.");
                }
                if (it + 1 >= end) {
                    throw std::runtime_error("Missing value for -max_area option (float value).");
                }
                std::string_view v = *(++it);
                double max_value = parseToDouble(v);
                if (max_value <= 0) {
                    throw std::runtime_error("Max value for area must be greater than zero.");
                }
                optArgs.insert({"-max_area", max_value});
                ++it;
            }
            else if ((*it) == "-not_delaunay") {
                optArgs.insert({"-not_delaunay", true});
                ++it;
            }
            else if ((*it) == "-check_delaunay") {
                optArgs.insert({"-check_delaunay", true});
                ++it;
            }
            else if ((*it) == "-centroid") {
                optArgs.insert({"-centroid", true});
                ++it;
            }
            else if ((*it) == "-polylla") {
                optArgs.insert({"-polylla", true});
                ++it;
            }
            else if ((*it) == "-dont_write") {
                optArgs.insert({"-dont_write", true});
                ++it;
            }
            else if ((*it) == "-print_region") {
                optArgs.insert({"-print_region", true});
                ++it;
            }
            else if ((*it) == "-circle") {
                if (optArgs.count("-region") > 0) {
                    throw std::runtime_error("Region already specified.");
                }
                if (it + 3 >= end) {
                    throw std::runtime_error("Must specify the center and radius with circle option. (-circle <center_x> <center_y> <radius>)");
                }
                double cx = parseToDouble(*(++it));
                double cy = parseToDouble(*(++it));
                double r = parseToDouble(*(++it));
                optArgs.insert({"-circle", std::vector<double>({cx, cy, r})});
                optArgs.insert({"-region", std::string("-circle")});
                ++it;
            }
            else if ((*it) == "-rectangle") {
                if (optArgs.count("-region") > 0) {
                    throw std::runtime_error("region already specified.");
                }
                if (it + 4 >= end) {
                    throw std::runtime_error("Must specify the center, width and height with rectangle option. (-rectangle <center_x> <center_y> <width> <height>)");
                }
                double cx = parseToDouble(*(++it));
                double cy = parseToDouble(*(++it));
                double width = parseToDouble(*(++it));
                double height = parseToDouble(*(++it));
                optArgs.insert({"-rectangle", std::vector<double>({cx, cy, width, height})});
                optArgs.insert({"-region", std::string("-rectangle")});
                ++it;
            }
            else if ((*it) == "-polygon") {
                if (optArgs.count("-region") > 0) {
                    throw std::runtime_error("region already specified.");
                }
                if (it + 1 >= end) {
                    throw std::runtime_error("Must specify the file in .off format with the polygon (-polygon <filename>).");
                }
                optArgs.insert({"-polygon", std::string((*(++it)))});
                optArgs.insert({"-region", std::string("-polygon")});
                ++it;
            }
            else if ((*it) == "-points_range") {
                if (optArgs.count("-region") > 0) {
                    throw std::runtime_error("region already specified.");
                }
                if (it + 2 >= end) {
                    throw std::runtime_error("Must specify a range of vertices. (-points_range <start> <end>)");
                }
                int start = parseToInt(*(++it));
                int end = parseToInt(*(++it));
                optArgs.insert({"-points_range", std::vector<int>({start, end})});
                optArgs.insert({"-region", std::string("-points_range")});
                ++it;
            }
            else if ((*it) == "-points_file") {
                if (optArgs.count("-region") > 0) {
                    throw std::runtime_error("region already specified.");
                }
                if (it + 1 >= end) {
                    throw std::runtime_error("Must specify a file indicating the vertices. (-points_file <filename>)");
                }
                optArgs.insert({"-points_file", std::string((*(++it)))});
                optArgs.insert({"-region", std::string("-points_file")});
                ++it;
            }
            else if ((*it) == "-print_iter") {
                if (it + 1 >= end) {
                    throw std::runtime_error("Must specify an already existing directory to save the iterations.");
                }
                optArgs.insert({"-print_iter", std::string((*(++it)))});
                ++it;
            }
            else if ((*it) == "-filename") {
                if (it + 1 >= end) {
                    throw std::runtime_error("Must specify a filename for the final mesh.");
                }
                optArgs.insert({"-filename", std::string((*(++it)))});
                ++it;
            }
            else {
                std::stringstream msg;
                msg << "Unknown value or option: " << *it; 
                throw std::runtime_error(msg.str());
            }
        }
        /*
        if (optArgs.count("-region") > 0 && optArgs.count("-factor") == 0) {
            throw std::runtime_error("-factor option must be specified if a region is given.");
        }
        */
    }


    std::vector<std::string_view> getArgs() {
        return args;
    }


    bool hasOption(std::string option) {
        auto it = optArgs.find(option);
        return it != optArgs.end();
    }


    std::any getOption(std::string option) {
        auto it = optArgs.find(option);
        if (it != optArgs.end()) {
            return it->second;
        }
        return false;
    }


    int getNumberFiles() {
        return args.size();
    }


    private:
    
    int parseToInt(std::string_view &str) {
        char *end;
        int x = std::strtol(str.data(), &end, 10);
        if (end != str.data() + str.size()) {
            std::stringstream msg;
            msg << str << " is not a valid integer value.";
            throw std::runtime_error(msg.str());
        }
        return x;
    }


    float parseToFloat(std::string_view &str) {
        char *end;
        float x = std::strtof(str.data(), &end);
        if (end != str.data() + str.size()) {
            std::stringstream msg;
            msg << str << " is not a valid float value.";
            throw std::runtime_error(msg.str());
        }
        return x;
    }


    double parseToDouble(std::string_view &str) {
        char *end;
        double x = std::strtod(str.data(), &end);
        if (end != str.data() + str.size()) {
            std::stringstream msg;
            msg << str << " is not a valid double value.";
            throw std::runtime_error(msg.str());
        }
        return x;
    }
};


#endif