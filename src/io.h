#ifndef CITYCFD_IO_H
#define CITYCFD_IO_H

#include "iostream"
#include "definitions.h"
#include "config.h"
#include "TopoFeature.h"
#include "Boundary.h"

namespace IO {
    //-- Input functions

    //-- Output functions
    void output_obj(std::vector<TopoFeature*>& allFeatures);
    void output_stl(std::vector<TopoFeature*>& allFeatures);
    void output_cityjson(std::vector<TopoFeature*>& allFeatures);

    void get_obj_pts(const Mesh& mesh, std::string& fs, std::string& bs, std::unordered_map<std::string, unsigned long>& dPts);
    void get_stl_pts(Mesh& mesh, std::string& fs);
    void get_cityjson_geom(const Mesh& mesh, nlohmann::json& g, std::unordered_map<std::string, unsigned long>& dPts, std::string primitive);

    //-- Templated function
    template <typename T>
    std::string gen_key_bucket(const T& p) {
        std::stringstream ss;
        ss << std::fixed << std::setprecision(3) << p.x() << " " << p.y() << " " << p.z();
        return ss.str();
    }

}

#endif //CITYCFD_IO_H