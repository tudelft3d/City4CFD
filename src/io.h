#ifndef CITYCFD_IO_H
#define CITYCFD_IO_H

#include "iostream"
#include "config.h"
#include "TopoFeature.h"
#include "Boundary.h"

namespace IO {
    //-- Input functions
    bool read_config(const char* file);
    bool read_point_cloud(const char* file, Point_set_3& pc);
    bool read_polygons(const char* file, nlohmann::json& j);

    //-- Output functions
    void output_obj(const OutputFeatures& allFeatures);
    void output_stl(const OutputFeatures& allFeatures);
    void output_cityjson(const OutputFeatures& allFeatures);

    void get_obj_pts(const Mesh& mesh, std::string& fs, std::string& bs, std::unordered_map<std::string, unsigned long>& dPts);
    void get_stl_pts(Mesh& mesh, std::string& fs);
    void get_cityjson_geom(const Mesh& mesh, nlohmann::json& g, std::unordered_map<std::string, unsigned long>& dPts, std::string primitive);

    //-- Templated function
    template <typename T> std::string gen_key_bucket(const T& p);
}

#endif //CITYCFD_IO_H