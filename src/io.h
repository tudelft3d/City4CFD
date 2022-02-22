#ifndef CITYCFD_IO_H
#define CITYCFD_IO_H

#include "types.h"
#include "CGALTypes.h"

namespace IO {
    //-- Input functions
    void read_config(std::string& config_path);
    bool read_point_cloud(std::string& file, Point_set_3& pc);
    void read_geojson_polygons(std::string& file, JsonVector& jsonPolygons);
    void read_explicit_geometries(std::string& file, JsonVector& importedBuildings, std::vector<Point_3>& importedBuildingPts);

    //-- Output functions
    void print_progress_bar(int percent);

    void output_obj(const OutputFeatures& allFeatures);
    void output_stl(const OutputFeatures& allFeatures);
    void output_cityjson(const OutputFeatures& allFeatures);

    void get_obj_pts(const Mesh& mesh, std::string& fs, std::string& bs, std::unordered_map<std::string, int>& dPts);
    void get_stl_pts(Mesh& mesh, std::string& fs);
    void get_cityjson_geom(const Mesh& mesh, nlohmann::json& g, std::unordered_map<std::string, int>& dPts, std::string primitive);

    bool not_small(std::vector<int> idxLst);

    void output_log();

    std::string gen_key_bucket(const Point_2 p);

    //-- Templated function
    template <typename T> std::string gen_key_bucket(const T& p);
}

#endif //CITYCFD_IO_H