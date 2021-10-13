#ifndef CITYCFD_CONFIG_H
#define CITYCFD_CONFIG_H

#include "definitions.h"

namespace config {
    //-- Input info
    extern std::string points_xyz;
    extern std::string gisdata;
    extern std::string buildings_xyz;
    extern std::vector<std::string> topoLayers;

    //-- Output info
    extern std::string outputFileName;

    //-- Domain size
    extern Point_2             pointOfInterest;
    extern double              dimOfDomain;
    extern double              topHeight;
    extern boost::variant<bool, double, std::string, Polygon_2> influenceRegion;

    //-- Reconstruction
    extern double       lod;
    extern double       buildingPercentile;

    //-- Polygons related
    extern double            edgeMaxLen;

    //-- Output flags
    extern fs::path     outputDir;
    extern OutputFormat outputFormat;
    extern bool         outputSeparately;
    // note: handle when radiusOfInterst is larger than dimOfDomain

    void set_config(nlohmann::json& j);
};

#endif //CITYCFD_CONFIG_H