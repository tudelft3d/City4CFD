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

    //-- Domain setup
    extern Point_2             pointOfInterest;
    extern double              topHeight;
    extern boost::variant<bool, double, std::string, Polygon_2> influRegionConfig;
    extern boost::variant<bool, double, std::string, Polygon_2> domainBndConfig;
    extern DomainType            bpgDomainType;
    extern Vector_2              flowDirection;
    extern std::vector<double>   bpgDomainSize;
    extern std::vector<Vector_2> enlargeDomainVec;
    extern double                domainBuffer;

    //-- Reconstruction
    extern double       lod;
    extern double       buildingPercentile;

    //-- Polygons related
    extern double            edgeMaxLen;

    //-- Output flags
    extern fs::path     outputDir;
    extern OutputFormat outputFormat;
    extern bool         outputSeparately;

    void set_config(nlohmann::json& j);
    void set_region(boost::variant<bool, double, std::string, Polygon_2>& regionType,
                    std::string& regionName,
                    nlohmann::json& j);
};

#endif //CITYCFD_CONFIG_H