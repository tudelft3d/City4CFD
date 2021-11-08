#ifndef CITYCFD_CONFIG_H
#define CITYCFD_CONFIG_H

#include "types.h"
#include "CGALTypes.h"
#include "boost/filesystem.hpp"

namespace fs = boost::filesystem;

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
    extern boost::variant<bool, double, Polygon_2> influRegionConfig;
    extern boost::variant<bool, double, Polygon_2> domainBndConfig;
    extern DomainType            bpgDomainType;
    extern bool                  bpgBlockageRatioFlag;
    extern double                bpgBlockageRatio;
    extern Vector_2              flowDirection;
    extern std::vector<double>   bpgDomainSize;
    extern double                domainBuffer;

    //-- Reconstruction
    extern std::string  lod;
    extern double       buildingPercentile;

    //-- Polygons related
    extern double            edgeMaxLen;

    //-- Output setup
    extern fs::path                  workDir;
    extern fs::path                  outputDir;
    extern OutputFormat              outputFormat;
    extern bool                      outputSeparately;
    extern std::vector<std::string>  outputSurfaces;
    extern int                       numSides;
    extern std::vector<int>          surfaceLayerIDs;

    //-- Data log
    extern bool               outputLog;
    extern std::string        logName;
    extern std::ostringstream log;
    extern std::ostringstream logSummary;
    extern std::vector<int>   failedBuildings;

    void validate(nlohmann::json& j);
    void set_config(nlohmann::json& j);
    void set_region(boost::variant<bool, double, Polygon_2>& regionType,
                    std::string& regionName,
                    nlohmann::json& j);
};

#endif //CITYCFD_CONFIG_H