#ifndef CITY4CFD_CONFIG_H
#define CITY4CFD_CONFIG_H

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
    extern std::string importedBuildings;

    //-- Domain setup
    extern Point_2  pointOfInterest;
    extern double   topHeight;
    //- Influ region and domain bnd
    extern boost::variant<bool, double, Polygon_2> influRegionConfig;
    extern boost::variant<bool, double, Polygon_2> domainBndConfig;
    extern DomainType            bpgDomainType;
    extern bool                  bpgBlockageRatioFlag;
    extern double                bpgBlockageRatio;
    extern Vector_2              flowDirection;
    extern std::vector<double>   bpgDomainSize;
    extern double                domainBuffer;

    //-- Reconstruction
    //- Terrain
    extern double    terrainSimplification;
    extern bool      smoothTerrain;
    //- Buildings
    extern std::string    lod;
    extern double         buildingPercentile;
    //- Imported Buildings
    extern bool        importAdvantage;
    extern bool        importTrueHeight;
    extern std::string importLoD;

    //-- Polygons related
    extern double                edgeMaxLen;
    extern std::map<int, double> averageSurfaces;

    //-- Output
    extern fs::path                  workDir;
    extern fs::path                  outputDir;
    extern std::string               outputFileName;
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

#endif //CITY4CFD_CONFIG_H