#ifndef CITYCFD_CONFIG_H
#define CITYCFD_CONFIG_H

#include "definitions.h"

namespace config {
    //-- Input info
    extern const char* points_xyz;
    extern const char* gisdata;
    extern const char* buildings_xyz;
    extern std::vector<const char*> topoLayers;

    //-- Output info
    extern std::string outputFileName;

    //-- Domain size
   extern Point_2      pointOfInterest;
   extern double       radiusOfInfluRegion;
   extern double       dimOfDomain;
   extern double       topHeight;

   //-- Reconstruction
   extern double       lod;
   extern double       buildingPercentile;

   //-- Polygon
   extern double       edgeMaxLen;

   //-- Output flags
   extern OutputFormat outputFormat;
   extern bool         outputSeparately;
    // note: handle when radiusOfInterst is larger than dimOfDomain

    bool read_config_file();
};

#endif //CITYCFD_CONFIG_H