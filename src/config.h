#ifndef CITYCFD_CONFIG_H
#define CITYCFD_CONFIG_H

#include "definitions.h"

namespace config {

   extern Point_2      pointOfInterest;
   extern double       radiusOfInfluRegion;
   extern double       dimOfDomain;
   extern double       topHeight;
   extern double       buildingPercentile;
   extern OutputFormat outputFormat;
   extern bool         outputSeparately;
    // note: handle when radiusOfInterst is larger than dimOfDomain

    bool read_config_file();
};

#endif //CITYCFD_CONFIG_H