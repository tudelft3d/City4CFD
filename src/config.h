/*
  Copyright (c) 2021-2022,
  Ivan Pađen <i.paden@tudelft.nl>
  3D Geoinformation,
  Delft University of Technology

  This file is part of City4CFD.

  City4CFD is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  City4CFD is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

  You should have received a copy of the GNU General Public License
  along with this program. If not, see <http://www.gnu.org/licenses/>
*/

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
    extern double    terrainThinning;
    extern bool      smoothTerrain;
    //- Buildings
    extern std::string    buildingUniqueId;
    extern std::string    lod;
    extern double         buildingPercentile;
    extern bool           clip;
    extern bool           handleSelfIntersections;
    // Attributes
    extern std::string    buildingHeightAttribute;
    extern std::string    floorAttribute;
    extern double         floorHeight;
    extern bool           buildingHeightAttributeAdvantage;
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
    extern bool                 outputLog;
    extern std::string          logName;
    extern std::ostringstream   log;
    extern std::ostringstream   logSummary;
    extern std::vector<int>     failedBuildings;

    void validate(nlohmann::json& j);
    void set_config(nlohmann::json& j);
    void set_region(boost::variant<bool, double, Polygon_2>& regionType,
                    std::string& regionName,
                    nlohmann::json& j);
};

#endif //CITY4CFD_CONFIG_H