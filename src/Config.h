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

class Config {
public:
    static Config& get() {
        static Config inst;
        return inst;
    }

protected:
    Config() = default;
    Config(const Config&) = delete;
    Config& operator=(const Config&) = delete;
    ~Config() = default;

public:
    void validate(nlohmann::json& j);
    void set_config(nlohmann::json& j);
    void set_region(boost::variant<bool, double, Polygon_2>& regionType,
                    std::string& regionName,
                    nlohmann::json& j);

    //-- Input info
    std::string              points_xyz;         // Ground
    std::string              buildings_xyz;      // Buildings
    std::string              gisdata;            // Building Polygons
    std::vector<std::string> topoLayers = {};         // Other polygons
    std::string              importedBuildings;  // Additional pre-reconstructed buildings

    //-- Domain setup
    Point_2     pointOfInterest;
    double      topHeight = 0;
    //- Influ region and domain bnd
    boost::variant<bool, double, Polygon_2> influRegionConfig;
    boost::variant<bool, double, Polygon_2> domainBndConfig;
    DomainType            bpgDomainType;
    bool                  bpgBlockageRatioFlag = false;
    double                bpgBlockageRatio     = 0.03;
    Vector_2              flowDirection        = Vector_2(1,0);
    std::vector<double>   bpgDomainSize;
    double                domainBuffer         = -g_largnum;

    //-- Reconstruction
    //- Terrain
    double    terrainThinning = 0.;
    bool      smoothTerrain   = false;
    //- Buildings
    std::string buildingUniqueId;
    std::string lod;
    double      buildingPercentile;
    bool        clip                     = false;
    bool        handleSelfIntersections  = false;
    // Height from attributes
    std::string buildingHeightAttribute;
    std::string floorAttribute;
    double      floorHeight;
    bool        buildingHeightAttributeAdvantage = false;
    //- Imported buildings
    bool        importAdvantage;
    bool        importTrueHeight;
    std::string importLoD;
    //- Boundary
    bool  reconstructBoundaries = false;

    //-- Polygons related
    double                edgeMaxLen;
    std::map<int, double> averageSurfaces;

    //-- Output
    fs::path                  workDir;
    fs::path                  outputDir        = fs::current_path();
    std::string               outputFileName;
    OutputFormat              outputFormat;
    bool                      outputSeparately = false;
    std::vector<std::string>  outputSurfaces   = {"Terrain", "Buildings"};
    int                       numSides         = 1;
    std::vector<int>          surfaceLayerIDs;

    //-- Data log
    bool               outputLog = false;
    std::string        logName   = std::string("log");
    std::ostringstream log;
    std::ostringstream logSummary;
    std::vector<int>   failedBuildings;
};

#endif //CITY4CFD_CONFIG_H