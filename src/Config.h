/*
  City4CFD
 
  Copyright (c) 2021-2023, 3D Geoinformation Research Group, TU Delft

  This file is part of City4CFD.

  City4CFD is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  City4CFD is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with City4CFD.  If not, see <http://www.gnu.org/licenses/>.

  For any information or further details about the use of City4CFD, contact
  Ivan Pađen
  <i.paden@tudelft.nl>
  3D Geoinformation Research Group
  Delft University of Technology
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
    static void write_to_log(const std::string& msg);

    //-- Input info
    std::string              ground_xyz;             // Ground points
    std::string              buildings_xyz;          // Building points
    std::string              gisdata;                // Building Polygons
    std::vector<std::string> topoLayers = {};        // Other polygons
    std::string              importedBuildingsPath;  // Additional pre-reconstructed buildings

    bool                     avoidBadPolys      = false;

    //-- Domain setup
    Point_2     pointOfInterest;
    double      topHeight                       = 0.;
    //- Influ region and domain bnd
    boost::variant<bool, double, Polygon_2> influRegionConfig;
    boost::variant<bool, double, Polygon_2> domainBndConfig;
    DomainType            bpgDomainType;
    bool                  bpgBlockageRatioFlag  = false;
    double                bpgBlockageRatio      = 0.03;
    Vector_2              flowDirection         = Vector_2(1,0);
    std::vector<double>   bpgDomainSize;
    double                domainBuffer          = -global::largnum;

    //-- Reconstruction
    //- Terrain
    double    terrainThinning                   = 0.;
    bool      smoothTerrain                     = false;
    int       nSmoothIterations                 = 0;
    int       maxSmoothPts                      = -9999;
    bool      flatTerrain                       = false;
    //- Buildings
    std::string buildingUniqueId;
    std::string lod;
    bool        refineReconstructedBuildings    = false;
    double      buildingPercentile              = -9999.; // Handled by schema
    double      minHeight                       = 2.;
    bool        reconstructFailed               = false;
    // Height from attributes
    std::string buildingHeightAttribute;
    std::string floorAttribute;
    double      floorHeight                     = 9999.; // Handled by schema
    bool        buildingHeightAttrAdv           = false;
    //- Imported buildings
    bool        importAdvantage                 = false;
    bool        importTrueHeight                = true;
    bool        refineImportedBuildings         = false;
    std::string importLoD                       = "9999";
    //- Boundary
    bool  reconstructBoundaries                 = false;

    //-- Polygons related
    double                edgeMaxLen            = -9999.; // Handled by schema
    std::map<int, double> flattenSurfaces;
    std::vector<int>     flattenVertBorder;

    //-- Output
    fs::path                  workDir;
    fs::path                  outputDir         = fs::current_path();
    std::string               outputFileName;
    GeomFormat                outputFormat;
    bool                      outputSeparately  = false;
    std::vector<std::string>  outputSurfaces    = {"Terrain", "Buildings"};
    int                       numSides          = 1;
    std::vector<int>          surfaceLayerIDs;

    //-- Data log
    bool               outputLog                = false;
    std::string        logName                  = std::string("log");
    std::ostringstream log;
    std::ostringstream logSummary;
    std::vector<int>   failedBuildings;

    //-- Experimental
    bool       clip                             = false;
    bool       handleSelfIntersect              = false;
    bool       alphaWrap                        = false;

    //-- Other settings
    const int searchtree_bucket_size = 100;
};

#endif //CITY4CFD_CONFIG_H