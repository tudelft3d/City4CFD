/*
  City4CFD
 
  Copyright (c) 2021-2025, 3D Geoinformation Research Group, TU Delft

  This file is part of City4CFD.

  City4CFD is free software: you can redistribute it and/or modify
  it under the terms of the GNU Affero General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  City4CFD is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Affero General Public License for more details.

  You should have received a copy of the GNU Affero General Public License
  along with City4CFD.  If not, see <http://www.gnu.org/licenses/>.

  For any information or further details about the use of City4CFD, contact
  Ivan Pađen
  <i.paden@tudelft.nl>
  3D Geoinformation Research Group
  Delft University of Technology
*/

#include "Config.h"

#include "valijson/adapters/nlohmann_json_adapter.hpp"
#include "valijson/schema.hpp"
#include "valijson/schema_parser.hpp"
#include "valijson/validator.hpp"

#include "io.h"
#include "geomutils.h"

#include "configSchema.inc"

#include  <boost/algorithm/string/predicate.hpp>

void Config::validate(nlohmann::json& j) {
    using valijson::Schema;
    using valijson::SchemaParser;
    using valijson::Validator;
    using valijson::adapters::NlohmannJsonAdapter;
    using valijson::ValidationResults;

    //- Load schema
    Schema configSchema;
    SchemaParser parser;
    NlohmannJsonAdapter configSchemaAdapter(jsonschema::schema);
    // debug
//    nlohmann::json schema = nlohmann::json::parse(std::ifstream("../../data/input/schema.json"), nullptr, true, true);
//    NlohmannJsonAdapter configSchemaAdapter(schema);

    parser.populateSchema(configSchemaAdapter, configSchema);

    //- Validate with schema
    Validator validator;
    NlohmannJsonAdapter configAdapter(j);
    ValidationResults results;
    if (!validator.validate(configSchema, configAdapter, &results)) {
        std::stringstream err_oss;
        err_oss << "Validation failed." << std::endl;
        ValidationResults::Error error;
        int error_num = 1;
        while (results.popError(error))
        {
            std::string context;
            std::vector<std::string>::iterator itr = error.context.begin();
            for (; itr != error.context.end(); itr++)
                context += *itr;

            err_oss << "Error #" << error_num << std::endl
                    << "  context: " << context << std::endl
                    << "  desc:    " << error.description << std::endl;
            ++error_num;
        }
        throw city4cfd_error(err_oss.str());
    }
}

void Config::set_config(nlohmann::json& j) {
    //-- Point cloud configuration
    if (j.contains("point_clouds")) {
        if (j["point_clouds"].contains("ground")) ground_xyz = j["point_clouds"]["ground"];
        if (j["point_clouds"].contains("buildings")) buildings_xyz = j["point_clouds"]["buildings"];
    }
    //-- Domain setup
    pointOfInterest = Point_2(j["point_of_interest"][0], j["point_of_interest"][1]);

    //- Reconstruction regions setup
    // need to define and add all reconstruction-specific parameters
    int buildingOutputLayerID = 1;
    for (auto regionJson : j["reconstruction_regions"]) {
        std::shared_ptr<ReconRegion> reconRegion = std::make_shared<ReconRegion>();

        //- Influence region
        Config::get().set_region(reconRegion->influRegionConfig, "influence_region", regionJson);

        reconRegion->lod = regionJson["lod"].front();
        if (regionJson.contains("import_advantage"))
            reconRegion->importAdvantage = regionJson["import_advantage"];
        if (regionJson.contains("bpg_influence_region_extra"))
            reconRegion->bpgInfluExtra = regionJson["bpg_influence_region_extra"];
        if (regionJson.contains("complexity_factor"))
            reconRegion->complexityFactor = regionJson["complexity_factor"];
        if (regionJson.contains("lod13_step_height"))
            reconRegion->lod13StepHeight = regionJson["lod13_step_height"];
        if (regionJson.contains("validate"))
            reconRegion->validate = regionJson["validate"];
        if (regionJson.contains("enforce_validity"))
            reconRegion->enforceValidity = regionJson["enforce_validity"];
        if (regionJson.contains("relative_alpha"))
            reconRegion->relativeAlpha = regionJson["relative_alpha"];
        if (regionJson.contains("relative_offset"))
            reconRegion->relativeOffset = regionJson["relative_offset"];
        if (regionJson.contains("skip_gap_closing"))
            reconRegion->skipGapClosing = regionJson["skip_gap_closing"];

        // validate regardless if enforce_validity is present
        if (!reconRegion->enforceValidity.empty()) reconRegion->validate = true;

        // set the building (belonging to recon region) output layer id
        reconRegion->outputLayerID = buildingOutputLayerID;
        ++buildingOutputLayerID;

        reconRegions.push_back(reconRegion);
    }

    // Handle output surface names for buildings
    if (j["reconstruction_regions"].size() < 2) {
        outputSurfaces.emplace_back("Buildings");
    } else {
        for (int i = 0; i < j["reconstruction_regions"].size(); ++i) {
            outputSurfaces.emplace_back("Buildings_" + std::to_string(i));
        }
    }

    //- Domain boundaries
    Config::get().set_region(domainBndConfig, "domain_bnd", j);
    // Define domain type if using BPG
    if (std::holds_alternative<bool>(domainBndConfig)) {
        if (j.contains("flow_direction"))
            flowDirection = Vector_2(j["flow_direction"][0], j["flow_direction"][1]);

        std::string bpgDomainConfig = j["bnd_type_bpg"];
        if (boost::iequals(bpgDomainConfig, "round")) {
            bpgDomainType = ROUND;
            if (j.contains("bpg_domain_size")) {
                bpgDomainSize = {j["bpg_domain_size"][0], j["bpg_domain_size"][1]};
            } else bpgDomainSize = {15, 6}; // BPG
        } else if (boost::iequals(bpgDomainConfig, "rectangle")) {
            bpgDomainType = RECTANGLE;
        } else if (boost::iequals(bpgDomainConfig, "oval")) {
            bpgDomainType = OVAL;
        }
    }
    // Sort out few specifics for domain types
    if (bpgDomainType == RECTANGLE || bpgDomainType == OVAL) {
        if (j.contains("bpg_domain_size")) {
            bpgDomainSize = {j["bpg_domain_size"][0], j["bpg_domain_size"][1], j["bpg_domain_size"][2], j["bpg_domain_size"][3]};
        } else bpgDomainSize = {5, 5, 15, 6}; // BPG
    }

    // Set domain side and top
    if (std::holds_alternative<Polygon_2>(domainBndConfig)) {
        Polygon_2 poly = std::get<Polygon_2>(domainBndConfig);
        numSides = poly.size();
        for (int i = 0; i < poly.size(); ++i) {
            outputSurfaces.emplace_back("Side_" + std::to_string(i % poly.size()));
        }
    } else if (std::holds_alternative<double>(domainBndConfig) || bpgDomainType != RECTANGLE) {
        outputSurfaces.emplace_back("Sides");
    } else if (bpgDomainType == RECTANGLE) { // Expand output surfaces with front and back
        numSides = 4;
        outputSurfaces.emplace_back("Side_1");
        outputSurfaces.emplace_back("Back");
        outputSurfaces.emplace_back("Side_2");
        outputSurfaces.emplace_back("Front");
    }
    outputSurfaces.emplace_back("Top");

    //-- Polygon configuration
    int i = 0;
    int surfLayerIdx = outputSurfaces.size(); // 0 terrain, 1-x buildings, x-y boundaries, y-z surface layers
    for (auto& poly : j["polygons"]) {
        if (poly["type"] == "Building") {
            gisdata = poly["path"];
            if (poly.contains("unique_id"))
                buildingUniqueId = poly["unique_id"];
            if (poly.contains("height_attribute"))
                buildingHeightAttribute = poly["height_attribute"];
            if (poly.contains("height_attribute_advantage"))
                buildingHeightAttrAdv = poly["height_attribute_advantage"];
            if (poly.contains("floor_attribute"))
                floorAttribute = poly["floor_attribute"];
            if (poly.contains("floor_height"))
                floorHeight = (double)poly["floor_height"];
            if (poly.contains("avoid_bad_polys"))
                avoidBadPolys = poly["avoid_bad_polys"];
            if (poly.contains("refine"))
                refineReconstructed= poly["refine"];
        }
        if (poly["type"] == "SurfaceLayer") {
            topoLayers.push_back(poly["path"]);
            if (poly.contains("layer_name")) {
                outputSurfaces.push_back(poly["layer_name"]);
            } else {
                outputSurfaces.push_back("SurfaceLayer" + std::to_string(++i));
            }
            if (poly.contains("flatten_surface")) {
                if (poly["flatten_surface"]) {
                    flattenSurfaces[surfLayerIdx] = poly["surface_percentile"];
                }
            }
            if (poly.contains("flatten_vertical_border")) {
                if (poly["flatten_vertical_border"]) {
                    flattenVertBorder.push_back(surfLayerIdx);
                }
            }
            ++surfLayerIdx;
        }
    }

    // Blockage ratio
    if (j.contains("bpg_blockage_ratio")) {
        if (j["bpg_blockage_ratio"].is_boolean()) {
            bpgBlockageRatioFlag = j["bpg_blockage_ratio"];
        } else if (j["bpg_blockage_ratio"].is_number()) {
            bpgBlockageRatioFlag = true;
            bpgBlockageRatio = (double)j["bpg_blockage_ratio"] / 100;
        }
    }

    // Top height
    if (j.contains("top_height"))
        topHeight = j["top_height"];

    // Buffer region
    if (j.contains("buffer_region"))
        domainBuffer = j["buffer_region"];

    //-- Reconstruction
    // Terrain
    if (j.contains("terrain_thinning"))
        terrainThinning = j["terrain_thinning"];
    if (j.contains("smooth_terrain")) {
        if (j["smooth_terrain"].contains("iterations")) {
            smoothTerrain = true;
            if ((int)j["smooth_terrain"]["iterations"] > 0) {
                nSmoothIterations = (int)j["smooth_terrain"]["iterations"];
            } else { smoothTerrain = false; }
        }
        if (j["smooth_terrain"].contains("max_pts")) {
            maxSmoothPts = (int)j["smooth_terrain"]["max_pts"];
        } else { smoothTerrain = false; }
    }
    if (j.contains("flat_terrain"))
        flatTerrain = j["flat_terrain"];

    // Buildings
    if (j.contains("building_percentile"))
        buildingPercentile = (double)j["building_percentile"].front() / 100.;
    if (j.contains("min_height"))
        minHeight = j["min_height"];
    if (j.contains("min_area"))
        minArea = j["min_area"];
    if (j.contains("reconstruct_failed"))
        reconstructFailed = j["reconstruct_failed"];
    if (j.contains("intersect_buildings_terrain"))
       intersectBuildingsTerrain = j["intersect_buildings_terrain"];

    // Imported buildings
    if (j.contains("import_geometries")) {
        importedBuildingsPath = j["import_geometries"]["path"];
        importTrueHeight      = j["import_geometries"]["true_height"];
        if (j["import_geometries"].contains("lod"))
            importLoD = j["import_geometries"]["lod"];
        if (j["import_geometries"].contains("refine"))
            refineImported = j["import_geometries"]["refine"];
    }

    // Boundary
    if (j.contains("reconstruct_boundaries"))
        reconstructBoundaries = j["reconstruct_boundaries"];

    //-- Polygons related
    edgeMaxLen = j["edge_max_len"];

    //-- Output
    //- File name
    outputFileName = j["output_file_name"];

    //- Output format
    std::string outputFormatConfig = j["output_format"];
    if (boost::iequals(outputFormatConfig, "obj")) {
        outputFormat = OBJ;
    } else if (boost::iequals(outputFormatConfig, "stl")) {
        outputFormat = STL;
    } else if (boost::iequals(outputFormatConfig, "cityjson")) {
        outputFormat = CityJSON;
    } else throw std::invalid_argument(std::string("'" + outputFormatConfig + "'" + " is unsupported file format!"));

    outputSeparately = j["output_separately"];

    //-- Data log
    if (j.contains("output_log")) {
        outputLog = true;
        if (j.contains("log_file")) logName = j["log_file"];
    }
    logSummary     <<"// ========================================CITY4CFD SUMMARY ======================================= //" << std::endl;
    log         << "\n// ========================================= CITY4CFD LOG ========================================= //" << std::endl;
    // Add point of interest info to log
    logSummary << "All coordinates are translated by -(" << pointOfInterest << ")"  << std::endl;

    //-- Experimental
    if (j.contains("experimental")) {
        if (j["experimental"].contains("clip"))
            clip = j["experimental"]["clip"];
        if (j["experimental"].contains("handle_self_intersections"))
            handleSelfIntersect = j["experimental"]["handle_self_intersections"];
        if (j["experimental"].contains("alpha_wrap_all"))
            alphaWrapAll = j["experimental"]["alpha_wrap_all"];
    }
}

//-- influRegion and domainBndConfig flow control
void Config::set_region(std::variant<bool, double, Polygon_2>& regionType,
                        const std::string regionName,
                        nlohmann::json& j) {
    if (j[regionName].is_string()) { // Search for GeoJSON polygon
        std::string polyFilePath = (std::string)j[regionName];
        if(!fs::exists(polyFilePath)) {
            throw std::invalid_argument(std::string("Cannot find polygon file '" +
                                                    polyFilePath + "' for " + regionName));
        }
        //-- Read poly
        Polygon_2 tempPoly;
        JsonVectorPtr influJsonPoly;
        IO::read_geojson_polygons(polyFilePath, influJsonPoly);
        for (auto& coords : influJsonPoly.front()->at("geometry").at("coordinates").front()) { // I know it should be only 1 polygon with 1 ring
            tempPoly.push_back(Point_2((double)coords[0] - Config::get().pointOfInterest.x(),
                                       (double)coords[1] - Config::get().pointOfInterest.y()));
        }
        //-- Prepare poly
        geomutils::pop_back_if_equal_to_front(tempPoly);
        if (tempPoly.is_clockwise_oriented()) tempPoly.reverse_orientation();

        regionType = tempPoly;
    } else if (j[regionName].size() > 2) { // Explicitly defined region polygon with points
        Polygon_2 tempPoly;
        for (auto& pt : j[regionName]) tempPoly.push_back(Point_2((double)pt[0] - Config::get().pointOfInterest.x(),
                                                                  (double)pt[1] - Config::get().pointOfInterest.y()));
        regionType = tempPoly;
    } else if (j[regionName].is_number() || j[regionName].is_array() && j[regionName][0].is_number()) { // Influ region radius
        regionType = (double)j[regionName].front();
    } else regionType = true; // Leave it to BPG
}

void Config::write_to_log(const std::string& msg) {
    #pragma omp critical
    Config::get().log << msg << std::endl;
}