/*
  City4CFD

  Copyright (c) 2021-2024, 3D Geoinformation Research Group, TU Delft

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

#include "ReconstructedBuilding.h"

#include "geomutils.h"
#include "LoD12.h"
#include "ImportedBuilding.h"

ReconstructedBuilding::ReconstructedBuilding()
        : Building(), m_attributeHeight(-global::largnum),
          m_attributeHeightAdvantage(Config::get().buildingHeightAttrAdv) {}

ReconstructedBuilding::ReconstructedBuilding(const Mesh& mesh)
        : ReconstructedBuilding() {
    m_mesh = mesh;
}

/*
ReconstructedBuilding::ReconstructedBuilding(const nlohmann::json& poly)
        : Building(poly), m_searchTree(nullptr),
          m_attributeHeight(-9999), m_attributeHeightAdvantage(Config::get().buildingHeightAttrAdv) {
    if (!Config::get().buildingUniqueId.empty()) {
        m_id = poly["properties"][Config::get().buildingUniqueId].dump();
    }
    if (poly["properties"].contains(Config::get().buildingHeightAttribute)) {
        if (poly["properties"][Config::get().buildingHeightAttribute].is_number()) {
            m_attributeHeight = poly["properties"][Config::get().buildingHeightAttribute];
        }
    } else if (poly["properties"].contains(Config::get().floorAttribute)) {
        m_attributeHeight = (double)poly["properties"][Config::get().floorAttribute] * Config::get().floorHeight;
    }
}
*/

ReconstructedBuilding::ReconstructedBuilding(const nlohmann::json& poly)
        : Building(poly), m_attributeHeight(-global::largnum),
          m_attributeHeightAdvantage(Config::get().buildingHeightAttrAdv) {
    if (!Config::get().buildingUniqueId.empty() && poly["properties"].contains(Config::get().buildingUniqueId)) {
        m_id = poly["properties"][Config::get().buildingUniqueId].dump();
    } else {
        m_id = std::to_string(m_polyInternalID);
    }
    if (poly["properties"].contains(Config::get().buildingHeightAttribute)) {
        if (poly["properties"][Config::get().buildingHeightAttribute].is_number()) {
            m_attributeHeight = poly["properties"][Config::get().buildingHeightAttribute];
        }
    } else if (poly["properties"].contains(Config::get().floorAttribute)) {
        if (poly["properties"][Config::get().floorAttribute].is_number()) {
            m_attributeHeight = (double) poly["properties"][Config::get().floorAttribute] * Config::get().floorHeight;
        }
    }
    if (!this->is_active()) { // It can only fail if the polygon is not simple
        this->mark_as_failed();
        Config::write_to_log("Failed to import building polygon ID:" + m_id
                             + ". Polygon is not simple.");
    }
}

ReconstructedBuilding::ReconstructedBuilding(const Polygon_with_attr& poly)
        : Building(poly), m_attributeHeight(-global::largnum),
          m_attributeHeightAdvantage(Config::get().buildingHeightAttrAdv) {
    // Check for the polygon ID attribute
    auto idIt = poly.attributes.find(Config::get().buildingUniqueId);
    if (idIt != poly.attributes.end()) {
        m_id = idIt->second;
     } else {
        m_id = std::to_string(m_polyInternalID);
    }
    // Check for the building height attribute
    auto buildingHeightAttrIt = poly.attributes.find(Config::get().buildingHeightAttribute);
    auto numFloorsAttrIT = poly.attributes.find(Config::get().floorAttribute);
    if (buildingHeightAttrIt != poly.attributes.end()) {
        m_attributeHeight = std::stod(buildingHeightAttrIt->second);
    } else if (numFloorsAttrIT != poly.attributes.end()) { // Check for the number of floors attribute
        m_attributeHeight = std::stod(buildingHeightAttrIt->second);
    }
    if (!this->is_active()) { // It can only fail if the polygon is not simple
        this->mark_as_failed();
        Config::write_to_log("Failed to import building polygon ID:" + m_id
                             + ". Polygon is not simple.");
    }
}

ReconstructedBuilding::ReconstructedBuilding(const std::shared_ptr<ImportedBuilding>& importedBuilding)
        : Building(importedBuilding->get_poly_w_attr()),
          m_attributeHeight(-global::largnum), m_attributeHeightAdvantage(Config::get().buildingHeightAttrAdv) {
    m_ptsPtr = importedBuilding->get_points();
    m_groundElevations = importedBuilding->get_ground_elevations();
    m_id = importedBuilding->get_id();
    m_reconSettings = importedBuilding->get_reconstruction_settings();
    m_outputLayerID = importedBuilding->get_output_layer_id();
}

ReconstructedBuilding::~ReconstructedBuilding() = default;

/*
 * Calculate building elevation without mesh reconstruction
 */
double ReconstructedBuilding::get_elevation() {
    if (m_elevation < -global::largnum + global::smallnum) { // calculate if not already available
        if (m_attributeHeightAdvantage && m_attributeHeight > 0) { // get height from attribute
            m_elevation = this->ground_elevation() + m_attributeHeight;
        } else if (m_ptsPtr->empty()) { // set height as minimum if not able to calculate
            m_elevation = this->ground_elevation();
            Config::write_to_log("Building ID: " + this->get_id() + " Missing points for elevation calculation."
                                 + "Using minimum of " + std::to_string(Config::get().minHeight) + "m");
        } else { // else calculate as a percentile from config
            // gather all building elevations
            std::vector<double> buildingElevations;
            for (auto& pt : m_ptsPtr->points()) {
                buildingElevations.push_back(pt.z());
            }
            // calculate percentile
            m_elevation = geomutils::percentile(buildingElevations, Config::get().buildingPercentile);
        }
    }
    return m_elevation;
}

void ReconstructedBuilding::reconstruct() {
    m_mesh.clear();
    if (m_clipBottom || Config::get().intersectBuildingsTerrain) {
        this->translate_footprint(-5);
    }
    //-- Check if reconstructing from height attribute takes precedence
    if (m_attributeHeightAdvantage) {
        this->reconstruct_from_attribute();
        return;
    }
    //-- Reconstruction fallbacks if there are no points belonging to the polygon
    if (m_ptsPtr->empty()) {
        // reconstruct using attribute
        if (this->reconstruct_again_from_attribute("Found no points belonging to the building")) {
            return;
        } else if (Config::get().reconstructFailed) { // fall back to minimum height if defined as an argument
            Config::write_to_log("Building ID:" + this->get_id()
                                 + " Found no points belonging to the building."
                                 + "Reconstructing with the minimum height of "
                                 + std::to_string(Config::get().minHeight) + " m");
            m_elevation = this->ground_elevation() + Config::get().minHeight;
        } else { // exception handling when cannot reconstruct
            this->deactivate();
            throw std::domain_error("Found no points belonging to the building");
        }
    }
    //-- LoD12 reconstruction
    if (this->get_height() < Config::get().minHeight) { // elevation calculated here
        Config::write_to_log("Building ID: " + this->get_id()
                             + " Height lower than minimum prescribed height of "
                             + std::to_string(Config::get().minHeight) + " m");
        m_elevation = this->ground_elevation() + Config::get().minHeight;
    }
    LoD12 lod12(m_poly, m_groundElevations, m_elevation);
    lod12.reconstruct(m_mesh);

    if (m_clipBottom || Config::get().intersectBuildingsTerrain) {
        this->translate_footprint(5);
    }
    if (Config::get().refineReconstructed) this->refine();
}

void ReconstructedBuilding::reconstruct_flat_terrain() {
    m_mesh.clear();
    if (m_clipBottom || Config::get().intersectBuildingsTerrain) {
        this->translate_footprint(-5);
    }
    // the new height was previously calculated
    LoD12 lod12HeightAttribute(m_poly, m_groundElevations, this->get_elevation());
    lod12HeightAttribute.reconstruct(m_mesh);

    if (m_clipBottom || Config::get().intersectBuildingsTerrain) {
        this->translate_footprint(5);
    }
}

void ReconstructedBuilding::get_cityjson_info(nlohmann::json& b) const {
    b["type"] = "Building";
//  b["attributes"];
//    get_cityjson_attributes(b, _attributes);
//    float hbase = z_to_float(this->get_height_base());
//    float h = z_to_float(this->get_height());
//    b["attributes"]["TerrainHeight"] = m_baseElevations.back(); // temp - will calculate avg for every footprint
    b["attributes"]["measuredHeight"] = m_elevation - geomutils::avg(m_groundElevations[0]);
}

void ReconstructedBuilding::get_cityjson_semantics(nlohmann::json& g) const { // Temp for checking CGAL mesh properties
    Face_property semantics;
    bool foundProperty;
    boost::tie(semantics, foundProperty) = m_mesh.property_map<face_descriptor, std::string>("f:semantics");
    //   auto semantics = m_mesh.property_map<face_descriptor, std::string>("f:semantics");
    if (!foundProperty) throw std::runtime_error("Semantic property map not found!");

    std::unordered_map<std::string, int> surfaceId;
    surfaceId["RoofSurface"]   = 0; g["semantics"]["surfaces"][0]["type"] = "RoofSurface";
    surfaceId["GroundSurface"] = 1; g["semantics"]["surfaces"][1]["type"] = "GroundSurface";
    surfaceId["WallSurface"]   = 2; g["semantics"]["surfaces"][2]["type"] = "WallSurface";

    for (auto& faceIdx : m_mesh.faces()) {
        auto it = surfaceId.find(semantics[faceIdx]);
        if (it == surfaceId.end()) throw std::runtime_error("Could not find semantic attribute!");

        g["semantics"]["values"][faceIdx.idx()] = it->second;
    }
}

void ReconstructedBuilding::reconstruct_from_attribute() {
    //-- Check attribute height
    if (m_attributeHeight <= 0) {
        this->deactivate();
        throw std::runtime_error("Attribute height from geojson file is invalid!");
    }
    //-- Set the height from attribute and reconstruct
    if (m_attributeHeight < Config::get().minHeight) { // in case of a small height
        Config::write_to_log("Building ID: " + this->get_id()
                             + " Height lower than minimum prescribed height of "
                             + std::to_string(Config::get().minHeight) + " m");
        m_elevation = this->ground_elevation() + Config::get().minHeight;
    } else {
        m_elevation = this->ground_elevation() + m_attributeHeight;
    }
    LoD12 lod12HeightAttribute(m_poly, m_groundElevations, m_elevation);
    lod12HeightAttribute.reconstruct(m_mesh);
}

bool ReconstructedBuilding::reconstruct_again_from_attribute(const std::string& reason) {
    if (m_attributeHeight > 0) {
        Config::write_to_log("Building ID: " + m_id + " Failed to reconstruct using point cloud. Reason: "
                             + reason + ". Reconstructing using height attribute from GeoJSON polygon.");
        this->reconstruct_from_attribute();
        return true;
    } else {
        return false;
    }
}