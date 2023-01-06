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

#include "ReconstructedBuilding.h"

#include "geomutils.h"
#include "LoD12.h"

ReconstructedBuilding::ReconstructedBuilding()
        : Building(), _attributeHeight(-global::largnum),
          _attributeHeightAdvantage(Config::get().buildingHeightAttrAdv) {}

ReconstructedBuilding::ReconstructedBuilding(const int internalID)
        : Building(internalID), _attributeHeight(-global::largnum),
          _attributeHeightAdvantage(Config::get().buildingHeightAttrAdv) {}

ReconstructedBuilding::ReconstructedBuilding(const Mesh& mesh)
        : ReconstructedBuilding() {
    _mesh = mesh;
}

/*
ReconstructedBuilding::ReconstructedBuilding(const nlohmann::json& poly)
        : Building(poly), _searchTree(nullptr),
          _attributeHeight(-9999), _attributeHeightAdvantage(Config::get().buildingHeightAttrAdv) {
    if (!Config::get().buildingUniqueId.empty()) {
        _id = poly["properties"][Config::get().buildingUniqueId].dump();
    }
    if (poly["properties"].contains(Config::get().buildingHeightAttribute)) {
        if (poly["properties"][Config::get().buildingHeightAttribute].is_number()) {
            _attributeHeight = poly["properties"][Config::get().buildingHeightAttribute];
        }
    } else if (poly["properties"].contains(Config::get().floorAttribute)) {
        _attributeHeight = (double)poly["properties"][Config::get().floorAttribute] * Config::get().floorHeight;
    }
}
*/

ReconstructedBuilding::ReconstructedBuilding(const nlohmann::json& poly, const int internalID)
        : Building(poly, internalID), _attributeHeight(-global::largnum),
          _attributeHeightAdvantage(Config::get().buildingHeightAttrAdv) {
    if (!Config::get().buildingUniqueId.empty() && poly["properties"].contains(Config::get().buildingUniqueId)) {
        _id = poly["properties"][Config::get().buildingUniqueId].dump();
    } else {
        _id = std::to_string(internalID);
    }
    if (poly["properties"].contains(Config::get().buildingHeightAttribute)) {
        if (poly["properties"][Config::get().buildingHeightAttribute].is_number()) {
            _attributeHeight = poly["properties"][Config::get().buildingHeightAttribute];
        }
    } else if (poly["properties"].contains(Config::get().floorAttribute)) {
        if (poly["properties"][Config::get().floorAttribute].is_number()) {
            _attributeHeight = (double) poly["properties"][Config::get().floorAttribute] * Config::get().floorHeight;
        }
    }
    if (!this->is_active()) { // It can only fail if the polygon is not simple
        Config::get().failedBuildings.push_back(internalID);
        Config::write_to_log("Failed to import building polygon ID:" + _id
                          + ". Polygon is not simple.");
    }
}

ReconstructedBuilding::~ReconstructedBuilding() = default;

/*
 * Calculate building elevation without mesh reconstruction
 */
double ReconstructedBuilding::get_elevation() {
    if (_elevation < -global::largnum + global::smallnum) { // calculate if not already available
        if (_attributeHeightAdvantage && _attributeHeight > 0) { // get height from attribute
            _elevation = this->ground_elevation() + _attributeHeight;
        } else if (_ptsPtr->empty()) { // set height as minimum if not able to calculate
            _elevation = this->ground_elevation() + this->slope_height() + Config::get().minHeight;
            Config::write_to_log("Building ID: " + this->get_id() + " Missing points for elevation calculation."
                                 + "Using minimum of " + std::to_string(Config::get().minHeight) + "m");
        } else { // else calculate as a percentile from config
            // gather all building elevations
            std::vector<double> buildingElevations;
            for (auto& pt : _ptsPtr->points()) {
                buildingElevations.push_back(pt.z());
            }
            // calculate percentile
            _elevation = geomutils::percentile(buildingElevations, Config::get().buildingPercentile);
        }
    }
    return _elevation;
}

void ReconstructedBuilding::reconstruct() {
    _mesh.clear();
    if (_clip_bottom) {
        this->translate_footprint(-5);
    }
    //-- Check if reconstructing from height attribute takes precedence
    if (_attributeHeightAdvantage) {
        this->reconstruct_from_attribute();
        return;
    }
    //-- Reconstruction fallbacks if there are no points belonging to the polygon
    if (_ptsPtr->empty()) {
        // reconstruct using attribute
        if (this->reconstruct_again_from_attribute("Found no points belonging to the building")) {
            return;
        } else if (Config::get().reconstructFailed) { // fall back to minimum height if defined as an argument
            Config::write_to_log("Building ID:" + this->get_id()
                                 + " Found no points belonging to the building."
                                 + "Reconstructing with the minimum height of "
                                 + std::to_string(Config::get().minHeight) + " m");
            _elevation = this->ground_elevation() + this->slope_height() + Config::get().minHeight;
        } else { // exception handling when cannot reconstruct
            this->deactivate();
            throw std::domain_error("Found no points belonging to the building");
        }
    }
    //-- LoD12 reconstruction
    if (this->get_height() < Config::get().minHeight) { // elevation calculated here
        Config::write_to_log("Building ID:" + this->get_id()
                             + " Height lower than minimum prescribed height of "
                             + std::to_string(Config::get().minHeight) + " m");
        _elevation = this->ground_elevation() + this->slope_height() + Config::get().minHeight;
    }
    LoD12 lod12(_poly, _groundElevations, _elevation);
    lod12.reconstruct(_mesh);

    if (_clip_bottom) {
        this->translate_footprint(5);
    }
    if (Config::get().refineReconstructedBuildings) this->refine();
}

void ReconstructedBuilding::reconstruct_flat_terrain() {
    _mesh.clear();
    if (_clip_bottom) {
        this->translate_footprint(-5);
    }
    // the new height was previously calculated
    LoD12 lod12HeightAttribute(_poly, _groundElevations, this->get_elevation());
    lod12HeightAttribute.reconstruct(_mesh);

    if (_clip_bottom) {
        this->translate_footprint(5);
    }
}

void ReconstructedBuilding::get_cityjson_info(nlohmann::json& b) const {
    b["type"] = "Building";
//  b["attributes"];
//    get_cityjson_attributes(b, _attributes);
//    float hbase = z_to_float(this->get_height_base());
//    float h = z_to_float(this->get_height());
//    b["attributes"]["TerrainHeight"] = _baseElevations.back(); // temp - will calculate avg for every footprint
    b["attributes"]["measuredHeight"] = _elevation - geomutils::avg(_groundElevations[0]);
}

void ReconstructedBuilding::get_cityjson_semantics(nlohmann::json& g) const { // Temp for checking CGAL mesh properties
    Face_property semantics;
    bool foundProperty;
    boost::tie(semantics, foundProperty) = _mesh.property_map<face_descriptor, std::string>("f:semantics");
    //   auto semantics = _mesh.property_map<face_descriptor, std::string>("f:semantics");
    if (!foundProperty) throw std::runtime_error("Semantic property map not found!");

    std::unordered_map<std::string, int> surfaceId;
    surfaceId["RoofSurface"]   = 0; g["semantics"]["surfaces"][0]["type"] = "RoofSurface";
    surfaceId["GroundSurface"] = 1; g["semantics"]["surfaces"][1]["type"] = "GroundSurface";
    surfaceId["WallSurface"]   = 2; g["semantics"]["surfaces"][2]["type"] = "WallSurface";

    for (auto& faceIdx : _mesh.faces()) {
        auto it = surfaceId.find(semantics[faceIdx]);
        if (it == surfaceId.end()) throw std::runtime_error("Could not find semantic attribute!");

        g["semantics"]["values"][faceIdx.idx()] = it->second;
    }
}

void ReconstructedBuilding::reconstruct_from_attribute() {
    //-- Check attribute height
    if (_attributeHeight <= 0) {
        this->deactivate();
        throw std::runtime_error("Attribute height from geojson file is invalid!");
    }
    //-- Set the height from attribute and reconstruct
    if (_attributeHeight < Config::get().minHeight) { // in case of a small height
        Config::write_to_log("Building ID:" + this->get_id()
                             + " Height lower than minimum prescribed height of "
                             + std::to_string(Config::get().minHeight) + " m");
        _elevation = this->ground_elevation() + this->slope_height() + Config::get().minHeight;
    } else {
        _elevation = this->ground_elevation() + _attributeHeight;
    }
    LoD12 lod12HeightAttribute(_poly, _groundElevations, _elevation);
    lod12HeightAttribute.reconstruct(_mesh);
}

bool ReconstructedBuilding::reconstruct_again_from_attribute(const std::string& reason) {
    if (_attributeHeight > 0) {
        Config::write_to_log("Building ID: " + _id + " Failed to reconstruct using point cloud. Reason: "
                             + reason + ". Reconstructing using height attribute from GeoJSON polygon.");
        this->reconstruct_from_attribute();
        return true;
    } else {
        return false;
    }
}