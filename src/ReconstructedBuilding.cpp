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

#include "ReconstructedBuilding.h"

#include "geomutils.h"
#include "LoD12.h"

ReconstructedBuilding::ReconstructedBuilding()
        : Building(), _searchTree(nullptr),
        _attributeHeight(-9999), _attributeHeightAdvantage(Config::get().buildingHeightAttrAdv) {}

ReconstructedBuilding::ReconstructedBuilding(const int internalID)
        : Building(internalID), _searchTree(nullptr),
          _attributeHeight(-9999), _attributeHeightAdvantage(Config::get().buildingHeightAttrAdv) {}

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

ReconstructedBuilding::ReconstructedBuilding(const nlohmann::json& poly, const int internalID)
        : Building(poly, internalID), _searchTree(nullptr),
          _attributeHeight(-9999), _attributeHeightAdvantage(Config::get().buildingHeightAttrAdv) {
    if (!Config::get().buildingUniqueId.empty()) {
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
}

ReconstructedBuilding::~ReconstructedBuilding() = default;

void ReconstructedBuilding::set_search_tree(const std::shared_ptr<SearchTree>& searchTree) {
    _searchTree = searchTree;
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

    //-- Take tree subset bounded by the polygon
    std::vector<Point_3> subsetPts;
    Point_3 bbox1(_poly.bbox().xmin(), _poly.bbox().ymin(), -g_largnum);
    Point_3 bbox2(_poly.bbox().xmax(), _poly.bbox().ymax(), g_largnum);
    Fuzzy_iso_box pts_range(bbox1, bbox2);
    _searchTree->search(std::back_inserter(subsetPts), pts_range);

    //-- Check if subset point lies inside the polygon
    std::vector<double> building_pts;
    for (auto& pt : subsetPts) {
        if (geomutils::point_in_poly(pt, _poly)) {
            building_pts.push_back(pt.z());
        }
    }

    //-- Don't reconstruct if there are no points belonging to the polygon
    if (building_pts.empty()) {
        if (this->reconstruct_again_from_attribute("Found no points belonging to the building")) {
            return;
        } else {
            this->deactivate();
            throw std::domain_error("Found no points belonging to the building");
        }
    }

    //-- LoD12 reconstruction
    LoD12 lod12(_poly, _base_heights, building_pts);
    lod12.lod12_calc_height(_height);
    lod12.lod12_reconstruct(_mesh);

    if (lod12.get_height() < _lowHeight) { // In case of a small height
        Config::get().log << "Building height lower than minimum prescribed height, ID: " << this->get_id()
                          << std::endl;
        _height = _lowHeight;
        lod12.lod12_reconstruct(_mesh, _height);
    }

    if (_clip_bottom) {
        this->translate_footprint(5);
    }
}

void ReconstructedBuilding::reconstruct_flat_terrain() {
    _mesh.clear();
    if (_clip_bottom) {
        this->translate_footprint(-5);
    }

    LoD12 lod12HeightAttribute(_poly, _base_heights, {}, _height);
    lod12HeightAttribute.lod12_reconstruct(_mesh);

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
//    b["attributes"]["TerrainHeight"] = _baseHeights.back(); // temp - will calculate avg for every footprint
    b["attributes"]["measuredHeight"] = _height - geomutils::avg(_base_heights[0]);
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

    //-- Average the fooprint height to get ground-zero height
    std::vector<double> footprintElevations;
    for (auto& rings : _base_heights) {
        for (auto& pt : rings) footprintElevations.push_back(pt);
    }
    double baseHeight = geomutils::avg(footprintElevations);

    //-- Set the height from attribute and reconstruct
    _height = baseHeight + _attributeHeight;
    LoD12 lod12HeightAttribute(_poly, _base_heights, {}, _height);
    lod12HeightAttribute.lod12_reconstruct(_mesh);

    //-- Low height check
    if (lod12HeightAttribute.get_height() < _lowHeight) { // In case of a small height
        Config::get().log << "Building height lower than minimum prescribed height, ID: " << this->get_id()
                          << std::endl;
        _height = _lowHeight;
        lod12HeightAttribute.lod12_reconstruct(_mesh, _height);
    }
}

bool ReconstructedBuilding::reconstruct_again_from_attribute(const std::string& reason) {
    if (_attributeHeight > 0) {
        Config::get().log << "Failed to reconstruct using point cloud building ID: " << _id
                    << " Reason: " << reason
                    << ". Reconstructing using height attribute from JSON polygon." << std::endl;
        this->reconstruct_from_attribute();
        return true;
    } else {
        return false;
    }
}