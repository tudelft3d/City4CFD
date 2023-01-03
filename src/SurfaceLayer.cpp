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

#include "SurfaceLayer.h"

#include "geomutils.h"

SurfaceLayer::SurfaceLayer()
        : PolyFeature() {}

SurfaceLayer::SurfaceLayer(const int outputLayerID)
        : PolyFeature(outputLayerID) {}

SurfaceLayer::SurfaceLayer(const nlohmann::json& poly)
        : PolyFeature(poly) {}

SurfaceLayer::SurfaceLayer(const nlohmann::json& poly, const int outputLayerID)
        : PolyFeature(poly, outputLayerID) {}

SurfaceLayer::~SurfaceLayer() = default;

void SurfaceLayer::check_feature_scope(const Polygon_2& bndPoly) {
    //-- Exclude all polygons that have at least one
    //-- vertex outside the domain
    for (auto& vert : _poly.outer_boundary()) {
        if (!geomutils::point_in_poly(vert, bndPoly)) {
            this->deactivate();
            return;
        }
    }
}

void SurfaceLayer::get_cityjson_info(nlohmann::json& b) const {

}

void SurfaceLayer::get_cityjson_semantics(nlohmann::json& g) const {

}

std::string SurfaceLayer::get_cityjson_primitive() const {
    return "Dunno yet";
}

TopoClass SurfaceLayer::get_class() const {
    return SURFACELAYER;
}

std::string SurfaceLayer::get_class_name() const {
    return "SurfaceLayer";
}