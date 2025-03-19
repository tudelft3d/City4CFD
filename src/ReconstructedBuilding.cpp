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

#include "ReconstructedBuilding.h"

#include "geomutils.h"
#include "LoD12.h"
#include "LoD22.h"
#include "ImportedBuilding.h"

#include "val3dity/src/val3dity.h"
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/Polygon_mesh_processing/polygon_mesh_to_polygon_soup.h>

ReconstructedBuilding::ReconstructedBuilding()
        : Building(), m_attributeHeight(-global::largnum),
          m_attributeHeightAdvantage(Config::get().buildingHeightAttrAdv) {}

ReconstructedBuilding::ReconstructedBuilding(const Mesh& mesh)
        : ReconstructedBuilding() {
    m_mesh = mesh;
}

ReconstructedBuilding::ReconstructedBuilding(const roofer::Mesh& rooferMesh, const ReconstructedBuildingPtr& other)
        : ReconstructedBuilding() {
    if (!other) throw city4cfd_error("Trying to pass a nullptr to ReconstructedBuilding!");
    // create new LoD22 object
    LoD22 lod22(rooferMesh);
    m_mesh = lod22.get_mesh();
    m_poly = lod22.get_footprint();
    m_groundElevations = lod22.get_base_elevations();

    // transfer other attributes from the original
    m_id = std::string(other->m_id + "_1");
    m_reconSettings = other->m_reconSettings;
    m_outputLayerID = other->m_outputLayerID;
    m_attributeHeight = other->m_attributeHeight;
    m_attributeHeightAdvantage = other->m_attributeHeightAdvantage;
    m_groundPtsPtr = other->m_groundPtsPtr;
    m_ptsPtr = other->m_ptsPtr;
    m_clipBottom = other->m_clipBottom;
}

ReconstructedBuilding::ReconstructedBuilding(const nlohmann::json& poly)
        : Building(poly), m_attributeHeight(-global::largnum),
          m_attributeHeightAdvantage(Config::get().buildingHeightAttrAdv),
          m_groundPtsPtr(std::make_shared<Point_set_3 >()) {
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
          m_attributeHeightAdvantage(Config::get().buildingHeightAttrAdv),
          m_groundPtsPtr(std::make_shared<Point_set_3 >()) {
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
          m_attributeHeight(-global::largnum), m_attributeHeightAdvantage(Config::get().buildingHeightAttrAdv),
          m_groundPtsPtr(std::make_shared<Point_set_3 >()) {
    m_ptsPtr = importedBuilding->get_points();
    m_groundElevations = importedBuilding->get_ground_elevations();
    m_id = importedBuilding->get_id();
    m_reconSettings = importedBuilding->get_reconstruction_settings();
    m_outputLayerID = importedBuilding->get_output_layer_id();
}

/*
 * Calculate building elevation without mesh reconstruction
 */
void ReconstructedBuilding::calc_elevation() {
    if (m_elevation < -global::largnum + global::smallnum) { // calculate if not already available
        if (m_attributeHeightAdvantage && m_attributeHeight > 0) { // get height from attribute
            m_elevation = this->ground_elevation() + m_attributeHeight;
        } else { // else calculate as a percentile from config
            // gather all building elevations
            std::vector<double> buildingElevations;
            for (auto& pt : m_ptsPtr->points()) {
                if (geomutils::point_in_poly_and_boundary(pt, m_poly))
                    buildingElevations.push_back(pt.z());
            }
            if (buildingElevations.empty()) { // set height as minimum if not able to calculate percentile
                m_elevation = this->ground_elevation() + Config::get().minHeight + 0.001;
                Config::write_to_log("Building ID: " + this->get_id() + " Missing points for elevation calculation."
                                     + " Using minimum of " + std::to_string(Config::get().minHeight) + "m");
            } else { // calculate percentile
                m_elevation = geomutils::percentile(buildingElevations, Config::get().buildingPercentile);
            }
        }
    }
}

void ReconstructedBuilding::reconstruct() {
    m_mesh.clear();
    assert(m_reconSettings);
    if (m_clipBottom || Config::get().intersectBuildingsTerrain) {
        this->translate_footprint(-5);
    }

    //-- Check if reconstructing from height attribute takes precedence
    if (m_attributeHeightAdvantage) {
        this->reconstruct_from_attribute();
        return;
    }

    //-- Narrow down check for points without the buffer
    bool innerPts = false;
    for (auto& pt : m_ptsPtr->points()) {
        if (geomutils::point_in_poly_and_boundary(pt, m_poly)) {
            innerPts = true;
            break;
        }
    }

    //-- Reconstruction fallbacks if there are no points belonging to the polygon
    if (!innerPts) {
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
            this->mark_as_failed();
            this->deactivate();
            throw city4cfd_error("Found no points belonging to the building");
        }
    }

#ifdef CITY4CFD_DEBUG
    std::cout << "Checking if polygon is valid and simple" << std::endl;
    for (auto& poly : m_poly.rings()) {
        if (!poly.is_simple() || poly.is_empty() || !poly.is_convex()) {
            std::cout << "Polygon is not simple, empty or convex" << std::endl;
            std::cout << "Poly id: " << m_id << std::endl;
        }
     }
    std::cout << "End polygon check" << std::endl;
#endif

    if (m_reconSettings->lod == "2.2" || m_reconSettings->lod == "1.3") {
        try {
            LoD22 lod22;
            if (m_reconSettings->lod == "2.2")
                lod22.reconstruct(m_ptsPtr, m_groundPtsPtr, m_poly, m_groundElevations,
                                  LoD22::ReconstructionConfig()
                                          .lod(22)
                                          .lambda(m_reconSettings->complexityFactor));
            else
                lod22.reconstruct(m_ptsPtr, m_groundPtsPtr, m_poly, m_groundElevations,
                                  LoD22::ReconstructionConfig()
                                  .lod(13)
                                  .lambda(m_reconSettings->complexityFactor)
                                  .lod13_step_height(m_reconSettings->lod13StepHeight));

            auto mesh = lod22.get_mesh();
            if (mesh.is_empty()) throw city4cfd_error("Unsuccessful mesh reconstruction.");

            // try easy hole plugging fix just in case
            if (!m_reconSettings->skipGapClosing && !CGAL::is_closed(mesh)) {
                // collect boundary halfedges
                typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;
                std::vector<halfedge_descriptor> borderCycles;
                PMP::extract_boundary_cycles(mesh, std::back_inserter(borderCycles));
                // fill using boundary halfedges
                for(halfedge_descriptor h : borderCycles)
                    PMP::triangulate_hole(mesh, h);
            }
            //-- Validity check
            if (m_reconSettings->validate) {
                std::vector<std::array<double, 3>> points;
                std::vector<std::vector<int>> polys;
                PMP::polygon_mesh_to_polygon_soup(mesh, points, polys);

                auto validity = val3dity::validate(points, polys);
                if (!validity["validity"]) {
                    std::string valdtyReport = "Invalid geometry, errors: " + nlohmann::to_string(validity["all_errors"]);
                    if (m_reconSettings->enforceValidity.empty()) {
                        Config::write_to_log("Building ID: " + this->get_id()
                                             + valdtyReport);
                    } else if (m_reconSettings->enforceValidity == "surface_wrap") {
                        // gather mesh information before sending to alpha wrap
                        m_groundElevations = lod22.get_base_elevations();
                        m_poly = lod22.get_footprint();
                        m_mesh = mesh;
                        throw city4cfd_error(valdtyReport);
                    } else {
                        // for lod1.2 don't use mesh, polygon, and elevations from
                        // higher lod reconstruction
                        throw city4cfd_error(valdtyReport);
                    }
                }
            }
            // get the new footprint and elevations
            m_groundElevations = lod22.get_base_elevations();
            m_poly = lod22.get_footprint();
            m_mesh = mesh;
        } catch (const std::exception& e) {
#ifdef CITY4CFD_DEBUG
            std::cout << "LoD2.2/1.3 reconstruction failed!" << std::endl;
            std::cout << "Reason: " << e.what() << std::endl;
#endif

            Config::write_to_log("Building ID: " + this->get_id()
                                 + " Failed to reconstruct at LoD2.2/1.3."
                                 + " Reason: " + e.what()
                                 + ". Falling back to LoD1.2 reconstruction/alpha wrapping");

            if (m_mesh.is_empty() ||
                m_reconSettings->enforceValidity.empty() ||
                m_reconSettings->enforceValidity == "lod1.2")
                this->reconstruct_lod12();
            else if (m_reconSettings->enforceValidity == "surface_wrap") // only alpha wrap if 1. mesh was reconstructed and 2. validity is enforced
                this->alpha_wrap(m_reconSettings->relativeAlpha, m_reconSettings->relativeOffset);
            else
                throw city4cfd_error("Error with reconstruction fallback!");
        }
    } else {
        // just reconstruct LoD22
        this->reconstruct_lod12();
    }

    if (m_clipBottom || Config::get().intersectBuildingsTerrain) {
        this->translate_footprint(5);
    }
    if (Config::get().refineReconstructed) this->refine();
}

void ReconstructedBuilding::reconstruct_lod12() {
        m_mesh.clear();
        if (this->get_height() < Config::get().minHeight) { // elevation calculated here
            Config::write_to_log("Building ID: " + this->get_id()
                                 + " Height lower than minimum prescribed height of "
                                 + std::to_string(Config::get().minHeight) + " m");
            m_elevation = this->ground_elevation() + Config::get().minHeight;
        }
        LoD12 lod12(m_poly, m_groundElevations, m_elevation);
        lod12.reconstruct(m_mesh);
}

void ReconstructedBuilding::insert_terrain_point(const Point_3& pt) {
    m_groundPtsPtr->insert(pt); //todo sort out terrain pts
}

const std::vector<roofer::Mesh>& ReconstructedBuilding::get_roofer_meshes() const {
    return m_roofer_meshes;
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

void ReconstructedBuilding::get_cityjson_semantics(nlohmann::json& g) const { // Temp for checking CGAL mesh properties
    // for now handle only LoD1.2 with semantics, LoD1.3 and LoD2.2 loses information
    // when making final repairs
    if (m_reconSettings->lod == "1.2") {
        Face_property semantics;
        auto semanticsMap = m_mesh.property_map<face_descriptor, std::string>("f:semantics");
        if (semanticsMap.has_value()) {
            semantics = semanticsMap.value();
        } else throw city4cfd_error("Semantic property map not found!");

        std::unordered_map<std::string, int> surfaceId;
        surfaceId["RoofSurface"] = 0;   g["surfaces"][0]["type"] = "RoofSurface";
        surfaceId["GroundSurface"] = 1; g["surfaces"][1]["type"] = "GroundSurface";
        surfaceId["WallSurface"] = 2;   g["surfaces"][2]["type"] = "WallSurface";

        for (auto faceIdx: m_mesh.faces()) {
            auto it = surfaceId.find(semantics[faceIdx]);
            if (it == surfaceId.end()) throw city4cfd_error("Could not find semantic attribute!");

            g["values"][faceIdx.idx()] = it->second;
        }
    }
}

void ReconstructedBuilding::reconstruct_from_attribute() {
    //-- Check attribute height
    if (m_attributeHeight <= 0) {
        this->mark_as_failed();
        this->deactivate();
        throw city4cfd_error("Attribute height from geojson file is invalid!");
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