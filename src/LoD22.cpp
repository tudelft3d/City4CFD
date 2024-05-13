/*
  City4CFD

  Copyright (c) 2021-2024, 3D Geoinformation Research Group, TU Delft

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

#include "LoD22.h"
#include "Config.h"
#include "misc/cgal_utils.hpp"

LoD22::LoD22(roofer::Mesh roofer_mesh) {
    // shorten long poly edges
    double sq_maxdist = Config::get().edgeMaxLen * Config::get().edgeMaxLen;
    this->shorten_mesh_edges(roofer_mesh, sq_maxdist);

    // sort out the footprint back
    this->get_footprint_from_mesh(roofer_mesh, m_footprint, m_baseElevations);

    // sort out the mesh
    m_mesh = roofer::Mesh2CGALSurfaceMesh<Point_3>(roofer_mesh);
}

void LoD22::reconstruct(const PointSet3Ptr& buildingPtsPtr,
                        const PointSet3Ptr& groundPtsPtr,
                        const Polygon_with_holes_2& footprint,
                        const std::vector<std::vector<double>>& base_elevations,
                        ReconstructionConfig config
                        ) {

    // prep building pts
    roofer::PointCollection buildingPts;
    for (const auto& pt : buildingPtsPtr->points()) {
        buildingPts.push_back({(float)pt.x(),
                               (float)pt.y(),
                               (float)pt.z()});
    }

    // prep ground pts
    roofer::PointCollection groundPts;
    for (const auto& pt : groundPtsPtr->points()) {
        groundPts.push_back({(float)pt.x(),
                             (float)pt.y(),
                             (float)pt.z()});
    }

    // prep footprints
    //todo  base_elevations are in for a rewrite
    roofer::LinearRing linearRing;
    int i = 0;
    for (auto& p : footprint.outer_boundary()) {
        float x = p.x();
        float y = p.y();
        float z = base_elevations.front()[i++];
        linearRing.push_back({x, y, z});
    }
    int j = 1;
    for (auto it = footprint.holes_begin(); it != footprint.holes_end(); ++it) {
        roofer::LinearRing iring;
        i = 0;
        for (auto& p : *it) {
            float x = p.x();
            float y = p.y();
            float z = base_elevations[j][i++];
            iring.push_back({x, y, z});
        }
        linearRing.interior_rings().push_back(iring);
        ++j;
    }

    //todo groundPts flag
    // reconstruct
//    m_roofer_meshes = roofer::reconstruct_single_instance(buildingPts, groundPts, linearRing,
//                                                         {.lambda = config.m_lambda,
//                                                          .lod = config.m_lod,
//                                                          .lod13_step_height = config.m_lod13_step_height});
    m_roofer_meshes = roofer::reconstruct_single_instance(buildingPts, linearRing,
                                                         {.lambda = config.m_lambda,
                                                          .lod = config.m_lod,
                                                          .lod13_step_height = config.m_lod13_step_height});

    auto roofer_mesh = m_roofer_meshes.front();

    // shorten long poly edges
    double sq_maxdist = Config::get().edgeMaxLen * Config::get().edgeMaxLen;
    this->shorten_mesh_edges(roofer_mesh, sq_maxdist);

    // sort out the footprint back
    this->get_footprint_from_mesh(roofer_mesh, m_footprint, m_baseElevations);

    // sort out the mesh
    m_mesh = roofer::Mesh2CGALSurfaceMesh<Point_3>(roofer_mesh);
}

void LoD22::shorten_mesh_edges(roofer::Mesh& roofer_mesh, const double sq_maxdist) const{
    for (auto& poly: roofer_mesh.get_polygons()) {
        int i = 0;
        while (i != poly.size()) {
            Point_3 pt1 = Point_3(poly[i][0], poly[i][1], poly[i][2]);
            Point_3 pt2 = Point_3(poly[(i + 1) % poly.size()][0], poly[(i + 1) % poly.size()][1],
                                  poly[(i + 1) % poly.size()][2]);
            auto sq_dist = CGAL::squared_distance(pt1, pt2);
            if (sq_dist > sq_maxdist) {
                auto midpt = CGAL::midpoint(pt1, pt2);
                std::array<float, 3> newpt{(float)midpt.x(), (float)midpt.y(), (float)midpt.z()};
                poly.insert(poly.begin() + i + 1, newpt);
            } else ++i;
        }
        for (auto& hole: poly.interior_rings()) {
            i = 0;
            while (i != hole.size()) {
                Point_3 pt1 = Point_3(hole[i][0], hole[i][1], hole[i][2]);
                Point_3 pt2 = Point_3(hole[(i + 1) % hole.size()][0], hole[(i + 1) % hole.size()][1],
                                      hole[(i + 1) % hole.size()][2]);
                auto sq_dist = CGAL::squared_distance(pt1, pt2);
                if (sq_dist > sq_maxdist) {
                    auto midpt = CGAL::midpoint(pt1, pt2);
                    std::array<float, 3> newpt{(float)midpt.x(), (float)midpt.y(), (float)midpt.z()};
                    hole.insert(hole.begin() + i + 1, newpt);
                } else ++i;
            }
        }
    }
}

void LoD22::get_footprint_from_mesh(const roofer::Mesh& roofer_mesh, Polygon_with_holes_2& footprint,
                                    std::vector<std::vector<double>>& baseElevations) {
    footprint.rings().clear();
    baseElevations.clear();
    auto labels = roofer_mesh.get_labels();
    roofer::LinearRing roofer_new_footprint;
    for (int i = 0; i != labels.size(); ++i) {
        if (labels[i] == 0) {
            roofer_new_footprint = roofer_mesh.get_polygons()[i];
            break;
        }
        throw city4cfd_error("No footprint found in the reconstructed mesh!");
    }
    Polygon_2 poly2;
    std::vector<double> outer_base_elevations;
    for (auto& p: roofer_new_footprint) {
        poly2.push_back(Point_2(p[0], p[1]));
        outer_base_elevations.push_back(p[2]);
    }
    baseElevations.push_back(outer_base_elevations);
    std::vector<Polygon_2> holes;
    for (auto& lr_hole: roofer_new_footprint.interior_rings()) {
        Polygon_2 hole;
        std::vector<double> hole_elevations;
        for (auto& p: lr_hole) {
            hole.push_back(Point_2(p[0], p[1]));
            hole_elevations.push_back(p[2]);
        }
        holes.push_back(hole);
        baseElevations.push_back(hole_elevations);
        hole_elevations.clear();
    }
    CGAL::Polygon_with_holes_2<EPICK> new_footprint = CGAL::Polygon_with_holes_2<EPICK>(poly2, holes.begin(),
                                                                                        holes.end());
    footprint = Polygon_with_holes_2(new_footprint);
}

Polygon_with_holes_2 LoD22::get_footprint() const {
    if (m_footprint.rings().empty()) throw city4cfd_error("No footprint found!");
    return m_footprint;
}

std::vector<std::vector<double>> LoD22::get_base_elevations() const {
    if (m_baseElevations.empty()) throw city4cfd_error("No base elevations found!");
    return m_baseElevations;
}

Mesh LoD22::get_mesh() const {
    if (m_mesh.is_empty()) throw city4cfd_error("No mesh found!");
    return m_mesh;
}