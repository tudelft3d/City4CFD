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

#include "LoD22.h"
#include "Config.h"
#include "misc/cgal_utils.hpp"

LoD22::LoD22(roofer::Mesh rooferMesh) {
    // shorten long poly edges
    const double sq_maxdist = Config::get().edgeMaxLen * Config::get().edgeMaxLen;
    this->shorten_mesh_edges(rooferMesh, sq_maxdist);

    // sort out the footprint back
    this->get_footprint_from_mesh(rooferMesh, m_footprint, m_baseElevations);

    // sort out the mesh
    m_mesh = roofer::Mesh2CGALSurfaceMesh<Point_3>(rooferMesh);
}

void LoD22::reconstruct(const PointSet3Ptr& buildingPtsPtr,
                        const PointSet3Ptr& groundPtsPtr,
                        const Polygon_with_holes_2& footprint,
                        const std::vector<std::vector<double>>& baseElevations,
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
    //todo  baseElevations are in for a rewrite
    roofer::LinearRing linearRing;
    int i = 0;
    for (auto& p : footprint.outer_boundary()) {
        float x = p.x();
        float y = p.y();
        float z = baseElevations.front()[i++];
        linearRing.push_back({x, y, z});
    }
    int j = 1;
    for (auto it = footprint.holes_begin(); it != footprint.holes_end(); ++it) {
        roofer::LinearRing iring;
        i = 0;
        for (auto& p : *it) {
            float x = p.x();
            float y = p.y();
            float z = baseElevations[j][i++];
            iring.push_back({x, y, z});
        }
        linearRing.interior_rings().push_back(iring);
        ++j;
    }

    // reconstruct
    m_rooferMeshes = roofer::reconstruct_single_instance(buildingPts, linearRing,
//                                                       groundPts, //todo groundPts flag
                                                         {.lambda = config.m_lambda,
                                                          .lod = config.m_lod,
                                                          .lod13_step_height = config.m_lod13_step_height});

    // store the first mesh from the vector, rest should be handled separately
    auto rooferMesh = m_rooferMeshes.front();

    // shorten long poly edges
    const double sq_maxdist = Config::get().edgeMaxLen * Config::get().edgeMaxLen;
    this->shorten_mesh_edges(rooferMesh, sq_maxdist);

    // sort out the footprint back
    this->get_footprint_from_mesh(rooferMesh, m_footprint, m_baseElevations);

    // sort out the mesh
    m_mesh = roofer::Mesh2CGALSurfaceMesh<Point_3>(rooferMesh);
}

void LoD22::shorten_mesh_edges(roofer::Mesh& roofer_mesh, const double sq_maxdist) const {
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

void LoD22::get_footprint_from_mesh(const roofer::Mesh& rooferMesh, Polygon_with_holes_2& footprint,
                                    std::vector<std::vector<double>>& baseElevations) const {
    footprint.rings().clear();
    baseElevations.clear();
    auto labels = rooferMesh.get_labels();
    roofer::LinearRing rooferNewFootprint;
    for (int i = 0; i != labels.size(); ++i) {
        if (labels[i] == 0) {
            rooferNewFootprint = rooferMesh.get_polygons()[i];
            break;
        }
        throw city4cfd_error("No footprint found in the reconstructed mesh!");
    }
    Polygon_2 poly2;
    std::vector<double> outerBaseElevations;
    for (auto& p: rooferNewFootprint) {
        poly2.push_back(Point_2(p[0], p[1]));
        outerBaseElevations.push_back(p[2]);
    }
    baseElevations.push_back(outerBaseElevations);
    std::vector<Polygon_2> holes;
    for (auto& lrHole: rooferNewFootprint.interior_rings()) {
        Polygon_2 hole;
        std::vector<double> holeElevations;
        for (auto& p: lrHole) {
            hole.push_back(Point_2(p[0], p[1]));
            holeElevations.push_back(p[2]);
        }
        holes.push_back(hole);
        baseElevations.push_back(holeElevations);
        holeElevations.clear();
    }
    CGAL::Polygon_with_holes_2<EPICK> newFootprint = CGAL::Polygon_with_holes_2<EPICK>(poly2, holes.begin(),
                                                                                       holes.end());
    footprint = Polygon_with_holes_2(newFootprint);
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
