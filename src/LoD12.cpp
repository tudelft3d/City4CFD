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

#include "LoD12.h"

#include "geomutils.h"

LoD12::LoD12(const Polygon_with_holes_2& poly,
             const std::vector<std::vector<double>>& baseElevations)
        : m_elevation(), m_poly(poly), m_baseElevations(baseElevations) {}

LoD12::LoD12(const Polygon_with_holes_2& poly,
             const std::vector<std::vector<double>>& baseElevations,
             const double elevation)
        : m_elevation(elevation), m_poly(poly), m_baseElevations(baseElevations) {}

void LoD12::set_elevation(const double& elevation) {
    m_elevation = elevation;
}

void LoD12::reconstruct(Mesh& mesh) {
    mesh.clear();
    // Add semantics with face properties to the property map
    auto surfaceType = mesh.add_property_map<face_descriptor , std::string>("f:semantics", "").first;
    face_descriptor fIdx;

    CDT cdtBuildings;

    //-- Map CDT and Mesh vertices
    std::map<CDT::Vertex_handle, Mesh::Vertex_index> cdtToMesh;

    int polyCount = 0;
    for (auto& poly : m_poly.rings()) { // Loop over polys
        std::vector<Vertex_handle> cdtHandle;
        std::vector<Mesh::Vertex_index> meshVertex;
        int count = 0;
        for (auto vert = poly.vertices_begin(); vert != poly.vertices_end(); ++vert) { // Loop over poly vertices
            cdtHandle.emplace_back(cdtBuildings.insert(ePoint_3(vert->x(),
                                                                vert->y(),
                                                                m_baseElevations[polyCount][count])));
            meshVertex.emplace_back(mesh.add_vertex(Point_3(vert->x(),
                                                            vert->y(),
                                                            m_baseElevations[polyCount][count++])));
            cdtToMesh[cdtHandle.back()] = meshVertex.back();
            meshVertex.emplace_back(mesh.add_vertex(Point_3(vert->x(),
                                                            vert->y(),
                                                            m_elevation)));
        }
        cdtHandle.emplace_back(cdtHandle.front());
        meshVertex.emplace_back(meshVertex.front());
        meshVertex.emplace_back(meshVertex.front() + 1);

        //- Add constraints and create mesh faces for sides
        for (auto i = 0; i < cdtHandle.size() - 1; ++i) {
            cdtBuildings.insert_constraint(cdtHandle[i], cdtHandle[i + 1]);

            auto it1 = mesh.vertices_begin();
            auto it2 = mesh.vertices_begin();

            auto v1 = cdtToMesh[cdtHandle[i]];
            auto v2 = cdtToMesh[cdtHandle[i + 1]];

            std::advance(it1, v1.idx() + 1);
            std::advance(it2, v2.idx() + 1);

            fIdx = mesh.add_face(v1, v2, *it1);
            if (fIdx != Mesh::null_face()) {
                surfaceType[fIdx] = "WallSurface";
            }
            fIdx = mesh.add_face(v2, *it2, *it1);
            if (fIdx != Mesh::null_face()) {
                surfaceType[fIdx] = "WallSurface";
            }
        }
        ++polyCount;
    }
    //- Handle top and bottom
    geomutils::mark_domains(cdtBuildings);
    for (auto& it : cdtBuildings.finite_face_handles()) {
        if (!it->info().in_domain()) continue;

        auto it1 = mesh.vertices_begin();
        auto it2 = mesh.vertices_begin();
        auto it3 = mesh.vertices_begin();

        std::advance(it1, cdtToMesh[it->vertex(0)]);
        std::advance(it2, cdtToMesh[it->vertex(1)]);
        std::advance(it3, cdtToMesh[it->vertex(2)]);

        fIdx = mesh.add_face(*it1, *it3, *it2); // Bottom face
        if (fIdx != Mesh::null_face()) {
            surfaceType[fIdx] = "GroundSurface";
        }
        fIdx = mesh.add_face(*std::next(it1), *std::next(it2), *std::next(it3));
        if (fIdx != Mesh::null_face()) {
            surfaceType[fIdx] = "RoofSurface";
        }
    }
}