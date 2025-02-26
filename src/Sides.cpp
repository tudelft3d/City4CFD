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

#include "Sides.h"

#include "Config.h"

Sides::Sides(const int outputLayerID)
        : Boundary(outputLayerID) {
}

void Sides::reconstruct() {
    std::vector<Mesh::vertex_index> meshVertexSide;

    //-- Add mesh vertices and store them in a vector
    for (auto it = m_sideOutputPts.begin(); it != m_sideOutputPts.end(); ++it) {
        meshVertexSide.emplace_back(m_mesh.add_vertex(*it));
        meshVertexSide.emplace_back(m_mesh.add_vertex(Point_3(it->x(), it->y(), Config::get().topHeight)));
    }

    //-- Add middle top point to mesh
//    Mesh::vertex_index topMiddlePoint = _meshTop.add_vertex(Point_3(bndInfo.xcent, bndInfo.ycent, bndInfo.height));

    //-- Add mesh faces for side
    for (auto i = 0; i < meshVertexSide.size() - 3; i= i + 2) {
        // -- i + 1 is i lifted up
        int v1 = i;
        int v2 = i + 2;

        m_mesh.add_face(meshVertexSide[v2], meshVertexSide[v1], meshVertexSide[v1 + 1]);
        m_mesh.add_face(meshVertexSide[v2 + 1], meshVertexSide[v2], meshVertexSide[v1 + 1]);
    }
}

TopoClass Sides::get_class() const {
    return SIDES;
}

std::string Sides::get_class_name() const {
    return "Sides";
}