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

#include "Top.h"

#include "geomutils.h"

Top::Top(const int outputLayerID)
        : Boundary(outputLayerID) {}

void Top::reconstruct() {
    std::vector<Mesh::vertex_index> meshVertexTop;

    //-- Top is done by making a CDT of outerPts
    CDT cdt_top;
    for (auto& pt : s_outerPts) {
        cdt_top.insert(ePoint_3(pt.x(), pt.y(), Config::get().topHeight));
    }
    //-- Add mesh faces for top
    geomutils::cdt_to_mesh(cdt_top, m_mesh);
}

TopoClass Top::get_class() const {
    return TOP;
}

std::string Top::get_class_name() const {
    return "Top";
}