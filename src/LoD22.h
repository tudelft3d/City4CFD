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
#ifndef CITY4CFD_LOD22_H
#define CITY4CFD_LOD22_H

#include "types.h"
#include "CGALTypes.h"
#include "roofer.h"

class LoD22 {
public:
    LoD22() = default;
    LoD22(roofer::Mesh rooferMesh);
    ~LoD22() = default;

    struct ReconstructionConfig {
        // control optimisation
        float m_lambda;
        // requested LoD
        int m_lod;
        // step height used for LoD13 generalisation
        float m_lod13_step_height;

        // Constructor with default values
        ReconstructionConfig() : m_lambda(1./9.), m_lod(22), m_lod13_step_height(3.) {}

        ReconstructionConfig& lambda(float lambda) { m_lambda = lambda; return *this; };
        ReconstructionConfig& lod(int lod) { m_lod = lod; return *this; };
        ReconstructionConfig& lod13_step_height(float lod13_step_height) { m_lod13_step_height = lod13_step_height; return *this; };
    };

    void reconstruct(const PointSet3Ptr& buildingPtsPtr,
                     const PointSet3Ptr& groundPtsPtr,
                     const Polygon_with_holes_2& footprint,
                     const std::vector<std::vector<double>>& baseElevations,
                     ReconstructionConfig config = ReconstructionConfig()
    );

    Polygon_with_holes_2              get_footprint() const;
    std::vector<std::vector<double>>  get_base_elevations() const;
    Mesh                              get_mesh() const;
    std::vector<roofer::Mesh>         get_roofer_meshes() const;


private:
    Mesh                                     m_mesh;
    std::vector<roofer::Mesh>                m_rooferMeshes;
    Polygon_with_holes_2                     m_footprint;
    std::vector<std::vector<double>>         m_baseElevations;

    void shorten_mesh_edges(roofer::Mesh& mesh, const double sq_maxdist) const;
    void get_footprint_from_mesh(const roofer::Mesh& rooferMesh, Polygon_with_holes_2& footprint, std::vector<std::vector<double>>& baseElevations) const;
};


#endif //CITY4CFD_LOD22_H
