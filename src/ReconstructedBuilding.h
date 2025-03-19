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

#ifndef CITY4CFD_RECONSTRUCTEDBUILDING_H
#define CITY4CFD_RECONSTRUCTEDBUILDING_H

#include "Building.h"
#include "roofer.h"

class ReconstructedBuilding : public Building {
public:
    ReconstructedBuilding();
    ReconstructedBuilding(const Mesh& mesh);
    ReconstructedBuilding(const roofer::Mesh& rooferMesh, const ReconstructedBuildingPtr& other);
    ReconstructedBuilding(const nlohmann::json& poly);
    ReconstructedBuilding(const Polygon_with_attr& poly);
    ReconstructedBuilding(const std::shared_ptr<ImportedBuilding>& importedBuilding);
    ~ReconstructedBuilding() = default;

    const std::vector<roofer::Mesh>& get_roofer_meshes() const;

    virtual void   calc_elevation() override;
    virtual void   reconstruct() override;
    virtual void   insert_terrain_point(const Point_3& pt) override;
    virtual void   reconstruct_flat_terrain() override;
    virtual void   get_cityjson_semantics(nlohmann::json& g) const override;

protected:
    double m_attributeHeight;
    bool   m_attributeHeightAdvantage;
    PointSet3Ptr m_groundPtsPtr;
    std::vector<roofer::Mesh> m_roofer_meshes;

    void reconstruct_from_attribute();
    bool reconstruct_again_from_attribute(const std::string& reason);
    void reconstruct_lod12();
};


#endif //CITY4CFD_RECONSTRUCTEDBUILDING_H