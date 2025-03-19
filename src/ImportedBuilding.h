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

#ifndef CITY4CFD_EXPLICITBUILDING_H
#define CITY4CFD_EXPLICITBUILDING_H

#include "Building.h"

#include <CGAL/Polygon_set_2.h>

class ImportedBuilding : public Building {
public:
    static int noBottom;

    ImportedBuilding() = delete;
    ImportedBuilding(std::unique_ptr<nlohmann::json>& buildingJson,
                     PointSet3Ptr& importedBuildingPts);
    ImportedBuilding(Mesh& mesh);
    ~ImportedBuilding() = default;

    virtual void   calc_elevation() override;
    virtual void   reconstruct() override;
    virtual void   reconstruct_flat_terrain() override;
    virtual void   insert_terrain_point(const Point_3& /* pt */) override;

    void   append_nonground_part(const std::shared_ptr<ImportedBuilding>& other);

    const nlohmann::json& get_building_json() const;
    const int             get_lod_idx() const;
    const bool            is_appending() const;

protected:
    std::unordered_map<int, Point_3>  m_ptMap;
    std::unique_ptr<nlohmann::json>   m_buildingJson;
    std::vector<int>                  m_footprintIdxList;
    std::vector<std::vector<int>>     m_footprintPtsIdxList;
    bool                              m_appendToBuilding;
    bool                              m_trueHeight;
    int                               m_lodIdx;

    void check_simplicity(Polygon_2& ring);
    void polyset_to_polygon(const CGAL::Polygon_set_2<CGAL::Epeck>& polySet);
    void set_footprint_mesh_connectivity(const std::unordered_map<std::string, int>& pointConnectivity);
};

#endif //CITY4CFD_EXPLICITBUILDING_H