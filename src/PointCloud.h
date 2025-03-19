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

#ifndef CITY4CFD_POINTCLOUD_H
#define CITY4CFD_POINTCLOUD_H

#include "types.h"
#include "CGALTypes.h"

#include "lasreader.hpp"
#include "CSF/src/CSF.h"

class PointCloud {
public:
    PointCloud() = default;
    ~PointCloud() = default;

    void random_thin_pts();
    void create_flat_terrain(const PolyFeaturesPtr& lsFeatures);
    void set_flat_terrain();
    void smooth_terrain();
    void terrain_points_in_polygon(BuildingsPtr& features);
    void flatten_polygon_pts(const PolyFeaturesPtr& lsFeatures, std::vector<EPECK::Segment_3>& constrainedEdges,
                             std::vector<std::pair<Polygon_with_holes_2, int>>& newPolys);
    void buffer_flat_edges(const PolyFeaturesPtr& avgFeatures, std::vector<EPECK::Segment_3>& constrainedEdges);
    void read_point_clouds();

    Point_set_3& get_terrain();
    Point_set_3& get_buildings();
    const Point_set_3& get_terrain() const;
    const Point_set_3& get_buildings() const;

private:
    Point_set_3 m_pointCloudTerrain;
    Point_set_3 m_pointCloudBuildings;
};


#endif //CITY4CFD_POINTCLOUD_H
