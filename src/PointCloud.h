/*
  Copyright (c) 2021-2022,
  Ivan Pađen <i.paden@tudelft.nl>
  3D Geoinformation,
  Delft University of Technology

  This file is part of City4CFD.

  City4CFD is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  City4CFD is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

  You should have received a copy of the GNU General Public License
  along with this program. If not, see <http://www.gnu.org/licenses/>
*/

#ifndef CITY4CFD_POINTCLOUD_H
#define CITY4CFD_POINTCLOUD_H

#include "types.h"
#include "CGALTypes.h"

typedef std::shared_ptr<SearchTree> SearchTreePtr;

class PointCloud {
public:
    PointCloud();
    ~PointCloud();

    void random_thin_pts();
    void average_polygon_pts(const PolyFeatures& lsFeatures);
    SearchTreePtr make_search_tree_buildings();
    void read_point_clouds();

    Point_set_3& get_terrain();
    Point_set_3& get_buildings();
    const Point_set_3& get_terrain() const;
    const Point_set_3& get_buildings() const;

private:
    Point_set_3 _pointCloudTerrain;
    Point_set_3 _pointCloudBuildings;
};


#endif //CITY4CFD_POINTCLOUD_H
