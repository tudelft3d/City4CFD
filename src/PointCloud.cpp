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

#include "PointCloud.h"

#include "io.h"
#include "geomutils.h"
#include "Config.h"
#include "PolyFeature.h"

PointCloud::PointCloud()  = default;
PointCloud::~PointCloud() = default;

void PointCloud::random_thin_pts() {
    if (Config::get().terrainThinning > 0 + g_smallnum) {
        std::cout <<"\nRandomly thinning terrain points" << std::endl;
        _pointCloudTerrain.remove(CGAL::random_simplify_point_set(_pointCloudTerrain,
                                                                  Config::get().terrainThinning),
                                  _pointCloudTerrain.end());
        _pointCloudTerrain.collect_garbage();
        std::cout << "    Terrain points after thinning: " << _pointCloudTerrain.size() << std::endl;
    }
}

void PointCloud::smooth_terrain() {
    std::cout << "\nSmoothing terrain" << std::endl;
    DT dt(_pointCloudTerrain.points().begin(), _pointCloudTerrain.points().end());
    geomutils::smooth_dt<DT, EPICK>(_pointCloudTerrain, dt);

    //-- Return new pts to the point cloud
    _pointCloudTerrain.clear();
    for (auto& pt : dt.points()) _pointCloudTerrain.insert(pt);
    _pointCloudTerrain.add_property_map<bool> ("is_building_point", false);
}

void PointCloud::average_polygon_pts(const PolyFeatures& lsFeatures) {
    std::cout << "\n    Averaging surfaces" << std::endl;
    std::map<int, Point_3> averagedPts;

    //-- Construct a connectivity map and remove duplicates along the way
    auto is_building_pt = _pointCloudTerrain.property_map<bool>("is_building_point").first;
    std::unordered_map<Point_3, int> pointCloudConnectivity;
    auto it = _pointCloudTerrain.points().begin();
    int count = 0;
    while (it != _pointCloudTerrain.points().end()) {
        auto itPC = pointCloudConnectivity.find(*it);
        if (itPC != pointCloudConnectivity.end()) {
            _pointCloudTerrain.remove(_pointCloudTerrain.begin() + count);
        } else {
            pointCloudConnectivity[*it] = count;
            ++it;
            ++count;
        }
    }
    _pointCloudTerrain.collect_garbage();

    //-- Construct search tree from ground points
    SearchTree searchTree(_pointCloudTerrain.points().begin(), _pointCloudTerrain.points().end());

    //-- Perform averaging
    for (auto& f : lsFeatures) {
        auto it = Config::get().averageSurfaces.find(f->get_output_layer_id());
        if (it != Config::get().averageSurfaces.end()) {
            f->average_polygon_inner_points(_pointCloudTerrain, averagedPts, searchTree, pointCloudConnectivity);
        }
    }

    //-- Change points with averaged values
    int pcOrigSize = _pointCloudTerrain.points().size();
    for (auto& it : averagedPts) {
        _pointCloudTerrain.insert(it.second);
    }
    for (int i = 0; i < pcOrigSize; ++i) {
        auto it = averagedPts.find(i);
        if (it != averagedPts.end()) {
            _pointCloudTerrain.remove(i);
            averagedPts.erase(i);
        }
    }
    _pointCloudTerrain.collect_garbage();
}

SearchTreePtr PointCloud::make_search_tree_buildings() {
    return std::make_shared<SearchTree>(_pointCloudBuildings.points().begin(),
                                        _pointCloudBuildings.points().end());
}

void PointCloud::read_point_clouds() {
    //-- Read ground points
    if (!Config::get().points_xyz.empty()) {
        std::cout << "Reading ground points" << std::endl;
        IO::read_point_cloud(Config::get().points_xyz, _pointCloudTerrain);
        _pointCloudTerrain.add_property_map<bool> ("is_building_point", false);

        std::cout << "    Points read: " << _pointCloudTerrain.size() << std::endl;
    } else {
        std::cout << "INFO: Did not find any ground points! Calculating ground as flat surface\n" << std::endl;
        //todo needs to be implemented - handled with the schema for now
    }

    //-- Read building points
    if (!Config::get().buildings_xyz.empty()) {
        std::cout << "Reading building points" << std::endl;
        IO::read_point_cloud(Config::get().buildings_xyz, _pointCloudBuildings);
        if (_pointCloudBuildings.empty()) throw std::invalid_argument("Didn't find any building points!");

        std::cout << "    Points read: " << _pointCloudBuildings.size() << std::endl;
    }
}

Point_set_3& PointCloud::get_terrain() {
    return _pointCloudTerrain;
}
Point_set_3& PointCloud::get_buildings() {
    return _pointCloudBuildings;
}

const Point_set_3& PointCloud::get_terrain() const {
    return _pointCloudTerrain;
}

const Point_set_3& PointCloud::get_buildings() const {
    return _pointCloudBuildings;
}