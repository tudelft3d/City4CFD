/*
  City4CFD
 
  Copyright (c) 2021-2026, 3D Geoinformation Research Group, TU Delft

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
*/


#ifndef CITY4CFD_IO_H
#define CITY4CFD_IO_H

#include "types.h"
#include "CGALTypes.h"
#include "BuildingFootprintFilter.h"

#include <set>

namespace IO {

    struct PointCloudReadOptions {
        std::set<int> terrain_classes  = {2};   // ASPRS codes routed to terrain bucket
        std::set<int> building_classes = {6};   // ASPRS codes routed to buildings bucket
        BuildingFootprintFilter* filter  = nullptr; // nullable; non-const for bucket accumulation
        std::size_t keep_every_nth = 0; // 0 = disabled; keep 1 in N terrain points, drop the rest
        std::size_t drop_every_nth = 0; // 0 = disabled; drop 1 in N terrain points, keep the rest
    };

    // Streaming LAS/LAZ reader with classification split.
    // Terrain points are accumulated into terrain.
    // Building points are routed into per-footprint buckets inside opts.filter (must be non-null
    // when building_classes is non-empty; stray points outside every footprint are dropped).
    void read_and_split_point_clouds(const std::vector<std::string>& files,
                                     Point_set_3& terrain,
                                     const PointCloudReadOptions& opts);

    // High-level point cloud reader: new unified path (point_cloud_files) or legacy two-file
    // path (ground_xyz + buildings_xyz).  Populates terrain and buildings in place.
    void read_point_clouds(Point_set_3& terrain, Point_set_3& buildings,
                           BuildingFootprintFilter& filter);

    //-- Input functions
    void read_config(std::string& config_path);
    bool read_point_cloud(std::string& file, Point_set_3& pc);
    void read_geojson_polygons(std::string& file, JsonVectorPtr& jsonPolygons);
    void read_polygons(std::string& file, PolyVecPtr& polygons, std::string* crsInfo);
    void read_cityjson_geometries(std::string& file, JsonVectorPtr& importedBuildings,
                                  PointSet3Ptr& importedBuildingPts);
    void read_other_geometries(std::string& file, std::vector<Mesh>& meshes);

    //-- Output functions
    void print_progress_bar(int percent);

    void output_obj(const OutputFeaturesPtr& allFeatures);
    void output_stl(const OutputFeaturesPtr& allFeatures);
    void output_cityjson(const OutputFeaturesPtr& allFeatures);

    void get_obj_pts(const Mesh& mesh, std::string& fs, std::string& bs, std::unordered_map<std::string, int>& dPts);
    void get_stl_pts(Mesh& mesh, std::string& fs);
    void get_cityjson_geom(const Mesh& mesh,
                           nlohmann::json& g,
                           std::unordered_map<std::string,int>& dPts,
                           std::vector<double>& minMaxZ);

    bool not_same(std::vector<int> idxLst);
    bool is_degen(const Mesh& mesh, Mesh::Face_index face);

    void output_log();
    void output_log(const BuildingsPtr& failedBuildings);
    bool has_substr(const std::string& strMain, const std::string& subStr);

    std::string gen_key_bucket(const Point_2 p);

    //-- Templated function
    template <typename T> std::string gen_key_bucket(const T& p);
    template <typename T> std::string gen_key_bucket_int(const T& p, const double inverseScale, const double inverseTranslate);
}

#endif //CITY4CFD_IO_H