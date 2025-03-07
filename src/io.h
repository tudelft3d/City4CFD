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

#ifndef CITY4CFD_IO_H
#define CITY4CFD_IO_H

#include "types.h"
#include "CGALTypes.h"

namespace IO {
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