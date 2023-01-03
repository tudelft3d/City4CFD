/*
  City4CFD
 
  Copyright (c) 2021-2023, 3D Geoinformation Research Group, TU Delft

  This file is part of City4CFD.

  City4CFD is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  City4CFD is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with City4CFD.  If not, see <http://www.gnu.org/licenses/>.

  For any information or further details about the use of City4CFD, contact
  Ivan Pađen
  <i.paden@tudelft.nl>
  3D Geoinformation Research Group
  Delft University of Technology
*/

#ifndef CITY4CFD_POLYFEATURE_H
#define CITY4CFD_POLYFEATURE_H

#include "TopoFeature.h"

class PolyFeature : public TopoFeature {
public:
    PolyFeature();
    PolyFeature(const int outputLayerID);
    PolyFeature(const nlohmann::json& poly, const bool checkSimplicity = false);
    PolyFeature(const int outputLayerID, const int internalID);
    PolyFeature(const nlohmann::json& poly, const bool checkSimplicity, const int outputLayerID);
    PolyFeature(const nlohmann::json& poly, const bool checkSimplicity, const int outputLayerID, const int internalID);
    PolyFeature(const nlohmann::json& poly, const int outputLayerID);
    PolyFeature(const nlohmann::json& poly, const int outputLayerID, const int internalID);
    virtual ~PolyFeature();

    void  calc_footprint_elevation_nni(const DT& dt);
#ifndef NDEBUG
    void  calc_footprint_elevation_linear(const DT& dt);
#endif
    double ground_elevation();
    double slope_height();
    bool   flatten_polygon_inner_points(const Point_set_3& pointCloud, std::map<int, Point_3>& flattenedPts,
                                       const SearchTree& searchTree, const std::unordered_map<Point_3,
                                       int>& pointCloudConnectivity) const;
    void  set_zero_borders();
    void  calc_min_bbox();
    void  clear_feature();

    virtual void        get_cityjson_info(nlohmann::json& b) const = 0;
    virtual void        get_cityjson_semantics(nlohmann::json& g) const = 0;
    virtual std::string get_cityjson_primitive() const = 0;
    virtual TopoClass   get_class() const = 0;
    virtual std::string get_class_name() const = 0;

    Polygon_with_holes_2&                    get_poly();
    const Polygon_with_holes_2&              get_poly() const;
    const std::vector<std::vector<double>>&  get_ground_elevations() const;
    const int                                get_internal_id() const;
    MinBbox&                                 get_min_bbox();

protected:
    int                               _polyInternalID;
    Polygon_with_holes_2              _poly;
    std::vector<std::vector<double>>  _groundElevations;
    double                            _groundElevation;
    MinBbox                           _minBbox;

    void parse_json_poly(const nlohmann::json& poly, const bool checkSimplicity);
};

#endif //CITY4CFD_POLYFEATURE_H