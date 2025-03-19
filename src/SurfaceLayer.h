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

#ifndef CITY4CFD_SURFACELAYER_H
#define CITY4CFD_SURFACELAYER_H

#include "PolyFeature.h"

class SurfaceLayer : public PolyFeature {
public:
    SurfaceLayer(const int outputLayerID);
    SurfaceLayer(const nlohmann::json& poly, const int outputLayerID);
    SurfaceLayer(const Polygon_with_attr& poly, const int outputLayerID);
    ~SurfaceLayer() = default;

    void check_feature_scope(const Polygon_2& bndPoly);

    virtual void        get_cityjson_cityobj_info(nlohmann::json& f) const override;
    virtual void        get_cityjson_geomobj_info(nlohmann::json& g) const override;
    virtual TopoClass   get_class() const override;
    virtual std::string get_class_name() const override;

private:

};

#endif //CITY4CFD_SURFACELAYER_H