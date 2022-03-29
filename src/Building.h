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

#ifndef CITY4CFD_BUILDING_H
#define CITY4CFD_BUILDING_H

#include "PolyFeature.h"

class Building : public PolyFeature {
public:
    Building();
    Building(const int internalID);
    Building(const nlohmann::json& poly);
    Building(const nlohmann::json& poly, const int internalID);
    ~Building();

    virtual void reconstruct() = 0;

    void   check_feature_scope(const Polygon_2& influRegion);
    double max_dim();

    double get_height() const;

    virtual void        get_cityjson_info(nlohmann::json& b) const override;
    virtual void        get_cityjson_semantics(nlohmann::json& g) const override;
    virtual std::string get_cityjson_primitive() const override;
    virtual TopoClass   get_class() const override;
    virtual std::string get_class_name() const override;

protected:
    double _height;
};

#endif //CITY4CFD_BUILDING_H