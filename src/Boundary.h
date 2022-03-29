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

#ifndef CITY4CFD_BOUNDARY_H
#define CITY4CFD_BOUNDARY_H

#include "TopoFeature.h"

class Boundary : public TopoFeature {
public:
    Boundary();
    Boundary(const int outputLayerID);
    virtual ~Boundary();

    static void set_bnd_poly(Polygon_2& bndPoly, Polygon_2& pcBndPoly, Polygon_2& startBufferPoly);
    static void set_bounds_to_pc(Point_set_3& pointCloud, const Polygon_2& pcBndPoly);
    static void set_bounds_to_terrain(Point_set_3& pointCloud, const Polygon_2& bndPoly,
                                      const Polygon_2& pcBndPoly, const Polygon_2& startBufferPoly);
    static std::vector<double> get_domain_bbox();

    virtual void reconstruct() = 0;

    void prep_output();
    void prep_output(Vector_2 edge);

    virtual TopoClass   get_class() const = 0;
    virtual std::string get_class_name() const = 0;
    virtual void        get_cityjson_info(nlohmann::json& b) const;
    virtual void        get_cityjson_semantics(nlohmann::json& g) const;
    virtual std::string get_cityjson_primitive() const;

protected:
    static std::vector<Point_3> _outerPts;
    static double               _outerBndHeight;
    std::vector<Point_3>        _sideOutputPts;
};
#endif //CITY4CFD_BOUNDARY_H