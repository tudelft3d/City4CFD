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

#ifndef CITY4CFD_BOUNDARY_H
#define CITY4CFD_BOUNDARY_H

#include "TopoFeature.h"

class Boundary : public TopoFeature {
public:
    Boundary();
    Boundary(const int outputLayerID);
    virtual ~Boundary() = default;

    static void set_bnd_poly(Polygon_2& bndPoly, Polygon_2& pcBndPoly, Polygon_2& startBufferPoly);
    static void set_bounds_to_buildings_pc(Point_set_3& pointCloud, const Polygon_2& pcBndPoly);
    static void set_bounds_to_terrain_pc(Point_set_3& pointCloud, const Polygon_2& bndPoly,
                                         const Polygon_2& pcBndPoly, const Polygon_2& startBufferPoly);
    static std::vector<double> get_outer_bnd_bbox();

    virtual void reconstruct() = 0;

    void prep_output();
    void prep_output(Vector_2 edge);

    virtual TopoClass   get_class() const = 0;
    virtual std::string get_class_name() const = 0;

protected:
    static std::vector<Point_3> s_outerPts;
    static double               s_outerBndHeight;
    std::vector<Point_3>        m_sideOutputPts;
};
#endif //CITY4CFD_BOUNDARY_H