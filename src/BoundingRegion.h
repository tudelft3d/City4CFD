/*
  City4CFD
 
  Copyright (c) 2021-2022, 3D Geoinformation Research Group, TU Delft  

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

#ifndef CITY4CFD_BOUNDINGREGION_H
#define CITY4CFD_BOUNDINGREGION_H

#include "types.h"
#include "CGALTypes.h"

class BoundingRegion {
public:
    BoundingRegion();
    ~BoundingRegion();

    void operator()(double radius);
    void operator()(Polygon_2& poly);

    void calc_influ_region_bpg(const DT& dt, const Point_set_3& pointCloudBuildings, Buildings& buildings);
    void calc_bnd_bpg(const Polygon_2& influRegionPoly, const Buildings& buildings);

    Polygon_2& get_bounding_region();
    const Polygon_2& get_bounding_region() const;

protected:
    Polygon_2  _boundingRegion;

    Polygon_2  calc_bnd_poly(const std::vector<Point_2>& candidatePts, const double hMax,
                             const double angle, const double enlargeRatio = 1) const;
    double     calc_blockage_ratio(const Buildings& buildings, const double angle, Polygon_2& localPoly) const;

};

#endif //CITY4CFD_BOUNDINGREGION_H