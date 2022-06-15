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

#ifndef CITY4CFD_LOD12_H
#define CITY4CFD_LOD12_H

#include "CGALTypes.h"

class LoD12 {
public:
    LoD12() = delete;
    LoD12(const Polygon_with_holes_2& poly, const std::vector<std::vector<double>>& base_heights,
          const std::vector<double>& building_pts);
    LoD12(const Polygon_with_holes_2& poly, const std::vector<std::vector<double>>& base_heights,
          const std::vector<double>& building_pts, const double height);
    ~LoD12() = default;

    void   lod12_calc_height(double& height);
    void   lod12_reconstruct(Mesh& mesh);
    void   lod12_reconstruct(Mesh& mesh, const double height);
    double get_height() const;

private:
    double _height;
    const Polygon_with_holes_2&              _poly;
    const std::vector<std::vector<double>>&  _baseHeights;
    const std::vector<double>&               _buildingPts;

    void create_mesh(Mesh& mesh);
};

#endif //CITY4CFD_LOD12_H