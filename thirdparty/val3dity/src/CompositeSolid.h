/*
  val3dity 

  Copyright (c) 2011-2024, 3D geoinformation research group, TU Delft  

  This file is part of val3dity.

  val3dity is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  val3dity is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with val3dity.  If not, see <http://www.gnu.org/licenses/>.

  For any information or further details about the use of val3dity, contact
  Hugo Ledoux
  <h.ledoux@tudelft.nl>
  Faculty of Architecture & the Built Environment
  Delft University of Technology
  Julianalaan 134, Delft 2628BL, the Netherlands
*/

#ifndef CompositeSolid_h
#define CompositeSolid_h

#include "Primitive.h"
#include "Solid.h"

#include <vector>
#include <string>

namespace val3dity
{

class CompositeSolid : public Primitive 
{
public:
                CompositeSolid(std::string id = ""); 
                ~CompositeSolid(); 

  bool          validate(double tol_planarity_d2p, double tol_planarity_normals, double tol_overlap = -1);
  int           is_valid();
  bool          is_empty();
  std::vector<json> get_errors();
  Primitive3D   get_type();
  std::set<int> get_unique_error_codes();

  void          get_min_bbox(double& x, double& y);
  void          translate_vertices();

  Nef_polyhedron* get_nef_polyhedron();

  bool          add_solid(Solid* s);
  int           number_of_solids();

protected:
  std::vector<Solid*> _lsSolids;
  Nef_polyhedron*     _nef;
};

} // namespace val3dity

#endif /* CompositeSolid_h */
