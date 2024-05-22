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

#ifndef Validate_prim_toporel_h
#define Validate_prim_toporel_h

#include "definitions.h"
#include "Solid.h"

namespace val3dity
{


class Primitive;

bool do_primitives_interior_overlap(std::vector<Primitive*>& lsPrimitives, 
                                    int errorcode_to_assign, 
                                    std::vector<Error>& lsErrors, 
                                    double tol_overlap);

bool are_cells_interior_disconnected_with_aabb(std::vector<std::tuple<std::string,Solid*>>& lsCells,
                                               int errorcode_to_assign, 
                                               std::vector<Error>& lsErrors, 
                                               double tol_overlap);

int are_primitives_adjacent(Primitive* p1, Primitive* p2, double tol_overlap);



} // namespace val3dity

#endif /* Validate_prim_toporel_h */
