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

#ifndef __val3dity__Solid__
#define __val3dity__Solid__

#include "Primitive.h"
#include "definitions.h"
#include "Surface.h"

#include <CGAL/Polygon_mesh_processing/bbox.h>

#include <vector>
#include <string>

namespace val3dity
{

class Solid : public Primitive
{
public:
                  Solid(std::string id = "");
                  ~Solid();
  
  Surface*        get_oshell();
  void            set_oshell(Surface* sh);
  Surface*        get_ishell(int i);
  void            add_ishell(Surface* sh);

  int             num_ishells();
  int             num_faces();
  int             num_vertices();

  std::vector<json> get_errors();
 
  bool            validate(double tol_planarity_d2p, double tol_planarity_normals, double tol_overlap = -1);
  Nef_polyhedron* get_nef_polyhedron();
  void            get_min_bbox(double& x, double& y);
  void            translate_vertices();
  std::string     get_poly_representation();
  std::string     get_off_representation(int shellno = 0);
  int             is_valid();
  bool            is_empty();
  Primitive3D     get_type();
  std::set<int>   get_unique_error_codes();

  int             is_point_in_solid(Point3& p);

  CGAL::Bbox_3    get_bbox();
  
  const std::vector<Surface*>&  get_shells();


protected:
  std::vector<Surface*>  _shells;
  Nef_polyhedron*        _nef;

  bool validate_solid_with_nef();
};

} // namespace val3dity

#endif /* defined(__val3dity__Solid__) */
