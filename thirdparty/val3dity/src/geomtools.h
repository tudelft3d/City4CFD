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

#ifndef __val3dity__geomtools__
#define __val3dity__geomtools__

#include "definitions.h"

namespace val3dity
{

//-- misc
#define PI 3.14159265

CgalPolyhedron::Plane_3  get_best_fitted_plane(const std::vector< Point3 > &lsPts);

bool    cmpPoint3(Point3 &p1, Point3 &p2, double tol);
void    create_cgal_polygon(const std::vector<Point3>& lsPts, const std::vector<int>& ids, const CgalPolyhedron::Plane_3 &plane, Polygon &outpgn);
bool    is_face_planar_distance2plane(const std::vector<Point3> &pts, const CgalPolyhedron::Plane_3 &plane, double& value, float tolerance);
bool    is_face_planar_normals(const std::vector<int*> &trs, const std::vector<Point3>& lsPts, double& value, float angleTolerance);

void mark_domains(CT& ct);
void mark_domains(CT& ct, CT::Face_handle start, int index, std::list<CT::Edge>& border);

Nef_polyhedron* dilate_nef_polyhedron(Nef_polyhedron* nef, float r);
Nef_polyhedron* erode_nef_polyhedron (Nef_polyhedron* nef, float r);
Nef_polyhedron* get_structuring_element_cube(float r);
Nef_polyhedron* get_structuring_element_dodecahedron(float r);
Nef_polyhedron* get_aabb(Nef_polyhedron* mynef) ;

} // namespace val3dity

#endif /* defined(__val3dity__geomtools__) */
