/*
  val3dity 

  Copyright (c) 2011-2022, 3D geoinformation research group, TU Delft  

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

#include "Surface.h"
#include <CGAL/intersections.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/intersections.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>

#include <set>
#include <list>

namespace val3dity
{

template <class HDS>
class ConstructShell : public CGAL::Modifier_base<HDS> {
  std::vector< std::vector<int*> > *faces;
  std::vector<Point3> *lsPts;
  int _width;
  Surface* sh;
public:
  bool isValid;
  ConstructShell(std::vector< std::vector<int*> > *faces, std::vector<Point3> *lsPts, Surface* sh)
    :faces(faces), lsPts(lsPts), sh(sh), isValid(true), _width(static_cast<int>(lsPts->size()))
  {
  }
  void operator()( HDS& hds);
  void construct_faces_order_given(CGAL::Polyhedron_incremental_builder_3<HDS>& B);
  int m2a(int m, int n);
  void construct_faces_flip_when_possible(CGAL::Polyhedron_incremental_builder_3<HDS>& B);
  bool try_to_add_face(CGAL::Polyhedron_incremental_builder_3<HDS>& B, std::list<int*>& trFaces, bool* halfedges, bool bMustBeConnected);
  bool is_connected(int* tr, bool* halfedges);
  void add_one_face(CGAL::Polyhedron_incremental_builder_3<HDS>& B, int i0, int i1, int i2, std::string faceID) ;
};


CgalPolyhedron*   construct_CgalPolyhedron_incremental(std::vector< std::vector<int*> > *lsTr, std::vector<Point3> *lsPts, Surface* sh);
CgalPolyhedron*   construct_CgalPolyhedron_batch(const std::vector< std::vector<int*> >&lsTr, const std::vector<Point3>& lsPts);
bool              check_global_orientation_normals(CgalPolyhedron* p, bool bOuter);

} // namespace val3dity
