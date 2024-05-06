/*
  val3dity 

  Copyright (c) 2011-2023, 3D geoinformation research group, TU Delft  

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

#ifndef Surface_h
#define Surface_h

#include "definitions.h"
#include "nlohmann-json/json.hpp"
#include <string>
#include <vector>
#include <set>
#include <unordered_map>

using json = nlohmann::json;

namespace val3dity
{

class Surface
{
public:
  Surface  (int id = -1, double tol_snap = 0.0);
  ~Surface ();
  
  bool validate_as_shell(double tol_planarity_d2p, double tol_planarity_normals);
  bool validate_as_multisurface(double tol_planarity_d2p, double tol_planarity_normals);
  bool validate_as_compositesurface(double tol_planarity_d2p, double tol_planarity_normals);
  
  bool is_shell(double tol_planarity_d2p, double tol_planarity_normals);

  CgalPolyhedron* get_cgal_polyhedron();

  int    number_vertices();
  int    number_faces();
  void   get_min_bbox(double& x, double& y);
  int    get_id();

  bool   does_self_intersect();
  bool   is_empty();
  int    add_point(Point3 p);
  void   add_face(std::vector< std::vector<int> > f, std::string id = "");

  json          get_report_json(std::string preid = "");
  void          add_error(int code, std::string faceid = "", std::string info = "");
  bool          has_errors();
  std::set<int> get_unique_error_codes();
  void          translate_vertices();
  static void   set_translation_min_values(double minx, double miny);
  std::string   get_poly_representation();
  std::string   get_off_representation();

  int           side_of_triangle_surface(Point3& p);

  bool          were_vertices_merged_during_parsing();
  int           get_number_parsed_vertices();
  
private:
  int                                         _id;
  std::vector<Point3>                         _lsPts;
  std::vector<std::vector<std::vector<int>>>  _lsFaces;
  std::vector<std::string>                    _lsFacesID;
  std::vector<std::vector<int*>>              _lsTr;
  CgalPolyhedron*                             _polyhedron;
  double                                      _tol_snap;
  int                                         _is_valid_2d; //-1: not done yet; 0: nope; 1: yes it's valid
  int                                         _vertices_added;
  static double                               _shiftx;
  static double                               _shifty;

  std::map<int, std::vector<std::tuple<std::string, std::string> > > _errors;
  
  bool validate_2d_primitives(double tol_planarity_d2p, double tol_planarity_normals);
  std::string get_coords_key(Point3* p);
  bool triangulate_shell();
  std::vector<int*> construct_ct_one_face(const std::vector<std::vector<int>>& pgnids);
  bool validate_polygon(std::vector<Polygon> &lsRings, std::string polygonid);
  bool validate_projected_ring(Polygon &pgn, std::string id);
  bool has_face_rings_toofewpoints(const std::vector< std::vector<int> >& theface);
  bool has_face_2_consecutive_repeated_pts(const std::vector< std::vector<int> >& theface);
  bool contains_nonmanifold_vertices();

};

} // namespace val3dity

#endif /* Surface_h */
