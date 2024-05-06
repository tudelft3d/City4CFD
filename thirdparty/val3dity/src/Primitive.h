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

#ifndef Primitive_h
#define Primitive_h

#include "definitions.h"
#include "nlohmann-json/json.hpp"
#include <map>
#include <vector>
#include <set>
#include <string>

using json = nlohmann::json;


namespace val3dity
{

class Primitive
{
public:
  Primitive  ();
  ~Primitive ();

  virtual bool          validate(double tol_planarity_d2p, double tol_planarity_normals, double tol_overlap = -1) = 0;
  virtual int           is_valid() = 0;
  virtual bool          is_empty() = 0;
  virtual json          get_report_json(std::string preid = "") = 0;
  virtual Primitive3D   get_type() = 0;

  virtual void          get_min_bbox(double& x, double& y) = 0;
  virtual void          translate_vertices() = 0;
  static void           set_translation_min_values(double minx, double miny);

  std::string           get_id();
  void                  set_id(std::string id);
  std::string           get_lod();
  void                  set_lod(std::string lod);  
  void                  add_error(int code, std::string whichgeoms, std::string info);
  virtual std::set<int> get_unique_error_codes();

protected:
  std::string           _id;
  int                   _is_valid; 
  std::string           _lod;
  static double         _shiftx;
  static double         _shifty;

  std::map<int, std::vector< std::tuple< std::string, std::string > > > _errors;

};


} // namespace val3dity

#endif /* Primitive_h */
