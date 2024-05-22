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

#ifndef Feature_h
#define Feature_h

#include "Primitive.h"
#include "definitions.h"
#include "nlohmann-json/json.hpp"
#include <vector>
#include <set>
#include <string>

using json = nlohmann::json;


namespace val3dity
{

class Feature
{
public:
                          Feature  ();
  virtual                 ~Feature ();

  virtual bool            validate(double tol_planarity_d2p, double tol_planarity_normals, double tol_overlap = -1) = 0;
  virtual bool            is_valid() = 0;
  virtual std::string     get_type() = 0;

  std::string             get_id();
  void                    set_id(std::string id);
  bool                    is_empty();
  
  void                    add_primitive(Primitive* p);

  const std::vector<Primitive*>&  get_primitives();

  int                     number_of_primitives();

  void                    add_error(int code, std::string info, std::string whichgeoms);
  json                    get_report_json();
  std::set<int>           get_unique_error_codes();

protected:
  int                     _is_valid; 
  std::string             _id;
  std::string             _type;
  std::vector<Primitive*> _lsPrimitives;
  
  bool                    validate_generic(double tol_planarity_d2p, double tol_planarity_normals, double tol_overlap = -1);  
  
  std::map<int, std::vector< std::tuple< std::string, std::string > > > _errors;

};


} // namespace val3dity

#endif /* Feature_h */
