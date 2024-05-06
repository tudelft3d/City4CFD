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

#include "GeometryTemplate.h"
#include "input.h"

namespace val3dity
{

GeometryTemplate::GeometryTemplate(std::string id) {
  _id = id;
  _is_valid = -1;
}

GeometryTemplate::~GeometryTemplate() {
}

Primitive3D GeometryTemplate::get_type() 
{
  return GEOMETRYTEMPLATE;
}

bool GeometryTemplate::validate(double tol_planarity_d2p, double tol_planarity_normals, double tol_overlap) 
{
  bool isValid = true;
  for (auto& p : _lsPrimitives)
  {
    if (p->validate(tol_planarity_d2p, tol_planarity_normals) == false)
      isValid = false;
  }
  _is_valid = isValid;
  return isValid;
}


int GeometryTemplate::is_valid() 
{
  if ( (_is_valid == 1) && (this->is_empty() == false) && (_errors.empty() == true) )
    return 1;
  else
    return _is_valid;
}


void GeometryTemplate::get_min_bbox(double& x, double& y)
{
  x = 0.0;
  y = 0.0;
}


void GeometryTemplate::translate_vertices()
{}

bool GeometryTemplate::is_empty() 
{
  return _lsPrimitives.empty();
}


json GeometryTemplate::get_report_json(std::string preid)
{
  json j;
  bool isValid = true;
  j["type"] = "GeometryTemplate";
  if (this->get_id() != "")
    j["id"] = this->_id;
  else
    j["id"] = "none";
  j["errors"] = json::array();
  for (auto& err : _errors)
  {
    for (auto& e : _errors[std::get<0>(err)])
    {
      json jj;
      jj["type"] = "Error";
      jj["code"] = std::get<0>(err);
      jj["description"] = ALL_ERRORS[std::get<0>(err)];
      jj["id"] = std::get<0>(e);
      jj["info"] = std::get<1>(e);
      j["errors"].push_back(jj);
      isValid = false;
    }
  }
  for (auto& p : _lsPrimitives) //-- there should be just one by definition
  {
    j["template"].push_back(p->get_report_json());
    if (p->is_valid() == false)
      isValid = false;
  }
  j["validity"] = isValid;
  return j;
}



bool GeometryTemplate::add_primitive(Primitive* s) {
  _lsPrimitives.push_back(s);
  return true;
}

std::set<int> GeometryTemplate::get_unique_error_codes() {
  std::set<int> errs = Primitive::get_unique_error_codes();
  for (auto& p : _lsPrimitives) {
    std::set<int> tmp = p->get_unique_error_codes();
    errs.insert(tmp.begin(), tmp.end());
  }
  return errs;
}

// int GeometryTemplate::number_of_solids() {
//   return _lsSolids.size();
// }

} // namespace val3dity
