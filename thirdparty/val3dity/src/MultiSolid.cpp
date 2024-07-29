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

#include "MultiSolid.h"
#include "input.h"

namespace val3dity
{

MultiSolid::MultiSolid(std::string id) {
  _id = id;
  _is_valid = -1;
}

MultiSolid::~MultiSolid() 
{
  for (auto& sol : _lsSolids) {
    delete sol;
  }
}

Primitive3D MultiSolid::get_type() 
{
  return MULTISOLID;
}

bool MultiSolid::validate(double tol_planarity_d2p, double tol_planarity_normals, double tol_overlap) 
{
  bool isValid = true;
  for (auto& s : _lsSolids)
  {
    if (s->validate(tol_planarity_d2p, tol_planarity_normals) == false)
      isValid = false;
  }
  _is_valid = isValid;
  return isValid;
}


int MultiSolid::is_valid() 
{
  if ( (_is_valid == 1) && (this->is_empty() == false) && (_errors.empty() == true) )
    return 1;
  else
    return _is_valid;
}


void MultiSolid::get_min_bbox(double& x, double& y)
{
  double tmpx, tmpy;
  double minx = 9e10;
  double miny = 9e10;
  for (auto& s : _lsSolids)
  {
    s->get_min_bbox(tmpx, tmpy);
    if (tmpx < minx)
      minx = tmpx;
    if (tmpy < miny)
      miny = tmpy;
  }
  x = minx;
  y = miny;
}


void MultiSolid::translate_vertices()
{
  for (auto& s : _lsSolids)
    s->translate_vertices();
}

bool MultiSolid::is_empty() 
{
  return _lsSolids.empty();
}

std::vector<json> MultiSolid::get_errors()
{
  std::vector<json> js;
  for (auto& err : _errors)
  {
    for (auto& e : _errors[std::get<0>(err)])
    {
      json j;
      j["code"] = std::get<0>(err);
      j["description"] = ALL_ERRORS[std::get<0>(err)];
      j["id"] = std::get<0>(e);
      j["info"] = std::get<1>(e);
      js.push_back(j);
    }
  }
  int solid = 0;
  for (auto& sol : _lsSolids)
  {
    auto e = sol->get_errors();
    js.insert(js.end(), e.begin(), e.end());    
    solid++;
  }
  return js;
}

bool MultiSolid::add_solid(Solid* s) {
  _lsSolids.push_back(s);
  return true;
}

std::set<int> MultiSolid::get_unique_error_codes() {
  std::set<int> errs = Primitive::get_unique_error_codes();
  for (auto& s : _lsSolids) {
    std::set<int> tmp = s->get_unique_error_codes();
    errs.insert(tmp.begin(), tmp.end());
  }
  return errs;
}

int MultiSolid::number_of_solids() {
  return _lsSolids.size();
}

} // namespace val3dity
