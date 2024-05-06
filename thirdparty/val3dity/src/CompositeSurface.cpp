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

#include "CompositeSurface.h"
#include "Primitive.h"
#include "input.h"

namespace val3dity
{

CompositeSurface::CompositeSurface(std::string id) {
  _id = id;
  _is_valid = -1;
}

CompositeSurface::~CompositeSurface() {
}

bool CompositeSurface::validate(double tol_planarity_d2p, double tol_planarity_normals, double tol_overlap)
{
  if (this->is_valid() == 0)
    return false;
  if (_surface->validate_as_compositesurface(tol_planarity_d2p, tol_planarity_normals) == true) 
  {
    _is_valid = 1;
    return true;
  }
  else
  {
    _is_valid = 0;
    return false;
  }
}


Primitive3D CompositeSurface::get_type() 
{
  return COMPOSITESURFACE;
}

int CompositeSurface::is_valid() {
  if (_surface->has_errors() == true)
  {
    _is_valid = 0;
    return 0;
  }
  if ( (_is_valid == 1) && (this->is_empty() == false) )
    return 1;
  else
    return _is_valid;
}


void CompositeSurface::get_min_bbox(double& x, double& y)
{
  _surface->get_min_bbox(x, y);
}


void CompositeSurface::translate_vertices()
{
  _surface->translate_vertices();
}

bool CompositeSurface::is_empty() {
  return _surface->is_empty();
}


json CompositeSurface::get_report_json(std::string preid)
{
  json j;
  bool isValid = true;
  j["type"] = "CompositeSurface";
  if (this->get_id() != "")
    j["id"] = this->_id;
  else
    j["id"] = "none";
  j["numberfaces"] = this->num_faces();
  j["numbervertices"] = this->num_vertices();  
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
  for (auto& each: _surface->get_report_json())
    j["errors"].push_back(each); 
  if (_surface->has_errors() == true)
    isValid = false;
  j["validity"] = isValid;  
  return j;
}


std::string CompositeSurface::get_off_representation()
{
  return _surface->get_off_representation();
}


int CompositeSurface::num_faces() 
{
  return _surface->number_faces();
}

int CompositeSurface::num_vertices() 
{
  return _surface->number_vertices();
}

bool CompositeSurface::set_surface(Surface* s) 
{
  _surface = s;
  return true;
}

std::set<int> CompositeSurface::get_unique_error_codes() {
  std::set<int> errs = Primitive::get_unique_error_codes();
  std::set<int> tmp = _surface->get_unique_error_codes();
  errs.insert(tmp.begin(), tmp.end());
  return errs;
}


Surface* CompositeSurface::get_surface() 
{
  return _surface;
}

} // namespace val3dity
