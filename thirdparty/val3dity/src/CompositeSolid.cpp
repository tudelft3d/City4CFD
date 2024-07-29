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

#include "CompositeSolid.h"
#include "input.h"
#include "geomtools.h"

namespace val3dity
{

CompositeSolid::CompositeSolid(std::string id)
{
  _id = id;
  _is_valid = -1;
}


CompositeSolid::~CompositeSolid()
{
  delete _nef;
  for (auto& sol : _lsSolids) {
    delete sol;
  }
}


Primitive3D CompositeSolid::get_type() 
{
  return COMPOSITESOLID;
}


Nef_polyhedron* CompositeSolid::get_nef_polyhedron()
{
  if (_nef != NULL)
    return _nef;
  Nef_polyhedron* unioned = new Nef_polyhedron(Nef_polyhedron::EMPTY);
  for (int i = 0; i < _lsSolids.size(); i++)
  {
    Nef_polyhedron* tmp = _lsSolids[i]->get_nef_polyhedron();
    *unioned = *unioned + *tmp;
    delete tmp;
  }
  _nef = unioned;
  return unioned;
}


void CompositeSolid::get_min_bbox(double& x, double& y)
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


void CompositeSolid::translate_vertices()
{
  for (auto& s : _lsSolids)
    s->translate_vertices();
}


bool CompositeSolid::validate(double tol_planarity_d2p, double tol_planarity_normals, double tol_overlap) 
{
  bool isValid = true;
  for (auto& s : _lsSolids)
  {
    if (s->validate(tol_planarity_d2p, tol_planarity_normals) == false)
      isValid = false;
  }
  if (isValid == true) 
  {
    std::vector<Nef_polyhedron*> lsNefs;
    for (int i = 0; i < _lsSolids.size(); i++)
      lsNefs.push_back(_lsSolids[i]->get_nef_polyhedron());
//-- 1. check if any 2 are the same? ERROR:502
    // std::clog << "-----Are two solids duplicated" << std::endl;
    for (int i = 0; i < (lsNefs.size() - 1); i++)
    {
      for (int j = i + 1; j < lsNefs.size(); j++) 
      {
        if (*lsNefs[i] == *lsNefs[j])
        {
          std::stringstream msg1, msg2;
          msg1 << "Geometry (CompositeSolid) #" << this->get_id();
          msg2 << "solid #" << _lsSolids[i]->get_id() << " and solid #" << _lsSolids[j]->get_id();
          this->add_error(502, msg1.str(), msg2.str());
          isValid = false;
        }
      }
    }
    if (isValid == true)
    {
//-- 2. check if their interior intersects ERROR:501
      // std::clog << "-----Intersections of solids" << std::endl;
      std::vector<Nef_polyhedron*> lsNefsEroded;
      if (tol_overlap > 0.0)
      {
        for (auto each : lsNefs)
          lsNefsEroded.push_back(erode_nef_polyhedron(each, tol_overlap));
      }
      else
      {
        lsNefsEroded = lsNefs;
      }
      Nef_polyhedron emptynef(Nef_polyhedron::EMPTY);
      for (int i = 0; i < (lsNefsEroded.size() - 1); i++)
      {
        Nef_polyhedron* a = lsNefsEroded[i];
        for (int j = i + 1; j < lsNefsEroded.size(); j++) 
        {
          Nef_polyhedron* b = lsNefsEroded[j];
          if (a->interior() * b->interior() != emptynef)
          {
            std::stringstream msg1, msg2;
            msg1 << "Geometry (CompositeSolid) #" << this->get_id();
            msg2 << "solid #" << _lsSolids[i]->get_id() << " and solid #" << _lsSolids[j]->get_id();
            this->add_error(501, msg1.str(), msg2.str());
            isValid = false;
          }
        }
      }
      if (tol_overlap > 0.0)
      {
        for (auto each : lsNefsEroded)
          delete each;
      }
    }
    if (isValid == true)
    {
//-- 3. check if their union yields one solid ERROR:503
      // std::clog << "-----Forming one solid (union)" << std::endl;
      std::vector<Nef_polyhedron*> lsNefsDilated;
      if (tol_overlap > 0.0)
      {
        for (auto each : lsNefs)
          lsNefsDilated.push_back(dilate_nef_polyhedron(each, tol_overlap));
      }
      else
      {
        lsNefsDilated = lsNefs;
      }
      Nef_polyhedron unioned(Nef_polyhedron::EMPTY);
      for (auto each : lsNefsDilated)
        unioned = unioned + *each;
      if (unioned.number_of_volumes() != 2)
      {
        std::stringstream msg1, msg2;
        msg1 << "Geometry (CompositeSolid) #" << this->get_id();
        msg2 << "CompositeSolid is formed of " << (unioned.number_of_volumes() - 1) << " parts";
        this->add_error(503, msg1.str(), msg2.str());
        isValid = false;
      }
      for (auto each : lsNefsDilated)
        delete each;
    } 
    //-- clear the temp Nefs
    if (tol_overlap > 0.0)
    {
      for (auto each: lsNefs)
        delete each;
    }
  }
  _is_valid = isValid;
  return isValid;
}


int CompositeSolid::is_valid()
{
  if ( (_is_valid == 1) && (this->is_empty() == false) && (_errors.empty() == true) )
    return 1;
  else
    return _is_valid;
}



bool CompositeSolid::is_empty() {
  return _lsSolids.empty();
}


std::vector<json> CompositeSolid::get_errors()
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


std::set<int> CompositeSolid::get_unique_error_codes() {
  std::set<int> errs = Primitive::get_unique_error_codes();
  for (auto& s : _lsSolids) {
    std::set<int> tmp = s->get_unique_error_codes();
    errs.insert(tmp.begin(), tmp.end());
  }
  return errs;
}


bool CompositeSolid::add_solid(Solid* s) {
  _lsSolids.push_back(s);
  return true;
}

int CompositeSolid::number_of_solids() {
  return _lsSolids.size();
}

} // namespace val3dity
