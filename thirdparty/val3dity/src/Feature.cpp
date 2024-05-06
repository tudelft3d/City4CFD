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

#include "Feature.h"
#include "input.h"
#include "val3dity_ostream.h"

namespace val3dity
{

std::string Feature::get_id()
{
  return _id;
}


void Feature::set_id(std::string id)
{
  _id = id;
}


bool Feature::is_empty()
{
  return _lsPrimitives.empty();
}


json Feature::get_report_json() 
{
  json j;
  j["id"] = _id;
  j["type"] = _type;
  j["validity"] = this->is_valid();
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
    }
  }
  j["primitives"];
  for (auto& p : _lsPrimitives)
    j["primitives"].push_back(p->get_report_json()); 
  return j;
}


void Feature::add_primitive(Primitive* p)
{
  _lsPrimitives.push_back(p);
}

int Feature::number_of_primitives() 
{
  return _lsPrimitives.size();
}

const std::vector<Primitive*>& Feature::get_primitives()
{
  return _lsPrimitives;
}


bool Feature::validate_generic(double tol_planarity_d2p, double tol_planarity_normals, double tol_overlap)
{
  *val3dityclog << std::endl << "######### Validating Feature #########" << std::endl;
  *val3dityclog << "id:   " << this->get_id() << std::endl;
  *val3dityclog << "type: " << this->get_type() << std::endl;
  *val3dityclog << "--" << std::endl;
  bool bValid = true;
  if (this->is_empty() == true) {
    this->add_error(906, "Feature has no geometry defined (or val3dity doesn't handle this type).", "");
    bValid = false;
  }
  if (_lsPrimitives.size() > 500) {
    *val3ditycout << "Validating " << _lsPrimitives.size() << " geometric primitives, this could be slow." << std::endl << std::flush;
  }
  for (auto& p : _lsPrimitives)
  {
    *val3dityclog << "======== Validating Primitive ========" << std::endl;
    switch(p->get_type())
    {
      case 0: *val3dityclog << "Solid"             << std::endl; break;
      case 1: *val3dityclog << "CompositeSolid"    << std::endl; break;
      case 2: *val3dityclog << "MultiSolid"        << std::endl; break;
      case 3: *val3dityclog << "CompositeSurface"  << std::endl; break;
      case 4: *val3dityclog << "MultiSurface"      << std::endl; break;
      case 5: *val3dityclog << "GeometryTemplate"  << std::endl; break;
      case 9: *val3dityclog << "ALL"               << std::endl; break;
    }
    *val3dityclog << "id: " << p->get_id() << std::endl;
    *val3dityclog << "--" << std::endl;
    if (p->validate(tol_planarity_d2p, tol_planarity_normals, tol_overlap) == false)
    {
      *val3dityclog << "======== INVALID ========" << std::endl;
      bValid = false;
    }
    else
      *val3dityclog << "========= VALID =========" << std::endl;
  }
  _is_valid = bValid;
  return bValid;
}  


void Feature::add_error(int code, std::string whichgeoms, std::string info)
{
  _is_valid = 0;
  std::tuple<std::string, std::string> a(whichgeoms, info);
  _errors[code].push_back(a);
  *val3dityclog << "\tERROR " << code << ": " << ALL_ERRORS[code];
  if (whichgeoms.empty() == false)
    *val3dityclog << " (id: " << whichgeoms << ")";
  *val3dityclog << std::endl;
  if (info.empty() == false)
    *val3dityclog << "\t[" << info << "]" << std::endl;
}

std::set<int> Feature::get_unique_error_codes()
{
  std::set<int> errs;
  //-- errors related to Features
  for (auto& err : _errors)
    errs.insert(std::get<0>(err));
  //-- errors for each primitive
  for (auto& p : _lsPrimitives)
    for (auto& each : p->get_unique_error_codes())
      errs.insert(each);
  return errs;
}

} // namespace val3dity
