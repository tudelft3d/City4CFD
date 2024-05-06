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


#include "CityObject.h"
#include "definitions.h"
#include "input.h"
#include "validate_prim_toporel.h"
#include "val3dity_ostream.h"

namespace val3dity
{

CityObject::CityObject(std::string theid, std::string thetype)
{
  _id = theid;
  _is_valid = -1;
  _type = thetype; // TODO : validate the type here?
}


CityObject::~CityObject()
{}


bool CityObject::validate(double tol_planarity_d2p, double tol_planarity_normals, double tol_overlap) 
{
  if (_is_valid != -1)
    return _is_valid;
  bool bValid = Feature::validate_generic(tol_planarity_d2p, tol_planarity_normals, tol_overlap);
  //-- Building
  if ( (bValid == true) && (this->_type == "Building") )
    bValid = validate_building(tol_overlap);
  _is_valid = bValid;
  return bValid;
}


bool CityObject::validate_building(double tol_overlap)
{
  bool bValid = true;
  if (_lsPrimitives.size() == 0)
    return bValid;
  *val3dityclog << "--- Interactions between BuildingParts ---" << std::endl;
  std::vector<Error> lsErrors;
  //-- process each LoDs separately
  std::set<std::string> uniquelods;
  for (auto& each: _lsPrimitives) {
    uniquelods.insert(each->get_lod());
  }
  for (auto& lod: uniquelods) {
    std::vector<Primitive*> g1lod;
    for (auto& p: _lsPrimitives) {
      if (p->get_lod() == lod) {
        g1lod.push_back(p);
      }
    }
    // *val3ditycout << g1lod.size() << std::endl;
    if (g1lod.size() > 1) {
      *val3dityclog << "LoD" << lod << " has " << g1lod.size()
        << " geometries, validating for overlaps." << std::endl;
      if (do_primitives_interior_overlap(g1lod, 601, lsErrors, tol_overlap) == true)
      {
        bValid = false;
        *val3dityclog << "Error: overlapping building parts (LoD" << lod << ")" << std::endl;
        for (auto& e : lsErrors) {
          std::stringstream msg;
          msg << "Geometries for LoD" << lod;
          this->add_error(e.errorcode, e.info1, msg.str());
        }
      }
    }
  }
  return bValid;
}


bool CityObject::is_valid() 
{
  if (_is_valid == 1) 
    return true;
  else
    return false;
}


std::string CityObject::get_type() {
  return _type;
}


} // namespace val3dity
