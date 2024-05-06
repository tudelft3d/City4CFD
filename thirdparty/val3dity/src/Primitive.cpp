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

#include "Primitive.h"
#include "input.h"
#include "val3dity_ostream.h"

namespace val3dity
{

double Primitive::_shiftx = 0.0;
double Primitive::_shifty = 0.0;


Primitive::Primitive() {
}

Primitive::~Primitive() {
}

void Primitive::set_translation_min_values(double minx, double miny)
{
  Primitive::_shiftx = minx;
  Primitive::_shifty = miny;
}


std::string  Primitive::get_id()
{
  return _id;
}

void Primitive::set_id(std::string id)
{
  _id = id;
}

std::string  Primitive::get_lod()
{
  return _lod;
}

void Primitive::set_lod(std::string lod)
{
  _lod = lod;
}


void Primitive::add_error(int code, std::string whichgeoms, std::string info)
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

std::set<int> Primitive::get_unique_error_codes()
{
  std::set<int> errs;
  for (auto& err : _errors)
  {
    errs.insert(std::get<0>(err));
  }
  return errs;
}

} // namespace val3dity
