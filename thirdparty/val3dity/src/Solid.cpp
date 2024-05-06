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


#include "Solid.h"
#include "definitions.h"

#include "input.h"
#include "validate_shell.h"
#include "val3dity_ostream.h"

namespace val3dity
{

Solid::Solid(std::string id)
{
  _id = id;
  _is_valid = -1;
  _nef = NULL;
}


Solid::~Solid()
{}

Surface* Solid::get_oshell()
{
  return _shells[0];
}

Primitive3D Solid::get_type() 
{
  return SOLID;
}

void Solid::set_oshell(Surface* sh)
{
  if (_shells.empty())
    _shells.push_back(sh);
  else
    _shells[0] = sh;
}


const std::vector<Surface*>& Solid::get_shells()
{
  return _shells;
}


void Solid::add_ishell(Surface* sh)
{
  _shells.push_back(sh);
}


int Solid::is_valid()
{
  for (auto& sh : _shells)
  {
    if (sh->has_errors() == true)
    {
      _is_valid = 0;
      return 0;
    }
  }
  if ( (_is_valid == 1) && (this->is_empty() == false) )
    return 1;
  else
    return _is_valid;
}


bool Solid::is_empty()
{
  if (_shells.size() == 0)
    return true;
  for (auto& sh : _shells)
  {
    if (sh->is_empty() == true)
      return true;
  }
  return false;
}


void Solid::get_min_bbox(double& x, double& y)
{
  double tmpx, tmpy;
  double minx = 9e10;
  double miny = 9e10;
  for (auto& sh : _shells)
  {
    sh->get_min_bbox(tmpx, tmpy);
    if (tmpx < minx)
      minx = tmpx;
    if (tmpy < miny)
      miny = tmpy;
  }
  x = minx;
  y = miny;
}


void Solid::translate_vertices()
{
  for (auto& sh : _shells)
    sh->translate_vertices();
}


bool Solid::validate(double tol_planarity_d2p, double tol_planarity_normals, double tol_overlap)
{
  if (this->is_valid() == 0)
  {
    std::string s;
    for (auto& e : this->get_unique_error_codes()) {
      s += std::to_string(e);
      s += "; ";
    }
    *val3dityclog << "Errors: " << s << std::endl;
    return false;
  }
  bool isValid = true;
  if (this->is_empty() == true)
  {
    this->add_error(902, "", "empty Solid, contains no points and/or surfaces");
    return false;
  }
  for (auto& sh : _shells)
  {
    if (sh->validate_as_shell(tol_planarity_d2p, tol_planarity_normals) == false) 
      isValid = false;
  }
  if (isValid == true) 
  {
    if (validate_solid_with_nef() == false)
      isValid = false;
  }
  _is_valid = isValid;
  return isValid;
}


std::string Solid::get_poly_representation()
{
  std::ostringstream ss;
  for (auto& sh : _shells)
  {
    ss << sh->get_off_representation() << std::endl;
  }
  return ss.str();
}

std::string Solid::get_off_representation(int shellno)
{
  return (_shells[shellno])->get_off_representation();
}


json Solid::get_report_json(std::string preid)
{
  json j;
  j["type"] = "Solid";
  if (this->get_id() != "")
    j["id"] = this->_id;
  else
    j["id"] = "none";
  j["numbershells"] = (this->num_ishells() + 1);
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
    }
  }
  int shid = 0;
  for (auto& sh : _shells) {
    std::string t = "shell:" + std::to_string(shid);
    for (auto& each: sh->get_report_json(t)) {
      j["errors"].push_back(each); 
    }
    shid++;
  }
  if (this->is_valid() == 1)
    j["validity"] = true;
  else 
    j["validity"] = false;  
  return j;
}


Nef_polyhedron* Solid::get_nef_polyhedron()
{
  if (_nef != NULL)
    return _nef;
  std::vector<Nef_polyhedron> nefs;
  for (auto& sh : this->get_shells())
  {
    //-- convert to an EPEC Polyhedron so that convertion to Nef is possible
    CgalPolyhedronE pe;
    Polyhedron_convert polyhedron_converter(*(sh->get_cgal_polyhedron()));
    pe.delegate(polyhedron_converter);
    Nef_polyhedron onef(pe);
    nefs.push_back(onef);
  }
  Nef_polyhedron* re = new Nef_polyhedron;
  *re += nefs[0];
  for (int i = 1; i < nefs.size(); i++) 
  {
    *re -= nefs[i];
  }
  _nef = re;
  return re;
}


/*
 -2 = not valid polyhedron
 -1 = outside
  0 = directly on the boundary of polyhedron
  1 = inside
 */
int Solid::is_point_in_solid(Point3& p)
{
  int re = -2;
  if (this->is_valid() == 1)
  {
    auto& s = _shells[0];
    re = s->side_of_triangle_surface(p);
    for (int i = 1; i <= this->num_ishells(); i++)
    {
      s = _shells[i];
      int re2 = s->side_of_triangle_surface(p);
      if ( (re2 == 0) || (re2 == 1) )
        re = -1;
    }
  }
  return re;
}


std::set<int> Solid::get_unique_error_codes() {
  std::set<int> errs = Primitive::get_unique_error_codes();
  for (auto& sh : _shells) {
    std::set<int> tmp = sh->get_unique_error_codes();
    errs.insert(tmp.begin(), tmp.end());
  }
  return errs;
}


int Solid::num_ishells()
{
  return (_shells.size() - 1);
}

int Solid::num_faces()
{
  int total = 0;
  for (auto& sh : _shells)
    total += sh->number_faces();
  return total;
}

int Solid::num_vertices()
{
  int total = 0;
  for (auto& sh : _shells)
    total += sh->number_vertices();
  return total;
}


bool Solid::validate_solid_with_nef()
{
  bool isValid = true;
  //-- check orientation of the normals is outwards or inwards
  *val3dityclog << "-----Global orientation of normals" << std::endl;
  int i = 0;
  for (auto& sh : this->get_shells())
  {
    if (check_global_orientation_normals(sh->get_cgal_polyhedron(), i == 0) == false) 
    {
      this->add_error(405, "i", "");
      isValid = false;
    }
    i++;
  }
  if (isValid == false)
    return false;

  if (this->num_ishells() == 0)
    return true;
    
  *val3dityclog << "---Inspection interactions between the " << (this->num_ishells() + 1) << " shells" << std::endl;
  std::vector<Nef_polyhedron> nefs;
  for (auto& sh : this->get_shells())
  {
    //-- convert to an EPEC Polyhedron so that convertion to Nef is possible
    CgalPolyhedronE pe;
    Polyhedron_convert polyhedron_converter(*(sh->get_cgal_polyhedron()));
    pe.delegate(polyhedron_converter);
    Nef_polyhedron onef(pe);
    nefs.push_back(onef);
  }

  //-- test axiom #1 from the paper, Sect 4.5:
  //-- https://3d.bk.tudelft.nl/hledoux/pdfs/13_cacaie.pdf
  Nef_polyhedron nef;
  for (int i = 1; i < nefs.size(); i++) 
  {
    nef = !nefs[0] * nefs[i];
    if (nef.is_empty() == false)
    {
      nef = nefs[0] * nefs[i];
      if (nef.is_empty() == true)
      {
        std::stringstream msg;
        msg << "Inner shell (#" << i << ") is completely outside the outer shell (#0))";
        this->add_error(403, "i", msg.str());
        isValid = false;
      }
      else
      {
        std::stringstream ss;
        ss << 0 << "--" << i;
        this->add_error(401, ss.str(), "");
        isValid = false; 
      }
    }
    else
    {
      nef = nefs[0] - nefs[i];
      if (nef.number_of_volumes() < 3)
      {
        std::stringstream ss;
        ss << 0 << "--" << i;
        this->add_error(401, ss.str(), "");
        isValid = false; 
      }
      else
      {
        nef = nefs[0].boundary() * nefs[i].boundary();
        if (nef.number_of_facets() > 0)
        {
          std::stringstream ss;
          ss << 0 << "--" << i;
          this->add_error(401, ss.str(), "");
          isValid = false;
        }
      }
    }
  }

  //-- test axiom #2 from the paper 
  nef.clear();
  for (int i = 1; i < nefs.size(); i++) 
  {
    for (int j = (i + 1); j < nefs.size(); j++) 
    {
      //-- 1. are they the same?
      if (nefs[i] == nefs[j])
      {
        std::stringstream ss;
        ss << i << "--" << j;
        this->add_error(402, ss.str(), "");
        isValid = false;
        continue;
      }
      //-- 2. intersection
      nef = nefs[i] * nefs[j];
      if (nef.number_of_volumes() > 1)
      {
        std::stringstream ss;
        ss << i << "--" << j;
        this->add_error(401, ss.str(), "");
        isValid = false;
        continue;
      }
      //-- 3. union
      nef = nefs[i] + nefs[j];
      if (nef.number_of_volumes() < 3)
      {
        std::stringstream ss;
        ss << 0 << "--" << i;
        this->add_error(401, ss.str(), "");
        isValid = false;
        continue;
      }
    }
  }
  //-- test axiom #3 from the paper
  if (isValid == true)
  {
    nef.clear();
    nef += nefs[0];
    int numvol = 2;
    for (int i = 1; i < nefs.size(); i++) 
    {
      nef = nef - nefs[i];
      numvol++;
      if (nef.number_of_volumes() != numvol)
      {
        this->add_error(404, "", "");
        isValid = false;
        break;
      }
    }
  }
  return isValid;
}


CGAL::Bbox_3 Solid::get_bbox() 
{
  CgalPolyhedron* poly = this->get_oshell()->get_cgal_polyhedron();
  return CGAL::Polygon_mesh_processing::bbox(*poly);
}


} // namespace val3dity
