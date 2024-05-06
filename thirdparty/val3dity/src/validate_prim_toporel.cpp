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

#include "validate_prim_toporel.h"
#include "Primitive.h"
#include "geomtools.h"
#include "input.h"
#include "Solid.h"
#include "CompositeSolid.h"
#include "val3dity_ostream.h"
#include <sstream>

#include <CGAL/box_intersection_d.h>

namespace val3dity
{

typedef std::vector<Nef_polyhedron*>                                    Nefs;
typedef Nefs::iterator                                                  Iterator;
typedef CGAL::Box_intersection_d::Box_with_handle_d<double,3,Iterator>  AABB;


struct Report_intersections {
  Nefs* nefs;
  std::vector<std::string>* lsCellIDs;
  std::vector<Error>* lsErrors; 
  int ecode;
  int* thecount;

  Report_intersections(Nefs& nefs, std::vector<std::string>& lsCell, std::vector<Error>& le, int code, int& count)
    : nefs(&nefs), lsCellIDs(&lsCell), lsErrors(&le), ecode(code), thecount(&count)
  {}

  // callback functor that reports when 2 AABBs intersect
  void operator()(const AABB& a, const AABB& b)  
  {
    int id1 = (a.handle() - nefs->begin());
    int id2 = (b.handle() - nefs->begin());
    (*thecount)++;
    if (*thecount % 25 == 0) {
      int r = rand() % 3;
      if (r == 0)
        *val3ditycout << "\b|" << std::flush;
      else if (r == 1)
        *val3ditycout << "\b/" << std::flush;
      else if (r == 2)
        *val3ditycout << "\b-" << std::flush;
      else
        *val3ditycout << "\b\\" << std::flush;
    }
    //-- we check here if the 2 (usually eroded) Nefs are intersecting, 
    //-- if yes then an error is added to the list lsErrors
    Nef_polyhedron* n1 = nefs->at(id1);
    Nef_polyhedron* n2 = nefs->at(id2);
    Nef_polyhedron emptynef(Nef_polyhedron::EMPTY);
    if (n1->interior() * n2->interior() != emptynef)
    {
      Error e;
      std::stringstream msg;
      msg << lsCellIDs->at(id1) << " and " << lsCellIDs->at(id2);
      e.errorcode = ecode;
      e.info1 = msg.str();
      e.info2 = "";
      lsErrors->push_back(e);
    }
  }
};



bool are_cells_interior_disconnected_with_aabb(std::vector<std::tuple<std::string,Solid*>>& lsCells,
                                               int errorcode_to_assign, 
                                               std::vector<Error>& lsErrors, 
                                               double tol_overlap)
{
  *val3dityclog << "--- Constructing Nef Polyhedra ---" << std::endl;
  std::vector<Nef_polyhedron*>                lsNefs;
  std::vector<std::tuple<std::string,Solid*>> subsetCells;
  for (auto& c : lsCells)
  {
    Solid* ts = std::get<1>(c);
    // TODO: only valid Solids are processed: what's the best way here?
    if (ts->is_valid() != 1)
      continue;
    Nef_polyhedron* tmpnef = ts->get_nef_polyhedron();
    if (tol_overlap > 0)
    {
      lsNefs.push_back(erode_nef_polyhedron(tmpnef, tol_overlap));
      // delete tmpnef;
    }
    else
      lsNefs.push_back(tmpnef);
    subsetCells.push_back(c);
  }
  *val3dityclog << "--- Constructing AABB tree ---" << std::endl;
  std::vector<AABB> aabbs;
  int counter = 0;
  for ( Iterator i = lsNefs.begin(); i != lsNefs.end(); ++i)
  {
    Solid* ts = std::get<1>(subsetCells[counter]);
    aabbs.push_back( AABB( ts->get_bbox(), i) );
    counter++;
  }
  std::vector<std::string> lsCellIDs;
  for (auto each : subsetCells)
    lsCellIDs.push_back(std::get<0>(each));
  *val3dityclog << "--- Testing intersections between Nefs ---" << std::endl;
  int n = lsErrors.size();
  int count = 0; 
  if (lsNefs.size() > 500) {
    *val3ditycout << "Testing intersections between " << lsNefs.size() << " cells, this could be slow." << std::endl << std::flush;
  }
  CGAL::box_self_intersection_d( aabbs.begin(), aabbs.end(), Report_intersections(lsNefs, lsCellIDs, lsErrors, errorcode_to_assign, count));
  *val3dityclog << "Total AABB tests: " << count << std::endl;
  if (lsErrors.size() > n)
    return false;
  else
    return true;
}


bool do_primitives_interior_overlap(std::vector<Primitive*>& lsPrimitives, 
                                    int errorcode_to_assign, 
                                    std::vector<Error>& lsErrors, 
                                    double tol_overlap)
{
  bool isValid = true;
  //-- 1. create Nef for all primitives and erode if necessary
  std::vector<Nef_polyhedron*> lsNefs;
  for (auto& p : lsPrimitives)
  {
    Nef_polyhedron* tmpnef;
    if (p->get_type() == SOLID)
      tmpnef = dynamic_cast<Solid*>(p)->get_nef_polyhedron();
    else if (p->get_type() == COMPOSITESOLID)
      tmpnef = dynamic_cast<CompositeSolid*>(p)->get_nef_polyhedron();
    else
      continue;
    if (tol_overlap > 0)
    {
      lsNefs.push_back(erode_nef_polyhedron(tmpnef, tol_overlap));
      // delete tmpnef;
    }
    else
    {
      lsNefs.push_back(tmpnef);
    }
  }
  //-- 2. check whether pairwise intersection of interiors is empty
  Nef_polyhedron emptynef(Nef_polyhedron::EMPTY);
  for (int i = 0; i < lsNefs.size(); i++)
  {
    Nef_polyhedron* a = lsNefs[i];
    for (int j = i + 1; j < lsNefs.size(); j++) 
    {
      Nef_polyhedron* b = lsNefs[j];
      if (a->interior() * b->interior() != emptynef)
      {
        Error e;
        std::stringstream msg;
        msg << lsPrimitives[i]->get_id() << " and " << lsPrimitives[j]->get_id();
        e.errorcode = errorcode_to_assign;
        e.info1 = msg.str();
        e.info2 = "";
        lsErrors.push_back(e);
        isValid = false;
      }
    }
  }
  if (tol_overlap > 0.0)
  {
    for (auto each : lsNefs)
      delete each;
  }
  return !isValid;
}


//-- adjacent 
//-- == 
//-- eroded.interior() not overlapping
//-- + 
//-- dilated.interior() overlapping
//-------------
// -1: input primitives not valid
//  0: not adjancent
//  1: adjacent
int are_primitives_adjacent(Primitive* p1,
                             Primitive* p2,
                             double tol_overlap)
{
  //-- only process valid primitives
  if ( (p1->is_valid() != 1) || (p2->is_valid() != 1) )
    return -1;

  Nef_polyhedron emptynef(Nef_polyhedron::EMPTY);
  Nef_polyhedron* n1;
  Nef_polyhedron* n2;
  if (p1->get_type() == SOLID)
    n1 = dynamic_cast<Solid*>(p1)->get_nef_polyhedron();
  else if (p1->get_type() == COMPOSITESOLID)
    n1 = dynamic_cast<CompositeSolid*>(p1)->get_nef_polyhedron();
  if (p2->get_type() == SOLID)
    n2 = dynamic_cast<Solid*>(p2)->get_nef_polyhedron();
  else if (p2->get_type() == COMPOSITESOLID)
    n2 = dynamic_cast<CompositeSolid*>(p2)->get_nef_polyhedron();

  if (tol_overlap < 0.0)
  {
    if (n1->boundary() * n2->boundary() != emptynef)
      return 1;
    else
      return 0;
  }
  else 
  {
    //-- 1. erode the nefs
    Nef_polyhedron* ne1 = erode_nef_polyhedron(n1, tol_overlap);
    Nef_polyhedron* ne2 = erode_nef_polyhedron(n2, tol_overlap);
    if (ne1->interior() * ne2->interior() != emptynef)
      return 0;
    //-- 2. dilate the Nefs
    Nef_polyhedron* nd1 = dilate_nef_polyhedron(n1, tol_overlap);
    Nef_polyhedron* nd2 = dilate_nef_polyhedron(n2, tol_overlap);
    if (nd1->interior() * nd2->interior() == emptynef)
      return 0;
    return 1;
  }
}



} // namespace val3dity
