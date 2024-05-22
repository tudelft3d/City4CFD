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


#include "IndoorModel.h"
#include "definitions.h"
#include "input.h"
#include "Solid.h"
#include "validate_prim_toporel.h"

namespace val3dity
{

IndoorModel::IndoorModel(std::string theid)
{
  _id = theid;
  _type = "IndoorModel"; 
  _is_valid = -1;
}


IndoorModel::~IndoorModel()
{}


bool IndoorModel::validate(double tol_planarity_d2p, double tol_planarity_normals, double tol_overlap) 
{
  // 
  // 1. each Cell is valid Solid
  // 2. do Cells overlap each others? https://doc.cgal.org/latest/Polygon_mesh_processing/group__PMP__predicates__grp.html#ga1ff63ec6e762d45ea5775bf7b49f9270
  // 3. is dual vertex of each Cell locate inside the Cell
  // 4. are graphs valid: c1->v4 && v4->c1 (faulty xlinks?)
  // 5. adj in dual == adj in primal
  // 

  // std::clog << "--- Validation of IndoorModel ---" << std::endl;
  if (_is_valid != -1)
    return _is_valid;

  bool bValid = true;
//-- 1. 4xx - ISO19107 check for Solid validity
//--    validate each IndoorCell geometry (Solids)
  bValid = Feature::validate_generic(tol_planarity_d2p, tol_planarity_normals, tol_overlap);

//-- 2. 701 - CELLS_OVERLAP
//--    overlapping test
  // std::clog << "--- Overlapping tests between Cells ---" << std::endl;
  std::vector<std::tuple<std::string,Solid*>> lsCells;
  for (auto& el : _cells)
    lsCells.push_back(std::make_tuple(el.first, (Solid*)_lsPrimitives[std::get<0>(el.second)]));
  std::vector<Error> lsErrors;  
  if (are_cells_interior_disconnected_with_aabb(lsCells, 701, lsErrors, tol_overlap) == false)
  {
    bValid = false;
    // std::clog << "Error: Cells have overlapping interior" << std::endl;
    for (auto& e : lsErrors)
      this->add_error(e.errorcode, e.info1, e.info2);
  }

//-- 3. 702 - DUAL_VERTEX_OUTSIDE_CELL
//--    is dual vertex of each cell located inside its Cell?
  // std::clog << "======== Validating Dual Vertex (Point-in-Solid tests) ========" << std::endl;
  for (auto& el : _cells)
  {
    // std::clog << "Cell (" << std::get<2>(el.second) << ") id=" << el.first;
    // TODO : what to do with many graphs here?
    //-- check if there's a dual, something there's not (and it's valid)
    if (std::get<1>(el.second) == "") 
    {
      // std::clog << " has no dual vertex" << std::endl;
    }
    else 
    {
      //-- does the dual graph contain that vertex ID?
      if (this->get_number_graphs() == 0) {
          std::stringstream msg;
          msg << "CellSpace id=" << el.first;
          this->add_error(702, msg.str(), "The dual graph is not in the file.");
          bValid = false;
          continue;
      }
      int gno = this->get_graph_containing_vertex(std::get<1>(el.second));
      if (gno == -1)
      {
        // std::stringstream msg;
        // msg << "Cell (" << std::get<2>(el.second) << ") id=" << el.first << " dual doesn't exist";
        // this->add_error(704, msg.str(), "");
        // std::clog << " ==> dual does not exist" << std::endl;
      } 
      else 
      {
        Point3 p = std::get<0>(_graphs[gno]->get_vertex(std::get<1>(el.second)));
        Solid* s = (Solid*)_lsPrimitives[std::get<0>(el.second)];
        int inside = s->is_point_in_solid(p);
        if (inside == -1)
        {
          std::stringstream msg;
          msg << "CellSpace id=" << el.first;
          this->add_error(702, msg.str(), "");
          bValid = false;
        }
        // else
          // std::clog << " --> ok" << std::endl;
      }
    }
  }

//-- 4. 703 - PRIMAL_DUAL_XLINKS_ERROR
//--    are primal-dual graphs valid
//--    this validates the XLinks basically, which is not done by XSD
  // std::clog << "======== Validating Primal-Dual links ========" << std::endl;
  for (auto& el : _cells)
  {
    std::string pdid = std::get<1>(el.second);
    //-- no dual to the cell
    if (pdid == "") 
    {
      //-- this is allowed, a cell does not need to have a dual vertex
      // std::clog << "Cell id=" << el.first << " has no dual vertex" << std::endl;
      continue;
    }
    if (this->get_number_graphs() == 0) {
      std::stringstream msg;
      msg << "CellSpace id=" << el.first;
      this->add_error(703, msg.str(), "The dual graph is not in the file.");
      bValid = false;
      continue;
    }
    int gno = this->get_graph_containing_vertex(pdid);
    if (gno == -1) {
        std::stringstream msg;
        msg << "Cell id=" << el.first << " dual vertex doesn't exist";
        bValid = false;
        this->add_error(703, msg.str(), "");
    }
    else
    {
      //-- check 'reciprocity' of the XLinks between primal and dual
      if (std::get<1>(_graphs[gno]->get_vertex(pdid)) != el.first)
      {
        std::stringstream msg;
        msg << "Cell id=" << el.first << " and its dual vertex id=" << pdid << " are not reciprocally linked";
        bValid = false;
        this->add_error(703, msg.str(), "");
      }
    }
  }
  //-- test if the neighbours of the dual vertices actually exist
  //-- an IndoorGML dual node (a "State") has a 'connects', these should exists
  for (auto& g: _graphs) {
    auto allids = g->get_vertices_ids();
    for (auto &vid : allids) {
      auto vadjs = g->get_vertex(vid);
      for (auto &vadj: std::get<2>(vadjs)) {
        if (vadj == "") {
          std::stringstream msg;
          msg << "Dual vertex id=" << vid << " reference a non-existing edge/Transition";
          bValid = false;
          this->add_error(703, msg.str(), "");
        }
        else if (g->has_vertex(vadj) == false) {
          std::stringstream msg;
          msg << "Dual vertex id=" << vadj << " does not exist (referenced from id=" << vid;
          bValid = false;
          this->add_error(703, msg.str(), "");
        }
      }
    }
  }

//-- 5. 704 - PRIMAL_DUAL_ADJACENCIES_INCONSISTENT  
//--    if 2 cells are adjacent in the primal, are they also in the dual?
  // std::clog << "======== Validating Primal-Dual links ========" << std::endl;
  for (auto& el : _cells)
  {
    std::string pdid = std::get<1>(el.second);
    if (pdid == "")
    {
      //-- this is allowed, a cell does not need to have a dual vertex
      // std::clog << "Cell id=" << el.first << " has no dual vertex" << std::endl;
      continue;
    }
    if (this->get_number_graphs() == 0) {
      std::stringstream msg;
      msg << "Cell id=" << el.first << " dual vertex doesn't exist";
      bValid = false;
      this->add_error(704, msg.str(), "");
      continue;
    }
    for (auto& g: _graphs) {
      for (auto &vadj: std::get<2>(g->get_vertex(pdid))) {
        std::string cadjid = std::get<1>(g->get_vertex(vadj));
        if (el.first > cadjid)
          continue;
        // std::clog << "Cells id=" << el.first << " id=" << cadjid;
        int re = are_primitives_adjacent(_lsPrimitives[std::get<0>(el.second)],
                                         _lsPrimitives[std::get<0>(_cells[cadjid])],
                                         tol_overlap);
        if (re == 0) {
          std::stringstream msg;
          bValid = false;
          msg << "Cells id=" << el.first << " & id=" << cadjid;
          this->add_error(704, msg.str(), "");
        }
      }
    }
  }
      
//-- bye-bye
  _is_valid = bValid;
  return bValid;
}


void IndoorModel::add_cell(std::string id, Primitive* p, std::string dual, std::string CellSpaceType)
{
  this->add_primitive(p);
  int pos = (_lsPrimitives.size() - 1);
  _cells[id] = std::make_tuple(pos, dual, CellSpaceType);
}


void IndoorModel::add_graph(IndoorGraph* g)
{
  _graphs.push_back(g);
}


bool IndoorModel::is_valid() {
  if (_is_valid == 1) 
    return true;
  else
    return false;
}


std::string IndoorModel::get_type() {
  return _type;
}


int IndoorModel::get_number_graphs() {
  return _graphs.size();
}

bool IndoorModel::graph_has_dual_vertex(std::string theid) {
  for (auto& g : _graphs)
    if (g->has_vertex(theid) == true)
      return true;
  return false;
}


int IndoorModel::get_graph_containing_vertex(std::string theid) {
  int i = 0;
  for (auto& g : _graphs) {
    if (g->has_vertex(theid) == true) {
      return i;
    }
    i++;
  }
  return -1;
}


} // namespace val3dity
