/*
  val3dity 

  Copyright (c) 2011-2020, 3D geoinformation research group, TU Delft

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

#ifndef __val3dity__IndoorModel__
#define __val3dity__IndoorModel__

#include "Feature.h"
#include "definitions.h"
#include "IndoorGraph.h"


namespace val3dity
{

class IndoorModel : public Feature
{
public:
  IndoorModel(std::string theid);
  ~IndoorModel();
  
  bool            validate(double tol_planarity_d2p, double tol_planarity_normals, double tol_overlap = -1);
  bool            is_valid();
  std::string     get_type();

  void            add_cell(std::string id, Primitive* p, std::string dual, std::string CellSpaceType);
  void            add_graph(IndoorGraph* g);

  int             get_number_graphs();
  bool            graph_has_dual_vertex(std::string theid);
  int             get_graph_containing_vertex(std::string theid);

protected:
  //-- <id> // <pos_lsPrimitives, dual, CellSpaceType> 
  std::map<std::string, std::tuple<int,std::string,std::string>> _cells;  
  std::vector<IndoorGraph*>                                      _graphs;

};

} // namespace val3dity

#endif /* defined(__val3dity__IndoorModel__) */
