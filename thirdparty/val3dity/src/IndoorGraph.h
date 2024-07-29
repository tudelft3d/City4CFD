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

#ifndef __val3dity__IndoorGraph__
#define __val3dity__IndoorGraph__

#include "Feature.h"
#include "definitions.h"


namespace val3dity
{

class IndoorGraph
{
public:
  IndoorGraph(std::string theid = "");
  ~IndoorGraph();

  void            add_vertex(std::string theid, 
                             double x, 
                             double y, 
                             double z, 
                             std::string dual, 
                             std::vector<std::string> vadj);
  
  bool            validate();
  bool            is_valid();
  std::string     get_type();

  bool            has_vertex(std::string);
  std::tuple<Point3,std::string,std::vector<std::string>>& 
                  get_vertex(std::string id);
  std::vector<std::string> 
                  get_vertices_ids();

protected:
  std::string     _id;
  int             _is_valid; 
  
  // _vertices[theid] = std::make_tuple(v, dual, vadj);
  std::map<std::string, std::tuple<Point3,std::string,std::vector<std::string>>> _vertices;

};

} // namespace val3dity

#endif /* defined(__val3dity__IndoorGraph__) */
