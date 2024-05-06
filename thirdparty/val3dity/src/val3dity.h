/*
  val3dity 

  Copyright (c) 2011-2023, 3D geoinformation research group, TU Delft

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
  Faculty of the Built Environment & Architecture
  Delft University of Technology
  Julianalaan 134, Delft 2628BL, the Netherlands
*/

#ifndef val3dity_h
#define val3dity_h

#include "nlohmann-json/json.hpp"

using json = nlohmann::json;

namespace val3dity
{

struct Parameters
{
  double _tol_snap = 0.001;
  double _planarity_d2p_tol = 0.01;
  double _planarity_n_tol = 20.0;
  double _overlap_tol = -1.0;
  bool _terminal_output = true;

  Parameters& tol_snap(double tol_snap)
  {
    _tol_snap = tol_snap;
    return *this;
  }
  Parameters& planarity_d2p_tol(double planarity_d2p_tol)
  {
    _planarity_d2p_tol = planarity_d2p_tol;
    return *this;
  }
  Parameters& planarity_n_tol(double planarity_n_tol)
  {
    _planarity_n_tol = planarity_n_tol;
    return *this;
  }
  Parameters& overlap_tol(double overlap_tol)
  {
    _overlap_tol = overlap_tol;
    return *this;
  }
  Parameters& terminal_output(bool terminal_output)
  {
    _terminal_output = terminal_output;
    return *this;
  }
};

bool
is_valid(json& j,
         Parameters params = Parameters());

json
validate(json& j,
         Parameters params = Parameters());

bool
is_valid(std::string& input,
         std::string format,
         Parameters params = Parameters());

json
validate(std::string& input,
         std::string format,
         Parameters params = Parameters());

json
validate(const std::vector<std::array<double, 3>>& vertices,
         const std::vector<std::vector<int>>& triangle_ids,
         Parameters params = Parameters());

bool
is_valid(const std::vector<std::array<double, 3>>& vertices,
         const std::vector<std::vector<int>>& triangle_ids,
         Parameters params = Parameters());

} // namespace val3dity

#endif /* val3dity_h */