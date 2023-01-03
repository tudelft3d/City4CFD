/*
  City4CFD
 
  Copyright (c) 2021-2023, 3D Geoinformation Research Group, TU Delft

  This file is part of City4CFD.

  City4CFD is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  City4CFD is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with City4CFD.  If not, see <http://www.gnu.org/licenses/>.

  For any information or further details about the use of City4CFD, contact
  Ivan Pađen
  <i.paden@tudelft.nl>
  3D Geoinformation Research Group
  Delft University of Technology
*/

#ifndef CITY4CFD_TYPES_H
#define CITY4CFD_TYPES_H

#include "nlohmann/json.hpp"

//-- Typedefs for smart pointers
class Building;    class Boundary; class TopoFeature;
class PolyFeature; class Terrain;  class SurfaceLayer;
class ReconstructedBuilding; class ImportedBuilding;
typedef std::shared_ptr<Terrain>                               TerrainPtr;
typedef std::vector<std::shared_ptr<Building>>                 BuildingsPtr;
typedef std::vector<std::shared_ptr<ReconstructedBuilding>>    ReconstructedBuildingsPtr;
typedef std::vector<std::shared_ptr<ImportedBuilding>>         ImportedBuildingsPtr;
typedef std::vector<std::shared_ptr<Boundary>>                 BoundariesPtr;
typedef std::vector<std::shared_ptr<TopoFeature>>              OutputFeaturesPtr;
typedef std::vector<std::shared_ptr<PolyFeature>>              PolyFeaturesPtr;
typedef std::vector<std::shared_ptr<SurfaceLayer>>             SurfaceLayersPtr;
typedef std::vector<std::unique_ptr<nlohmann::json>>           JsonVectorPtr;

//-- TopoClasses
typedef enum {
    TERRAIN          = 0,
    BUILDING         = 1,
    SIDES            = 2,
    TOP              = 3,
    SURFACELAYER     = 4,
} TopoClass;

//-- Domain types
typedef enum {
    ROUND        = 0,
    RECTANGLE    = 1,
    OVAL         = 2
} DomainType;

//-- Output Formats
typedef enum {
    OBJ          = 0,
    CityJSON     = 1,
    STL          = 2,
} GeomFormat;

//-- Global Constants
namespace global {
    const double largnum  = 1e7;
    const double smallnum = 1e-7;
}

#endif //CITY4CFD_TYPES_H