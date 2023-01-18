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

#ifndef CITY4CFD_MAP3D_H
#define CITY4CFD_MAP3D_H

#include "types.h"
#include "CGALTypes.h"
#include "BoundingRegion.h"
#include "PointCloud.h"

class Map3d {
public:
    Map3d();
    ~Map3d();

    void reconstruct();

    void read_data();
    void output();

private:
    PointCloud                  _pointCloud;
    JsonVectorPtr               _polygonsBuildings;
    JsonVectorPtr               _importedBuildingsJSON;
    std::vector<JsonVectorPtr>  _polygonsSurfaceLayers;
    PointSet3Ptr                _importedBuildingsPts;
    std::vector<Mesh>           _importedBuildingsOther;

    TerrainPtr                  _terrainPtr;
    BuildingsPtr                _buildingsPtr;
    ReconstructedBuildingsPtr   _reconstructedBuildingsPtr;
    ImportedBuildingsPtr        _importedBuildingsPtr;
    SurfaceLayersPtr            _surfaceLayersPtr;
    BoundariesPtr               _boundariesPtr;
    PolyFeaturesPtr             _allFeaturesPtr;
    OutputFeaturesPtr           _outputFeaturesPtr;

    BoundingRegion              _influRegion;
    BoundingRegion              _domainBnd;
    DT                          _dt;

    bool                        _influRegionBPG = false;
    bool                        _bndBPG         = false;
    bool                        _cityjsonInput  = false;

    void set_features();
    void set_influ_region();
    void set_bnd();
    void bnd_sanity_check();
    void add_building_pts();
    void remove_extra_terrain_pts();
    void reconstruct_terrain();
    void reconstruct_buildings();
    void reconstruct_boundaries();
    void reconstruct_with_flat_terrain();
    void solve_building_conflicts();
    void clip_buildings();
    void wrap();

    void prep_feature_output();
    void prep_cityjson_output();

    void clear_inactives();

    //-- Templated functions
    template<typename T> void shorten_polygons(T& features);
    template<typename T> void set_footprint_elevation(T& features);
};

#endif //CITY4CFD_MAP3D_H