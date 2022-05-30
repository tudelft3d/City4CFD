/*
  Copyright (c) 2021-2022,
  Ivan Pađen <i.paden@tudelft.nl>
  3D Geoinformation,
  Delft University of Technology

  This file is part of City4CFD.

  City4CFD is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  City4CFD is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

  You should have received a copy of the GNU General Public License
  along with this program. If not, see <http://www.gnu.org/licenses/>
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
    JsonVector                  _polygonsBuildings;
    JsonVector                  _importedBuildingsJson;
    std::vector<JsonVector>     _polygonsSurfaceLayers;
    std::vector<Point_3>        _importedBuildingsPts;

    Terrainptr                  _terrain;
    Buildings                   _buildings;
    ReconstructedBuildings      _reconstructedBuildings;
    ImportedBuildings           _importedBuildings;
    SurfaceLayers               _surfaceLayers;
    Boundaries                  _boundaries;
    PolyFeatures                _lsFeatures;
    OutputFeatures              _outputFeatures;

    BoundingRegion              _influRegion;
    BoundingRegion              _domainBnd;
    DT                          _dt;

    bool                        _influRegionBPG = false;
    bool                        _bndBPG = false;

    void set_features();
    void set_influ_region();
    void set_bnd();
    void bnd_sanity_check();
    void reconstruct_terrain();
    void reconstruct_buildings();
    void reconstruct_boundaries();
    void reconstruct_with_flat_terrain();
    void solve_building_conflicts();
    void clip_buildings();

    void prep_feature_output();
    void prep_cityjson_output();

    void clear_inactives();

    //-- Templated functions
    template<typename T> void shorten_polygons(T& features);
    template<typename T> void set_footprint_elevation(T& features);

    // Testing
    void one_mesh();
};

#endif //CITY4CFD_MAP3D_H