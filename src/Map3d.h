/*
  City4CFD
 
  Copyright (c) 2021-2025, 3D Geoinformation Research Group, TU Delft

  This file is part of City4CFD.

  City4CFD is free software: you can redistribute it and/or modify
  it under the terms of the GNU Affero General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  City4CFD is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Affero General Public License for more details.

  You should have received a copy of the GNU Affero General Public License
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
    Map3d() = default;
    ~Map3d() = default;

    void reconstruct();

    const BuildingsPtr& get_failed_buildings() const;

    void read_data();
    void output();

private:
    PointCloud                  m_pointCloud;
//    JsonVectorPtr               m_polygonsBuildings;     // for json-only reader
    PolyVecPtr                  m_polygonsBuildings;
    JsonVectorPtr               m_importedBuildingsJSON;
//    std::vector<JsonVectorPtr>  m_polygonsSurfaceLayers; // for json-only reader
    std::vector<PolyVecPtr>     m_polygonsSurfaceLayers;
    PointSet3Ptr                m_importedBuildingsPts;
    std::vector<Mesh>           m_importedBuildingsOther;

    TerrainPtr                  m_terrainPtr;
    BuildingsPtr                m_buildingsPtr;
    BuildingsPtr                m_failedBuildingsPtr;
    ReconstructedBuildingsPtr   m_reconstructedBuildingsPtr;
    ImportedBuildingsPtr        m_importedBuildingsPtr;
    SurfaceLayersPtr            m_surfaceLayersPtr;
    BoundariesPtr               m_boundariesPtr;
    PolyFeaturesPtr             m_allFeaturesPtr;
    OutputFeaturesPtr           m_outputFeaturesPtr;

    std::vector<BoundingRegion> m_reconRegions; // one influ region -> vector of reconstruction regions
    BoundingRegion              m_domainBnd;
    DT                          m_dt;

    bool                        m_bndBPG         = false;
    bool                        m_cityjsonInput  = false;

    void set_features();
    void set_influ_region();
    void set_bnd();
    void bnd_sanity_check();
    void add_building_pts();
    void remove_extra_terrain_pts();
    void reconstruct_terrain();
    void reconstruct_buildings();
    void reconstruct_one_building(std::shared_ptr<Building>& building);
    void reconstruct_boundaries();
    void reconstruct_with_flat_terrain();
    void solve_building_conflicts();
    void clip_buildings();
    void wrap();

    void prep_cityjson_output();

    void clear_inactives();

    //-- Templated functions
    template<typename T> void shorten_polygons(T& features);
    template<typename T> void set_footprint_elevation(T& features);
};

#endif //CITY4CFD_MAP3D_H