#ifndef CITYCFD_MAP3D_H
#define CITYCFD_MAP3D_H

#include "config.h"
#include "io.h"
#include "TopoFeature.h"
#include "Terrain.h"
#include "Building.h"
#include "Boundary.h"
#include "SurfaceLayer.h"

class Map3d {
public:
    Map3d();
    ~Map3d();

    void reconstruct();

    void read_data();
    void output();

private:
    Point_set_3                _pointCloud;
    Point_set_3                _pointCloudBuildings;
    JsonPolygons               _polygonsBuildings;
    std::vector<JsonPolygons>  _polygonsSurfaceLayers;
    Terrainptr                 _terrain;
    Buildings                  _buildings;
    SurfaceLayers              _surfaceLayers;
    Boundaries                 _boundaries;
    PolyFeatures               _lsFeatures;
    OutputFeatures             _outputFeatures;
    BoundingRegion             _influRegion;
    BoundingRegion             _domainBnd;
    DT                         _dt;
    bool                       _influRegionBPG = false;
    bool                       _bndBPG = false;

    void set_features();
    void set_influ_region();
    void set_bnd();
    void bnd_sanity_check();
    void reconstruct_terrain();
    void reconstruct_buildings();
    void reconstruct_boundaries();

    void prep_feature_output();
    void prep_cityjson_output();

    void clear_inactives();

    //-- Templated functions
    template<typename T> void shorten_polygons(T& features);
    template<typename T> void set_footprint_elevation(T& features);
};

#endif //CITYCFD_MAP3D_H