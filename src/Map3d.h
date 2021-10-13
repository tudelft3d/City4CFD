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

    void set_features();
    void set_boundaries();
    void triangulate_terrain();
    void polygon_processing();
    void set_footprint_elevation();
    void threeDfy();
    void prep_feature_output();
    void prep_cityjson_output();
    void collect_garbage();
};

#endif //CITYCFD_MAP3D_H