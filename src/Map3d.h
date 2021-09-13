#ifndef CITYCFD_MAP3D_H
#define CITYCFD_MAP3D_H

#include "definitions.h"
#include "config.h"
#include "io.h"
#include "TopoFeature.h"
#include "Terrain.h"
#include "Building.h"
#include "Boundary.h"

class Map3d {
public:
    Map3d()  = default;
    ~Map3d();

    void reconstruct();

    bool read_config(const char* points_xyz);
    bool read_point_cloud(const char* points_xyz);
    bool read_point_cloud_buildings(const char* points_xyz);
    bool read_polygons(const char* gisdata);

    void output();

private:
     Point_set_3                _pointCloud;
     Point_set_3                _pointCloudBuildings;
     json                       _polygons;
     Terrain*                   _terrain;
     std::vector<Boundary*>     _boundaries;
     std::vector<PolyFeature*>  _lsFeatures;
     std::vector<TopoFeature*>  _allFeatures;

    void set_features();
    void set_boundaries();
    void set_footprint_elevation();
    void threeDfy();
    void prep_feature_output();
    void collect_garbage();
    void clear_features();
};

#endif //CITYCFD_MAP3D_H