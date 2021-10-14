#ifndef CITYCFD_TERRAIN_H
#define CITYCFD_TERRAIN_H

#include "config.h"
#include "geomtools.h"
#include "io.h"
#include "TopoFeature.h"
#include "SurfaceLayer.h"

class Terrain : public TopoFeature {
public:
    using TopoFeature::TopoFeature;
    Terrain();
    Terrain(int pid);
    ~Terrain();

    void set_cdt(const Point_set_3 &pointCloud);
    void create_mesh(const PolyFeatures& features);

    CDT&         get_cdt();
    const CDT&   get_cdt() const;

    void         get_cityjson_info(nlohmann::json& b) const override;
    std::string  get_cityjson_primitive() const override;
    TopoClass    get_class() const override;
    std::string  get_class_name() const override;

    const SurfaceLayers& get_surface_layers() const;

    //-- Templated functions
    template<typename T> void constrain_features(const T& features);

protected:
    CDT  _cdt;
    SurfaceLayers _surfaceLayersTerrain;
};

#endif //CITYCFD_TERRAIN_H