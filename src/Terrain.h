#ifndef CITYCFD_TERRAIN_H
#define CITYCFD_TERRAIN_H

#include "TopoFeature.h"

class Terrain : public TopoFeature {
public:
    using TopoFeature::TopoFeature;
    Terrain();
    Terrain(int pid);
    ~Terrain();

    void set_cdt(const Point_set_3 &pointCloud);
    void prep_constraints(const PolyFeatures& features, Point_set_3& pointCloud);
    void constrain_features();
    void create_mesh(const PolyFeatures& features);

    CDT&         get_cdt();
    const CDT&   get_cdt() const;

    void         get_cityjson_info(nlohmann::json& b) const override;
    std::string  get_cityjson_primitive() const override;
    TopoClass    get_class() const override;
    std::string  get_class_name() const override;

    const SurfaceLayers& get_surface_layers() const;

protected:
    CDT                    _cdt;
    SurfaceLayers          _surfaceLayersTerrain;
    std::vector<Polygon_3> _constrainedPolys;
};

#endif //CITYCFD_TERRAIN_H