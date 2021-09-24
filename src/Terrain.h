#ifndef CITYCFD_TERRAIN_H
#define CITYCFD_TERRAIN_H

#include "definitions.h"
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

    void threeDfy(const Point_set_3& pointCloud, const std::vector<PolyFeature*>& features);

    CDT&         get_cdt();
    const CDT&   get_cdt() const;

    void         get_cityjson_info(nlohmann::json& b) const override;
    std::string  get_cityjson_primitive() const override;
    TopoClass    get_class() const override;
    std::string  get_class_name() const override;
    void         constrain_footprint(const Polygon_with_holes_2& poly, const std::vector<double>& heights); //testing

    const std::vector<SurfaceLayer*>& get_surface_layers() const;

protected:
    CDT  _cdt;
    std::vector<SurfaceLayer*> _surfaceLayersTerrain; // Vector or list for the mesh?

    void set_cdt(const Point_set_3 &pointCloud);
    void smooth(const Point_set_3& pointCloud);
    void create_mesh();
};

#endif //CITYCFD_TERRAIN_H
