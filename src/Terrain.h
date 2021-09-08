#ifndef CITYCFD_TERRAIN_H
#define CITYCFD_TERRAIN_H

#include "definitions.h"
#include "geomtools.h"
#include "io.h"
#include "TopoFeature.h"

class Terrain : public TopoFeature {
public:
    using TopoFeature::TopoFeature;
    Terrain();
    Terrain(int pid);
    ~Terrain() = default;

    void threeDfy(const Point_set_3& pointCloud, const std::vector<PolyFeature*>& features);

    const CDT&   get_cdt() const;
    std::string  get_class_name() const;

    TopoClass    get_class() const override;

protected:
    CDT _cdt;

    void set_cdt(const Point_set_3 &pointCloud);
    void constrain_footprint(const Polygon_with_holes_2& poly, const std::vector<double>& heights);
    void smooth(const Point_set_3& pointCloud);
    void create_mesh();
};

#endif //CITYCFD_TERRAIN_H
