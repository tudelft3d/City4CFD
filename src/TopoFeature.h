#ifndef CITYCFD_TOPOFEATURE_H
#define CITYCFD_TOPOFEATURE_H

#include "definitions.h"
#include "config.h"
#include "geomtools.h"

class TopoFeature {
public:
    TopoFeature()                    = default;
    TopoFeature(const int pid) ;
    ~TopoFeature()                   = default;

    virtual TopoClass get_class() const = 0;

    Mesh&       get_mesh();
    const Mesh& get_mesh() const;
    int         get_id() const;
    bool        is_active() const;
    void        deactivate();

protected:
    Mesh       _mesh;
    int        _id;
    bool       _f_active;
};

//-- TopoFeature derived from polygons
class PolyFeature : public TopoFeature {
public:
    using TopoFeature::TopoFeature;
    PolyFeature() = default;
    PolyFeature(const int pid);
    PolyFeature(const json& poly, const int pid);
    ~PolyFeature() = default;

    virtual void      calc_footprint_elevation(const SearchTree& searchTree);
    virtual void      threeDfy(const SearchTree& searchTree) = 0;

    void check_influ_region();

    const Polygon_with_holes_2& get_poly() const;
    const std::vector<double>&  get_base_heights() const;

protected:
    Polygon_with_holes_2               _poly;
    std::vector<double>                _base_heights;

};

#endif //CITYCFD_TOPOFEATURE_H