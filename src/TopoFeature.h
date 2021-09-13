#ifndef CITYCFD_TOPOFEATURE_H
#define CITYCFD_TOPOFEATURE_H

#include "definitions.h"
#include "config.h"
#include "geomtools.h"

class TopoFeature {
public:
    TopoFeature();
    TopoFeature(const int pid);
    ~TopoFeature()                   = default;

    virtual TopoClass   get_class() const = 0;
    virtual std::string get_class_name() const = 0;
    virtual void get_cityjson_info(nlohmann::json& b);
    virtual std::string get_cityjson_primitive() const;

    Mesh&       get_mesh();
    const Mesh& get_mesh() const;
    void        set_id(unsigned long id);
    std::string get_id() const;
    bool        is_active() const;
    void        deactivate();

protected:
    Mesh           _mesh;
    std::string    _id;
    bool           _f_active;
};

//-- TopoFeature derived from polygons
class PolyFeature : public TopoFeature {
public:
    using TopoFeature::TopoFeature;
    PolyFeature() = default;
    PolyFeature(const json& poly);
    ~PolyFeature() = default;

    virtual TopoClass   get_class() const = 0;
    virtual std::string get_class_name() const = 0;
    virtual void        calc_footprint_elevation(const SearchTree& searchTree);
    virtual void        threeDfy(const SearchTree& searchTree) = 0;
    virtual void        get_cityjson_info(nlohmann::json& b) = 0;
    virtual std::string get_cityjson_primitive() const = 0;

    void check_influ_region();

    const Polygon_with_holes_2& get_poly() const;
    const std::vector<double>&  get_base_heights() const;

protected:
    Polygon_with_holes_2               _poly;
    std::vector<double>                _base_heights;

};

#endif //CITYCFD_TOPOFEATURE_H