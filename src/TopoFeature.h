#ifndef CITYCFD_TOPOFEATURE_H
#define CITYCFD_TOPOFEATURE_H

#include "config.h"
#include "geomtools.h"

class TopoFeature {
public:
    TopoFeature();
    TopoFeature(std::string pid);
    TopoFeature(int outputLayerID);
    virtual ~TopoFeature();

    virtual TopoClass    get_class() const = 0;
    virtual std::string  get_class_name() const = 0;
    virtual void         get_cityjson_info(nlohmann::json& b) const;
    virtual void         get_cityjson_semantics(nlohmann::json& g) const;
    virtual std::string  get_cityjson_primitive() const;

    static int           get_num_output_layers();

    Mesh&       get_mesh();
    const Mesh& get_mesh() const;
    void        set_id(unsigned long id);
    std::string get_id() const;
    const int   get_output_layer_id() const;
    bool        is_active() const;
    void        deactivate();

protected:
    static int _numOfOutputLayers;

    Mesh           _mesh;
    std::string    _id;
    bool           _f_active;
    int            _outputLayerID; // 0- Terrain
    // 1- Buildings
    // 2- Sides
    // 3- Top
    // 4 Onwards - surface layers
};

//-- TopoFeature derived from polygons
class PolyFeature : public TopoFeature {
public:
    PolyFeature();
    PolyFeature(const int outputLayerID);
    PolyFeature(const nlohmann::json& poly);
    PolyFeature(const nlohmann::json& poly, const int outputLayerID);
    virtual ~PolyFeature();

    virtual void        check_feature_scope();
    virtual void        get_cityjson_info(nlohmann::json& b) const = 0;
    virtual void        get_cityjson_semantics(nlohmann::json& g) const = 0;
    virtual std::string get_cityjson_primitive() const = 0;

    virtual TopoClass   get_class() const = 0;
    virtual std::string get_class_name() const = 0;

    void  calc_footprint_elevation_nni(const DT& dt);
    void  calc_footprint_elevation_linear(const DT& dt);
    void  calc_footprint_elevation_from_pc(const SearchTree& searchTree);
    void  clear_base_heights();

    Polygon_with_holes_2&                    get_poly();
    const Polygon_with_holes_2&              get_poly() const;
    const std::vector<std::vector<double>>&  get_base_heights() const;

protected:
    Polygon_with_holes_2              _poly;
    std::vector<std::vector<double>>  _base_heights;
};

#endif //CITYCFD_TOPOFEATURE_H