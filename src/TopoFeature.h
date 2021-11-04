#ifndef CITYCFD_TOPOFEATURE_H
#define CITYCFD_TOPOFEATURE_H

#include "types.h"
#include "CGALTypes.h"

class TopoFeature {
public:
    TopoFeature();
    TopoFeature(std::string pid);
    TopoFeature(int outputLayerID);
    virtual ~TopoFeature();

    static int           get_num_output_layers();

    virtual void         get_cityjson_info(nlohmann::json& b) const;
    virtual void         get_cityjson_semantics(nlohmann::json& g) const;
    virtual std::string  get_cityjson_primitive() const;
    virtual TopoClass    get_class() const = 0;
    virtual std::string  get_class_name() const = 0;

    Mesh&       get_mesh();
    const Mesh& get_mesh() const;
    void        set_id(unsigned long id);
    std::string get_id() const;
    const int   get_output_layer_id() const;
    bool        is_active() const;
    void        deactivate();
    void test();

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

#endif //CITYCFD_TOPOFEATURE_H