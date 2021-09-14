#ifndef CITYCFD_BUILDING_H
#define CITYCFD_BUILDING_H

#include "definitions.h"
#include "geomtools.h"
#include "io.h"
#include "TopoFeature.h"
#include "LoD12.h"

class Building : public PolyFeature {
public:
    using PolyFeature::PolyFeature;
    Building() = default;
    Building(const json& poly);
    ~Building() = default;

    void        calc_footprint_elevation(const SearchTree& searchTree) override;
    void        threeDfy(const SearchTree& searchTree) override;
    void        get_cityjson_info(nlohmann::json& b) const override;
    void        get_cityjson_semantics(nlohmann::json& g) const override;
    std::string get_cityjson_primitive() const override;

    TopoClass   get_class() const override;
    std::string get_class_name() const override;

private:
    double _height;
};

#endif //CITYCFD_BUILDING_H