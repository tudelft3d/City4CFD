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
    Building(const int pid);
    Building(const json& poly, const int pid);
    ~Building() = default;

    void      calc_footprint_elevation(const SearchTree& searchTree) override;
    void      threeDfy(const SearchTree& searchTree) override;

    TopoClass   get_class() const override;
    std::string get_class_name() const override;

private:
    double _height;
};

#endif //CITYCFD_BUILDING_H