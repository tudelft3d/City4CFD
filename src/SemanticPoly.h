#ifndef CITYCFD_SEMANTICPOLY_H
#define CITYCFD_SEMANTICPOLY_H

#include "definitions.h"
#include "config.h"
#include "TopoFeature.h"

class SemanticPoly : public PolyFeature {
public:
    using PolyFeature::PolyFeature;
    SemanticPoly();
    SemanticPoly(const nlohmann::json& poly, int semanticLayerID);
    ~SemanticPoly();

    virtual void        calc_footprint_elevation(const SearchTree& searchTree) override;
    virtual void        threeDfy(const SearchTree& searchTree) override;
    virtual void        get_cityjson_info(nlohmann::json& b) const override;
    virtual void        get_cityjson_semantics(nlohmann::json& g) const override;
    virtual std::string get_cityjson_primitive() const override;

    virtual TopoClass   get_class() const override;
    virtual std::string get_class_name() const override;

private:
    int    _semanticLayerID;

};


#endif //CITYCFD_SEMANTICPOLY_H