#ifndef CITY4CFD_SURFACELAYER_H
#define CITY4CFD_SURFACELAYER_H

#include "PolyFeature.h"

class SurfaceLayer : public PolyFeature {
public:
    SurfaceLayer();
    SurfaceLayer(const int outputLayerID);
    SurfaceLayer(const nlohmann::json& poly);
    SurfaceLayer(const nlohmann::json& poly, const int outputLayerID);
    ~SurfaceLayer();

    void check_feature_scope(const Polygon_2& bndPoly);

    virtual void        get_cityjson_info(nlohmann::json& b) const override;
    virtual void        get_cityjson_semantics(nlohmann::json& g) const override;
    virtual std::string get_cityjson_primitive() const override;
    virtual TopoClass   get_class() const override;
    virtual std::string get_class_name() const override;

private:

};

#endif //CITY4CFD_SURFACELAYER_H