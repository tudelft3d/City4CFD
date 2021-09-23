#ifndef CITYCFD_SURFACELAYER_H
#define CITYCFD_SURFACELAYER_H

#include "definitions.h"
#include "config.h"
#include "TopoFeature.h"

class SurfaceLayer : public PolyFeature {
public:
    using PolyFeature::PolyFeature;
    SurfaceLayer();
    SurfaceLayer(const nlohmann::json& poly, int semanticLayerID);
    ~SurfaceLayer();

    virtual void        check_feature_scope() override;
    virtual void        get_cityjson_info(nlohmann::json& b) const override;
    virtual void        get_cityjson_semantics(nlohmann::json& g) const override;
    virtual std::string get_cityjson_primitive() const override;

    virtual TopoClass   get_class() const override;
    virtual std::string get_class_name() const override;

    void threeDfy(CDT& cdt) ;

private:

};


#endif //CITYCFD_SURFACELAYER_H