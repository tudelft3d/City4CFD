#include "SurfaceLayer.h"

SurfaceLayer::SurfaceLayer()
    : PolyFeature() {}

SurfaceLayer::SurfaceLayer(const int outputLayerID)
    : PolyFeature(outputLayerID) {}

SurfaceLayer::SurfaceLayer(const nlohmann::json& poly, const int outputLayerID)
    : PolyFeature(poly, outputLayerID) {}

SurfaceLayer::~SurfaceLayer() = default;

void SurfaceLayer::check_feature_scope() { //-- Overrides default, as it doesn't depend on influence region
    //-- Exclude all polygons that have at least one vertex outside the domain
    for (auto& poly : _poly.rings()) {
        for (auto& vert : poly) {
            if (pow(config::pointOfInterest.x() - vert.x(), 2)
                + pow(config::pointOfInterest.y() - vert.y(), 2)
                > pow(0.8 * config::dimOfDomain,2)) {
                this->deactivate();
                return;
            }
        }
    }
}

void SurfaceLayer::get_cityjson_info(nlohmann::json& b) const {

}

void SurfaceLayer::get_cityjson_semantics(nlohmann::json& g) const {

}

std::string SurfaceLayer::get_cityjson_primitive() const {
    return "Dunno yet";
}

TopoClass SurfaceLayer::get_class() const {
    return SURFACELAYER;
}

std::string SurfaceLayer::get_class_name() const {
    return "SurfaceLayer";
}