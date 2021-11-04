#include "SurfaceLayer.h"

#include "geomutils.h"

SurfaceLayer::SurfaceLayer()
    : PolyFeature() {}

SurfaceLayer::SurfaceLayer(const int outputLayerID)
    : PolyFeature(outputLayerID) {}

SurfaceLayer::SurfaceLayer(const nlohmann::json& poly)
        : PolyFeature(poly) {}

SurfaceLayer::SurfaceLayer(const nlohmann::json& poly, const int outputLayerID)
    : PolyFeature(poly, outputLayerID) {}

SurfaceLayer::~SurfaceLayer() = default;

void SurfaceLayer::check_feature_scope(const Polygon_2& bndPoly) {
    //-- Exclude all polygons that have at least one
    //-- vertex outside the domain
        for (auto& vert : _poly.outer_boundary()) {
            if (!geomutils::point_in_poly(vert, bndPoly)) {
                this->deactivate();
                return;
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