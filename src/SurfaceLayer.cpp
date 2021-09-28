#include "SurfaceLayer.h"

SurfaceLayer::SurfaceLayer()
    : PolyFeature() {}

SurfaceLayer::SurfaceLayer(const int outputLayerID)
    : PolyFeature(outputLayerID) {}

SurfaceLayer::SurfaceLayer(const nlohmann::json& poly, const int outputLayerID)
    : PolyFeature(poly, outputLayerID) {}

SurfaceLayer::~SurfaceLayer() = default;

void SurfaceLayer::check_feature_scope() {
    // TODO: really gotta rewrite those polygons, cgal implementation is awful
    //-- Temporary - will rewrite polygons
    std::vector<Polygon_2> rings;
    //- Add outer poly and holes into one data structure
    rings.push_back(_poly.outer_boundary());
    for (auto& hole : _poly.holes()) {
        rings.push_back(hole);
    }

    //-- Exclude all polygons that have at least one vertex outside the domain
    for (auto& poly : rings) {
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

void SurfaceLayer::threeDfy(CDT& cdt) {
//    geomtools::cdt_to_mesh(cdt, _mesh, this->get_layer_id());
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
    return "SemanticLayer";
}