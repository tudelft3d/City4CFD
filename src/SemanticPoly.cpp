#include "SemanticPoly.h"

SemanticPoly::SemanticPoly() = default;

SemanticPoly::SemanticPoly(const nlohmann::json& poly, int semanticLayerID)
    : PolyFeature(poly), _semanticLayerID(semanticLayerID) {}

SemanticPoly::~SemanticPoly() = default;

void SemanticPoly::calc_footprint_elevation(const SearchTree& searchTree) {

}

void SemanticPoly::threeDfy(const SearchTree& searchTree) {

}

void SemanticPoly::get_cityjson_info(nlohmann::json& b) const {

}

void SemanticPoly::get_cityjson_semantics(nlohmann::json& g) const {

}

std::string SemanticPoly::get_cityjson_primitive() const {
    return "Dunno yet";
}

TopoClass SemanticPoly::get_class() const {
    return SEMANTICLAYER;
}

std::string SemanticPoly::get_class_name() const {
    return "SemanticLayer";
}