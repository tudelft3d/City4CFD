#include "SemanticPoly.h"

SemanticPoly::SemanticPoly() = default;

SemanticPoly::SemanticPoly(const nlohmann::json& poly, int semanticLayerID)
    : PolyFeature(poly), _semanticLayerID(semanticLayerID) {}

SemanticPoly::~SemanticPoly() = default;

void SemanticPoly::calc_footprint_elevation(const SearchTree& searchTree) {
    //-- Calculate elevation of polygon outer boundary
    //-- Point elevation is the average of 5 nearest neighbors from the PC
    for (auto& polypt : _poly.outer_boundary()) {
        Point_3 query(polypt.x() , polypt.y(), 0);
        Neighbor_search search(searchTree, query, 5);
        // TODO: radius search instead of NN?
//        Fuzzy_sphere search_radius(query, 5);
//        std::list<Point_3> result;
//        searchTree.search(std::back_inserter(result), search_radius);

        std::vector<double> poly_height;
        for (Neighbor_search::iterator it = search.begin(); it != search.end(); ++it) {
            poly_height.push_back(it->first.z());
        }
//        for (auto& pt : result) {
//            poly_height.push_back(pt.z());
//        }
        _base_heights.emplace_back(geomtools::avg(poly_height));
    }

    //-- In case of inner rings, set inner points as average of outer points, as the last element in _baseHeights
    if (_poly.has_holes()) {
        _base_heights.emplace_back(geomtools::avg(_base_heights));
    }
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