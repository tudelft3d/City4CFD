#include "TopoFeature.h"

//-- TopoFeature class
TopoFeature::TopoFeature()
        : _mesh(), _id(), _f_active(true) {}

TopoFeature::TopoFeature(const int pid)
    : _mesh(), _id(std::to_string(pid)), _f_active(true) {}

TopoFeature::~TopoFeature() = default;

Mesh& TopoFeature::get_mesh() {
    return _mesh;
}

const Mesh& TopoFeature::get_mesh() const {
    return _mesh;
}

void TopoFeature::set_id(unsigned long id) {
    _id = std::to_string(id);
}

std::string TopoFeature::get_id() const {
    return _id;
}

bool TopoFeature::is_active() const {
    return _f_active;
}

void TopoFeature::deactivate() {
    _f_active = false;
}

void TopoFeature::get_cityjson_info(nlohmann::json& j) const {
    //TEMP UNTIL ALL FUNCTIONS ARE IMPLEMENTED
}

void TopoFeature::get_cityjson_semantics(nlohmann::json& g) const {
    // TEMP until I figure what to do with this
}


std::string TopoFeature::get_cityjson_primitive() const {
    //TEMP UNTIL ALL FUNCTIONS ARE IMPLEMENTED
    return "Nope";
}

//-- PolyFeature class
PolyFeature::PolyFeature() = default;

PolyFeature::PolyFeature(const nlohmann::json& poly)
    : TopoFeature() {
    //-- Store the polygon
    nlohmann::json polygonStart;
    if (poly["geometry"]["type"] == "Polygon") {
        polygonStart = poly["geometry"]["coordinates"];
    } else if (poly["geometry"]["type"] == "MultiPolygon") {
        polygonStart = poly["geometry"]["coordinates"];
//        if (polygonStart.size() > 1) throw std::runtime_error(poly["geometry"]["type"]);

        //-- GOTTA SEE WHAT TO DO HERE
        polygonStart = polygonStart[0];
    } else {
        throw std::runtime_error(poly["geometry"]["type"]);
    }
//    for (auto& polyEdges : poly["geometry"]["coordinates"]) {
    for (auto& polyEdges : polygonStart) {
        Polygon_2 tempPoly;
        for (auto& coords : polyEdges) {
           tempPoly.push_back(Point_2(coords[0], coords[1]));
        }

        if (_poly.is_unbounded()) {
            _poly = Polygon_with_holes_2(tempPoly);
        } else {
            _poly.add_hole(tempPoly);
        }
    }
}

PolyFeature::~PolyFeature() = default;

void PolyFeature::calc_footprint_elevation(const SearchTree& searchTree) {
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
//            poly_height.push_back(0); // test
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

void PolyFeature::check_feature_scope() {
    // TODO: really gotta rewrite those polygons, cgal implementation is awful
    //-- Temporary - will rewrite polygons
    std::vector<Polygon_2> rings;
    //- Add outer poly and holes into one data structure
    rings.push_back(_poly.outer_boundary());
    for (auto& hole : _poly.holes()) {
        rings.push_back(hole);
    }

    //-- Include all polygons that have at least one vertex in the influence region
    for (auto& poly : rings) {
        for (auto& vert : poly) {
            if (pow(config::pointOfInterest.x() - vert.x(), 2)
              + pow(config::pointOfInterest.y() - vert.y(), 2)
              < pow(config::radiusOfInfluRegion,2)) {
                return;
            }
        }
    }
//    std::cout << "Poly ID " << this->get_id() << " is outside the influ region. Deactivating." << std::endl;
    this->deactivate();
}

const Polygon_with_holes_2& PolyFeature::get_poly() const {
    return _poly;
}

const std::vector<double>& PolyFeature::get_base_heights() const {
    return _base_heights;
}

void PolyFeature::threeDfy(const SearchTree& searchTree) {
    //TEMP
}