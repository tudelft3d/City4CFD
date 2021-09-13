#include "TopoFeature.h"

//-- TopoFeature class
TopoFeature::TopoFeature()
        : _mesh(), _id(), _f_active(true) {}

TopoFeature::TopoFeature(const int pid)
    : _mesh(), _id(pid), _f_active(true) {}

Mesh& TopoFeature::get_mesh() {
    return _mesh;
}

const Mesh& TopoFeature::get_mesh() const {
    return _mesh;
}

int TopoFeature::get_id() const {
    return _id;
}

bool TopoFeature::is_active() const {
    return _f_active;
}

void TopoFeature::deactivate() {
    _f_active = false;
}

void TopoFeature::get_obj_pts(std::string &fs,
                              std::string &bs,
                              std::unordered_map<std::string, unsigned long> &dPts) {
    IO::get_obj_pts(_mesh, fs, bs, dPts);
}

void TopoFeature::get_stl_pts(std::string &fs) {
    IO::get_stl_pts(_mesh, fs);
}

//-- PolyFeature class
PolyFeature::PolyFeature(const int pid)
    : TopoFeature(pid), _poly(), _base_heights() {}


PolyFeature::PolyFeature(const json& poly, const int pid)
    : PolyFeature(pid) {
    //-- Store the polygon
    for (auto& polyEdges : poly["geometry"]["coordinates"]) {
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

void PolyFeature::calc_footprint_elevation(const SearchTree& searchTree) {}

void PolyFeature::check_influ_region() {
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