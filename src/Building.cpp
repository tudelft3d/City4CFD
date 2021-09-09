#include "Building.h"

Building::Building(const int pid)
    : PolyFeature(pid), _height(-9999.0) {}

Building::Building(const json& poly, const int pid)
    : PolyFeature(poly, pid), _height(-9999.0) {}

void Building::calc_footprint_elevation(const SearchTree& searchTree) {
    //-- Calculate elevation of polygon outer boundary
    //-- Point elevation is the average of 5 nearest neighbors from the PC
    for (auto& polypt : _poly.outer_boundary()) {
        Point_3 query(polypt.x() , polypt.y(), 0);
        Neighbor_search search(searchTree, query, 5);
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
        _base_heights.emplace_back(avg(poly_height));
    }

    //-- In case of inner rings, set inner points as average of outer points, as the last element in _base_heights
    if (_poly.has_holes()) {
        _base_heights.emplace_back(avg(_base_heights));
    }
}

void Building::threeDfy(const SearchTree& searchTree) {
    //-- Take tree subset bounded by the polygon
    std::vector<Point_3> subsetPts;
    Point_3 bbox1(_poly.bbox().xmin(), _poly.bbox().ymin(), -infty);
    Point_3 bbox2(_poly.bbox().xmax(), _poly.bbox().ymax(), infty);
    Fuzzy_iso_box pts_range(bbox1, bbox2);
    searchTree.search(std::back_inserter(subsetPts), pts_range);

    //-- Check if subset point lies inside the polygon
    std::vector<double> building_pts;
    for (auto& pt : subsetPts) {
        if (check_inside(pt, _poly)) {
            building_pts.push_back(pt.z());
        }
    }

    //-- Don't reconstruct if there are no points belonging to the polygon
    // TODO: exception/warning handling
    if (building_pts.empty()) {
        this->deactivate();
        return;
    }

    //-- LoD12 reconstruction
    LoD12 lod12(_poly, _base_heights, building_pts);
    lod12.lod12reconstruct(_mesh);

    double lowHeight = 4.0; // Hardcoded low height here
    // TODO: exception/warning handling
    if (lod12.get_height() < lowHeight) { // in case of a small height
//        std::cerr << "Building with a low height, building ID: " << this->get_id() << std::endl;
        this->deactivate();
        return;
    }
}

TopoClass Building::get_class() const {
    return BUILDING;
}