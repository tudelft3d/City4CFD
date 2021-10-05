#include "Building.h"

Building::Building()
    : PolyFeature(), _height(-infty) {}

Building::Building(const nlohmann::json& poly)
    : PolyFeature(poly, 1), _height(-infty) {}

Building::~Building() = default;

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
        if (geomtools::check_inside(pt, _poly)) {
            building_pts.push_back(pt.z());
        }
    }

    //-- Don't reconstruct if there are no points belonging to the polygon
    // TODO: exception/warning handling - add this information to log
    if (building_pts.empty()) {
        this->deactivate();
        throw std::domain_error("Found no points belonging to the building.");
    }

    //-- LoD12 reconstruction
    LoD12 lod12(_poly, _base_heights, building_pts);
    lod12.lod12reconstruct(_mesh, _height);

    double lowHeight = 3.0; // Hardcoded low height here
    if (lod12.get_height() < lowHeight) { // in case of a small height
        this->deactivate();
        throw std::domain_error("Building height lower than minimum prescribed height, building ID: " + std::string(this->get_id()));
    }
}

void Building::get_cityjson_info(nlohmann::json& b) const {
    b["type"] = "Building";
//  b["attributes"];
//    get_cityjson_attributes(b, _attributes);
//    float hbase = z_to_float(this->get_height_base());
//    float h = z_to_float(this->get_height());
//    b["attributes"]["TerrainHeight"] = _baseHeights.back(); // temp - will calculate avg for every footprint
    b["attributes"]["measuredHeight"] = _height - geomtools::avg(_base_heights[0]);
}

void Building::get_cityjson_semantics(nlohmann::json& g) const { // Temp for checking CGAL mesh properties
    Face_property semantics;
    bool foundProperty;
    boost::tie(semantics, foundProperty) = _mesh.property_map<face_descriptor, std::string>("f:semantics");
    //   auto semantics = _mesh.property_map<face_descriptor, std::string>("f:semantics");
    if (!foundProperty) throw std::runtime_error("Semantic property map not found!");

    std::unordered_map<std::string, int> surfaceId;
    surfaceId["RoofSurface"]   = 0; g["semantics"]["surfaces"][0]["type"] = "RoofSurface";
    surfaceId["GroundSurface"] = 1; g["semantics"]["surfaces"][1]["type"] = "GroundSurface";
    surfaceId["WallSurface"]   = 2; g["semantics"]["surfaces"][2]["type"] = "WallSurface";

    for (auto& faceIdx : _mesh.faces()) {
        auto it = surfaceId.find(semantics[faceIdx]);
        if (it == surfaceId.end()) throw std::runtime_error("Could not find semantic attribute!");

        g["semantics"]["values"][faceIdx.idx()] = it->second;
    }
}

std::string Building::get_cityjson_primitive() const {
    return "MultiSurface";
};

TopoClass Building::get_class() const {
    return BUILDING;
}

std::string Building::get_class_name() const {
    return "Building";
}