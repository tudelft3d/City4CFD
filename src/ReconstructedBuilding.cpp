#include "ReconstructedBuilding.h"

#include "geomutils.h"
#include "LoD12.h"

ReconstructedBuilding::ReconstructedBuilding()
        : Building(), _searchTree(nullptr) {}

ReconstructedBuilding::ReconstructedBuilding(const int internalID)
        : Building(internalID), _searchTree(nullptr) {}

ReconstructedBuilding::ReconstructedBuilding(const nlohmann::json& poly)
        : Building(poly), _searchTree(nullptr) {}

ReconstructedBuilding::ReconstructedBuilding(const nlohmann::json& poly, const int internalID)
        : Building(poly, internalID), _searchTree(nullptr) {}

ReconstructedBuilding::~ReconstructedBuilding() = default;

void ReconstructedBuilding::set_search_tree(const std::shared_ptr<SearchTree>& searchTree) {
    _searchTree = searchTree;
}

void ReconstructedBuilding::reconstruct() {
    //-- Take tree subset bounded by the polygon
    std::vector<Point_3> subsetPts;
    Point_3 bbox1(_poly.bbox().xmin(), _poly.bbox().ymin(), -g_largnum);
    Point_3 bbox2(_poly.bbox().xmax(), _poly.bbox().ymax(), g_largnum);
    Fuzzy_iso_box pts_range(bbox1, bbox2);
    _searchTree->search(std::back_inserter(subsetPts), pts_range);

    //-- Check if subset point lies inside the polygon
    std::vector<double> building_pts;
    for (auto& pt : subsetPts) {
        if (geomutils::point_in_poly(pt, _poly)) {
            building_pts.push_back(pt.z());
        }
    }

    //-- Don't reconstruct if there are no points belonging to the polygon
    if (building_pts.empty()) {
        this->deactivate();
        throw std::domain_error("Found no points belonging to the building");
    }

    //-- LoD12 reconstruction
    LoD12 lod12(_poly, _base_heights, building_pts);
    lod12.lod12reconstruct(_mesh, _height);

    double lowHeight = 2.; // Hardcoded low height here
    if (lod12.get_height() < lowHeight) { // In case of a small height
        this->deactivate();
        throw std::domain_error("Building height lower than minimum prescribed height");
    }
}

void ReconstructedBuilding::get_cityjson_info(nlohmann::json& b) const {
    b["type"] = "Building";
//  b["attributes"];
//    get_cityjson_attributes(b, _attributes);
//    float hbase = z_to_float(this->get_height_base());
//    float h = z_to_float(this->get_height());
//    b["attributes"]["TerrainHeight"] = _baseHeights.back(); // temp - will calculate avg for every footprint
    b["attributes"]["measuredHeight"] = _height - geomutils::avg(_base_heights[0]);
}

void ReconstructedBuilding::get_cityjson_semantics(nlohmann::json& g) const { // Temp for checking CGAL mesh properties
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