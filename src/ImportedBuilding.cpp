#include "ImportedBuilding.h"

#include "geomutils.h"

ImportedBuilding::ImportedBuilding(nlohmann::json buildingJson, const std::vector<Point_3>& importedBuildingPts, const int internalID)
        : Building(internalID), _buildingJson(std::move(buildingJson)), _dPts(importedBuildingPts),
          _avgFootprintHeight(-9999), _footprintIdxList(), _parentBuildingID(),
          _appendToBuilding(false), _lodIdx(-1) {

    //-- Get parent building ID
    _parentBuildingID = _buildingJson["parents"].front();

    //-- Define LoD
    int idx = 0;
    for (auto& lodGeom : _buildingJson["geometry"]) {
//                    if (building["lod"] == "2.2") {
        if (lodGeom["lod"] == 2.2) { //todo make it an option and bump to CityJSON 1.1
            _lodIdx = idx;
            break; // todo make it search for highest LoD as option
        }
        ++idx;
    }
    if (_lodIdx == -1) throw std::runtime_error(std::string("Didn't find LoD in question for buildingPart ID:" + std::to_string(internalID)));

    nlohmann::json& geometry = _buildingJson["geometry"][_lodIdx];

    //-- Get the footprint polygon
    //- Find GroundSurf semantic index
    nlohmann::json& semanticSurfaces = geometry["semantics"]["surfaces"];
    int groundSemanticIdx = -9999;
    for (int i = 0; i < semanticSurfaces.size(); ++i) {
        if (semanticSurfaces[i]["type"] == "GroundSurface") {
            groundSemanticIdx = i;
            break;
        }
    }
    if (groundSemanticIdx == -9999) throw std::runtime_error("Cannot find 'GroundSurface' in imported buildings CityJSON file!");

    //- Find boundary ID of the footprint
    nlohmann::json& semanticValues = geometry["semantics"]["values"].front();
    for (int i = 0; i < semanticValues.size(); ++i) {
        if (semanticValues[i] == groundSemanticIdx) {
            _footprintIdxList.push_back(i);
        }
    }

    //- Handle building part in case it is not a ground part
    if (_footprintIdxList.empty()) {
        _appendToBuilding = true;
        return;
    }

    //- Construct footprint polygon from ground surface
    std::vector<double> footprintElevations;
    for (int & footprintIdx : _footprintIdxList) {
        nlohmann::json coordBnd = geometry["boundaries"].front()[footprintIdx];
        for (auto& ring: coordBnd) {
            Polygon_2 tempPoly;
            //todo need to map those points to _dPts so I can edit them later
            for (auto& ptIdx: ring) {
                tempPoly.push_back(Point_2(_dPts[ptIdx].x(), _dPts[ptIdx].y()));
                footprintElevations.push_back(_dPts[ptIdx].z());
            }
            if (!tempPoly.is_simple()) {
                config::log << "Failed to import building: " << this->get_parent_building_id()
                            << " Reason: " << "Footprint polygon is not valid." << std::endl;
                this->deactivate();
                return;
            }
            CGAL::internal::pop_back_if_equal_to_front(tempPoly);
            if (_poly._rings.empty()) {
                if (tempPoly.is_clockwise_oriented()) tempPoly.reverse_orientation();
            } else {
                if (tempPoly.is_counterclockwise_oriented()) tempPoly.reverse_orientation();
            }
            _poly._rings.push_back(tempPoly);
        }
    }

    /*
    //todo need to rewrite this in case of triangulation! Or not if I can avoid handling triangulated footprints
//    CGAL::Polygon_2<EPECK> inftyPoly;
//    inftyPoly.push_back(ePoint_2(0,0));
//    inftyPoly.push_back(ePoint_2(g_largnum,0));
//    inftyPoly.push_back(ePoint_2(g_largnum,g_largnum));
//    inftyPoly.push_back(ePoint_2(0,g_largnum));

    CGAL::Polygon_set_2<EPECK> polySet;
//    polySet.insert(inftyPoly);
    std::vector<double> footprintElevations;
//    int polyNo = 0;
    for (int i = 0; i < _footprintIdxList.size(); ++i) {
        auto& footprintIdx = _footprintIdxList[i];
        //-- Construct footprint polygon from ground surface
        nlohmann::json coordBnd = _buildingJson["boundaries"].front()[footprintIdx].front();
        CGAL::Polygon_2<EPECK> facePoly;
        for (auto& ptIdx: coordBnd) {
            facePoly.push_back(ePoint_2(_dPts[ptIdx].x(), _dPts[ptIdx].y()));
            footprintElevations.push_back(_dPts[ptIdx].z());
        }
        CGAL::internal::pop_back_if_equal_to_front(facePoly);
        if (facePoly.is_clockwise_oriented()) facePoly.reverse_orientation();

//        std::cout << facePoly << std::endl;
//        std::cout << "NOW HANDLING POLY NO: " << polyNo++ << std::endl;
        polySet.join(facePoly);
    }
//    polySet.remove_redundant_edges();
//    std::cout << "NUMBER OF POLYGONS IN THIS COMPOUND:" << polySet.number_of_polygons_with_holes() << std::endl;
    std::list<CGAL::Polygon_with_holes_2<EPECK>> res;
    polySet.polygons_with_holes(std::back_inserter (res));

    Converter<EPECK, EPICK> to_inexact;
    Polygon_2 transferKernelPoly;
    for (auto& outerPt : res.front().outer_boundary()) {
        Point_2 polyPt = to_inexact(outerPt);
        transferKernelPoly.push_back(polyPt);
    }
    _poly._rings.push_back(transferKernelPoly);
     */

    _avgFootprintHeight = geomutils::avg(footprintElevations);
}

ImportedBuilding::~ImportedBuilding() = default;

void ImportedBuilding::reconstruct() {
    typedef EPICK::FT                 FT;
    typedef std::array<FT, 3>         Custom_point;
    typedef std::vector<std::size_t>  CGAL_Polygon;
    nlohmann::json& geometry = _buildingJson["geometry"][_lodIdx];

    std::vector<std::array<FT, 3>> points;
    std::vector<CGAL_Polygon> polygons;
    int surfIdx = -1;
    for (auto& faces : geometry["boundaries"].front()) {
        //-- Remove bottom surface
        ++surfIdx;
        if (std::find(_footprintIdxList.begin(), _footprintIdxList.end(), surfIdx) != _footprintIdxList.end())
            continue;

        for (auto& faceLst : faces) {
            CGAL_Polygon p;
            for (auto& face : faceLst) {
                points.push_back(CGAL::make_array<FT>(_dPts[face].x(), _dPts[face].y(), _dPts[face].z()));
                p.push_back(points.size() - 1);
            }
            polygons.push_back(p);
        }
    }
    //-- Do mumbo-jumbo to make everything work
    PMP::repair_polygon_soup(points, polygons, CGAL::parameters::geom_traits(geomutils::Array_traits()));
    PMP::orient_polygon_soup(points, polygons);
    PMP::polygon_soup_to_polygon_mesh(points, polygons, _mesh);
    PMP::triangulate_faces(_mesh);

        /*
        //-- Add other surfaces
        std::vector<Mesh::Vertex_index> faceVertices;
        int facid = -1; //todo temp
        for (auto& faceLst : faces) {
            if (++facid > 0) std::cout << "YO THERE's A FACE NO: " << facid << std::endl;
//            if (!(facid > 0)) continue;
            for (auto& face: faceLst) {
                faceVertices.emplace_back(_mesh.add_vertex(_dPts[face]));
            }
            bool isReconstruct = _mesh.add_face(faceVertices);
            //todo temp
            if (!isReconstruct) {
                std::cout << "I have a failed surface!!!" << std::endl;
                CGAL::Polygon_with_holes_2<EPICK> tempPoly;
                for (auto& face : faceLst) {
                    tempPoly.outer_boundary().push_back(Point_2(_dPts[face].x(), _dPts[face].y()));
                }
                std::cout << "IS THAT FAILED SURFACE VALID POLY?? : " << tempPoly.outer_boundary().is_simple() << std::endl;
//                CGAL::draw(tempPoly);
            }
            PMP::duplicate_non_manifold_vertices(_mesh);
        }
    }
    PMP::stitch_borders(_mesh);
    PMP::triangulate_faces(_mesh);
    */
}

void ImportedBuilding::append_nonground_part(const std::shared_ptr<ImportedBuilding>& other) {
    this->_buildingJson["geometry"][this->_lodIdx]["boundaries"].front()
          .push_back(other->get_building_json()["geometry"][other->get_lod_idx()]["boundaries"].front());
}

const nlohmann::json& ImportedBuilding::get_building_json() const {
    return _buildingJson;
}

const std::string& ImportedBuilding::get_parent_building_id() const {
    return _parentBuildingID;
}

const int ImportedBuilding::get_lod_idx() const {
    return _lodIdx;
}

const bool ImportedBuilding::is_appending() const {
    return _appendToBuilding;
}