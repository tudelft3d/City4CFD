#include "ImportedBuilding.h"

#include "geomutils.h"

ImportedBuilding::ImportedBuilding(nlohmann::json buildingJson, std::vector<Point_3>& importedBuildingPts, const int internalID)
        : Building(internalID), _buildingJson(std::move(buildingJson)), _dPts(importedBuildingPts),
          _avgFootprintHeight(-9999), _footprintIdxList(), _parentBuildingID(),
          _appendToBuilding(false), _lodIdx(-1), _footprintPtsIdxList() {

    _f_imported = true;
    //-- Get parent building ID
    _parentBuildingID = _buildingJson["parents"].front();

    //-- Define LoD
    std::map<std::string, int> lodGeomLst;
    int idx = 0;
    for (auto& lodGeom : _buildingJson["geometry"]) {
        lodGeomLst[lodGeom["lod"]] = idx++;
    }
    auto it = lodGeomLst.find(config::importLoD);
    if (it == lodGeomLst.end()) {
        --it;
    }
    _lodIdx = it->second;

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
    //todo exception handling
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
    for (auto& footprintIdx : _footprintIdxList) {
        nlohmann::json coordBnd = geometry["boundaries"].front()[footprintIdx];
        for (auto& ring: coordBnd) {
            std::vector<int> ringFootprintIdxs;
            Polygon_2 tempPoly;
            for (auto& ptIdx: ring) {
                tempPoly.push_back(Point_2(_dPts[ptIdx].x(), _dPts[ptIdx].y()));
                ringFootprintIdxs.push_back(ptIdx);
                footprintElevations.push_back(_dPts[ptIdx].z());
            }
            if (!tempPoly.is_simple()) {
                config::log << "Failed to import building: " << this->get_parent_building_id()
                            << " Reason: " << "Footprint polygon is not valid." << std::endl;
                this->deactivate();
                return;
            }
//            CGAL::internal::pop_back_if_equal_to_front(tempPoly);
            if (_poly._rings.empty()) {
                if (tempPoly.is_clockwise_oriented()) {
                    tempPoly.reverse_orientation();
                    std::reverse(ringFootprintIdxs.begin() + 1, ringFootprintIdxs.end()); // Following the reverse of CGAL
//                    std::reverse(footprintElevations.begin(), footprintElevations.end());
                }
            } else {
                if (tempPoly.is_counterclockwise_oriented()) {
                    tempPoly.reverse_orientation();
                    std::reverse(ringFootprintIdxs.begin() + 1, ringFootprintIdxs.end());
//                    std::reverse(footprintElevations.begin(), footprintElevations.end());
                }
            }
            _footprintPtsIdxList.push_back(ringFootprintIdxs);
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

    //-- Adjust building height points
    if (!config::importTrueHeight) {
        Vector_3 movePt(0, 0, _avgFootprintHeight);
        std::vector<int> checkedPt;
        for (auto& faces : geometry["boundaries"].front()) {
            for (auto& faceLst : faces) {
                for (auto& facePt: faceLst) {
                    if (std::find(checkedPt.begin(), checkedPt.end(), facePt) == checkedPt.end()) {
                        _dPts[facePt] += movePt;
                        checkedPt.push_back(facePt);
                    } else continue;
                }
            }
        }
    }

    //-- Adjust footprints to terrain
    int count = 0;
    for (int i = 0; i < _footprintPtsIdxList.size(); ++i) {
        for (int j = 0; j < _footprintPtsIdxList[i].size(); ++j) {
            _dPts[_footprintPtsIdxList[i][j]] = Point_3(_dPts[_footprintPtsIdxList[i][j]].x(),
                                                         _dPts[_footprintPtsIdxList[i][j]].y(),
                                                         _base_heights[i][j]);
        }
    }

    //-- Add points to mesh
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
            for (auto& facePt : faceLst) {
                points.push_back(CGAL::make_array<FT>(_dPts[facePt].x(), _dPts[facePt].y(), _dPts[facePt].z()));
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