/*
  City4CFD
 
  Copyright (c) 2021-2023, 3D Geoinformation Research Group, TU Delft

  This file is part of City4CFD.

  City4CFD is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  City4CFD is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with City4CFD.  If not, see <http://www.gnu.org/licenses/>.

  For any information or further details about the use of City4CFD, contact
  Ivan Pađen
  <i.paden@tudelft.nl>
  3D Geoinformation Research Group
  Delft University of Technology
*/

#include "ImportedBuilding.h"

#include <utility>

#include "geomutils.h"
#include "io.h"

#include <CGAL/alpha_wrap_3.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

int ImportedBuilding::noBottom = 0;

ImportedBuilding::ImportedBuilding(std::unique_ptr<nlohmann::json>& buildingJson, PointSet3Ptr& importedBuildingPts, const int internalID)
        : Building(internalID), _buildingJson(std::move(buildingJson)),
          _footprintIdxList(), _parentBuildingID(), _ptMap(),
          _appendToBuilding(false), _lodIdx(-1), _footprintPtsIdxList(), _trueHeight(Config::get().importTrueHeight) {

    _f_imported = true; // the flag is here to avoid shorten polygons later. todo to fix
    //-- Get parent building ID
    _parentBuildingID = (*_buildingJson)["parents"].front();

    //-- Define LoD
    std::map<std::string, int> lodGeomLst;
    int idx = 0;
    for (auto& lodGeom : (*_buildingJson)["geometry"]) {
        lodGeomLst[lodGeom["lod"]] = idx++;
    }
    auto it = lodGeomLst.find(Config::get().importLoD);
    if (it == lodGeomLst.end()) {
        --it;
    }
    _lodIdx = it->second;

    nlohmann::json& geometry = (*_buildingJson)["geometry"][_lodIdx];

    //-- Collect points belonging to a building and create a point map
    double highestPt = -global::largnum;
    for (auto& faces : geometry["boundaries"].front()) {
        for (auto& faceLst : faces) {
            for (const int& facePt : faceLst) {
                _ptMap[facePt]=(importedBuildingPts->point(facePt));
                //-- Store max height to calculate building height later
                if (importedBuildingPts->point(facePt).z() > highestPt)
                    highestPt = importedBuildingPts->point(facePt).z();
            }
        }
    }
    // sort out the highest point depending on whether buildings with height or elevation were imported
    if (Config::get().importTrueHeight) {
        _elevation = highestPt;
    } else {
        _height = highestPt;
    }

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
    // if CityJSON is turned to mesh here, I can use the footprint extraction algorithm from CAD as a
    // workaround
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

    std::unordered_map<std::string, int> pointConnectivity;
    CGAL::Polygon_set_2<EPECK> polySet;
    std::vector<double> footprintElevations;
    //    int polyNo = 0;
    for (int& footprintIdx : _footprintIdxList) {
        //-- Construct footprint polygon from ground surface
        nlohmann::json coordBnd = geometry["boundaries"].front()[footprintIdx].front();
        CGAL::Polygon_2<EPECK> facePoly;
        for (const int& ptIdx: coordBnd) {
            facePoly.push_back(ePoint_2(_ptMap.at(ptIdx).x(), _ptMap.at(ptIdx).y()));
            footprintElevations.push_back(_ptMap.at(ptIdx).z());
            pointConnectivity[IO::gen_key_bucket(Point_2(_ptMap.at(ptIdx).x(), _ptMap.at(ptIdx).y()))] = ptIdx;
        }
        if (!facePoly.is_simple()) {
            Config::write_to_log("Failed to import building: " + this->get_parent_building_id()
                                       + " Reason: Footprint polygon is not simple.");
            this->deactivate();
            return;
        }
        geomutils::pop_back_if_equal_to_front(facePoly);
        if (facePoly.is_clockwise_oriented()) facePoly.reverse_orientation();

        polySet.join(facePoly);
    }
    //-- Polyset to polygon data structure
    this->polyset_to_polygon(polySet);

    //-- Connect footprint polygons to JSON object
    this->set_footprint_mesh_connectivity(pointConnectivity);
}

ImportedBuilding::ImportedBuilding(Mesh& mesh, const int internalID)
    : Building(internalID), _buildingJson(std::make_unique<nlohmann::json>()),
      _footprintIdxList(), _parentBuildingID(), _appendToBuilding(false), _ptMap(),
      _lodIdx(-1), _footprintPtsIdxList(), _trueHeight(Config::get().importTrueHeight) {

    _f_imported = true; // the flag is here to avoid shorten polygons later. todo to fix
    //-- Get the polygon from the building bottom
    // group faces pointing down
    CGAL::Vector_3<EPICK> downVector(0, 0, -1);
    std::vector<int> visited(mesh.faces().size(), false);

    std::vector<std::vector<Mesh::face_index>> downwardFaceGroups;
    auto vnormals = mesh.add_property_map<vertex_descriptor, Vector_3>("v:normals", CGAL::NULL_VECTOR).first;
    auto fnormals = mesh.add_property_map<face_descriptor, Vector_3>("f:normals", CGAL::NULL_VECTOR).first;
    PMP::compute_normals(mesh, vnormals, fnormals);
    const double angle = 40;
    //-- Group adjacent downward-facing faces
    for (int i = 0; i < visited.size(); ++i) {
        if (visited[i]) continue;
        visited[i] = true;
        std::vector<Mesh::face_index> downwardFaceGroup;
        // Check if a face is pointing downwards
        auto faceIt = mesh.faces_begin();
        std::advance(faceIt, i);
        if (CGAL::approximate_angle(fnormals[*faceIt], downVector) < angle) {
            std::list<Mesh::face_index> checkList; checkList.push_front(*faceIt);
            // DFS
            while (!checkList.empty()) {
                auto currEdge = mesh.halfedge(checkList.front());
                downwardFaceGroup.push_back(mesh.face(currEdge));
                checkList.pop_front();
                // Loop over adjacent faces
                for (auto adjHalfedge : mesh.halfedges_around_face(currEdge)) {
                    auto nbrFace = mesh.face(mesh.opposite(adjHalfedge));
                    if (nbrFace.idx() > visited.size()) continue;
                    if (visited[nbrFace.idx()])
                        continue;
                    else
                        visited[nbrFace.idx()] = true;
                    if (CGAL::approximate_angle(fnormals[nbrFace], downVector) < angle) {
                        checkList.push_back(nbrFace);
                    }
                }
            }
        }
        if (!downwardFaceGroup.empty()) downwardFaceGroups.push_back(downwardFaceGroup);
    }
    if (downwardFaceGroups.empty()) {
//        std::cout << "I have a building without found bottom surface!" << std::endl;
        ++noBottom;
        this->deactivate();
        return;
    }
    //-- Get the lowest group
    //todo what if there are more?
    int lowestGroup = 0;
    double minHeight = mesh.vertex(mesh.edge(mesh.halfedge(downwardFaceGroups.front().front())), 0);
    for (int i = 0; i < downwardFaceGroups.size(); ++i) {
        double groupPt = mesh.point(mesh.vertex(mesh.edge(mesh.halfedge(downwardFaceGroups[i].front())), 0)).z();
        if (groupPt < minHeight) {
            lowestGroup = i;
            minHeight = groupPt;
        }
    }
    //-- Get all pts for reconstruction from JSON
    for (auto vert : mesh.vertices()) {
        _ptMap[vert.idx()] = mesh.point(vert);
    }

    //-- Stitch faces of the lowest group into a polygon
    std::unordered_map<std::string, int> pointConectivity;
    CGAL::Polygon_set_2<EPECK> polySet;
    std::vector<double> footprintElevations;
    for (auto& face : downwardFaceGroups[lowestGroup]) {
        CGAL::Polygon_2<EPECK> facePoly;
        for (auto pt : mesh.vertices_around_face(mesh.halfedge(face))) {
            facePoly.push_back(ePoint_2(_ptMap.at(pt.idx()).x(), _ptMap.at(pt.idx()).y()));
            footprintElevations.push_back(_ptMap.at(pt.idx()).z());
            pointConectivity[IO::gen_key_bucket(Point_2(_ptMap.at(pt.idx()).x(), _ptMap.at(pt.idx()).y()))] = pt.idx();
        }
        if (!facePoly.is_simple()) {
            Config::write_to_log("Failed to import building: " + std::to_string(this->get_internal_id())
                              + " Reason: Footprint polygon is not simple.");
            this->deactivate();
            return;
        }
        geomutils::pop_back_if_equal_to_front(facePoly);
        if (facePoly.is_clockwise_oriented()) facePoly.reverse_orientation();

        polySet.join(facePoly);
    }
    //-- Polyset to polygon data structure
    this->polyset_to_polygon(polySet);

    //-- Connect footprint polygons to JSON object
    this->set_footprint_mesh_connectivity(pointConectivity);

    //-- Shove it back to JSON for now
    _lodIdx = 0;
    nlohmann::json geom;
    for (auto& face : mesh.faces()) {
        nlohmann::json faceJson;
        for (auto vert : mesh.vertices_around_face(mesh.halfedge(face))) {
            faceJson.push_back(vert.idx());
        }
        nlohmann::json faceLstJson; faceLstJson.push_back(faceJson);
        geom.push_back(faceLstJson);
    }
    (*_buildingJson)["geometry"].push_back({});
    (*_buildingJson)["geometry"].front()["boundaries"].push_back(geom);
    mesh.clear();
}

ImportedBuilding::~ImportedBuilding() = default;

/*
 * Calculate building elevation without reconstruction.
 * Defined as the highest point
 */
double ImportedBuilding::get_elevation() {
    if (_elevation < -global::largnum + global::smallnum) {
        if (_height > 0) {
            _elevation = _height - this->ground_elevation();
        } else {
            if (_ptMap.empty()) throw std::runtime_error("Building missing points!");
            // loop over all points and find the highest one
            for (auto& pt: _ptMap) {
                if (pt.second.z() > _elevation) _elevation = pt.second.z();
            }
        }
    }
    return _elevation;
}

void ImportedBuilding::reconstruct() {
    typedef EPICK::FT                 FT;
    typedef std::vector<std::size_t>  CGAL_Polygon;

    nlohmann::json& geometry = (*_buildingJson)["geometry"][_lodIdx];

    _mesh.clear();
    if (_clip_bottom) {
        this->translate_footprint(-5);
    }
    //-- Adjust building height points
    if (!_trueHeight || Config::get().ground_xyz.empty()) {
        Vector_3 movePt(0, 0, this->ground_elevation());
        std::vector<int> checkedPt;
        for (auto& faces : geometry["boundaries"].front()) {
            for (auto& faceLst : faces) {
                for (const int facePt: faceLst) {
                    if (std::find(checkedPt.begin(), checkedPt.end(), facePt) == checkedPt.end()) {
                        _ptMap.at(facePt) += movePt;
                        checkedPt.push_back(facePt);
                    } else continue;
                }
            }
        }
    }
    //-- Adjust footprints to terrain
    for (int i = 0; i < _footprintPtsIdxList.size(); ++i) {
        for (int j = 0; j < _footprintPtsIdxList[i].size(); ++j) {
            _ptMap.at(_footprintPtsIdxList[i][j]) = Point_3(_ptMap.at(_footprintPtsIdxList[i][j]).x(),
                                                            _ptMap.at(_footprintPtsIdxList[i][j]).y(),
                                                            _groundElevations[i][j]);
        }
    }
    //-- Add points to mesh
    std::vector<std::array<FT, 3>> points;
    std::vector<CGAL_Polygon> polygons;
//    int surfIdx = -1;
    for (auto& faces : geometry["boundaries"].front()) {
        /*
        //-- Remove bottom surface
        ++surfIdx;
        if (std::find(_footprintIdxList.begin(), _footprintIdxList.end(), surfIdx) != _footprintIdxList.end())
            continue;
        */

        for (auto& faceLst : faces) {
            CGAL_Polygon p;
            for (const int facePt : faceLst) {
                points.push_back(CGAL::make_array<FT>(_ptMap.at(facePt).x(),
                                                      _ptMap.at(facePt).y(),
                                                      _ptMap.at(facePt).z()));
                p.push_back(points.size() - 1);
                //-- Store max elevation
                if (_ptMap.at(facePt).z() > _elevation) _elevation = _ptMap.at(facePt).z();
            }
            polygons.push_back(p);
        }
    }
    //-- Polygon soup to polygon mesh
    PMP::repair_polygon_soup(points, polygons, CGAL::parameters::geom_traits(geomutils::Array_traits()));
    PMP::orient_polygon_soup(points, polygons);
    PMP::polygon_soup_to_polygon_mesh(points, polygons, _mesh);
    PMP::triangulate_faces(_mesh);

    /*
    Mesh wrap;
    const double relative_alpha = 300.;
    const double relative_offset = 5000.;
    CGAL::Bbox_3 bbox = CGAL::Polygon_mesh_processing::bbox(_mesh);
    const double diag_length = std::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                                         CGAL::square(bbox.ymax() - bbox.ymin()) +
                                         CGAL::square(bbox.zmax() - bbox.zmin()));
    const double alpha = diag_length / relative_alpha;
    const double offset = diag_length / relative_offset;
//    CGAL::alpha_wrap_3(mesh, alpha, offset, wrap);
    CGAL::alpha_wrap_3(_mesh, alpha, offset, wrap);
//    CGAL::alpha_wrap_3(_mesh, 10, 10, wrap);
    _mesh = wrap;
     */

    if (_clip_bottom) {
        this->translate_footprint(5);
    }
        /*
        //-- Add other surfaces
        std::vector<Mesh::Vertex_index> faceVertices;
        int facid = -1; //todo temp
        for (auto& faceLst : faces) {
            if (++facid > 0) std::cout << "YO THERE's A FACE NO: " << facid << std::endl;
//            if (!(facid > 0)) continue;
            for (auto& face: faceLst) {
                faceVertices.emplace_back(_mesh.add_vertex(_ptsPtr[face]));
            }
            bool isReconstruct = _mesh.add_face(faceVertices);
            //todo temp
            if (!isReconstruct) {
                std::cout << "I have a failed surface!!!" << std::endl;
                CGAL::Polygon_with_holes_2<EPICK> tempPoly;
                for (auto& face : faceLst) {
                    tempPoly.outer_boundary().push_back(Point_2(_ptsPtr[face].x(), _ptsPtr[face].y()));
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
    if (Config::get().refineImportedBuildings) this->refine();
}

void ImportedBuilding::reconstruct_flat_terrain() {
    _trueHeight = false;
    _groundElevation = 0;
    this->reconstruct();
}

void ImportedBuilding::append_nonground_part(const std::shared_ptr<ImportedBuilding>& other) {
    // append json
    (*_buildingJson)["geometry"][this->_lodIdx]["boundaries"].front()
          .push_back(other->get_building_json()["geometry"][other->get_lod_idx()]["boundaries"].front());
    // append point map
    _ptMap.insert(other->_ptMap.begin(), other->_ptMap.end());
    // reset elevation and height
    _elevation = -global::largnum;
    _height    = -global::largnum;
}

const nlohmann::json& ImportedBuilding::get_building_json() const {
    return *_buildingJson;
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

void ImportedBuilding::check_simplicity(Polygon_2& ring) {
    if (!ring.is_simple()) {
        Config::write_to_log("Failed to import building: " + this->get_parent_building_id()
                             + " Reason: Footprint polygon is not simple.");
        this->deactivate();
        return;
    }
}

void ImportedBuilding::polyset_to_polygon(const CGAL::Polygon_set_2<CGAL::Epeck>& polySet) {
//    polySet.remove_redundant_edges();
    std::list<CGAL::Polygon_with_holes_2<EPECK>> res;
    polySet.polygons_with_holes(std::back_inserter(res));

    Converter<EPECK, EPICK> to_inexact;
    Polygon_2 transferKernelPoly;
    for (auto& outerPt : res.front().outer_boundary()) {
        Point_2 polyPt = to_inexact(outerPt);
        transferKernelPoly.push_back(polyPt);
    }
    this->check_simplicity(transferKernelPoly);
    _poly._rings.push_back(transferKernelPoly);

    for (auto& hole : res.front().holes()) {
        transferKernelPoly.clear();
        for (auto& pt : hole) {
            Point_2 polyPt = to_inexact(pt);
            transferKernelPoly.push_back(polyPt);
        }
        this->check_simplicity(transferKernelPoly);
        _poly._rings.push_back(transferKernelPoly);
    }
}

void ImportedBuilding::set_footprint_mesh_connectivity(const std::unordered_map<std::string, int>& pointConnectivity) {
    for (auto& ring: _poly._rings) {
        std::vector<int> ringFootprintIdxs;
        for (auto& pt : ring) {
            auto itPt = pointConnectivity.find(IO::gen_key_bucket(pt));
            assert(itPt != pointConnectivity.end());
            ringFootprintIdxs.push_back(itPt->second);
        }
        _footprintPtsIdxList.push_back(ringFootprintIdxs);
    }
}