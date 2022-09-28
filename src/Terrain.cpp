/*
  City4CFD
 
  Copyright (c) 2021-2022, 3D Geoinformation Research Group, TU Delft  

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

#include "Terrain.h"

#include "geomutils.h"
#include "io.h"
#include "SurfaceLayer.h"

Terrain::Terrain()
        : TopoFeature(0), _cdt(), _surfaceLayersTerrain(),
          _constrainedPolys(), _vertexFaceMap(), _searchTree() {}

Terrain::Terrain(int pid)
        : TopoFeature(pid), _cdt(), _surfaceLayersTerrain(),
          _constrainedPolys(), _vertexFaceMap(), _searchTree() {}

Terrain::~Terrain() = default;

void Terrain::set_cdt(const Point_set_3& pointCloud) {
    Converter<EPICK, EPECK> to_exact;

    std::cout << "\n    Preparing triangulation" << std::endl;
    int count = 0;
    std::vector<ePoint_3> pts;
    for (auto& pt : pointCloud.points()) {
        pts.push_back(to_exact(pt));

        if ((count % 5000) == 0) IO::print_progress_bar(100 * count / pointCloud.size());
        ++count;
    }
    IO::print_progress_bar(100); std::clog << std::endl;
    std::cout << "    Triangulating..." << std::flush;
    _cdt.insert(pts.begin(), pts.end());
    std::cout << "\r    Triangulating...done" << std::endl;

    /*
    //-- Smoothing
    if (Config::get().smoothTerrain) {
        std::cout << "\n    Smoothing" << std::endl;
        geomutils::smooth_dt<CDT, EPECK>(pointCloud, _cdt);
    }
   */
}

void Terrain::prep_constraints(const PolyFeatures& features, Point_set_3& pointCloud) {
    std::cout << "    Lifting polygon edges to terrain height" << std::endl;
    int countFeatures = 0;
    auto is_building_pt = pointCloud.property_map<bool>("is_building_point").first;
    for (auto& f : features) {
        if (!f->is_active()) continue;
        bool is_building = false;
        if (f->get_class() == BUILDING) is_building = true;
        int polyCount = 0;
        for (auto& ring : f->get_poly().rings()) {
            auto& heights = f->get_base_heights();
            //-- Add ring points
            int i = 0;
            Polygon_3 pts;
            for (auto& polyVertex : ring) {
                pts.push_back(ePoint_3(polyVertex.x(), polyVertex.y(), heights[polyCount][i]));
                auto it = pointCloud.insert(Point_3(polyVertex.x(), polyVertex.y(), heights[polyCount][i++]));
                if (is_building) is_building_pt[*it] = true;
            }
            _constrainedPolys.push_back(pts);
            ++polyCount;
        }
        ++countFeatures;
    }
    std::clog << "\n    Num of polygons to constrain: " << countFeatures << std::endl;
}

void Terrain::constrain_features() {
    std::cout << "\n    Imprinting polygons" << std::endl;

    int count = 0;
    for (auto& ring : _constrainedPolys) {
        //-- Set added points as constraints
        _cdt.insert_constraint(ring.begin(), ring.end(), true);

        if ((count % 100) == 0) IO::print_progress_bar(100 * count / _constrainedPolys.size());
        ++count;
    }
    IO::print_progress_bar(100); std::clog << std::endl;
}

void Terrain::create_mesh(const PolyFeatures& features) {
    _mesh.clear();
    //-- Mark surface layer
    geomutils::mark_domains(_cdt, features);

    //-- Create the mesh for the terrain
    geomutils::cdt_to_mesh(_cdt, _mesh);

    // -- Surface layer meshes are stored here
    for (int i : Config::get().surfaceLayerIDs) {
        auto layer = std::make_shared<SurfaceLayer>(i);
        geomutils::cdt_to_mesh(_cdt, layer->get_mesh(), i); // Create mesh for surface layers
        _surfaceLayersTerrain.push_back(layer);
    }
}

void Terrain::prepare_subset() {
    typedef std::vector<std::size_t>  CGAL_Polygon;

    //-- Make terrain mesh without surface layers marked
    geomutils::cdt_to_mesh(_cdt, _mesh);
    if (_mesh.is_empty()) throw std::runtime_error("Cannot create vertex-face map of empty mesh!");

    //-- Construct vertex-face map used for search
    _vertexFaceMap.clear();
    for (auto& face : _mesh.faces()) {
        for (auto vertex : CGAL::vertices_around_face(_mesh.halfedge(face), _mesh)) {
            auto pt = _mesh.point(vertex);
            auto it = _vertexFaceMap.find(pt);
            if (it == _vertexFaceMap.end()) {
                _vertexFaceMap[pt].push_back(face);
            } else {
                it->second.push_back(face);
            }
        }
    }
    //-- Construct search tree
    _searchTree.clear();
    for (auto& pt : _mesh.points()) {
        _searchTree.insert(pt);
    }
}

Mesh Terrain::mesh_subset(const Polygon_with_holes_2& poly) const {
    //-- Get mesh vertices that are bounded by the polygon
    std::vector<Point_3> subsetPts;
    double expandSearch = 1;
    while (subsetPts.size() <= poly.outer_boundary().size()) {
        subsetPts.clear();
        Point_3 bbox1(poly.bbox().xmin() - expandSearch, poly.bbox().ymin() - expandSearch, -global::largnum);
        Point_3 bbox2(poly.bbox().xmax() + expandSearch, poly.bbox().ymax() + expandSearch, global::largnum);
        Fuzzy_iso_box pts_range(bbox1, bbox2);
        _searchTree.search(std::back_inserter(subsetPts), pts_range);
        expandSearch *= 2;
    }

    //-- Generate mesh subset from vertices
    Mesh subsetMesh;
    std::unordered_map<Point_3, vertex_descriptor> dPts;
    std::vector<face_descriptor> faceLst;
    for (auto& pt : subsetPts) {
//        auto it = _vertexFaceMap.find(IO::gen_key_bucket(pt));
        auto it = _vertexFaceMap.find(pt);
        assert(it != _vertexFaceMap.end());
        for (auto& face : it->second) {
            //- Make sure the same face isn't added multiple times
            if (std::find(faceLst.begin(), faceLst.end(), face) == faceLst.end()) {
                faceLst.push_back(face);
            } else {
                continue;
            }
            //- Take subset
            std::vector<vertex_descriptor> facePts; facePts.reserve(3);
            for (auto vertex : CGAL::vertices_around_face(_mesh.halfedge(face), _mesh)) {
                auto meshPt = _mesh.point(vertex);
                auto vertexIt = dPts.find(meshPt);
                if (vertexIt == dPts.end()) {
                    auto meshVertId = subsetMesh.add_vertex(meshPt);
                    dPts[meshPt] = meshVertId;
                    facePts.push_back(meshVertId);
                } else {
                    facePts.push_back(vertexIt->second);
                }
            }
            subsetMesh.add_face(facePts[0], facePts[1], facePts[2]);
        }
    }
    return subsetMesh;
}

void Terrain::clear_subset() {
    _vertexFaceMap.clear();
    _searchTree.clear();
}

void Terrain::get_cityjson_info(nlohmann::json& b) const {
    b["type"] = "TINRelief";
//    b["attributes"]; // commented out until I have attributes to add
}

std::string Terrain::get_cityjson_primitive() const {
    return "CompositeSurface";
}

CDT& Terrain::get_cdt() {
    return _cdt;
}

const CDT& Terrain::get_cdt() const {
    return _cdt;
}

const vertex_face_map& Terrain::get_vertex_face_map() const {
    return _vertexFaceMap;
}

const SearchTree& Terrain::get_mesh_search_tree() const {
    return _searchTree;
}

TopoClass Terrain::get_class() const {
    return TERRAIN;
}

std::string Terrain::get_class_name() const {
    return "Terrain";
}

const SurfaceLayers& Terrain::get_surface_layers() const {
    return _surfaceLayersTerrain;
}