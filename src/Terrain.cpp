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

#include "Terrain.h"

#include "geomutils.h"
#include "io.h"
#include "SurfaceLayer.h"

Terrain::Terrain()
        : TopoFeature(0), _cdt(), _surfaceLayersTerrain(),
          _constrainedPolys(), _vertexFaceMap(), _extraConstrainedEdges(),
          _searchTree(Config::get().searchtree_bucket_size) {}

Terrain::Terrain(int pid)
        : TopoFeature(pid), _cdt(), _surfaceLayersTerrain(),
          _constrainedPolys(), _vertexFaceMap(), _extraConstrainedEdges(),
          _searchTree(Config::get().searchtree_bucket_size) {}

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
}

void Terrain::prep_constraints(const PolyFeaturesPtr& features, Point_set_3& pointCloud) {
    std::cout << "    Lifting polygon edges to terrain elevation" << std::endl;
    int countFeatures = 0;
    auto is_building_pt = pointCloud.property_map<bool>("is_building_point").first;
    for (auto& f : features) {
        if (!f->is_active()) continue;
        bool is_building = false;
        if (f->get_class() == BUILDING) is_building = true;
        int polyCount = 0;
        for (auto& ring : f->get_poly().rings()) {
            auto& elevations = f->get_ground_elevations();
            //-- Add ring points
            int i = 0;
            Polygon_3 pts;
            for (auto& polyVertex : ring) {
                pts.push_back(ePoint_3(polyVertex.x(), polyVertex.y(), elevations[polyCount][i]));
                auto it = pointCloud.insert(Point_3(polyVertex.x(), polyVertex.y(), elevations[polyCount][i++]));
                if (is_building) is_building_pt[*it] = true;
            }
            _constrainedPolys.push_back(pts);
            ++polyCount;
        }
        ++countFeatures;
    }
    std::clog << "\n    Number of polygons to constrain: " << countFeatures << std::endl;
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
    // extra edges to constrain when whole polygons couldn't be added
    if (!_extraConstrainedEdges.empty()) std::cout << "\n    Adding extra constrained edges" << std::endl;
    for (auto& extraEdge : _extraConstrainedEdges) {
        _cdt.insert_constraint(extraEdge.source(), extraEdge.target());
        ++count;
    }
}

void Terrain::create_mesh(const PolyFeaturesPtr& features) {
    _mesh.clear();
    //-- Mark surface layers
    this->tag_layers(features);

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
    _searchTree.insert(_mesh.points().begin(), _mesh.points().end());
}

Mesh Terrain::mesh_subset(const Polygon_with_holes_2& poly) const {
    //-- Get mesh vertices that are bounded by the polygon
    std::vector<Point_3> subsetPts;
    double expandSearch = 1;
    while (subsetPts.size() <= poly.outer_boundary().size()) {
        subsetPts.clear();
        Point_2 bbox1(poly.bbox().xmin() - expandSearch, poly.bbox().ymin() - expandSearch);
        Point_2 bbox2(poly.bbox().xmax() + expandSearch, poly.bbox().ymax() + expandSearch);
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

/*
 * Domain marker enriched with surface layer tagging
 */
void Terrain::tag_layers(const Face_handle& start,
                         int index,
                         std::list<CDT::Edge>& border,
                         const PolyFeaturesPtr& features)
{
    if (start->info().nesting_level != -1) {
        return;
    }

    //-- Check which polygon contains the constrained (i.e. non-terrain) point
    Point_3 chkPoint;
    Converter<EPECK, EPICK> to_inexact;
    if (!features.empty()) {
        chkPoint = CGAL::centroid(to_inexact(start->vertex(0)->point()),
                                  to_inexact(start->vertex(1)->point()),
                                  to_inexact(start->vertex(2)->point()));
    }
    int surfaceLayer = -1; //-- Default value is unmarked triangle, i.e. general terrain
    if (index != 0) {
        for (auto& feature : features) {
            if (!feature->is_active()) continue;
            //- Polygons are already ordered according to importance - find first polygon
            if (geomutils::point_in_poly(chkPoint, feature->get_poly())) {
                if (feature->get_class() == BUILDING) {
//                    surfaceLayer = 9999; //- Remove building footprints from terrain
                    surfaceLayer = -1; //- Leave building footprints as part of terrain
                    break;
                } else {
                    surfaceLayer = feature->get_output_layer_id();
                    break;
                }
            }
        }
    }
    std::list<Face_handle> queue;
    queue.push_back(start);
    while (! queue.empty()) {
        Face_handle fh = queue.front();
        queue.pop_front();
        if (fh->info().nesting_level == -1) {
            fh->info().nesting_level = index;
            if (surfaceLayer != -1) {
                fh->info().surfaceLayer = surfaceLayer;
//                check_layer(fh, surfaceLayer);
            }
            for (int i = 0; i < 3; i++) {
                CDT::Edge e(fh,i);
                Face_handle n = fh->neighbor(i);
                if (n->info().nesting_level == -1) {
                    if (_cdt.is_constrained(e)) {
                    #pragma omp critical
                        border.push_back(e);
                    } else queue.push_back(n);
                }
            }
        }
    }
}

void Terrain::tag_layers(const PolyFeaturesPtr& features) {
    for (CDT::Face_handle f : _cdt.all_face_handles()) {
        f->info().nesting_level = -1;
    }
    std::list<CDT::Edge> border;
    tag_layers(_cdt.infinite_face(), 0, border, features);
    #pragma omp parallel
    while (!border.empty()) {
        bool sstop = false;
        CDT::Edge e;
        #pragma omp critical
        {
            if (!border.empty()) {
                e = border.front();
                border.pop_front();
            } else {
                sstop = true;
            }
        }
        if (sstop) continue;

        Face_handle n = e.first->neighbor(e.second);
        if (n->info().nesting_level == -1) {
            tag_layers(n, e.first->info().nesting_level + 1, border, features);
        }
    }
    for (CDT::Face_handle f : _cdt.all_face_handles()) {
        if (f->info().surfaceLayer == -2) {
            f->info().surfaceLayer = -9999;
        }
    }
}

/*
void Terrain::check_layer(const Face_handle& fh, int surfaceLayer) {
    if (fh->info().surfaceLayer == 9999) return;
    auto it = Config::get().flattenSurfaces.find(surfaceLayer);
    if (it != Config::get().flattenSurfaces.end()) {
        Converter<EPECK, EPICK> to_inexact;
        Vector_3 vertical(0, 0, 1);
        Vector_3 norm = CGAL::normal(to_inexact(fh->vertex(0)->point()),
                                     to_inexact(fh->vertex(1)->point()),
                                     to_inexact(fh->vertex(2)->point()));

        if (CGAL::approximate_angle(norm, vertical) < 60.0) {
            fh->info().surfaceLayer = surfaceLayer;
        } else {
            fh->info().surfaceLayer = -2;
        }
    } else {
        fh->info().surfaceLayer = surfaceLayer;
    }
}
*/

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

std::vector<Polygon_3>& Terrain::get_constrained_polys() {
    return _constrainedPolys;
}

const vertex_face_map& Terrain::get_vertex_face_map() const {
    return _vertexFaceMap;
}

const SearchTree& Terrain::get_mesh_search_tree() const {
    return _searchTree;
}

std::vector<Polygon_3::Segment_2>& Terrain::get_extra_constrained_edges() {
    return _extraConstrainedEdges;
}

TopoClass Terrain::get_class() const {
    return TERRAIN;
}

std::string Terrain::get_class_name() const {
    return "Terrain";
}

const SurfaceLayersPtr& Terrain::get_surface_layers() const {
    return _surfaceLayersTerrain;
}