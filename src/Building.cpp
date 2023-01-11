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

#include "Building.h"

#include "geomutils.h"
#include "LoD12.h"
#include "Terrain.h"

#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/clip.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Mesh_3/dihedral_angle_3.h>
#include <CGAL/alpha_wrap_3.h>

Building::Building()
        : PolyFeature(1), _elevation(-global::largnum), _height(-global::largnum),
          _ptsPtr(std::make_shared<Point_set_3>()) {}

Building::Building(const int internalID)
        : PolyFeature(1, internalID), _elevation(-global::largnum), _height(-global::largnum),
          _ptsPtr(std::make_shared<Point_set_3>()) {}

Building::Building(const nlohmann::json& poly)
        : PolyFeature(poly, true, 1), _elevation(-global::largnum), _height(-global::largnum),
          _ptsPtr(std::make_shared<Point_set_3>()) {}
        // 'true' here to check for polygon simplicity

Building::Building(const nlohmann::json& poly, const int internalID)
        : PolyFeature(poly, true, 1, internalID), _elevation(-global::largnum), _height(-global::largnum),
          _ptsPtr(std::make_shared<Point_set_3>()) {}
        // 'true' here to check for polygon simplicity

Building::~Building() = default;

void Building::insert_point(const Point_3& pt) {
    _ptsPtr->insert(pt);
}

double Building::get_height() {
    if (_height < -global::largnum + global::smallnum) {
        _height = this->get_elevation() - this->ground_elevation();
    }
    return _height;
}

void Building::alpha_wrap(const BuildingsPtr& buildings, Mesh& newMesh) {
    typedef EPICK::FT                 FT;
    typedef std::vector<std::size_t>  CGAL_Polygon;

    //-- Make a single mesh out of all individual buildings
    std::vector<std::array<FT, 3>> points;
    std::vector<CGAL_Polygon> polygons;
    for (auto& b : buildings) {
        auto& mesh = b->get_mesh();
        for (auto& face: mesh.faces()) {
            CGAL_Polygon p;
            auto vertices = mesh.vertices_around_face(mesh.halfedge(face));
            for (auto vertex = vertices.begin(); vertex != vertices.end(); ++vertex) {
                points.push_back(CGAL::make_array<FT>(mesh.point(*vertex).x(),
                                                      mesh.point(*vertex).y(),
                                                      mesh.point(*vertex).z()));
                p.push_back(points.size() - 1);
            }
            polygons.push_back(p);
        }
    }
    PMP::repair_polygon_soup(points, polygons, CGAL::parameters::geom_traits(geomutils::Array_traits()));
    PMP::orient_polygon_soup(points, polygons);
    PMP::polygon_soup_to_polygon_mesh(points, polygons, newMesh);
    PMP::triangulate_faces(newMesh);

    /*
    typedef Mesh::Halfedge_index           halfedge_descriptor;
    typedef Mesh::Edge_index               edge_descriptor;
    //-- Set the property map for constrained edges
    Mesh::Property_map<edge_descriptor,bool> is_constrained =
            newMesh.add_property_map<edge_descriptor,bool>("e:is_constrained",false).first;

    //-- Detect sharp features
    for (auto& e : edges(newMesh)) {
        halfedge_descriptor hd = halfedge(e,newMesh);
        if (!is_border(e,newMesh)) {
            double angle = CGAL::Mesh_3::dihedral_angle(newMesh.point(source(hd,newMesh)),
                                                        newMesh.point(target(hd,newMesh)),
                                                        newMesh.point(target(next(hd,newMesh),newMesh)),
                                                        newMesh.point(target(next(opposite(hd,newMesh),newMesh),newMesh)));
            if (CGAL::abs(angle)<179.5)
                is_constrained[e]=true;
        }
    }

    Mesh wrap;
    CGAL::alpha_wrap_3(newMesh, 0.1, 0.01, wrap,
                             CGAL::parameters::edge_is_constrained_map(is_constrained));
    */

    //-- Perform CGAL's alpha wrapping
//    const double relative_alpha = 900.;
//    const double relative_offset = 15000.;
//    CGAL::Bbox_3 bbox = CGAL::Polygon_mesh_processing::bbox(newMesh);
//    const double diag_length = std::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
//                                         CGAL::square(bbox.ymax() - bbox.ymin()) +
//                                         CGAL::square(bbox.zmax() - bbox.zmin()));
//    const double alpha = diag_length / relative_alpha;
//    const double offset = diag_length / relative_offset;
//    CGAL::alpha_wrap_3(newMesh, alpha, offset, wrap);

//    CGAL::alpha_wrap_3(newMesh, 2, 0.01, newMesh);   // 'coarse'
    CGAL::alpha_wrap_3(newMesh, 1.5, 0.03, newMesh); // 'medium'
//    CGAL::alpha_wrap_3(newMesh, 0.7, 0.03, newMesh); // 'fine'

//    CGAL::alpha_wrap_3(newMesh, 0.3, 0.03, newMesh); // that one takes long time
//    CGAL::alpha_wrap_3(points, polygons, 0.1, 0.001, newMesh);
//    newMesh = wrap;
}

void Building::clip_bottom(const TerrainPtr& terrain) {
    if (!_clip_bottom) return;
    if (this->has_self_intersections() && !Config::get().handleSelfIntersect) throw
                std::runtime_error(std::string("Clip error in building ID " + this->get_id() +
                                               + ". Cannot clip if there are self intersections!"));
    //-- Get terrain subset
    Mesh terrainSubsetMesh = terrain->mesh_subset(_poly);
    PMP::reverse_face_orientations(terrainSubsetMesh);

    //-- Set exact point maps
    Exact_point_map mesh1_exact_points =
            terrainSubsetMesh.add_property_map<vertex_descriptor,EK::Point_3>("v:exact_point").first;
    Exact_point_map mesh2_exact_points =
            _mesh.add_property_map<vertex_descriptor,EK::Point_3>("v:exact_point").first;
    Exact_vertex_point_map mesh1_vpm(mesh1_exact_points, terrainSubsetMesh);
    Exact_vertex_point_map mesh2_vpm(mesh2_exact_points, _mesh);

    //-- Mesh processing and clip
    PMP::remove_degenerate_faces(_mesh);
    PMP::remove_degenerate_edges(_mesh);
    if (Config::get().handleSelfIntersect) geomutils::remove_self_intersections(_mesh);
    PMP::clip(_mesh, terrainSubsetMesh, params::vertex_point_map(mesh2_vpm), params::vertex_point_map(mesh1_vpm));
}

void Building::refine() {
    typedef Mesh::Halfedge_index           halfedge_descriptor;
    typedef Mesh::Edge_index               edge_descriptor;

    const double target_edge_length = 5; //5;
    const unsigned int nb_iter =  30;   //30;

    PMP::remove_degenerate_faces(_mesh);
    /*
    if (PMP::does_self_intersect(_mesh)) {
        ++config::selfIntersecting;
        PMP::remove_self_intersections(_mesh);
    }
    */

    //-- Set the property map for constrained edges
    Mesh::Property_map<edge_descriptor,bool> is_constrained =
            _mesh.add_property_map<edge_descriptor,bool>("e:is_constrained",false).first;

    //-- Detect sharp features
    for (auto& e : edges(_mesh)) {
        halfedge_descriptor hd = halfedge(e,_mesh);
        if (!is_border(e,_mesh)) {
            double angle = CGAL::Mesh_3::dihedral_angle(_mesh.point(source(hd,_mesh)),
                                                        _mesh.point(target(hd,_mesh)),
                                                        _mesh.point(target(next(hd,_mesh),_mesh)),
                                                        _mesh.point(target(next(opposite(hd,_mesh),_mesh),_mesh)));
            if (CGAL::abs(angle)<179.5)
                is_constrained[e]=true;
        }
    }

    PMP::isotropic_remeshing(faces(_mesh), target_edge_length, _mesh,
                             PMP::parameters::number_of_iterations(nb_iter)
                                     .edge_is_constrained_map(is_constrained));

//    PMP::remove_self_intersections(_mesh);
}

void Building::translate_footprint(const double h) {
    for (auto& ring : _groundElevations) {
        for (auto& pt : ring) {
            pt += h;
        }
    }
}

void Building::check_feature_scope(const Polygon_2& otherPoly) {
    for (auto& ring: _poly.rings()) {
        for (auto& vert : ring) {
            if (geomutils::point_in_poly(vert, otherPoly))
                return;
        }
    }
//    std::cout << "Poly ID " << this->get_id() << " is outside the influ region. Deactivating." << std::endl;
    this->deactivate();
}

void Building::set_clip_flag(const bool flag) {
    _clip_bottom = flag;
}

bool Building::has_self_intersections() const {
    return PMP::does_self_intersect(_mesh);
}

void Building::set_to_zero_terrain() {
    //-- Get average footprint height
    std::vector<double> avgRings;
    for (auto& ring : _groundElevations) {
        avgRings.emplace_back(geomutils::avg(ring));
    }
    _elevation = this->get_elevation() - geomutils::avg(avgRings);
    this->set_zero_borders();
    this->reconstruct_flat_terrain();
}

double Building::sq_max_dim() {
    std::vector<double> dims;
    /*
    EPICK::Vector_2 diag(_poly.bbox().xmax() - _poly.bbox().xmin(), _poly.bbox().ymax() - _poly.bbox().ymin());
    dims.emplace_back(diag.squared_length() * pow(cos(M_PI_4), 2));
    dims.emplace_back(diag.squared_length() * pow(sin(M_PI_4), 2));
    */
    MinBbox& minBbox = this->get_min_bbox();
    dims.emplace_back(minBbox.vec1.squared_length());
    dims.emplace_back(minBbox.vec2.squared_length());
    dims.emplace_back(_elevation * _elevation);

    return *(std::max_element(dims.begin(), dims.end()));
}

void Building::get_cityjson_info(nlohmann::json& b) const {
    b["type"] = "Building";
//  b["attributes"];
//    get_cityjson_attributes(b, _attributes);
//    float hbase = z_to_float(this->get_height_base());
//    float h = z_to_float(this->get_height());
//    b["attributes"]["TerrainHeight"] = _baseElevations.back(); // temp - will calculate avg for every footprint
    b["attributes"]["measuredHeight"] = _elevation - geomutils::avg(_groundElevations[0]);
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