/*
  City4CFD
 
  Copyright (c) 2021-2025, 3D Geoinformation Research Group, TU Delft

  This file is part of City4CFD.

  City4CFD is free software: you can redistribute it and/or modify
  it under the terms of the GNU Affero General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  City4CFD is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Affero General Public License for more details.

  You should have received a copy of the GNU Affero General Public License
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
        : PolyFeature(-1), m_elevation(-global::largnum), m_height(-global::largnum),
          m_ptsPtr(std::make_shared<Point_set_3>()), m_hasFailed(false), m_reconSettings(nullptr) {}

Building::Building(const nlohmann::json& poly)
        : PolyFeature(poly, true, -1), m_elevation(-global::largnum), m_height(-global::largnum),
          m_ptsPtr(std::make_shared<Point_set_3>()), m_hasFailed(false), m_reconSettings(nullptr) {}
        // 'true' here to check for polygon simplicity

Building::Building(const Polygon_with_attr& poly)
        : PolyFeature(poly, true, -1), m_elevation(-global::largnum), m_height(-global::largnum),
          m_ptsPtr(std::make_shared<Point_set_3>()), m_hasFailed(false), m_reconSettings(nullptr) {}
        // 'true' here to check for polygon simplicity

void Building::insert_point(const Point_3& pt) {
    m_ptsPtr->insert(pt);
}

double Building::get_elevation() {
    if (m_elevation < -global::largnum + global::smallnum) // calculate if not set
        this->calc_elevation();

    return m_elevation;
}

double Building::get_height() {
    if (m_height < -global::largnum + global::smallnum) {
        m_height = this->get_elevation() - this->ground_elevation();
    }
    return m_height;
}

void Building::alpha_wrap_all(const BuildingsPtr& buildings, Mesh& newMesh) {
    typedef EPICK::FT                 FT;
    typedef std::vector<std::size_t>  CGAL_Polygon;

    //-- Make a single mesh out of all individual buildings
    std::vector<std::array<FT, 3>> points;
    std::vector<CGAL_Polygon> polygons;
    for (auto& b : buildings) {
        if (!b->is_active()) continue; // skip failed reconstructions
        auto& mesh = b->get_mesh();
        for (auto face: mesh.faces()) {
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

    //-- Perform CGAL's alpha wrapping
    const double relativeAlpha = 2000.; //1000.
    const double relativeOffset = 7000.; // 12000.
    CGAL::Bbox_3 bbox = CGAL::Polygon_mesh_processing::bbox(newMesh);
    const double diagLength = std::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                                        CGAL::square(bbox.ymax() - bbox.ymin()) +
                                        CGAL::square(bbox.zmax() - bbox.zmin()));
    const double alpha = diagLength / relativeAlpha;
    const double offset = diagLength / relativeOffset;
    CGAL::alpha_wrap_3(newMesh, alpha, offset, newMesh);

//    CGAL::alpha_wrap_3(newMesh, 2, 0.01, newMesh);   // 'coarse'
//    CGAL::alpha_wrap_3(newMesh, 1.5, 0.03, newMesh); // 'medium'
//    CGAL::alpha_wrap_3(newMesh, 0.7, 0.03, newMesh); // 'fine'

//    CGAL::alpha_wrap_3(newMesh, 0.3, 0.03, newMesh); // that one takes long time
//    CGAL::alpha_wrap_3(points, polygons, 0.1, 0.001, newMesh);
//    newMesh = wrap;
}

void Building::clip_bottom(const TerrainPtr& terrain) {
    if (!m_clipBottom) return;
    if (this->has_self_intersections() && !Config::get().handleSelfIntersect) throw
                city4cfd_error(std::string("Clip error in building ID " + this->get_id() +
                                               + ". Cannot clip if there are self intersections!"));
    //-- Get terrain subset
    Mesh terrainSubsetMesh = terrain->mesh_subset(m_poly);
    PMP::reverse_face_orientations(terrainSubsetMesh);

    //-- Set exact point maps
    Exact_point_map mesh1_exact_points =
            terrainSubsetMesh.add_property_map<vertex_descriptor,EK::Point_3>("v:exact_point").first;
    Exact_point_map mesh2_exact_points =
            m_mesh.add_property_map<vertex_descriptor,EK::Point_3>("v:exact_point").first;
    Exact_vertex_point_map mesh1_vpm(mesh1_exact_points, terrainSubsetMesh);
    Exact_vertex_point_map mesh2_vpm(mesh2_exact_points, m_mesh);

    //-- Mesh processing and clip
    PMP::remove_degenerate_faces(m_mesh);
    PMP::remove_degenerate_edges(m_mesh);
    if (Config::get().handleSelfIntersect) geomutils::remove_self_intersections(m_mesh);
    PMP::clip(m_mesh, terrainSubsetMesh, params::vertex_point_map(mesh2_vpm), params::vertex_point_map(mesh1_vpm));
}

void Building::refine() {
    typedef Mesh::Halfedge_index           halfedge_descriptor;
    typedef Mesh::Edge_index               edge_descriptor;

    const double targetEdgeLength = 5; //5;
    const unsigned int nbIter =  30;   //30;

    PMP::remove_degenerate_faces(m_mesh);
    /*
    if (PMP::does_self_intersect(m_mesh)) {
        ++config::selfIntersecting;
        PMP::remove_self_intersections(m_mesh);
    }
    */

    //-- Set the property map for constrained edges
    Mesh::Property_map<edge_descriptor,bool> isConstrained =
            m_mesh.add_property_map<edge_descriptor,bool>("e:isConstrained", false).first;

    //-- Detect sharp features
    for (auto e : edges(m_mesh)) {
        halfedge_descriptor hd = halfedge(e, m_mesh);
        if (!is_border(e, m_mesh)) {
            double angle = CGAL::Mesh_3::dihedral_angle(m_mesh.point(source(hd, m_mesh)),
                                                        m_mesh.point(target(hd, m_mesh)),
                                                        m_mesh.point(target(next(hd, m_mesh), m_mesh)),
                                                        m_mesh.point(target(next(opposite(hd, m_mesh), m_mesh), m_mesh)));
            if (CGAL::abs(angle)<179.5)
                isConstrained[e]=true;
        }
    }

    PMP::isotropic_remeshing(faces(m_mesh), targetEdgeLength, m_mesh,
                             PMP::parameters::number_of_iterations(nbIter)
                                     .edge_is_constrained_map(isConstrained));

//    PMP::remove_self_intersections(m_mesh);
}

void Building::alpha_wrap(double relativeAlpha, double relativeOffset) {
    typedef EPICK::FT                 FT;

    //-- Perform CGAL's alpha wrapping
    CGAL::Bbox_3 bbox = CGAL::Polygon_mesh_processing::bbox(m_mesh);
    const double diagLength = std::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                                        CGAL::square(bbox.ymax() - bbox.ymin()) +
                                        CGAL::square(bbox.zmax() - bbox.zmin()));
    const double alpha = diagLength / relativeAlpha;
    const double offset = diagLength / relativeOffset;
    CGAL::alpha_wrap_3(m_mesh, alpha, offset, m_mesh);
}

void Building::translate_footprint(const double h) {
    for (auto& ring : m_groundElevations) {
        for (auto& pt : ring) {
            pt += h;
        }
    }
}

// Store the reconstruction settings for a defined reconstruction (influence) region
void Building::set_reconstruction_rules(const BoundingRegion& reconRegion) {
    m_reconSettings = reconRegion.m_reconSettings;
    // set the output layer ID
    m_outputLayerID = m_reconSettings->outputLayerID;
}

// Remove stored reconstruction settings
void Building::remove_reconstruction_rules() {
    m_reconSettings.reset();
    m_outputLayerID = 0;
}

std::shared_ptr<const Config::ReconRegion> Building::get_reconstruction_settings() const {
    if (!m_reconSettings)
        throw city4cfd_error("Building " + m_id + " missing reconstruction settings."
                                 " Imported building? " + std::to_string(m_f_imported));
    return m_reconSettings;
}

bool Building::is_part_of(const Polygon_2& otherPoly) const {
    for (auto& ring: m_poly.rings()) {
        for (auto& vert : ring) {
            if (geomutils::point_in_poly(vert, otherPoly))
                return true;
        }
    }
//    std::cout << "Poly ID " << this->get_id() << " is outside the influ region." << std::endl;
    return false;
}

bool Building::has_reconstruction_region() const {
    return m_outputLayerID > 0;
}

void Building::set_clip_flag(const bool flag) {
    m_clipBottom = flag;
}

bool Building::has_self_intersections() const {
    return PMP::does_self_intersect(m_mesh);
}

void Building::mark_as_failed() {
    m_hasFailed = true;
}

bool Building::has_failed_to_reconstruct() const {
    return m_hasFailed;
}

void Building::set_to_zero_terrain() {
    //-- Get average footprint height
    std::vector<double> avgRings;
    for (auto& ring : m_groundElevations) {
        avgRings.emplace_back(geomutils::avg(ring));
    }
    m_elevation = this->get_elevation() - geomutils::avg(avgRings);
    this->set_zero_borders();
    this->reconstruct_flat_terrain();
}

double Building::sq_max_dim() {
    std::vector<double> dims;
    /*
    EPICK::Vector_2 diag(m_poly.bbox().xmax() - m_poly.bbox().xmin(), m_poly.bbox().ymax() - m_poly.bbox().ymin());
    dims.emplace_back(diag.squared_length() * pow(cos(M_PI_4), 2));
    dims.emplace_back(diag.squared_length() * pow(sin(M_PI_4), 2));
    */
    MinBbox& minBbox = this->get_min_bbox();
    dims.emplace_back(minBbox.vec1.squared_length());
    dims.emplace_back(minBbox.vec2.squared_length());
    dims.emplace_back(m_elevation * m_elevation);

    return *(std::max_element(dims.begin(), dims.end()));
}

PointSet3Ptr Building::get_points() const {
    return m_ptsPtr;
}

std::string Building::get_lod() const {
    return m_reconSettings->lod;
}

void Building::get_cityjson_cityobj_info(nlohmann::json& f) const {
    f["type"] = "Building";
//  b["attributes"];
//    get_cityjson_attributes(b, _attributes);
//    float hbase = z_to_float(this->get_height_base());
//    float h = z_to_float(this->get_height());
//    b["attributes"]["TerrainHeight"] = m_baseElevations.back(); // temp - will calculate avg for every footprint
    f["attributes"]["measuredHeight"] = m_elevation - geomutils::avg(m_groundElevations[0]);
}

void Building::get_cityjson_geomobj_info(nlohmann::json& g) const {
    g["type"] = "MultiSurface";
    g["lod"] = this->get_lod();
}

void Building::get_cityjson_semantics(nlohmann::json& g) const {
    Face_property semantics;
    auto semanticsMap = m_mesh.property_map<face_descriptor, std::string>("f:semantics");
    if (semanticsMap.has_value()) {
        semantics = semanticsMap.value();
    } else throw city4cfd_error("Semantic property map not found!");

    std::unordered_map<std::string, int> surfaceId;
    surfaceId["RoofSurface"]   = 0; g["surfaces"][0]["type"] = "RoofSurface";
    surfaceId["GroundSurface"] = 1; g["surfaces"][1]["type"] = "GroundSurface";
    surfaceId["WallSurface"]   = 2; g["surfaces"][2]["type"] = "WallSurface";

    for (auto faceIdx : m_mesh.faces()) {
        auto it = surfaceId.find(semantics[faceIdx]);
        if (it == surfaceId.end()) throw city4cfd_error("Could not find semantic attribute!");

        g["values"][faceIdx.idx()] = it->second;
    }
}

TopoClass Building::get_class() const {
    return BUILDING;
}

std::string Building::get_class_name() const {
    return "Building";
}
