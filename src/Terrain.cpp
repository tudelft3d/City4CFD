#include "Terrain.h"

Terrain::Terrain()
    : TopoFeature() {}

Terrain::Terrain(int pid)
    : TopoFeature(pid) {}

void Terrain::threeDfy(const Point_set_3& pointCloud, const std::vector<PolyFeature*>& features) {
    //-- TODO: smoothing

    //-- Add buildings as constraints to the terrain
    for (auto &feature : features) {
        if (feature->is_active()) {
            this->constrain_footprint(feature->get_poly(), feature->get_base_heights());
        }
    }

    //-- Add ground points from the point cloud to terrain
    this->set_cdt(pointCloud);

    //-- Store the CGAL terrain in the triangle-vertex data structure
    this->triangulate_mesh();
}

void Terrain::constrain_footprint(const Polygon_with_holes_2& poly, const std::vector<double>& heights) {
    std::vector<Vertex_handle> vh;
    //-- Add outer ring points
    int count = 0;
    for (auto& polyVertex : poly.outer_boundary()) {
        vh.emplace_back(_cdt.insert(Point_3(polyVertex.x(), polyVertex.y(), heights[count++])));
    }

    //-- Set added points as constraints
    for (auto i = 0; i < vh.size() - 1; ++i) {
        _cdt.insert_constraint(vh[i], vh[i + 1]);
    }

    //-- Add inner ring points
    if (poly.has_holes()) {
        for (auto& polyHoles : poly.holes()) {
            vh.clear();
            for (auto& holeVertex : polyHoles) {
                vh.emplace_back(_cdt.insert(Point_3(holeVertex.x(), holeVertex.y(), heights.back())));
            }
            //-- Set added points as constraints
            for (auto i = 0; i < vh.size() - 1; ++i) {
                _cdt.insert_constraint(vh[i], vh[i + 1]);
            }
        }
    }

}

void Terrain::set_cdt(const Point_set_3& pointCloud) {
    _cdt.insert(pointCloud.points().begin(), pointCloud.points().end());
}

void Terrain::triangulate_mesh() { // Don't know how to handle it for now
    std::map<CDT::Vertex_handle, int> indices;
    std::vector<Mesh::vertex_index> mesh_vertex;
    std::vector<Mesh::face_index> face_index;
    mesh_vertex.reserve(_cdt.dimension());

    int counter = 0;
    for (CDT::Finite_vertices_iterator it = _cdt.finite_vertices_begin(); it != _cdt.finite_vertices_end(); ++it) {
        mesh_vertex.emplace_back(_mesh.add_vertex(it->point()));
        //        outstream << it->point() << std::endl;
        indices.insert(std::pair<CDT::Vertex_handle, int>(it, counter++));
    }

    for (CDT::Finite_faces_iterator it = _cdt.finite_faces_begin(); it != _cdt.finite_faces_end(); ++it) {

        int v1 = indices[it->vertex(0)];
        int v2 = indices[it->vertex(1)];
        int v3 = indices[it->vertex(2)];
        _mesh.add_face(mesh_vertex[v1], mesh_vertex[v2], mesh_vertex[v3]);
    }
}

void Terrain::output_feature(std::string& fs,
                             std::string& bs,
                             std::unordered_map<std::string, unsigned long>& dPts) const {
    bs += "\ng Terrain";
    get_obj(_mesh, fs, bs, dPts);
};

const CDT& Terrain::get_cdt() const {
    return _cdt;
}

TopoClass Terrain::get_class() const {
    return TERRAIN;
}
std::string Terrain::get_class_name() const {
    return "Terrain";
}
