#include "LoD12.h"

LoD12::LoD12(const Polygon_with_holes_2& poly,
             const std::vector<double>& base_heights,
             const std::vector<double>& building_pts)
    : _height(),  _poly(poly), _base_heights(base_heights), _building_pts(building_pts) {}

void LoD12::lod12reconstruct(Mesh& mesh) {
    //-- Reconstruction is just simple average/median/percentile
    _height = avg(_building_pts);

    this->get_mesh(mesh);
}

void LoD12::get_mesh(Mesh& mesh) {
    CDT cdt_buildings;

    //-- Map CDT and Mesh vertices
    std::map<CDT::Vertex_handle, Mesh::Vertex_index> cdtToMesh;

    //-- Temporary - will rewrite polygons
    std::vector<Polygon_2> rings;
    //- Add outer poly and holes into one data structure
    rings.push_back(_poly.outer_boundary());
    for (auto& hole : _poly.holes()) {
        rings.push_back(hole);
    }

    int polyCount = 0;
    for (auto& poly : rings) { // Loop over polys
        std::vector<Vertex_handle> cdt_handle;
        std::vector<Mesh::Vertex_index> mesh_vertex;
        int count = 0;
        ++polyCount;
        for (auto vert = poly.vertices_begin(); vert != poly.vertices_end(); ++vert) { // Loop over poly vertices
            if (polyCount == 1) {
                cdt_handle.emplace_back(cdt_buildings.insert(Point_3(vert->x(), vert->y(), _base_heights[count++])));
            } else {
                cdt_handle.emplace_back(cdt_buildings.insert(Point_3(vert->x(), vert->y(), _base_heights.back())));
            }
            mesh_vertex.emplace_back(mesh.add_vertex(cdt_handle.back()->point()));
            cdtToMesh[cdt_handle.back()] = mesh_vertex.back();

            mesh_vertex.emplace_back(mesh.add_vertex(Point_3(vert->x(), vert->y(), _height)));
        }

        //- Add constraints and create mesh faces for sides
        for (auto i = 0; i < cdt_handle.size() - 1; ++i) {
            cdt_buildings.insert_constraint(cdt_handle[i], cdt_handle[i + 1]);

            auto it1 = mesh.vertices_begin();
            auto it2 = mesh.vertices_begin();

            auto v1 = cdtToMesh[cdt_handle[i]];
            auto v2 = cdtToMesh[cdt_handle[i + 1]];

            std::advance(it1, v1.idx() + 1);
            std::advance(it2, v2.idx() + 1);

            mesh.add_face(v1, v2, *it1);
            mesh.add_face(v2, *it2, *it1);
        }
    }

    //- Handle top
    mark_domains(cdt_buildings);
    for (auto& it : cdt_buildings.finite_face_handles()) {
        if (!it->info().in_domain()) continue;

        auto it1 = mesh.vertices_begin();
        auto it2 = mesh.vertices_begin();
        auto it3 = mesh.vertices_begin();

        std::advance(it1, cdtToMesh[it->vertex(0)]);
        std::advance(it2, cdtToMesh[it->vertex(1)]);
        std::advance(it3, cdtToMesh[it->vertex(2)]);

//        mesh.add_face(*it1, *it2, *it3); // Don't need bottom face
        mesh.add_face(*std::next(it1), *std::next(it2), *std::next(it3));
    }

}

double LoD12::get_height() {
    return _height;
}


