#include "LoD12.h"

LoD12::LoD12(const Polygon_with_holes_2& poly,
             const std::vector<std::vector<double>>& base_heights,
             const std::vector<double>& building_pts)
    : _height(), _poly(poly), _baseHeights(base_heights), _buildingPts(building_pts) {}

void LoD12::lod12reconstruct(Mesh& mesh, double& height) {
    _height = geomtools::percentile(_buildingPts, config::buildingPercentile);
    height = _height;

    this->create_mesh(mesh);
}

void LoD12::create_mesh(Mesh& mesh) {
    // Test - add semantics with face properties
    auto surfaceType = mesh.add_property_map<face_descriptor , std::string>("f:semantics", "").first;
    face_descriptor fIdx;

    CDT cdt_buildings;

    //-- Map CDT and Mesh vertices
    std::map<CDT::Vertex_handle, Mesh::Vertex_index> cdtToMesh;

    int polyCount = 0;
    for (auto& poly : _poly.rings()) { // Loop over polys
        std::vector<Vertex_handle> cdt_handle;
        std::vector<Mesh::Vertex_index> mesh_vertex;
        int count = 0;
        for (auto vert = poly.vertices_begin(); vert != poly.vertices_end(); ++vert) { // Loop over poly vertices
            cdt_handle.emplace_back(cdt_buildings.insert(ePoint_3(vert->x(), vert->y(), _baseHeights[polyCount][count])));
            mesh_vertex.emplace_back(mesh.add_vertex(Point_3(vert->x(), vert->y(), _baseHeights[polyCount][count++])));
            cdtToMesh[cdt_handle.back()] = mesh_vertex.back();
            mesh_vertex.emplace_back(mesh.add_vertex(Point_3(vert->x(), vert->y(), _height)));
        }
        cdt_handle.emplace_back(cdt_handle.front());
        mesh_vertex.emplace_back(mesh_vertex.front());
        mesh_vertex.emplace_back(mesh_vertex.front() + 1);

        //- Add constraints and create mesh faces for sides
        for (auto i = 0; i < cdt_handle.size() - 1; ++i) {
            cdt_buildings.insert_constraint(cdt_handle[i], cdt_handle[i + 1]);

            auto it1 = mesh.vertices_begin();
            auto it2 = mesh.vertices_begin();

            auto v1 = cdtToMesh[cdt_handle[i]];
            auto v2 = cdtToMesh[cdt_handle[i + 1]];

            std::advance(it1, v1.idx() + 1);
            std::advance(it2, v2.idx() + 1);

            fIdx = mesh.add_face(v1, v2, *it1);   surfaceType[fIdx] = "WallSurface";
            fIdx = mesh.add_face(v2, *it2, *it1); surfaceType[fIdx] = "WallSurface";
        }
        ++polyCount;
    }

    //- Handle top
    geomtools::mark_domains(cdt_buildings);
    for (auto& it : cdt_buildings.finite_face_handles()) {
        if (!it->info().in_domain()) continue;

        auto it1 = mesh.vertices_begin();
        auto it2 = mesh.vertices_begin();
        auto it3 = mesh.vertices_begin();

        std::advance(it1, cdtToMesh[it->vertex(0)]);
        std::advance(it2, cdtToMesh[it->vertex(1)]);
        std::advance(it3, cdtToMesh[it->vertex(2)]);

//        mesh.add_face(*it1, *it3, *it2); // Bottom face
//        surfaceType[fIdx] = "GroundSurface";

        fIdx = mesh.add_face(*std::next(it1), *std::next(it2), *std::next(it3));
        surfaceType[fIdx] = "RoofSurface";
    }
}

double LoD12::get_height() const {
    return _height;
}