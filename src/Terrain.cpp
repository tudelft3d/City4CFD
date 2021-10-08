#include "Terrain.h"

Terrain::Terrain()
    : TopoFeature(0) {}

Terrain::Terrain(int pid)
    : TopoFeature(pid) {}

Terrain::~Terrain() = default;

void Terrain::set_cdt(const Point_set_3& pointCloud) {
    Converter<EPICK, EPECK> to_exact;
    for (auto& pt : pointCloud.points()) _cdt.insert(to_exact(pt));

    //-- Smoothing
#ifdef SMOOTH
    geomtools::smooth_dt<CDT, EPECK>(pointCloud, _cdt);
#endif
}

void Terrain::threeDfy(const Point_set_3& pointCloud, const PolyFeatures& features) {
    int count = 0;
    //-- Constrain surface layers
    for (auto& f : features) {
        if (!f->is_active() || f->get_class() != SURFACELAYER) continue;
        std::cout << "Constraining feature " << count++ << " of class" << f->get_class_name() << std::endl;
        this->constrain_footprint(f->get_poly(), f->get_base_heights());
    }
    std::cout << "Done constraining" << std::endl;

    //-- Constrain buildings
    for (auto& feature : features) {
        if (!feature->is_active() || feature->get_class() != BUILDING) continue;
        std::cout << "Constrained feature " << count++ << " of class" << feature->get_class_name() << std::endl;
        this->constrain_footprint(feature->get_poly(), feature->get_base_heights());
    }

    std::cout << "Done constructing CDT" << std::endl;

    //-- Mark surface layer
    geomtools::mark_domains(this->get_cdt(), features);

    //-- Store the CGAL terrain and surface layers in separate meshes
    this->create_mesh();

    std::cout << "Done making mesh" << std::endl;
}

void Terrain::constrain_footprint(const Polygon_with_holes_2& poly,
                                  const std::vector<std::vector<double>>& heights) {
    int polyCount = 0;
    for (auto& ring : poly.rings()) {
        //-- Add ring points
        int count = 0;
        Polygon_3 pts;
        std::vector<Vertex_handle> vertex;
        for (auto& polyVertex : ring) {
            pts.push_back(ePoint_3(polyVertex.x(), polyVertex.y(), heights[polyCount][count++]));
        }
        //-- Set added points as constraints
        _cdt.insert_constraint(pts.begin(), pts.end(), true);

        ++polyCount;
    }
}

void Terrain::create_mesh() {
    geomtools::cdt_to_mesh(_cdt, _mesh);

    // -- Create placeholders for surface layer mesh
    int layerNum = 2; // Config or some other way
    for (int i = 4; i < layerNum + 4; ++i) { // Surface layers start with output ID 4
        auto layer = std::make_shared<SurfaceLayer>(i);
        geomtools::cdt_to_mesh(_cdt, layer->get_mesh(), i);
        _surfaceLayersTerrain.push_back(layer);
    }
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

TopoClass Terrain::get_class() const {
    return TERRAIN;
}

std::string Terrain::get_class_name() const {
    return "Terrain";
}

const SurfaceLayers& Terrain::get_surface_layers() const {
    return _surfaceLayersTerrain;
}