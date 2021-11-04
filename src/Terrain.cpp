#include "Terrain.h"

#include "geomutils.h"
#include "io.h"
#include "SurfaceLayer.h"

Terrain::Terrain()
    : TopoFeature(0) {}

Terrain::Terrain(int pid)
    : TopoFeature(pid) {}

Terrain::~Terrain() = default;

void Terrain::set_cdt(const Point_set_3& pointCloud) {
    Converter<EPICK, EPECK> to_exact;

    std::cout << "--- Triangulating terrain ---" << std::endl;
    int count = 0;
    int pcSize = pointCloud.size();
    IO::print_progress_bar(0);
    for (auto& pt : pointCloud.points()) {
        _cdt.insert(to_exact(pt));

        IO::print_progress_bar(100 * count++ / pcSize);
    }
    IO::print_progress_bar(100); std::clog << std::endl;

    //-- Smoothing
#ifdef SMOOTH
    geomutils::smooth_dt<CDT, EPECK>(pointCloud, _cdt);
#endif
}

void Terrain::constrain_features(const PolyFeatures& features) {
    int count = 0;
    int numFeatures = features.size();
    //-- Constrain polygon features;

    std::cout << "--- Constraining polygons ---" << std::endl;
    IO::print_progress_bar(0);
    for (auto& f : features) {
        int polyCount = 0;
        if (!f->is_active()) continue;
        for (auto& ring : f->get_poly().rings()) {
            auto& heights = f->get_base_heights();
            //-- Add ring points
            int i = 0;
            Polygon_3 pts;
            std::vector<Vertex_handle> vertex;
            for (auto& polyVertex : ring) {
                pts.push_back(ePoint_3(polyVertex.x(), polyVertex.y(), heights[polyCount][i++]));
            }
            //-- Set added points as constraints
            _cdt.insert_constraint(pts.begin(), pts.end(), true);

            ++polyCount;
        }

        IO::print_progress_bar(100 * count++ / numFeatures);
    }
    IO::print_progress_bar(100); std::clog << std::endl;
    std::clog << "--> Num of constrained polygons: " << count << " <--" << std::endl;
    std::cout << "--- Done constraining ---" << std::endl;
}

void Terrain::create_mesh(const PolyFeatures& features) {
    //-- Mark surface layer
    geomutils::mark_domains(this->get_cdt(), features);

    //-- Create the mesh for the terrain
    geomutils::cdt_to_mesh(_cdt, _mesh);

    // -- Surface layer meshes are stored here
    for (int i : config::surfaceLayerIDs) {
        auto layer = std::make_shared<SurfaceLayer>(i);
        geomutils::cdt_to_mesh(_cdt, layer->get_mesh(), i); // Create mesh for surface layers
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