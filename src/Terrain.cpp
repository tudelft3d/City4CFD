#include "Terrain.h"

#include "geomutils.h"
#include "io.h"
#include "SurfaceLayer.h"

Terrain::Terrain()
        : TopoFeature(0), _cdt(), _surfaceLayersTerrain(), _constrainedPolys() {}

Terrain::Terrain(int pid)
        : TopoFeature(pid), _cdt(), _surfaceLayersTerrain(), _constrainedPolys() {}

Terrain::~Terrain() = default;

void Terrain::set_cdt(const Point_set_3& pointCloud) {
    Converter<EPICK, EPECK> to_exact;

    std::cout << "\n    Triangulating" << std::endl;
    int count = 0;
    int pcSize = pointCloud.size();
    std::vector<ePoint_3> pts; //todo check progress bar
    IO::print_progress_bar(0);
    for (auto& pt : pointCloud.points()) {
        pts.emplace_back(to_exact(pt));

        IO::print_progress_bar(99 * count++ / pcSize);
    }
    IO::print_progress_bar(99);
    _cdt.insert(pts.begin(), pts.end());
    IO::print_progress_bar(100); std::clog << std::endl;

    //-- Smoothing
    if (config::smoothTerrain) {
        std::cout << "\n    Smoothing" << std::endl;
        geomutils::smooth_dt<CDT, EPECK>(pointCloud, _cdt);
    }
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
    int count = 0;
    int numFeatures = _constrainedPolys.size();

    std::cout << "\n    Constraining polygons" << std::endl;
    IO::print_progress_bar(0);
    for (auto& ring : _constrainedPolys) {
        //-- Set added points as constraints
        _cdt.insert_constraint(ring.begin(), ring.end(), true);

        IO::print_progress_bar(100 * count++ / numFeatures);
    }
    IO::print_progress_bar(100); std::clog << std::endl;
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