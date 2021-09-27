#include "Terrain.h"
#include "SurfaceLayer.h"

Terrain::Terrain()
    : TopoFeature() {}

Terrain::Terrain(int pid)
    : TopoFeature(pid) {}

Terrain::~Terrain() {
    for (auto& layer : _surfaceLayersTerrain) {
        layer = nullptr; delete layer;
    }
}

void Terrain::threeDfy(const Point_set_3& pointCloud, const std::vector<PolyFeature*>& features) {
    int count = 0;
    //-- Constrain surface layers
    for (auto& f : features) {
        if (!f->is_active() || f->get_class() != SURFACELAYER) continue;
        std::cout << "Constraining feature " << count++ << " of class" << f->get_class_name() << std::endl;
        this->constrain_footprint(f->get_poly(), f->get_base_heights());
    }
    std::cout << "Done constraining" << std::endl;

    //-- Add ground points from the point cloud to terrain
    this->set_cdt(pointCloud); // CDT's got to go first if performing smoothing

    //-- Smoothing
    this->smooth(pointCloud);

    //-- Constrain buildings
    for (auto& feature : features) {
        //debug
        if (!feature->is_active() || feature->get_class() != BUILDING) continue;
        std::cout << "Constrained feature " << count++ << " of class" << feature->get_class_name() << std::endl;
        this->constrain_footprint(feature->get_poly(), feature->get_base_heights());
    }

    std::cout << "Done constructing CDT" << std::endl;

    //-- Mark surface layer
    geomtools::mark_surface_layers(this->get_cdt(), features);

    //-- Store the CGAL terrain and surface layers in separate meshes
    this->create_mesh();

    std::cout << "Done making mesh" << std::endl;
}

void Terrain::constrain_footprint(const Polygon_with_holes_2& poly,
                                  const std::vector<std::vector<double>>& heights) {
    // necessary reminder to rewrite polygons to something more comfortable
    //-- Temporary - will rewrite polygons
    std::vector<Polygon_2> rings;
    //- Add outer poly and holes into one data structure
    rings.push_back(poly.outer_boundary());
    for (auto& hole: poly.holes()) {
        rings.push_back(hole);
    }

    int polyCount = 0;
    for (auto& ring : rings) {
        //-- Add ring points
        int count = 0;
        Polygon_3 pts;
        for (auto& polyVertex : ring) {
            pts.push_back(Point_3(polyVertex.x(), polyVertex.y(), heights[polyCount][count++]));
        }

        //-- Set added points as constraints
        _cdt.insert_constraint(pts.begin(), pts.end());

        ++polyCount;
    }
}

void Terrain::set_cdt(const Point_set_3& pointCloud) {
    _cdt.insert(pointCloud.points().begin(), pointCloud.points().end());
}

//-- Taken from CGAL's example
void Terrain::smooth(const Point_set_3& pointCloud) {
    // Smooth heights with 5 successive Gaussian filters
#ifdef CGAL_LINKED_WITH_TBB
    using Concurrency_tag = CGAL::Parallel_tag;
#else
    using Concurrency_tag = CGAL::Sequential_tag;
#endif
    double spacing = CGAL::compute_average_spacing<Concurrency_tag>(pointCloud, 6);
    spacing *= 2;
    double gaussian_variance = 4 * spacing * spacing;
    for (CDT::Vertex_handle vh : _cdt.finite_vertex_handles())
    {
        double z = vh->point().z();
        double total_weight = 1;
        CDT::Vertex_circulator circ = _cdt.incident_vertices (vh),
                start = circ;
        do
        {
            if (!_cdt.is_infinite(circ))
            {
                double sq_dist = CGAL::squared_distance (vh->point(), circ->point());
                double weight = std::exp(- sq_dist / gaussian_variance);
                z += weight * circ->point().z();
                total_weight += weight;
            }
        }
        while (++ circ != start);
        z /= total_weight;
        vh->point() = Point_3 (vh->point().x(), vh->point().y(), z);
    }
}

void Terrain::create_mesh() {
    geomtools::cdt_to_mesh(_cdt, _mesh);

    // -- Create placeholders for surface layer mesh
    int layerNum = 1; // Config or some other way
    for (int i = 1; i <= layerNum; ++i) {
        SurfaceLayer* layer = new SurfaceLayer(i);
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

const std::vector<SurfaceLayer*>& Terrain::get_surface_layers() const {
    return _surfaceLayersTerrain;
}