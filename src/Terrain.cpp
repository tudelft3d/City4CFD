#include "Terrain.h"

Terrain::Terrain()
    : TopoFeature() {}

Terrain::Terrain(int pid)
    : TopoFeature(pid) {}

void Terrain::threeDfy(const Point_set_3& pointCloud, const std::vector<PolyFeature*>& features) {
    //-- Add ground points from the point cloud to terrain
    this->set_cdt(pointCloud); // CDT's got to go first if performing smoothing

    //-- Smoothing
    this->smooth(pointCloud);

    //-- Add buildings as constraints to the terrain
    for (auto& feature : features) {
        if (feature->is_active()) {
            this->constrain_footprint(feature->get_poly(), feature->get_base_heights());
        }
    }

    //-- Add ground points from the point cloud to terrain
//    this->set_cdt(pointCloud);

    //-- Store the CGAL terrain in the triangle-vertex data structure
    this->create_mesh();
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

    //-- Add inner ring points - TODO:necessary reminder to rewrite polygons into something more comfortable
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

//-- Taken from CGAL's example
void Terrain::smooth(const Point_set_3& pointCloud) {
    // Smooth terrain
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
    cdt_to_mesh(_cdt, _mesh);
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