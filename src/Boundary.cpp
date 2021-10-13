#include "Boundary.h"

Boundary::Boundary()
    : PolyFeature() {}

Boundary::Boundary(const int outputLayerID)
    : PolyFeature(outputLayerID) {}

Boundary::~Boundary() = default;

//-- Static member definition
std::vector<Point_3> Boundary::_outerPts;
double               Boundary::_outerBndHeight;
Polygon_2            Boundary::_bndPoly;

//-- Deactivate point cloud points that are out of bounds
void Boundary::set_bounds_to_pc(Point_set_3& pointCloud) { // Will try to template it to include CDT
    //- 80% of the total domain size. The rest is left for the buffer zone
    auto it = pointCloud.points().begin();
    int count = 0;
    while (it != pointCloud.points().end()) {
        if (!geomtools::point_in_circle(*it, config::pointOfInterest, 0.8 * config::dimOfDomain)) {
            pointCloud.remove(pointCloud.begin() + count);
        } else {
            ++it;
            ++count;
        }
    }
    pointCloud.collect_garbage(); // Free removed points from the memory
}

//-- Add buffer between the terrain and domain edge
void Boundary::add_buffer(Point_set_3& pointCloud) {
    const int nPts      = 360; // Hardcoded
    const double angInt = 2 * M_PI / (double)nPts;
    double ang = 0;
    for (auto i = 0; i < nPts; ++i) {
        double xPt = config::pointOfInterest.x() + config::dimOfDomain * cos(ang + angInt);
        double yPt = config::pointOfInterest.y() + config::dimOfDomain * sin(ang + angInt);
        ang = ang + angInt;
        Point_3 pt(xPt, yPt, 0.0); // Height set at 0 for now. Could be average of the edge or average of the whole terrain
        _outerPts.push_back(pt); // for top and sides
        pointCloud.insert(pt);
    }
    _outerPts.push_back(_outerPts[0]); // Put first point at the end to close the loop
}

std::vector<double> Boundary::get_domain_bbox() {
    //todo: proper bbox calculation
    double maxx(-infty), maxy(-infty), maxz(-infty);
    double minx(infty),  miny(infty),  minz(infty);

    for (auto& pt : _outerPts) {
        if (pt.x() > maxx) maxx = pt.x();
        else if (pt.x() < minx) minx = pt.x();

        if (pt.y() > maxy) maxy = pt.y();
        else if (pt.y() < miny) miny = pt.y();

//        if (pt.z() > maxz) maxz = pt.z();
//        else if (pt.z() < minz) minz = pt.z();
    }
    minz = -5;
    maxz = 100;

    return std::vector<double> {minx, miny, minz, maxx, maxy, maxz};
}

//-- TEMP
void Boundary::get_cityjson_info(nlohmann::json& b) const {
    //temp
}

void Boundary::get_cityjson_semantics(nlohmann::json& g) const {
    //temp
}

std::string Boundary::get_cityjson_primitive() const {
    return "";
};

//-- SIDES CLASS --//
Sides::Sides()
    : Boundary(2) {
}

Sides::~Sides() = default;

void Sides::threeDfy() {
    std::vector<Mesh::vertex_index> mesh_vertex_side;

    //-- Add mesh vertices and store them in a vector
    for (auto it = _outerPts.begin(); it != _outerPts.end(); ++it) {
        mesh_vertex_side.emplace_back(_mesh.add_vertex(*it));
        mesh_vertex_side.emplace_back(_mesh.add_vertex(Point_3(it->x(), it->y(), config::topHeight)));
    }

    //-- Add middle top point to mesh
//    Mesh::vertex_index topMiddlePoint = _meshTop.add_vertex(Point_3(bndInfo.xcent, bndInfo.ycent, bndInfo.height));

    //-- Add mesh faces for side
    for (auto i = 0; i < mesh_vertex_side.size() - 3; i= i + 2) {
        // -- i + 1 is i lifted up
        int v1 = i;
        int v2 = i + 2;

        _mesh.add_face(mesh_vertex_side[v2], mesh_vertex_side[v1], mesh_vertex_side[v1 + 1]);
        _mesh.add_face(mesh_vertex_side[v2 + 1], mesh_vertex_side[v2], mesh_vertex_side[v1 + 1]);
    }
}

// TBD
void Sides::set_bnd_poly(double radius, SearchTree& searchTree) {
    // todo either gonna make buffer the extrusion of domain, or shrink the domain to fit buffer
    geomtools::make_round_poly(config::pointOfInterest, radius, _poly);

    //-- Height for buffer - average of outer pts
    this->calc_footprint_elevation_from_pc(searchTree);
    _outerBndHeight = geomtools::avg(_base_heights.front());
}

TopoClass Sides::get_class() const {
    return SIDES;
}

std::string Sides::get_class_name() const {
    return "Sides";
}

//-- TOP CLASS --//
Top::Top()
    : Boundary(3) {
}

Top::~Top() = default;

void Top::threeDfy() {
    std::vector<Mesh::vertex_index> mesh_vertex_top;

    //-- Top is done by making a CDT of outerPts
    CDT cdt_top;
    for (auto& pt : _outerPts) {
       cdt_top.insert(ePoint_3(pt.x(), pt.y(), config::topHeight));
    }

    //-- Add mesh faces for top
    geomtools::cdt_to_mesh(cdt_top, _mesh);
}

TopoClass Top::get_class() const {
    return TOP;
}

std::string Top::get_class_name() const {
    return "Top";
}

//-- INFLU REGION CLASS--//
InfluRegion::InfluRegion() = default;
InfluRegion::~InfluRegion() = default;

void InfluRegion::operator()(double radius) {
    geomtools::make_round_poly(config::pointOfInterest, radius, _influRegion);
}

void InfluRegion::operator()(Polygon_2& poly) {
    _influRegion = poly;
}

void InfluRegion::operator()(std::string& polyPath) {
    JsonPolygons influJsonPoly;
    IO::read_geojson_polygons(polyPath, influJsonPoly);

    for (auto& coords : influJsonPoly.front()->front()) { // I know it should be only 1 polygon with 1 ring
            _influRegion.push_back(Point_2(coords[0], coords[1]));
    }
    CGAL::internal::pop_back_if_equal_to_front(_influRegion);
    if (_influRegion.is_clockwise_oriented()) _influRegion.reverse_orientation();
}

void InfluRegion::operator()(Point_set_3& pointCloud, Point_set_3& pointCloudBuildings, Buildings& buildings) {
#ifndef NDEBUG
    assert(boost::get<bool>(config::influenceRegion));
#endif
    double influRegionRadius = this->calc_influ_region_radius_bpg(pointCloud,
                                                                  pointCloudBuildings,
                                                                  buildings);
    geomtools::make_round_poly(config::pointOfInterest, influRegionRadius, _influRegion);
}

const Polygon_2& InfluRegion::get_influ_region() const {
    return _influRegion;
}

double
InfluRegion::calc_influ_region_radius_bpg(Point_set_3& pointCloud, Point_set_3& pointCloudBuildings, Buildings& buildings) {
    double influRegionRadius;
    //-- Find building where the point of interest lies in and define radius of interest with BPG
    SearchTree searchTree, searchTreeBuildings;
    searchTree.insert(pointCloud.points().begin(), pointCloud.points().end());
    searchTreeBuildings.insert(pointCloudBuildings.points().begin(), pointCloudBuildings.points().end()); //TODO Maybe save as member variable

    bool foundBuilding = false;
    for (auto& f : buildings) {
        if (geomtools::point_in_poly(config::pointOfInterest, f->get_poly())) {
            f->calc_footprint_elevation_from_pc(searchTree); //- Quick calculation from point cloud before CDT is constructed
            try {
                f->threeDfy(searchTreeBuildings);
            } catch (std::exception& e) {
                std::cerr << std::endl << "Error: " << e.what() << std::endl;
                throw std::runtime_error("Impossible to automatically determine influence region");
            }
            influRegionRadius = f->max_dim() * 3.; //- BPG by Liu
            f->clear_feature(); //- They will be properly interpolated from DT later

            foundBuilding = true;
            break;
        }
    }
    if (!foundBuilding) throw std::invalid_argument("Point of interest does not belong to any building! Impossible to determine influence region.");

    return influRegionRadius;
}