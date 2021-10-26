#include "Boundary.h"

//todo maybe merge top and sides into one
Boundary::Boundary()
    : PolyFeature(), _sideOutputPts() {}

Boundary::Boundary(const int outputLayerID)
    : PolyFeature(outputLayerID), _sideOutputPts() {}

Boundary::~Boundary() = default;

//-- Static member definition
std::vector<Point_3> Boundary::_outerPts;
double               Boundary::_outerBndHeight;

void Boundary::set_bnd_poly(Polygon_2& bndPoly, Polygon_2& pcBndPoly, Polygon_2& startBufferPoly) {
    Polygon_2 bufferPoly;
    //-- Set bndPoly and pcBndPoly depending on the buffer setup
    if (config::domainBuffer > -g_largnum) {
        double bufferLen = config::domainBuffer / 100.;
        Point_2 center = CGAL::centroid(bndPoly.begin(), bndPoly.end(),
                                        CGAL::Dimension_tag<0>());
        for (auto& pt: bndPoly)
            bufferPoly.push_back(pt + (pt - center) * bufferLen);

        if (config::domainBuffer < 0) {
            startBufferPoly = bufferPoly;
        } else {
            startBufferPoly = bndPoly;
            bndPoly = bufferPoly;
        }

        for (auto& pt : startBufferPoly)
            pcBndPoly.push_back(pt - (pt - center) * 0.05 * std::abs(bufferLen));
    } else {
        startBufferPoly = bndPoly;
        pcBndPoly = startBufferPoly;
    }
}

void Boundary::set_bounds_to_pc(Point_set_3& pointCloud, const Polygon_2& pcBndPoly) {
    //-- Remove points out of the boundary region
    auto it = pointCloud.points().begin();
    int count = 0;
    while (it != pointCloud.points().end()) {
        if (!geomtools::point_in_poly(*it, pcBndPoly)) {
            pointCloud.remove(pointCloud.begin() + count);
        } else {
            ++it;
            ++count;
        }
    }
    pointCloud.collect_garbage(); // Free removed points from the memory
}

//-- Deactivate point cloud points that are out of bounds
void Boundary::set_bounds_to_terrain(Point_set_3& pointCloud, const Polygon_2& bndPoly,
                                     const Polygon_2& pcBndPoly, const Polygon_2& startBufferPoly) {
    //-- Remove points out of the boundary region
    Boundary::set_bounds_to_pc(pointCloud, pcBndPoly);

    //-- Add points to match the domain size to prescribed one
    SearchTree searchTree(pointCloud.points().begin(),pointCloud.points().end());

    std::vector<double> bndHeights;
    geomtools::interpolate_poly_from_pc(bndPoly, bndHeights, pointCloud);
    _outerBndHeight = geomtools::avg(bndHeights); // Height for buffer (for now) - average of outer pts

    if (config::domainBuffer > -g_largnum + g_smallnum) {
        for (auto& pt : startBufferPoly) {
            pointCloud.insert(Point_3(pt.x(), pt.y(), _outerBndHeight));
        }
        for (auto& pt : bndPoly) {
            _outerPts.emplace_back(Point_3(pt.x(), pt.y(), _outerBndHeight));
            pointCloud.insert(_outerPts.back());
        }
    } else {
        int i = 0;
        for (auto& pt : bndPoly) {
            _outerPts.emplace_back(Point_3(pt.x(), pt.y(), bndHeights[i++]));
            pointCloud.insert(_outerPts.back());
        }
    }
    _outerPts.push_back(_outerPts.front());
}

void Boundary::prep_output() {
    _sideOutputPts = _outerPts;
}

//-- Find all outerPts along the defined edge
void Boundary::prep_output(Vector_2 edge) {
    edge /= sqrt(edge.squared_length());

    //-- Search the outerPts for the same vector
    for (int i = 0; i < _outerPts.size() - 1; ++i) {
        Vector_2 checkEdge(_outerPts[i + 1].x() - _outerPts[i].x(),
                           _outerPts[i + 1].y() - _outerPts[i].y());
        checkEdge /= sqrt(checkEdge.squared_length());

        if (edge * checkEdge > 1 - g_smallnum && edge * checkEdge < 1 + g_smallnum) {
            _sideOutputPts.push_back(_outerPts[i]);
            _sideOutputPts.push_back(_outerPts[i + 1]);

            bool collinear = true;
            while (collinear) {
                int j = i + 2;
                if (_outerPts.begin() + j == _outerPts.end()) break;

                Vector_2 nextEdge(_outerPts[j].x() - _outerPts[i + 1].x(),
                                  _outerPts[j].y() - _outerPts[i + 1].y());
                nextEdge /= sqrt(nextEdge.squared_length());

                if (nextEdge * edge > 1 - g_smallnum && nextEdge * edge < 1 + g_smallnum) {
                    _sideOutputPts.push_back(_outerPts[j]);
                } else {
                    collinear = false;
                }
                ++i;
            }
            return;
        }
    }
    throw std::runtime_error("Cannot find side for output!");
}

std::vector<double> Boundary::get_domain_bbox() {
    //todo: proper bbox calculation
    double maxx(-g_largnum), maxy(-g_largnum), maxz(-g_largnum);
    double minx(g_largnum),  miny(g_largnum),  minz(g_largnum);

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
Sides::Sides(const int outputLayerID)
    : Boundary(outputLayerID) {
}

Sides::~Sides() = default;

void Sides::reconstruct() {
    std::vector<Mesh::vertex_index> mesh_vertex_side;

    //-- Add mesh vertices and store them in a vector
    for (auto it = _sideOutputPts.begin(); it != _sideOutputPts.end(); ++it) {
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

TopoClass Sides::get_class() const {
    return SIDES;
}

std::string Sides::get_class_name() const {
    return "Sides";
}

//-- TOP CLASS --//
Top::Top(const int outputLayerID)
    : Boundary(outputLayerID) {
}

Top::~Top() = default;

void Top::reconstruct() {
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

//-- BOUNDED REGION CLASS--//
BoundingRegion::BoundingRegion() = default;
BoundingRegion::~BoundingRegion() = default;

//-- Operators to read bounded region if explicitly defined in config
void BoundingRegion::operator()(double radius) {
    geomtools::make_round_poly(config::pointOfInterest, radius, _boundingRegion);
}

void BoundingRegion::operator()(Polygon_2& poly) {
    _boundingRegion = poly;
}

void BoundingRegion::operator()(std::string& polyPath) {
    JsonPolygons influJsonPoly;
    IO::read_geojson_polygons(polyPath, influJsonPoly);

    for (auto& coords : influJsonPoly.front()->front()) { // I know it should be only 1 polygon with 1 ring
            _boundingRegion.push_back(Point_2(coords[0], coords[1]));
    }
    CGAL::internal::pop_back_if_equal_to_front(_boundingRegion);
    if (_boundingRegion.is_clockwise_oriented()) _boundingRegion.reverse_orientation();
}

void
BoundingRegion::calc_influ_region_bpg(const DT& dt, const Point_set_3& pointCloudBuildings, Buildings& buildings) {
#ifndef NDEBUG
    assert(boost::get<bool>(config::influRegionConfig));
#endif
    double influRegionRadius;
    //-- Find building where the point of interest lies in and define radius of interest with BPG
    SearchTree searchTreeBuildings;
    searchTreeBuildings.insert(pointCloudBuildings.points().begin(), pointCloudBuildings.points().end()); //TODO Maybe save as member variable

    bool foundBuilding = false;
    for (auto& f : buildings) {
        if (geomtools::point_in_poly(config::pointOfInterest, f->get_poly())) {
            f->calc_footprint_elevation_nni(dt);
            try {
                f->reconstruct(searchTreeBuildings);
            } catch (std::exception& e) {
                std::cerr << std::endl << "Error: " << e.what() << std::endl;
                throw std::runtime_error("Impossible to automatically determine influence region");
            }
            influRegionRadius = f->max_dim() * 3.; //- BPG by Liu

            f->clear_feature();
            foundBuilding = true;
            break;
        }
    }
    if (!foundBuilding)
        throw std::invalid_argument("Point of interest does not belong to any building! "
                                    "Impossible to determine influence region.");

    geomtools::make_round_poly(config::pointOfInterest, influRegionRadius, _boundingRegion);
}

void BoundingRegion::calc_bnd_bpg(double hMax,
                                  const Polygon_2& influRegionPoly,
                                  const Buildings& buildings) {
    //-- Set relations for coordinate transformation
    double sinPhi = config::flowDirection.y() / sqrt(config::flowDirection.squared_length());
    double cosPhi = config::flowDirection.x() / sqrt(config::flowDirection.squared_length());
    CGAL::Aff_transformation_2<EPICK> rotate(CGAL::ROTATION, -sinPhi, cosPhi);
    CGAL::Aff_transformation_2<EPICK> rotate_back(CGAL::ROTATION, sinPhi, cosPhi);

    //-- Find candidate points for the AABB
    std::vector<Point_2> candidatePts;
    for (auto& pt : influRegionPoly) {
        candidatePts.push_back(pt);
    }
    //- Add buildings points that may end up outside the influ
    //  region poly due to the way influ region is calculated
    for (auto& b : buildings) {
        if (!b->is_active()) continue;
        for (auto& pt : b->get_poly().outer_boundary()) {
            if (!geomtools::point_in_poly(pt, influRegionPoly)) {
                candidatePts.push_back(pt);
            }
        }
    }

    //-- Create domain polygons
    Polygon_2 localPoly;
    if (config::bpgDomainType == ROUND) {
        double bndRadius = 0;
        for (auto& pt: candidatePts) {
            double sqDist = CGAL::squared_distance(config::pointOfInterest, pt);
            if (sqDist > bndRadius * bndRadius) bndRadius = sqrt(sqDist);
        }
        bndRadius += hMax * config::bpgDomainSize.front();
        geomtools::make_round_poly(config::pointOfInterest, bndRadius, localPoly);

        std::cout << "--- INFO: Defined boundary radius is: "
                  << bndRadius << " ---" << std::endl;

    } else {
        //-- Axes aligning transformation
        for (auto& pt : candidatePts) {
            pt = rotate(pt);
        }

        //-- Get bbox and translation vector
        Polygon_2 bbox = geomtools::calc_bbox_poly(candidatePts);
        auto& translateBoundary = config::enlargeDomainVec;

        if (config::bpgDomainType == RECTANGLE) {

            int i = 0;
            for (auto& pt: bbox) {
                pt += hMax * (translateBoundary[i] + translateBoundary[(i + 1) % 4]);
                localPoly.push_back(rotate_back(pt));
                ++i;
            }
        } else if (config::bpgDomainType == ELLIPSE) {
            // well it's todo eh
        }
    }
    //-- Set the top
    config::topHeight = hMax * config::bpgDomainSize.back();

    //todo blockage ratio calc here
    _boundingRegion = localPoly;
}

Polygon_2& BoundingRegion::get_bounding_region() {
    return _boundingRegion;
}
const Polygon_2& BoundingRegion::get_bounding_region() const {
    return _boundingRegion;
}