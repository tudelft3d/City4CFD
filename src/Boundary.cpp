#include "Boundary.h"

//todo maybe merge top and sides into one
//todo bump down to topoFeature?
//todo break apart this code into something more maintainable
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

//-- Deactivate point cloud points that are out of bounds
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

void Boundary::set_bounds_to_terrain(Point_set_3& pointCloud, const Polygon_2& bndPoly,
                                     const Polygon_2& pcBndPoly, const Polygon_2& startBufferPoly) {
    //-- Remove points out of the boundary region
    Boundary::set_bounds_to_pc(pointCloud, pcBndPoly);

    //-- Add outer points to match the domain size to prescribed one
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

void BoundingRegion::calc_bnd_bpg(const Polygon_2& influRegionPoly,
                                  const Buildings& buildings) {
    double angle = std::atan2(config::flowDirection.y(), config::flowDirection.x());

    //-- Find candidate points for the AABB
    std::vector<Point_2> candidatePts;
    for (auto& pt : influRegionPoly) {
        candidatePts.push_back(pt);
    }
    //- Add buildings points that may end up outside the influ
    //  region poly due to the way influ region is calculated
    //- While looping through buildings, also find the highest one
    double hMax = 0;
    for (auto& b : buildings) {
        if (!b->is_active()) continue;
        if (b->get_height() > hMax) hMax = b->get_height();

        for (auto& pt : b->get_poly().outer_boundary()) {
            if (!geomtools::point_in_poly(pt, influRegionPoly)) {
                candidatePts.push_back(pt);
            }
        }
    }

    //-- Axes aligning transformation
    for (auto& pt : candidatePts)
        pt = geomtools::rotate_pt(pt, -angle);

    //-- Get the bnd poly
    Polygon_2 localPoly = this->calc_bnd_poly(candidatePts, hMax, angle);

    //-- Set the top
    config::topHeight = hMax * config::bpgDomainSize.back();

    //-- Blockage ratio handling
    std::cout << "Calculating blockage ratio" << std::endl;
    double blockRatio = this->calc_blockage_ratio(buildings, angle, localPoly);
    std::cout << "Blockage ratio is: " << blockRatio << std::endl;
    if (blockRatio > 0.03) {
        std::cout << "Blockage ratio is more than 3%. Expanding domain cross section to meet the guideline" << std::endl;
        double expRatio = std::sqrt(blockRatio / 0.03);
        //-- Recalculate the bnd poly and height with new values
        localPoly.clear();
        localPoly = this->calc_bnd_poly(candidatePts, hMax, angle, expRatio);
        config::topHeight = hMax * config::bpgDomainSize.back() * expRatio;
    }
    //-- Return the points back to global coordinates
    for (auto& pt : localPoly) _boundingRegion.push_back(geomtools::rotate_pt(pt, angle));
}

Polygon_2& BoundingRegion::get_bounding_region() {
    return _boundingRegion;
}
const Polygon_2& BoundingRegion::get_bounding_region() const {
    return _boundingRegion;
}

Polygon_2 BoundingRegion::calc_bnd_poly(const std::vector<Point_2>& candidatePts,
                                        const double hMax, const double angle, const double enlargeRatio) const {
    Polygon_2 localPoly;
    if (config::bpgDomainType == ROUND) {
        Point_2 localPOI = geomtools::rotate_pt(config::pointOfInterest, -angle);

        double bndRadius = 0;
        for (auto& pt: candidatePts) {
            double sqDist = CGAL::squared_distance(localPOI, pt);
            if (sqDist > bndRadius * bndRadius) bndRadius = sqrt(sqDist);
        }
        bndRadius = (bndRadius + hMax * config::bpgDomainSize.front()) * enlargeRatio;
        geomtools::make_round_poly(localPOI, bndRadius, localPoly);

        std::cout << "--- INFO: Defined boundary radius is: "
                  << bndRadius << " ---" << std::endl;

    } else {
        //-- Get bbox
        Polygon_2 bbox = geomtools::calc_bbox_poly(candidatePts);

        if (config::bpgDomainType == RECTANGLE) {
            //-- Construct enlargement vector
            std::vector<Vector_2> translateBoundary;
            translateBoundary.emplace_back(Vector_2(-config::bpgDomainSize[0], 0));
            translateBoundary.emplace_back(Vector_2(0, -config::bpgDomainSize[1] * enlargeRatio));
            translateBoundary.emplace_back(Vector_2(config::bpgDomainSize[2], 0));
            translateBoundary.emplace_back(Vector_2(0, config::bpgDomainSize[1] * enlargeRatio));

            //-- Additional enlargement for the bbox in case of large blockage ratio
            std::vector<Vector_2> addEnlargementVector;
            addEnlargementVector.emplace_back((bbox[0] - CGAL::midpoint(bbox[0], bbox[3])) * (enlargeRatio - 1)); //front, -front, -front
            addEnlargementVector.emplace_back(addEnlargementVector.front());
            addEnlargementVector.emplace_back(-addEnlargementVector.front());
            addEnlargementVector.emplace_back(-addEnlargementVector.front());

            int i = 0;
            for (auto& pt: bbox) {
                pt += hMax * (translateBoundary[i] + translateBoundary[(i + 1) % 4]) + addEnlargementVector[i];
//                pt += hMax * (translateBoundary[i] + translateBoundary[(i + 1) % 4]);
                localPoly.push_back(pt);
                ++i;
            }
        } else if (config::bpgDomainType == OVAL) {
            std::vector<double> bpgDomainDist;
            bpgDomainDist.push_back(config::bpgDomainSize[1]); // Down
            bpgDomainDist.push_back(config::bpgDomainSize[2]); // Right (Back)
            bpgDomainDist.push_back(config::bpgDomainSize[1]); // Up
            bpgDomainDist.push_back(config::bpgDomainSize[0]); // Left (Front)

            Point_2 centerPt = CGAL::centroid(bbox.begin(), bbox.end());
            std::vector<double> distances;
            for (int i = 0; i < 4; ++i) {
                Point_2 pt = CGAL::midpoint(bbox[i], bbox[(i + 1) % 4]);
                distances.emplace_back(sqrt(CGAL::squared_distance(pt, centerPt)) + bpgDomainDist[i] * hMax);
            }
            std::vector<double> radiuses;
            if (distances[0] > distances[2])
                radiuses.push_back(distances[0]); // Sides
            else
                radiuses.push_back(distances[2]);
            radiuses.push_back(distances[3]);    // Front
            radiuses.push_back(distances[1]);    // Back

            radiuses.front() *= enlargeRatio;

            //-- Make front half of the oval domain
            geomtools::make_round_poly(centerPt, radiuses[1], radiuses[0],
                                       180, M_PI/180, M_PI_2, localPoly);
            //-- Make the back half of the oval domain
            geomtools::make_round_poly(centerPt, radiuses[2], radiuses[0],
                                       180, M_PI/180, 3*M_PI_2, localPoly);
        }
    }
    return localPoly;
}

//-- Blockage ratio based on the flow direction
double BoundingRegion::calc_blockage_ratio(const Buildings& buildings, const double angle, Polygon_2& localPoly) const {
    //-- We're working with a local coordinate system, normal to yz plane
    CDT projCDT;
    double blockArea = 0;
    for (auto& b : buildings) {
        if (!b->is_active()) continue;
        double height = b->get_height();
        std::vector<std::vector<double>> baseHeights = b->get_base_heights();

        //-- Project building pts onto 2d plane
        std::vector<Point_2> buildingPts;
        int i = 0;
        for (auto pt : b->get_poly().outer_boundary()) {
            pt = geomtools::rotate_pt(pt, -angle); // Get points to local system
            buildingPts.emplace_back(Point_2(pt.y(), baseHeights.front()[i++]));
            buildingPts.emplace_back(Point_2(pt.y(), height));
        }

        //-- Approximate the blockArea of the projection with the convex hull
        Polygon_2 projConvHull;
        CGAL::convex_hull_2(buildingPts.begin(), buildingPts.end(), std::back_inserter(projConvHull));

        //-- Add convex hull to CDT
        Polygon_3 projConvHullCDT;
        for (auto& pt : projConvHull) {
            projConvHullCDT.push_back(ePoint_3(pt.x(), pt.y(), 0));
        }
        projCDT.insert_constraint(projConvHullCDT.begin(), projConvHullCDT.end(), true);
    }
    //-- Mark constrained regions
    geomtools::mark_domains(projCDT);

    //-- Calculate blockArea of constrained region
    for (auto& face : projCDT.finite_face_handles()) {
        if (!face->info().in_domain_noholes()) continue;
        std::vector<Point_2> pts;
        for (int i = 0; i < 3; ++i)
            pts.emplace_back(Point_2(CGAL::to_double(face->vertex(i)->point().x()),
                                        CGAL::to_double(face->vertex(i)->point().y())));

        blockArea += CGAL::area(pts[0], pts[1], pts[2]);
    }
    //-- Get blockArea of the domain cross section at the influence region
    Polygon_2 bbox = geomtools::calc_bbox_poly(localPoly);
    double domainCrossArea = std::sqrt(Vector_2(bbox.vertex(3) - bbox.vertex(0)).squared_length()) * config::topHeight;

    //-- Return the blockage ration
    return blockArea / domainCrossArea;
}