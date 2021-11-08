#include "BoundingRegion.h"

#include "geomutils.h"
#include "Building.h"

BoundingRegion::BoundingRegion() = default;
BoundingRegion::~BoundingRegion() = default;

//-- Operators to read bounded region if explicitly defined in config
void BoundingRegion::operator()(double radius) {
    geomutils::make_round_poly(config::pointOfInterest, radius, _boundingRegion);
}

void BoundingRegion::operator()(Polygon_2& poly) {
    _boundingRegion = poly;
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
        if (geomutils::point_in_poly(config::pointOfInterest, f->get_poly())) {
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
                                    "Impossible to determine influence region");

    geomutils::make_round_poly(config::pointOfInterest, influRegionRadius, _boundingRegion);
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
            if (!geomutils::point_in_poly(pt, influRegionPoly)) {
                candidatePts.push_back(pt);
            }
        }
    }

    //-- Axes aligning transformation
    for (auto& pt : candidatePts)
        pt = geomutils::rotate_pt(pt, -angle);

    //-- Get the bnd poly
    Polygon_2 localPoly = this->calc_bnd_poly(candidatePts, hMax, angle);

    //-- Set the top
    config::topHeight = hMax * config::bpgDomainSize.back();

    //-- Blockage ratio handling
    std::cout << "\nCalculating blockage ratio for flow direction (" << config::flowDirection
              << ")" << std::endl;
    double blockRatio = this->calc_blockage_ratio(buildings, angle, localPoly);
    std::cout << "    Blockage ratio is: " << blockRatio << std::endl;
    if (config::bpgBlockageRatioFlag && blockRatio > config::bpgBlockageRatio) {
        std::cout << "INFO: Blockage ratio is more than " << config::bpgBlockageRatio * 100
                  << "%. Expanding domain cross section to meet the guideline"
                  << std::endl;
        double expRatio = std::sqrt(blockRatio / 0.03);
        //-- Recalculate the bnd poly and height with new values
        localPoly.clear();
        localPoly = this->calc_bnd_poly(candidatePts, hMax, angle, expRatio);
        config::topHeight = hMax * config::bpgDomainSize.back() * expRatio;
    }
    //-- Return the points back to global coordinates
    for (auto& pt : localPoly) _boundingRegion.push_back(geomutils::rotate_pt(pt, angle));
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
        Point_2 localPOI = geomutils::rotate_pt(config::pointOfInterest, -angle);

        double bndRadius = 0;
        for (auto& pt: candidatePts) {
            double sqDist = CGAL::squared_distance(localPOI, pt);
            if (sqDist > bndRadius * bndRadius) bndRadius = sqrt(sqDist);
        }
        bndRadius = (bndRadius + hMax * config::bpgDomainSize.front()) * enlargeRatio;
        geomutils::make_round_poly(localPOI, bndRadius, localPoly);

        std::cout << "Calculated boundary radius is: "
                  << bndRadius << std::endl;

    } else {
        //-- Get bbox
        Polygon_2 bbox = geomutils::calc_bbox_poly(candidatePts);

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
            geomutils::make_round_poly(centerPt, radiuses[1], radiuses[0],
                                       180, M_PI/180, M_PI_2, localPoly);
            //-- Make the back half of the oval domain
            geomutils::make_round_poly(centerPt, radiuses[2], radiuses[0],
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
            pt = geomutils::rotate_pt(pt, -angle); // Get points to local system
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
    geomutils::mark_domains(projCDT);

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
    Polygon_2 bbox = geomutils::calc_bbox_poly(localPoly);
    double domainCrossArea = std::sqrt(Vector_2(bbox.vertex(3) - bbox.vertex(0)).squared_length())
                           * config::topHeight;

    //-- Return the blockage ration
    return blockArea / domainCrossArea;
}