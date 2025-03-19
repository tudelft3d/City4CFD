/*
  City4CFD
 
  Copyright (c) 2021-2025, 3D Geoinformation Research Group, TU Delft

  This file is part of City4CFD.

  City4CFD is free software: you can redistribute it and/or modify
  it under the terms of the GNU Affero General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  City4CFD is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Affero General Public License for more details.

  You should have received a copy of the GNU Affero General Public License
  along with City4CFD.  If not, see <http://www.gnu.org/licenses/>.

  For any information or further details about the use of City4CFD, contact
  Ivan Pađen
  <i.paden@tudelft.nl>
  3D Geoinformation Research Group
  Delft University of Technology
*/

#include "BoundingRegion.h"

#include "geomutils.h"
#include "Building.h"

#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Alpha_shape_vertex_base_2.h>
#include <CGAL/Alpha_shape_face_base_2.h>
#include <CGAL/centroid.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/Delaunay_triangulation_2.h>

BoundingRegion::BoundingRegion(std::shared_ptr<Config::ReconRegion> reconRegion)
: m_reconSettings(reconRegion)
{}

//-- Operators to read bounded region if explicitly defined in config
//   called through boost:apply_visitor
void BoundingRegion::operator()(double radius) {
    geomutils::make_round_poly(global::nullPt, radius, m_boundingRegion);
}

void BoundingRegion::operator()(Polygon_2& poly) {
    m_boundingRegion = poly;
}

//-- Reconstruction (influ) regions related functions
double
BoundingRegion::calc_influ_region_bpg(const DT& dt, BuildingsPtr& buildings) {
    double influRegionRadius;
    //-- Find building where the point of interest lies in and define radius of interest with BPG
    for (auto& f : buildings) {
        if (geomutils::point_in_poly(global::nullPt, f->get_poly())) {
            f->calc_footprint_elevation_nni(dt);
            try {
                f->set_reconstruction_rules(*this); // temporary add settings for this reconstruction
                f->reconstruct();
                f->remove_reconstruction_rules();
            } catch (std::exception& e) {
                std::cerr << std::endl << "Error: " << e.what() << std::endl;
                throw city4cfd_error("Impossible to automatically determine influence region");
            }
            double maxDim = sqrt(f->sq_max_dim());
            influRegionRadius = maxDim * (3. + m_reconSettings->bpgInfluExtra); //- BPG by Liu hardcoded

            f->clear_feature();

            geomutils::make_round_poly(global::nullPt, influRegionRadius, m_boundingRegion);
            return maxDim;
        }
    }
    throw std::invalid_argument("Point of interest does not belong to any building! "
                                "Impossible to determine influence region");
}

void BoundingRegion::calc_influ_region_bpg(const double maxDim) {
    double influRegionRadius = maxDim * (3. + m_reconSettings->bpgInfluExtra); //- BPG by Liu hardcoded
    geomutils::make_round_poly(global::nullPt, influRegionRadius, m_boundingRegion);
}

//-- Boundary related functions
void BoundingRegion::calc_bnd_bpg(const Polygon_2& influRegionPoly,
                                  const BuildingsPtr& buildings) {
    double angle = std::atan2(Config::get().flowDirection.y(), Config::get().flowDirection.x());

    //-- Find candidate points for the AABB
    std::vector<Point_2> candidatePts;
    for (auto& pt : influRegionPoly) {
        candidatePts.push_back(pt);
    }
    //- Add buildings points that may end up outside the influ
    //  region poly due to the way influ region is calculated
    //- While looping through buildings, also find the highest one
    double hMax = 0;
    double elevationMax = 0;
    for (auto& b : buildings) {
        if (!b->is_active()) continue;
        if (b->get_height() > hMax) {
            hMax = b->get_height();
            elevationMax = b->get_elevation();
        }
        for (auto& pt : b->get_poly().outer_boundary()) {
            if (!geomutils::point_in_poly(pt, influRegionPoly)) {
                candidatePts.push_back(pt);
            }
        }
    }
    //-- Axes aligning transformation - get points
    //   aligned with the flow direction
    for (auto& pt : candidatePts)
        pt = geomutils::rotate_pt(pt, -angle);

    //-- Get the bnd poly
    Polygon_2 localPoly = this->calc_bnd_poly(candidatePts, hMax, angle);

    //-- Set the top
    Config::get().topHeight = elevationMax + hMax * (Config::get().bpgDomainSize.back() - 1);

    //-- Blockage ratio handling
    if (Config::get().bpgBlockageRatioFlag) {
        std::cout << "\nCalculating blockage ratio for flow direction (" << Config::get().flowDirection
                  << ")" << std::endl;

//    double blockRatio = this->calc_blockage_ratio_from_chull(buildings, angle, localPoly);
//    double blockRatio = this->calc_blockage_ratio_from_ashape(buildings, angle, localPoly);
        double blockRatio = this->calc_blockage_ratio_comb(buildings, angle, localPoly);
//    double blockRatio = this->calc_blockage_ratio_from_ashape_alt(buildings, angle, localPoly);
//    double blockRatio = this->calc_blockage_ratio_from_edges(buildings, angle, localPoly);

        std::cout << "    Blockage ratio is: " << blockRatio << std::endl;
        if (blockRatio > Config::get().bpgBlockageRatio) {
            std::cout << "INFO: Blockage ratio is more than " << Config::get().bpgBlockageRatio * 100
                      << "%. Expanding domain cross section to meet the guideline"
                      << std::endl;
            double expRatio = std::sqrt(blockRatio / 0.03);
            //-- Recalculate the bnd poly and height with new values
            localPoly.clear();
            localPoly = this->calc_bnd_poly(candidatePts, hMax, angle, expRatio);
            Config::get().topHeight = hMax * Config::get().bpgDomainSize.back() * expRatio;
        }
    }
    //-- Return the points back to global coordinates
    for (auto& pt : localPoly) m_boundingRegion.push_back(geomutils::rotate_pt(pt, angle));
}

// Check if all points of this region fall within otherRegion
bool BoundingRegion::is_subset_of(const BoundingRegion& otherRegion) const {
    auto& poly = otherRegion.get_bounding_region();
    for (auto& vert : m_boundingRegion) {
        if (!geomutils::point_in_poly(vert, poly))
            return false;
    }
    return true;
}

Polygon_2& BoundingRegion::get_bounding_region() {
    return m_boundingRegion;
}
const Polygon_2& BoundingRegion::get_bounding_region() const {
    return m_boundingRegion;
}

Polygon_2 BoundingRegion::calc_bnd_poly(const std::vector<Point_2>& candidatePts,
                                        const double hMax, const double angle, const double enlargeRatio) {
    Polygon_2 localPoly;
    if (Config::get().bpgDomainType == ROUND) {
        Point_2 localPOI = geomutils::rotate_pt(global::nullPt, -angle);

        double bndRadius = 0;
        for (auto& pt: candidatePts) {
            double sqDist = CGAL::squared_distance(localPOI, pt);
            if (sqDist > bndRadius * bndRadius) bndRadius = sqrt(sqDist);
        }
        bndRadius = (bndRadius + hMax * Config::get().bpgDomainSize.front()) * enlargeRatio;
        geomutils::make_round_poly(localPOI, bndRadius, localPoly);

        std::cout << "Calculated boundary radius is: "
                  << bndRadius << std::endl;
    } else {
        //-- Get bbox
        Polygon_2 bbox = geomutils::calc_bbox_poly(candidatePts);

        if (Config::get().bpgDomainType == RECTANGLE) {
            //-- Construct enlargement vector
            std::vector<Vector_2> translateBoundary;
            translateBoundary.emplace_back(Vector_2(-Config::get().bpgDomainSize[0], 0));
            translateBoundary.emplace_back(Vector_2(0, -Config::get().bpgDomainSize[1] * enlargeRatio));
            translateBoundary.emplace_back(Vector_2(Config::get().bpgDomainSize[2], 0));
            translateBoundary.emplace_back(Vector_2(0, Config::get().bpgDomainSize[1] * enlargeRatio));

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
        } else if (Config::get().bpgDomainType == OVAL) {
            std::vector<double> bpgDomainDist;
            bpgDomainDist.push_back(Config::get().bpgDomainSize[1]); // Down
            bpgDomainDist.push_back(Config::get().bpgDomainSize[2]); // Right (Back)
            bpgDomainDist.push_back(Config::get().bpgDomainSize[1]); // Up
            bpgDomainDist.push_back(Config::get().bpgDomainSize[0]); // Left (Front)

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

/*
 * Calculate blockage ratio from a convex hull of building points projected to a plane normal to the flow direction.
 *
 * Points are transformed and projected to yz plane for blockage ratio calculation. Coordinates are then flipped to the xy plane
 * to work with CDT which is in xy projection plane.
 *
 * Exact for LoD1.2, approximation for others; fast
 */
/*
double BoundingRegion::calc_blockage_ratio_from_chull(const BuildingsPtr& buildings, const double angle,
                                                      Polygon_2& localPoly) const {
    CDT projCDT;
    for (auto& b : buildings) {
        if (!b->is_active()) continue;
        std::vector<Point_2> buildingPts;

        //-- Project building mesh points onto yz plane
        this->project_mesh_pts(b->get_mesh(), angle, buildingPts);

        //-- Approximate the blockArea of the projection with the convex hull
        this->chull_to_cdt(buildingPts, projCDT);
    }
    //-- Get the area of constrained regions
    double blockArea       = 0;
    double domainCrossArea = 0;
    this->calc_cross_sec_areas(projCDT, localPoly, buildings, blockArea, domainCrossArea);

    //-- Return the blockage ration
    return blockArea / domainCrossArea;
}
*/

/*
 * Calculate blockage ratio from an alpha shape of building points projected to a plane normal to the flow direction.
 *
 * Points are transformed and projected to yz plane for blockage ratio calculation. Coordinates are then flipped to the xy plane
 * to work with CDT which is in xy projection plane.
 *
 * Can be slow in case of many edges
 */
/*
double BoundingRegion::calc_blockage_ratio_from_ashape(const BuildingsPtr& buildings, const double angle,
                                                       Polygon_2& localPoly) const {
    CDT projCDT;
    for (auto& b : buildings) {
        std::vector<Point_2> buildingPts;
        if (!b->is_active()) continue;

        //-- Project building mesh points onto yz plane
        this->project_mesh_pts(b->get_mesh(), angle, buildingPts);

        const double aVal = b->get_height() * 12;
//        const double aVal = 1100;
        //-- Calculate the alpha shape
        this->ashape_to_cdt(buildingPts, projCDT, aVal);
    }
    //-- Get the area of constrained regions
    double blockArea       = 0;
    double domainCrossArea = 0;
    this->calc_cross_sec_areas(projCDT, localPoly, buildings, blockArea, domainCrossArea);

    //-- Return the blockage ration
    return blockArea / domainCrossArea;
}
*/

/*
 * Calculate blockage ratio from a combination of chull (for LoD1.2) and ashape to a plane normal to the flow direction.
 *
 * Points are transformed and projected to yz plane for blockage ratio calculation. Coordinates are then flipped to the xy plane
 * to work with CDT which is in xy projection plane.
 */
double BoundingRegion::calc_blockage_ratio_comb(const BuildingsPtr& buildings, const double angle,
                                                Polygon_2& localPoly) const {
    CDT projCDT;
    for (auto& b : buildings) {
        if (!b->is_active()) continue;
        std::vector<Point_2> buildingPts;

        //-- Project building mesh points onto yz plane
        this->project_mesh_pts(b->get_mesh(), angle, buildingPts);

        //-- Approximate the blockArea of the projection with the convex hull
        if (b->is_imported()) {
//            const double aVal = 1100;
            const double aVal = b->get_height() * 12;
            this->ashape_to_cdt(buildingPts, projCDT, aVal);
        } else {
            this->chull_to_cdt(buildingPts, projCDT);
        }
    }
    //-- Get the area of constrained regions
    double blockArea       = 0;
    double domainCrossArea = 0;
    this->calc_cross_sec_areas(projCDT, localPoly, buildings, blockArea, domainCrossArea);

    //-- Return the blockage ration
    return blockArea / domainCrossArea;
}

/*
 * Calculate blockage ratio from an alpha shape of *ALL* points projected to a plane normal to the flow direction.
 *
 * Points are transformed and projected to yz plane for blockage ratio calculation. Coordinates are then flipped to the xy plane
 * to work with CDT which is in xy projection plane.
 */
/*
double BoundingRegion::calc_blockage_ratio_from_ashape_alt(const BuildingsPtr& buildings, const double angle,
                                                           Polygon_2& localPoly) const {
    //-- Project building pts onto 2d plane
    CDT projCDT;
    std::vector<Point_2> buildingPts;
    for (auto& b : buildings) {
        if (!b->is_active()) continue;

        //-- Project building mesh points onto yz plane
        this->project_mesh_pts(b->get_mesh(), angle, buildingPts);
    }
    const double aVal = 1100;
    //-- Calculate the alpha shape
    this->ashape_to_cdt(buildingPts, projCDT, aVal);

    //-- Get the area of constrained regions
    double blockArea       = 0;
    double domainCrossArea = 0;
    this->calc_cross_sec_areas(projCDT, localPoly, buildings, blockArea, domainCrossArea);

    //-- Return the blockage ration
    return blockArea / domainCrossArea;
}
*/

/*
 * Calculate blockage ratio by constraining building edges.
 *
 * Points are transformed and projected to yz plane for blockage ratio calculation. CDT is projected onto xy plane.
 *
 * It is easier to project all coordinates to yz and then flip to xy, instead of
 * writing a separate triangulation structures with yz projection traits.
 *
 * Takes a long time to calculate as there are many building edges.
 */
/*
double BoundingRegion::calc_blockage_ratio_from_edges(const BuildingsPtr& buildings, const double angle,
                                                      Polygon_2& localPoly) const {
    Converter<EPICK, EPECK> to_exact;
    CDT projCDT;
    for (auto& b : buildings) {
        CDT buildingCDT;
        if (!b->is_active()) continue;
        auto& mesh = b->get_mesh();
        //-- Construct projected CDT from all building pts - time-consuming part
        for (auto& edge : mesh.edges()) {
            //-- Orient points normal to the flow direction
            ePoint_3 pt1 = to_exact(geomutils::rotate_pt_xy(mesh.point(mesh.vertex(edge, 0)), -angle));
            ePoint_3 pt2 = to_exact(geomutils::rotate_pt_xy(mesh.point(mesh.vertex(edge, 1)), -angle));
            //-- Flip xyz to yz
            pt1 = ePoint_3(pt1.y(), pt1.z(), 0);
            pt2 = ePoint_3(pt2.y(), pt2.z(), 0);

            buildingCDT.insert_constraint(pt1, pt2);
        }
        //-- Add only outer shell of the building triangulation to overall projected CDT
        geomutils::tag_layers(buildingCDT);
        for (auto it = buildingCDT.edges_begin(); it != buildingCDT.edges_end(); ++it) {
            if (!buildingCDT.is_constrained(*it)) continue;
            auto fh1 = it->first;
            auto fh2 = fh1->neighbor(it->second);
            if (fh1->info().in_domain_noholes() && !fh2->info().in_domain_noholes() ||
               !fh1->info().in_domain_noholes() &&  fh2->info().in_domain_noholes()) {
                auto vs = fh1->vertex(fh1->cw(it->second));
                auto vt = fh1->vertex(fh1->ccw(it->second));
                projCDT.insert_constraint(buildingCDT.point(vs), buildingCDT.point(vt));
            }
        }
    }
    //-- Get the area of constrained regions
    double blockArea       = 0;
    double domainCrossArea = 0;
    this->calc_cross_sec_areas(projCDT, localPoly, buildings, blockArea, domainCrossArea);

    //-- Return the blockage ration
    return blockArea / domainCrossArea;
}
*/

void BoundingRegion::project_mesh_pts(const Mesh& mesh, const double angle,
                                      std::vector<Point_2>& buildingPts) const {
    for (auto buildingPt : mesh.points()) {
        buildingPt = geomutils::rotate_pt_xy(buildingPt, -angle);
        Point_2 pt(buildingPt.y(), buildingPt.z());
        buildingPts.push_back(pt);
    }
}

void BoundingRegion::chull_to_cdt(const std::vector<Point_2>& buildingPts, CDT& projCDT) const {
    Polygon_2 projConvHull;
    CGAL::convex_hull_2(buildingPts.begin(), buildingPts.end(), std::back_inserter(projConvHull));

    //-- Add convex hull to CDT
    Polygon_3 projConvHullCDT;
    for (auto& pt : projConvHull) {
        projConvHullCDT.push_back(ePoint_3(pt.x(), pt.y(), 0));
    }
    projCDT.insert_constraint(projConvHullCDT.begin(), projConvHullCDT.end(), true);
}

void BoundingRegion::ashape_to_cdt(const std::vector<Point_2>& buildingPts, CDT& projCDT, const double aVal) const {
    typedef CGAL::Alpha_shape_vertex_base_2<EPICK>                            Vba;
    typedef CGAL::Alpha_shape_face_base_2<EPICK>                              Fba;
    typedef CGAL::Triangulation_data_structure_2<Vba,Fba>                     Tds;
    typedef CGAL::Delaunay_triangulation_2<EPICK ,Tds>                        Triangulation_2;
    typedef CGAL::Alpha_shape_2<Triangulation_2>                              Alpha_shape_2;

    Alpha_shape_2 A(buildingPts.begin(), buildingPts.end(), aVal, Alpha_shape_2::GENERAL);

    for (auto it = A.finite_edges_begin(); it != A.finite_edges_end(); ++it) {
        auto fh1 = it->first;
        auto fh2 = it->first->neighbor(it->second);
        if (A.classify(fh1) == Alpha_shape_2::INTERIOR && A.classify(fh2) != Alpha_shape_2::INTERIOR ||
            A.classify(fh1) != Alpha_shape_2::INTERIOR && A.classify(fh2) == Alpha_shape_2::INTERIOR) {
            auto vs = fh1->vertex(fh1->cw(it->second));
            auto vt = fh1->vertex(fh1->ccw(it->second));
            ePoint_3 pt1(A.point(vs).x(), A.point(vs).y(), 0);
            ePoint_3 pt2(A.point(vt).x(), A.point(vt).y(), 0);
            projCDT.insert_constraint(pt1, pt2);
        }
    }
}

/*
 * Get the cross-section areas for the blocked domain, and for the full domain.
 */
void BoundingRegion::calc_cross_sec_areas(CDT& projCDT, Polygon_2& localPoly, const BuildingsPtr& buildings,
                                          double& blockArea, double& domainCrossArea) const {
    //-- Mark constrained regions
    geomutils::mark_domains(projCDT);

    //-- Get surface area of constrained regions
    for (auto& face : projCDT.finite_face_handles()) {
        if (!face->info().in_domain_noholes()) continue;
        std::vector<Point_2> pts;
        for (int i = 0; i < 3; ++i)
            pts.emplace_back(CGAL::to_double(face->vertex(i)->point().x()),
                             CGAL::to_double(face->vertex(i)->point().y()));

        blockArea += CGAL::area(pts[0], pts[1], pts[2]);
    }
    //-- Get the area of the domain cross section at the influence region
    Polygon_2 bbox = geomutils::calc_bbox_poly(localPoly);
    //- Get the mean height of building bases
    std::vector<double> footprintElevations;
    for (auto& b : buildings) {
        //- One elevation per building as there might be different num of points per polygon
        footprintElevations.push_back(b->ground_elevation());
    }
    double medianFootprintHeight = geomutils::percentile(footprintElevations, 0.5);

    //- Calculate domain cross section area
    domainCrossArea = std::sqrt(Vector_2(bbox.vertex(3) - bbox.vertex(0)).squared_length())
                      * (Config::get().topHeight - medianFootprintHeight);
}