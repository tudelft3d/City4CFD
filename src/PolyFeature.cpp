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

#ifdef CITY4CFD_VERBOSE
  #define CITY4CFD_POLYFEATURE_VERBOSE
#endif

#include "PolyFeature.h"

#include "geomutils.h"
#include "Building.h"

#include <CGAL/convex_hull_2.h>
#include <CGAL/min_quadrilateral_2.h>
#include <CGAL/natural_neighbor_coordinates_2.h>
#include <CGAL/Straight_skeleton_builder_2.h>
#include <CGAL/compute_outer_frame_margin.h>
#include <CGAL/create_offset_polygons_from_polygon_with_holes_2.h>
#include <CGAL/Polygon_set_2.h>

#ifndef NDEBUG
#include <CGAL/Barycentric_coordinates_2/triangle_coordinates_2.h>
#endif

PolyFeature::PolyFeature()
        : TopoFeature(), m_poly(), m_groundElevations(), m_polyInternalID(new_internal_id()),
          m_groundElevation(-global::largnum), m_minBbox() {}

PolyFeature::PolyFeature(const int outputLayerID)
        : TopoFeature(outputLayerID), m_poly(), m_groundElevations(), m_polyInternalID(new_internal_id()),
          m_groundElevation(-global::largnum), m_minBbox() {}

PolyFeature::PolyFeature(const nlohmann::json& poly, const bool checkSimplicity)
        : TopoFeature(), m_groundElevations(), m_polyInternalID(new_internal_id()),
          m_groundElevation(-global::largnum), m_minBbox() {
    this->parse_json_poly(poly, checkSimplicity);
}

PolyFeature::PolyFeature(const nlohmann::json& poly, const bool checkSimplicity, const int outputLayerID)
        : PolyFeature(poly, checkSimplicity) {
    m_outputLayerID = outputLayerID;
    if (m_outputLayerID >= s_numOfOutputLayers) s_numOfOutputLayers = m_outputLayerID + 1;
}

PolyFeature::PolyFeature(const nlohmann::json& poly, const int outputLayerID)
        : PolyFeature(poly, false, outputLayerID) {}

PolyFeature::PolyFeature(const Polygon_with_attr& poly, const bool checkSimplicity)
        : TopoFeature(), m_groundElevations(), m_polyInternalID(new_internal_id()),
          m_groundElevation(-global::largnum), m_minBbox() {
    bool isOuterRing = true;
    for (auto& ring : poly.polygon.rings()) {
        Polygon_2 tempPoly;
        const double minDist = 0.0001; // hardcoded
        Point_2 prev(global::largnum, global::largnum);
        for (auto& pt: ring.vertices()) {
            if (CGAL::squared_distance(pt, prev) > minDist) {
                tempPoly.push_back(pt);
                prev = pt;
            }
        }
        geomutils::pop_back_if_equal_to_front(tempPoly);
        if (tempPoly.size() < 3) { // Sanity check if it is even a polygon
            std::cout << "WARNING: Skipping import of a zero-area polygon" << std::endl;
            this->deactivate();
            return;
        }
        if (checkSimplicity) {
            if (!tempPoly.is_simple()) {
                if (Config::get().avoidBadPolys) {
                    this->deactivate();
                    return;
                } else {
                    std::cout << "WARNING: Bad building polygon found! This might effect reconstruction quality! "
                                 "If you end up having problems, try to fix the dataset with GIS software or 'pprepair'."
                              << std::endl;
                    std::cout << "    Alternatively, you can use the 'avoid_bad_polys' flag to skip"
                                 " the import of problematic polygons.\n" << std::endl;
                }
            }
        }
        if (isOuterRing) {
            if (tempPoly.is_clockwise_oriented()) tempPoly.reverse_orientation();
            isOuterRing = false;
        } else {
            if (tempPoly.is_counterclockwise_oriented()) tempPoly.reverse_orientation();
        }
        m_poly.m_rings.push_back(tempPoly);
    }
}

PolyFeature::PolyFeature(const Polygon_with_attr& poly, const bool checkSimplicity, const int outputLayerID)
        : PolyFeature(poly, checkSimplicity) {
    m_outputLayerID = outputLayerID;
    if (m_outputLayerID >= s_numOfOutputLayers) s_numOfOutputLayers = m_outputLayerID + 1;
}

PolyFeature::PolyFeature(const Polygon_with_attr& poly, const int outputLayerID)
        : PolyFeature(poly, false, outputLayerID) {}

int PolyFeature::s_numOfPolyFeatures = 0;

void PolyFeature::calc_footprint_elevation_nni(const DT& dt) {
    typedef std::vector<std::pair<DT::Point, double>> Point_coordinate_vector;

    DT::Face_handle fh = nullptr;
    for (auto& ring: m_poly.rings()) {
        std::vector<double> ringElevations;
        for (auto& polypt: ring) {
            Point_coordinate_vector coords;
            DT::Point pt(polypt.x(), polypt.y(), 0);
            fh = dt.locate(pt, fh);
            CGAL::Triple<std::back_insert_iterator<Point_coordinate_vector>, double, bool> result =
                    CGAL::natural_neighbor_coordinates_2(dt, pt, std::back_inserter(coords), fh);

            if (!result.third) {
//                throw city4cfd_error("Trying to interpolate the point that lies outside the convex hull!");
                this->deactivate();
                return;
            }

            double elevation = 0;
            for (auto& coord : coords) {
                elevation += coord.first.z() * coord.second / result.second;
            }
            ringElevations.push_back(elevation);
        }
        m_groundElevations.push_back(ringElevations);
    }
}

#ifndef NDEBUG
void PolyFeature::calc_footprint_elevation_linear(const DT& dt) {
    DT::Face_handle fh = nullptr;
    for (auto& ring : m_poly.rings()) {
        std::vector<double> ringHeights;
        for (auto& polypt : ring) {
            DT::Point pt(polypt.x(), polypt.y(), 0);

            DT::Locate_type lt;
            int li;
            fh = dt.locate(pt, lt, li, fh);
            if (lt == DT::OUTSIDE_CONVEX_HULL) {
                this->deactivate();
                return;
            }

            std::vector<double> coords;
            CGAL::Barycentric_coordinates::triangle_coordinates_2(
                Point_2(fh->vertex(0)->point().x(), fh->vertex(0)->point().y()),
                Point_2(fh->vertex(1)->point().x(), fh->vertex(1)->point().y()),
                Point_2(fh->vertex(2)->point().x(), fh->vertex(2)->point().y()),
                Point_2(pt.x(), pt.y()),
                std::back_inserter(coords)
            );

            double h = 0;
            for (int i = 0; i < 3; ++i) {
                h += fh->vertex(i)->point().z() * coords[i];
            }
            ringHeights.push_back(h);
        }
        m_groundElevations.push_back(ringHeights);
    }
}
#endif

double PolyFeature::ground_elevation() {
    if (m_groundElevation < -global::largnum + global::smallnum) {
        if (m_groundElevations.empty()) throw city4cfd_error("Polygon elevations missing!"
                                                                " Cannot calculate average");
        // calculating base elevation as 95 percentile of outer ring
        m_groundElevation = geomutils::percentile(m_groundElevations.front(), 0.95);
    }
    return m_groundElevation;
}

///*
// * Difference between high and low elevations in the polygon
// * Calculated on the fly as rarely used
// */
//double PolyFeature::slope_height() {
//    if (m_groundElevations.empty())throw city4cfd_error("Polygon elevations missing!"
//                                                           " Cannot perform calculations");
//    return this->ground_elevation() - geomutils::percentile(m_groundElevations.front(), 0.2);
//}

bool PolyFeature::flatten_polygon_inner_points(const Point_set_3& pointCloud,
                                               std::map<int, Point_3>& flattenedPts,
                                               const SearchTree& searchTree,
                                               const std::unordered_map<Point_3, int>& pointCloudConnectivity,
                                               std::vector<EPECK::Segment_3>& constrainedEdges,
                                               std::vector<std::pair<Polygon_with_holes_2, int>>& newPolys,
                                               bool& isNextToBuilding) {

    typedef CGAL::Straight_skeleton_2<EPICK>           Ss;
    typedef std::shared_ptr<CGAL::Polygon_with_holes_2<EPICK>> PolygonPtrWH;
    typedef std::vector<PolygonPtrWH> PolygonPtrVectorWH;

#ifdef CITY4CFD_POLYFEATURE_VERBOSE
    std::cout << "\nINFO: Flattening a polygon..." << std::endl;
    int nonSimpleRings = 0;
    for (auto& ring : m_poly.rings()) {
        if (!ring.is_simple()) ++nonSimpleRings;
    }
    std::cout << "Non-simple rings: " << nonSimpleRings << std::endl;
#endif

    std::vector<int>    indices;
    std::vector<double> originalHeights;
    auto buildingPt = pointCloud.property_map<std::shared_ptr<Building>>("building_point").value();
    //-- Take tree subset bounded by the polygon
    std::vector<Point_3> subsetPts;
    Polygon_2 bbox = geomutils::calc_bbox_poly(m_poly.rings().front());
    Point_2 bbox1(bbox[0].x(), bbox[0].y());
    Point_2 bbox2(bbox[2].x(), bbox[2].y());
    Fuzzy_iso_box pts_range(bbox1, bbox2);
    searchTree.search(std::back_inserter(subsetPts), pts_range);

    //-- Check if the polygon is overlapping with a building
    std::map<int, std::shared_ptr<Building>> overlappingBuildings; //id-building map
    for (auto& pt3 : subsetPts) {
        Point_2 pt(pt3.x(), pt3.y());
        if (geomutils::point_in_poly_and_boundary(pt, m_poly)) {
            auto itIdx = pointCloudConnectivity.find(pt3);
            auto pointSetIt = pointCloud.begin();
            std::advance(pointSetIt, itIdx->second);

            auto currBuilding = buildingPt[*pointSetIt];
            if (currBuilding != nullptr) {
                int buildingId = currBuilding->get_internal_id();
                auto it = overlappingBuildings.find(buildingId);
                if (it == overlappingBuildings.end()) {
                    overlappingBuildings[buildingId] = currBuilding;
                }
            }
        }
    }
    if (!overlappingBuildings.empty()) isNextToBuilding = true;
    //-- If next intersecting a building, clip poly with building
    std::vector<Polygon_with_holes_2> flattenCandidatePolys;
    if (!overlappingBuildings.empty()) {
#ifdef CITY4CFD_POLYFEATURE_VERBOSE
        std::cout << "Overlapping buildings: " << overlappingBuildings.size() << std::endl;
#endif
        // Add clipping polys to set
        CGAL::Polygon_set_2<EPECK> polySet;
        Converter<EPICK, EPECK> to_exact;
        for (auto& intersectBuilding : overlappingBuildings) {
            polySet.join(intersectBuilding.second->get_poly().get_exact_outer_boundary());
        }
#ifdef CITY4CFD_POLYFEATURE_VERBOSE
        std::cout << "Polyset size? " << polySet.number_of_polygons_with_holes() << std::endl;
#endif
        CGAL_assertion(polySet.is_valid());
        // Clip with this polygon
        polySet.complement();
        polySet.intersection(m_poly.get_exact());
        // Store this as the new polygon
        std::vector<CGAL::Polygon_with_holes_2<EPECK>> resPolys;
        polySet.polygons_with_holes(std::back_inserter(resPolys));
        for (auto& flattenBndPoly : resPolys) {
            flattenCandidatePolys.emplace_back(geomutils::exact_poly_to_poly(flattenBndPoly));
        }
    }
#ifdef CITY4CFD_POLYFEATURE_VERBOSE
    std::cout << "After flatten polygon clipping, total new polygons: " << flattenCandidatePolys.size() << std::endl;
#endif
    if (!flattenCandidatePolys.empty()) {
        bool isFirst = true;
        for (auto& poly : flattenCandidatePolys) {
            const double offsetVal = 0.1; // offset value hardcoded
            PolygonPtrVectorWH offsetPoly =
                    CGAL::create_interior_skeleton_and_offset_polygons_with_holes_2(offsetVal,
                                                                                    poly.get_cgal_type());
#ifdef CITY4CFD_POLYFEATURE_VERBOSE
            std::cout << "Offset polygons created: " << offsetPoly.size() << std::endl;
#endif
            if (!offsetPoly.empty()) { // make a check whether the offset is successfully created
                poly = *(offsetPoly.front());
            } else {
                std::cout << "WARNING: Polygon ID" << m_id << " Skeleton construction failed! Some polygons will not be flattened." << std::endl;
                return false;
            }
            if (isFirst) {
                m_poly = Polygon_with_holes_2(*(offsetPoly.front()));
                isFirst = false;
            } else {
                // Some polys are cut -- save the new resulting polys and add them as new features in Map3D
                newPolys.emplace_back(Polygon_with_holes_2(*(offsetPoly.front())), m_outputLayerID);
            }
        }
    } else {
        if (!overlappingBuildings.empty()) {
            std::cout << "WARNING: Polygon ID" << m_id << " flattening failed due to unresolved"
                                                         " intersections with buildings!" << std::endl;
            return false;
        }
        flattenCandidatePolys.push_back(m_poly);
    }
    //-- Collect points that have not been already flattened
    for (auto& pt3 : subsetPts) {
        for (auto& flattenBndPoly : flattenCandidatePolys) {
            Point_2 pt(pt3.x(), pt3.y());
            if (geomutils::point_in_poly_and_boundary(pt, flattenBndPoly)) {
                auto itIdx = pointCloudConnectivity.find(pt3);

                auto it = flattenedPts.find(itIdx->second);
                if (it == flattenedPts.end()) {
                    indices.push_back(itIdx->second);
                    originalHeights.push_back(pointCloud.point(itIdx->second).z());
                }
            }
        }
    }
    //-- Flatten surface points
    if (indices.empty()) {
        return false;
    }
    double avgHeight = geomutils::percentile(originalHeights,
                                             Config::get().flattenSurfaces[this->get_output_layer_id()] / 100);

    //-- Add new points to the temp map
    for (auto& i : indices) {
        flattenedPts[i] = Point_3(pointCloud.point(i).x(), pointCloud.point(i).y(), avgHeight);
    }
    //-- Add additional new segments to constrain
    if (!flattenCandidatePolys.empty()) {
        for (auto& flattenBndPoly : flattenCandidatePolys) {
            for (auto& ring: flattenBndPoly.rings()) {
                for (const auto& newEdge: ring.edges()) {
                    ePoint_3 pt1(newEdge.point(0).x(), newEdge.point(0).y(), avgHeight);
                    ePoint_3 pt2(newEdge.point(1).x(), newEdge.point(1).y(), avgHeight);
                    constrainedEdges.emplace_back(pt1, pt2);
                }
            }
        }
    }
    return true;
}

void PolyFeature::set_zero_borders() {
    for (auto& ring : m_groundElevations) {
        for (auto& pt : ring) {
            pt = 0.;
        }
    }
}

void PolyFeature::calc_min_bbox() {
    if (m_poly.rings().front().is_empty()) throw city4cfd_error("Missing polygon!");
    //-- Point set needs to be convex for the rotating caliper algorithm
    std::vector<Point_2> chull;
    CGAL::convex_hull_2(m_poly.rings().front().vertices_begin(),
                        m_poly.rings().front().vertices_end(),
                        std::back_inserter(chull));

    std::vector<Point_2> obbPts; obbPts.reserve(4);
    CGAL::min_rectangle_2(chull.begin(),
                          chull.end(),
                          std::back_inserter(obbPts));
    assert(obbPts.size() == 4);

    m_minBbox.vec1 = obbPts[1] - obbPts[0];
    m_minBbox.vec2 = obbPts[3] - obbPts[0];

    m_minBbox.calc();
}

void PolyFeature::clear_feature() {
    m_groundElevations.clear();
    m_mesh.clear();
}

int PolyFeature::new_internal_id() {
    return ++s_numOfPolyFeatures;
}

Polygon_with_holes_2& PolyFeature::get_poly() {
    return m_poly;
}

const Polygon_with_holes_2& PolyFeature::get_poly() const {
    return m_poly;
}

Polygon_with_attr PolyFeature::get_poly_w_attr() const {
    Polygon_with_attr poly;
    poly.polygon = m_poly;
    return poly;
}

const std::vector<std::vector<double>>& PolyFeature::get_ground_elevations() const {
    return m_groundElevations;
}

const int PolyFeature::get_internal_id() const {
    return m_polyInternalID;
}

MinBbox& PolyFeature::get_min_bbox() {
    if (m_minBbox.empty()) {
        this->calc_min_bbox();
    }
    return m_minBbox;
}

void PolyFeature::parse_json_poly(const nlohmann::json& poly, const bool checkSimplicity) {
    for (auto& polyEdges : poly["geometry"]["coordinates"]) {
        Polygon_2 tempPoly;

        const double minDist = 0.0001;
        Point_2 prev(global::largnum, global::largnum);
        for (auto& coords: polyEdges) {
            Point_2 pt2((double)coords[0] - Config::get().pointOfInterest.x(),
                        (double)coords[1] - Config::get().pointOfInterest.y());
            if (CGAL::squared_distance(pt2, prev) > minDist) {
                tempPoly.push_back(pt2);
                prev = pt2;
            }
        }
        geomutils::pop_back_if_equal_to_front(tempPoly);
        if (tempPoly.size() < 3) { // Sanity check if it is even a polygon
            std::cout << "WARNING: Skipping import of a zero-area polygon" << std::endl;
            this->deactivate();
            return;
        }
        if (checkSimplicity) {
            if (!tempPoly.is_simple()) {
                if (Config::get().avoidBadPolys) {
                    this->deactivate();
                    return;
                } else {
                    std::cout << "WARNING: Bad building polygon found! This might effect reconstruction quality! "
                                 "If you end up having problems, try to fix the dataset with GIS software or 'pprepair'."
                              << std::endl;
                    std::cout << "    Alternatively, you can use the 'avoid_bad_polys' flag to skip"
                                 " the import of problematic polygons.\n" << std::endl;
                }
            }
        }
        if (m_poly.m_rings.empty()) {
            if (tempPoly.is_clockwise_oriented()) tempPoly.reverse_orientation();
        } else {
            if (tempPoly.is_counterclockwise_oriented()) tempPoly.reverse_orientation();
        }
        m_poly.m_rings.push_back(tempPoly);
    }
}