/*
  City4CFD
 
  Copyright (c) 2021-2023, 3D Geoinformation Research Group, TU Delft

  This file is part of City4CFD.

  City4CFD is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  City4CFD is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with City4CFD.  If not, see <http://www.gnu.org/licenses/>.

  For any information or further details about the use of City4CFD, contact
  Ivan Pađen
  <i.paden@tudelft.nl>
  3D Geoinformation Research Group
  Delft University of Technology
*/

#include "PolyFeature.h"

#include "geomutils.h"

#include <CGAL/convex_hull_2.h>
#include <CGAL/min_quadrilateral_2.h>
#include <CGAL/natural_neighbor_coordinates_2.h>
#ifndef NDEBUG
#include <CGAL/Barycentric_coordinates_2/Triangle_coordinates_2.h>
#endif

PolyFeature::PolyFeature()
    : TopoFeature(), _poly(), _groundElevations(), _polyInternalID(),
      _groundElevation(-global::largnum), _minBbox() {}

PolyFeature::PolyFeature(const int outputLayerID)
    : TopoFeature(outputLayerID), _poly(), _groundElevations(), _polyInternalID(),
      _groundElevation(-global::largnum), _minBbox() {}

PolyFeature::PolyFeature(const nlohmann::json& poly, const bool checkSimplicity)
    : TopoFeature(), _groundElevations(), _polyInternalID(),
      _groundElevation(-global::largnum), _minBbox() {
    this->parse_json_poly(poly, checkSimplicity);
}

PolyFeature::PolyFeature(const int outputLayerID, const int internalID)
    : TopoFeature(outputLayerID), _groundElevations(), _polyInternalID(internalID),
      _groundElevation(-global::largnum), _minBbox() {}

PolyFeature::PolyFeature(const nlohmann::json& poly, const bool checkSimplicity, const int outputLayerID)
    : PolyFeature(poly, checkSimplicity) {
    _outputLayerID = outputLayerID;
    if (_outputLayerID  >= _numOfOutputLayers) _numOfOutputLayers = _outputLayerID + 1;
}

PolyFeature::PolyFeature(const nlohmann::json& poly, const int outputLayerID)
        : PolyFeature(poly, false, outputLayerID) {}

PolyFeature::PolyFeature(const nlohmann::json& poly, const bool checkSimplicity,
                         const int outputLayerID, const int internalID)
    : PolyFeature(poly, checkSimplicity) {
    _polyInternalID = internalID;
    _outputLayerID  = outputLayerID;
    if (_outputLayerID  >= _numOfOutputLayers) _numOfOutputLayers = _outputLayerID + 1;
}

PolyFeature::PolyFeature(const nlohmann::json& poly, const int outputLayerID, const int internalID)
        : PolyFeature(poly, false, outputLayerID, internalID) {}

PolyFeature::~PolyFeature() = default;

void PolyFeature::calc_footprint_elevation_nni(const DT& dt) {
    typedef std::vector<std::pair<DT::Point, double>> Point_coordinate_vector;
    DT::Face_handle fh = nullptr;
    for (auto& ring: _poly.rings()) {
        std::vector<double> ringElevations;
        for (auto& polypt: ring) {
            Point_coordinate_vector coords;
            DT::Point pt(polypt.x(), polypt.y(), 0);
            fh = dt.locate(pt, fh);
            CGAL::Triple<std::back_insert_iterator<Point_coordinate_vector>, double, bool> result =
                    CGAL::natural_neighbor_coordinates_2(dt, pt, std::back_inserter(coords), fh);

            if (!result.third) {
//                throw std::runtime_error("Trying to interpolate the point that lies outside the convex hull!");
                this->deactivate();
                return;
            }

            double elevation = 0;
            for (auto& coord : coords) {
                elevation += coord.first.z() * coord.second / result.second;
            }
            ringElevations.push_back(elevation);
        }
        _groundElevations.push_back(ringElevations);
    }
}

#ifndef NDEBUG
void PolyFeature::calc_footprint_elevation_linear(const DT& dt) {
    typedef CGAL::Barycentric_coordinates::Triangle_coordinates_2<iProjection_traits>   Triangle_coordinates;
    DT::Face_handle fh = nullptr;
    for (auto& ring : _poly.rings()) {
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

            Triangle_coordinates triangle_coordinates(fh->vertex(0)->point(),
                                                      fh->vertex(1)->point(),
                                                      fh->vertex(2)->point());
            std::vector<double> coords;
            triangle_coordinates(pt, std::back_inserter(coords));

            double h = 0;
            for (int i = 0; i < 3; ++i) {
                h += fh->vertex(i)->point().z() * coords[i];
            }
            ringHeights.push_back(h);
        }
        _groundElevations.push_back(ringHeights);
    }
}
#endif

double PolyFeature::ground_elevation() {
    if (_groundElevation < -global::largnum + global::smallnum) {
        if (_groundElevations.empty())throw std::runtime_error("Polygon elevations missing!"
                                                           " Cannot calculate average");
        // calculating base elevation as 90 percentile of outer ring
        _groundElevation = geomutils::percentile(_groundElevations.front(), 0.2);
    }
    return _groundElevation;
}

double PolyFeature::slope_height() {
    if (_groundElevations.empty())throw std::runtime_error("Polygon elevations missing!"
                                                           " Cannot perform calculations");
    // calculating slope height as difference between high and low elevations
    // calculated on the fly as rarely used
    return geomutils::percentile(_groundElevations.front(), 0.9) - this->ground_elevation();
}

bool PolyFeature::flatten_polygon_inner_points(const Point_set_3& pointCloud,
                                               std::map<int, Point_3>& flattenedPts,
                                               const SearchTree& searchTree,
                                               const std::unordered_map<Point_3, int>& pointCloudConnectivity) const {
    std::vector<int>    indices;
    std::vector<double> originalHeights;
    auto is_building_pt = pointCloud.property_map<bool>("is_building_point").first;
    //-- Take tree subset bounded by the polygon
    std::vector<Point_3> subsetPts;
    Polygon_2 bbox = geomutils::calc_bbox_poly(_poly.rings().front());
    Point_2 bbox1(bbox[0].x(), bbox[0].y());
    Point_2 bbox2(bbox[2].x(), bbox[2].y());
    Fuzzy_iso_box pts_range(bbox1, bbox2);
    searchTree.search(std::back_inserter(subsetPts), pts_range);

    //-- Collect points that have not been already flattened
    for (auto& pt3 : subsetPts) {
        Point_2 pt(pt3.x(), pt3.y());
        if (CGAL::bounded_side_2(_poly._rings.front().begin(),
                                 _poly._rings.front().end(),
                                 pt) != CGAL::ON_UNBOUNDED_SIDE) {
            auto itIdx = pointCloudConnectivity.find(pt3);

            auto pointSetIt = pointCloud.begin();
            std::advance(pointSetIt, itIdx->second);
            if (is_building_pt[*pointSetIt])
                return false; //todo temp solution when having adjacent buildings

            auto it = flattenedPts.find(itIdx->second);
            if (it == flattenedPts.end()) {
                indices.push_back(itIdx->second);
                originalHeights.push_back(pointCloud.point(itIdx->second).z());
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
    return true;
}

void PolyFeature::set_zero_borders() {
    for (auto& ring : _groundElevations) {
        for (auto& pt : ring) {
            pt = 0.;
        }
    }
}

void PolyFeature::calc_min_bbox() {
    if (_poly.rings().front().is_empty()) throw std::runtime_error("Missing polygon!");
    //-- Point set needs to be convex for the rotating caliper algorithm
    std::vector<Point_2> chull;
    CGAL::convex_hull_2(_poly.rings().front().vertices_begin(),
                        _poly.rings().front().vertices_end(),
                        std::back_inserter(chull));

    std::vector<Point_2> obbPts; obbPts.reserve(4);
    CGAL::min_rectangle_2(chull.begin(),
                          chull.end(),
                          std::back_inserter(obbPts));
    assert(obbPts.size() == 4);

    _minBbox.vec1 = obbPts[1] - obbPts[0];
    _minBbox.vec2 = obbPts[3] - obbPts[0];

    _minBbox.calc();
}

void PolyFeature::clear_feature() {
    _groundElevations.clear();
    _mesh.clear();
}

Polygon_with_holes_2& PolyFeature::get_poly() {
    return _poly;
}

const Polygon_with_holes_2& PolyFeature::get_poly() const {
    return _poly;
}

const std::vector<std::vector<double>>& PolyFeature::get_ground_elevations() const {
    return _groundElevations;
}

const int PolyFeature::get_internal_id() const {
    return _polyInternalID;
}

MinBbox& PolyFeature::get_min_bbox() {
    if (_minBbox.empty()) {
        this->calc_min_bbox();
    }
    return _minBbox;
}

void PolyFeature::parse_json_poly(const nlohmann::json& poly, const bool checkSimplicity) {
    for (auto& polyEdges : poly["geometry"]["coordinates"]) {
        Polygon_2 tempPoly;

        const double minDist = 0.0001;
        Point_2 prev((double)polyEdges.front()[0] - Config::get().pointOfInterest.x(),
                     (double)polyEdges.front()[1] - Config::get().pointOfInterest.y());
        prev += (Point_2(global::largnum, global::largnum) - prev);
        for (auto& coords: polyEdges) {
            Point_2 pt2((double)coords[0] - Config::get().pointOfInterest.x(),
                        (double)coords[1] - Config::get().pointOfInterest.y());
            if (CGAL::squared_distance(pt2, prev) > minDist) {
                tempPoly.push_back(pt2);
                prev = pt2;
            }
        }
        geomutils::pop_back_if_equal_to_front(tempPoly);

        if (_poly._rings.empty()) {
            if (tempPoly.is_clockwise_oriented()) tempPoly.reverse_orientation();
        } else {
            if (tempPoly.is_counterclockwise_oriented()) tempPoly.reverse_orientation();
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
        _poly._rings.push_back(tempPoly);
    }
}