#include "PolyFeature.h"

#include "geomutils.h"

PolyFeature::PolyFeature()
    : TopoFeature(), _poly(), _base_heights(), _polyInternalID() {}

PolyFeature::PolyFeature(const int outputLayerID)
    : TopoFeature(outputLayerID), _poly(), _base_heights(), _polyInternalID() {}

PolyFeature::PolyFeature(const nlohmann::json& poly)
    : TopoFeature(), _base_heights(), _polyInternalID() {
    this->parse_json_poly(poly);
}

PolyFeature::PolyFeature(const int outputLayerID, const int internalID)
    : TopoFeature(outputLayerID), _polyInternalID(internalID) {}

PolyFeature::PolyFeature(const nlohmann::json& poly, const int outputLayerID)
    : PolyFeature(poly) {
    _outputLayerID = outputLayerID;
    if (_outputLayerID  >= _numOfOutputLayers) _numOfOutputLayers = _outputLayerID + 1;
}

PolyFeature::PolyFeature(const nlohmann::json& poly, const int outputLayerID, const int internalID)
    : PolyFeature(poly) {
    _polyInternalID = internalID;
    _outputLayerID    = outputLayerID;
    if (_outputLayerID  >= _numOfOutputLayers) _numOfOutputLayers = _outputLayerID + 1;
}

PolyFeature::~PolyFeature() = default;

void PolyFeature::calc_footprint_elevation_nni(const DT& dt) {
    typedef std::vector<std::pair<DT::Point, double>> Point_coordinate_vector;
    DT::Face_handle fh = nullptr;
    for (auto& ring: _poly.rings()) {
        std::vector<double> ringHeights;
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

            double height = 0;
            for (auto& coord : coords) {
                height += coord.first.z() * coord.second / result.second;
            }
            ringHeights.push_back(height);
        }
        _base_heights.push_back(ringHeights);
    }
}

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
        _base_heights.push_back(ringHeights);
    }
}

void PolyFeature::average_polygon_inner_points(const Point_set_3& pointCloud,
                                               std::map<int, Point_3>& averagedPts,
                                               const SearchTree& searchTree,
                                               const std::unordered_map<Point_3, int>& pointCloudConnectivity) const {
    std::vector<int>    indices;
    std::vector<double> originalHeights;
    //-- Take tree subset bounded by the polygon
    std::vector<Point_3> subsetPts;
    Polygon_2 bbox = geomutils::calc_bbox_poly(_poly.rings().front());
    Point_3 bbox1(bbox[0].x(), bbox[0].y(), -g_largnum);
    Point_3 bbox2(bbox[2].x(), bbox[2].y(), g_largnum);
    Fuzzy_iso_box pts_range(bbox1, bbox2);
    searchTree.search(std::back_inserter(subsetPts), pts_range);

    //-- Collect points that have not been already averaged
    for (auto& pt3 : subsetPts) {
        Point_2 pt(pt3.x(), pt3.y());
        if (CGAL::bounded_side_2(_poly._rings.front().begin(), _poly._rings.front().end(), pt) != CGAL::ON_UNBOUNDED_SIDE) {
            auto itIdx = pointCloudConnectivity.find(pt3);
            auto it = averagedPts.find(itIdx->second);
            if (it == averagedPts.end()) {
                indices.push_back(itIdx->second);
                originalHeights.push_back(pointCloud.point(itIdx->second).z());
            }
        }
    }
    //-- Average points
    if (indices.empty()) {
        return;
    }
    double avgHeight = geomutils::percentile(originalHeights, config::averageSurfaces[this->get_output_layer_id()] / 100);

    //-- Add new points to the temp map
    for (auto& i : indices) {
        averagedPts[i] = Point_3(pointCloud.point(i).x(), pointCloud.point(i).y(), avgHeight);
    }
}

void PolyFeature::clear_feature() {
    _base_heights.clear();
    _mesh.clear();
}

Polygon_with_holes_2& PolyFeature::get_poly() {
    return _poly;
}

const Polygon_with_holes_2& PolyFeature::get_poly() const {
    return _poly;
}

const std::vector<std::vector<double>>& PolyFeature::get_base_heights() const {
    return _base_heights;
}

const int PolyFeature::get_internal_id() const {
    return _polyInternalID;
}

void PolyFeature::parse_json_poly(const nlohmann::json& poly) {
    for (auto& polyEdges : poly) {
        Polygon_2 tempPoly;
        for (auto& coords : polyEdges) {
            tempPoly.push_back(Point_2(coords[0], coords[1]));
        }
        CGAL::internal::pop_back_if_equal_to_front(tempPoly);

        if (_poly._rings.empty()) {
            if (tempPoly.is_clockwise_oriented()) tempPoly.reverse_orientation();
        } else {
            if (tempPoly.is_counterclockwise_oriented()) tempPoly.reverse_orientation();
        }
        _poly._rings.push_back(tempPoly);
    }
}