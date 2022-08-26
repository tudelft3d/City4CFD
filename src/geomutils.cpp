/*
  City4CFD
 
  Copyright (c) 2021-2022, 3D Geoinformation Research Group, TU Delft  

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

#include "geomutils.h"

#include "PolyFeature.h"

#include <CGAL/Polygon_mesh_processing/repair.h>

double geomutils::avg(const std::vector<double>& values) {
    if (values.empty()) throw std::length_error("Can't calculate average of a zero-sized vector!");
    double average = 0;
    for (auto& value : values) {
        average += value;
    }
    return (average / (double)values.size());
}

double geomutils::percentile(std::vector<double> values, const double percentile) {
    assert(percentile >= 0 && percentile <= 1);
    if (values.empty()) throw std::length_error("Can't calculate percentile of a zero-sized vector!");
    std::sort(values.begin(), values.end());
    int i = values.size() * percentile;
    return values[i];
}

bool geomutils::point_in_circle(const Point_3& pt, const Point_2& center, const double& radius) {
    if (pow(pt.x() - center.x(), 2)
      + pow(pt.y() - center.y(), 2)
      < pow(radius, 2)) {
        return true;
    }
    return false;
}

void geomutils::cdt_to_mesh(CDT& cdt, Mesh& mesh, const int surfaceLayerID) {
    std::map<CDT::Vertex_handle, int> indices;
    std::vector<Mesh::vertex_index> mesh_vertex;
    std::vector<Mesh::face_index> face_index;
    mesh_vertex.reserve(cdt.number_of_vertices());

    int counter = 0;
    for (const auto& it : cdt.finite_vertex_handles()) {
        mesh_vertex.emplace_back(mesh.add_vertex(Converter<EPECK, EPICK>()(it->point())));
        //        outstream << it->point() << std::endl;
        indices.insert(std::pair<CDT::Vertex_handle, int>(it, counter++));
    }

    for (const auto& it : cdt.finite_face_handles()) {
        if (it->info().surfaceLayer != surfaceLayerID) continue;

        int v1 = indices[it->vertex(0)];
        int v2 = indices[it->vertex(1)];
        int v3 = indices[it->vertex(2)];
        mesh.add_face(mesh_vertex[v1], mesh_vertex[v2], mesh_vertex[v3]);
    }
}
//-- CGAL's constrained domain marker expanded to mark different polygon types
void geomutils::mark_domains(CDT& ct,
                             const Face_handle& start,
                             int index,
                             std::list<CDT::Edge>& border,
                             PolyFeatures& features)
{
    if (start->info().nesting_level != -1) {
        return;
    }

    //-- Check which polygon contains the constrained (i.e. non-terrain) point
    Point_3 chkPoint;
    Converter<EPECK, EPICK> to_inexact;
    if (!features.empty()) {
        chkPoint = CGAL::centroid(to_inexact(start->vertex(0)->point()),
                                  to_inexact(start->vertex(1)->point()),
                                  to_inexact(start->vertex(2)->point()));
    }
    int surfaceLayer = -1; //-- Default value is unmarked triangle, i.e. general terrain
    if (index != 0) {
        for (auto& feature : features) {
            if (!feature->is_active()) continue;
            //- Polygons are already ordered according to importance - find first polygon
            if (geomutils::point_in_poly(chkPoint, feature->get_poly())) {
                if (feature->get_class() == BUILDING) {
//                    surfaceLayer = 9999; //- Remove building footprints from terrain
                    surfaceLayer = -1; //- Leave building footprints as part of terrain
                    break;
                } else {
                    surfaceLayer = feature->get_output_layer_id();
                    break;
                }
            }
        }
    }

    std::list<Face_handle> queue;
    queue.push_back(start);
    while (! queue.empty()) {
        Face_handle fh = queue.front();
        queue.pop_front();
        if (fh->info().nesting_level == -1) {
            fh->info().nesting_level = index;
            if (surfaceLayer != -1) {
//                fh->info().surfaceLayer = surfaceLayer;
                check_layer(fh, surfaceLayer);
            }
            for (int i = 0; i < 3; i++) {
                CDT::Edge e(fh,i);
                Face_handle n = fh->neighbor(i);
                if (n->info().nesting_level == -1) {
                    if (ct.is_constrained(e)) border.push_back(e);
                    else queue.push_back(n);
                }
            }
        }
    }
}

void geomutils::mark_domains(CDT& cdt, PolyFeatures features) {
    for (CDT::Face_handle f : cdt.all_face_handles()) {
        f->info().nesting_level = -1;
    }
    std::list<CDT::Edge> border;
    mark_domains(cdt, cdt.infinite_face(), 0, border, features);
    while (!border.empty()) {
        CDT::Edge e = border.front();
        border.pop_front();
        Face_handle n = e.first->neighbor(e.second);
        if (n->info().nesting_level == -1) {
            mark_domains(cdt, n, e.first->info().nesting_level + 1, border, features);
        }
    }
    for (CDT::Face_handle f : cdt.all_face_handles()) {
        if (f->info().surfaceLayer == -2) {
            f->info().surfaceLayer = -9999;
        }
    }
}

void geomutils::check_layer(const Face_handle& fh, int surfaceLayer) {
    if (fh->info().surfaceLayer == 9999) return;
    auto it = Config::get().flattenSurfaces.find(surfaceLayer);
    if (it != Config::get().flattenSurfaces.end()) {
        Converter<EPECK, EPICK> to_inexact;
        Vector_3 vertical(0, 0, 1);
        Vector_3 norm = CGAL::normal(to_inexact(fh->vertex(0)->point()),
                                     to_inexact(fh->vertex(1)->point()),
                                     to_inexact(fh->vertex(2)->point()));

        if (CGAL::approximate_angle(norm, vertical) < 60.0) {
            fh->info().surfaceLayer = surfaceLayer;
        } else {
            fh->info().surfaceLayer = -2;
        }
    } else {
        fh->info().surfaceLayer = surfaceLayer;
    }
}

void geomutils::shorten_long_poly_edges(Polygon_2& poly, double maxLen) {
    auto& polyVec = poly.container();
    int i = 0;
    while (i != polyVec.size()) {
        auto edge = polyVec[(i + 1) % polyVec.size()] - polyVec[i];
        double edgeSqLen = edge.squared_length();
        if (edgeSqLen > maxLen*maxLen) {
            int  numSeg  = ceil(sqrt(edgeSqLen) / maxLen);
            auto segVec = edge / numSeg;
            ++i;
            for (auto j = 0; j < numSeg - 1; ++j) {
                Point_2 pt = polyVec[i - 1] + segVec;
                polyVec.insert(polyVec.begin() + i, pt);
                ++i;
            }
        } else ++i;
    }
}

Point_2 geomutils::rotate_pt(const Point_2& pt, const double angle, Point_2 centerPt) {
    return {cos(angle) * (pt.x() - centerPt.x()) - sin(angle) * (pt.y() - centerPt.y()) + centerPt.x(),
                   sin(angle) * (pt.x() - centerPt.x()) + cos(angle) * (pt.y() - centerPt.y()) + centerPt.y()};
}

Point_3 geomutils::rotate_pt_xy(const Point_3& pt, const double angle, Point_2 centerPt) {
    return {cos(angle) * (pt.x() - centerPt.x()) - sin(angle) * (pt.y() - centerPt.y()) + centerPt.x(),
            sin(angle) * (pt.x() - centerPt.x()) + cos(angle) * (pt.y() - centerPt.y()) + centerPt.y(),
            pt.z()};
}

void geomutils::interpolate_poly_from_pc(const Polygon_2& poly, std::vector<double>& heights,
                                         const Point_set_3& pointCloud) {
    SearchTree searchTree(pointCloud.points().begin(), pointCloud.points().end());
    //-- Calculate elevation of polygon outer boundary
    //-- Point elevation is the average of 5 nearest neighbors from the PC
    std::vector<double> ringHeights;
    for (auto& polypt : poly) {
        Point_3 query(polypt.x(), polypt.y(), 0);
        Neighbor_search search(searchTree, query, 5);
//        Fuzzy_sphere search_radius(query, 5);
//        std::list<Point_3> result;
//        searchTree.search(std::back_inserter(result), search_radius);

        std::vector<double> poly_height;
        for (Neighbor_search::iterator it = search.begin(); it != search.end(); ++it) {
            poly_height.push_back(it->first.z());
        }
//        for (auto& pt : result) {
//            poly_height.push_back(pt.z());
//        }
        heights.emplace_back(geomutils::avg(poly_height));
    }
}

void geomutils::remove_self_intersections(Mesh& mesh) {
#if CGAL_VERSION_NR >= 1050101000 // 5.1.0
#ifndef __clang__
    PMP::experimental::remove_self_intersections(mesh);
#else
    throw std::runtime_error(std::string("Function remove_self_intersections() does not work with clang compiler!"
                                         " Set 'handle_self_intersections' to false"
                                         " or recompile the program again using gcc"));
#endif
#else
    PMP::remove_self_intersections(mesh);
#endif
}

bool geomutils::polygons_in_contact(const Polygon_with_holes_2& firstPoly, const Polygon_with_holes_2& secondPoly) {
    if (CGAL::do_overlap(firstPoly.bbox(), secondPoly.bbox())) {
        if (CGAL::do_intersect(firstPoly.get_exact(),
                               secondPoly.get_exact())) {
            return true;
        }
    }
    return false;
}

//-- Templated functions
//-- Check if the point is inside a polygon on a 2D projection
template <typename T>
bool geomutils::point_in_poly(const T& pt2, const Polygon_2& polygon) {
    Point_2 pt(pt2.x(), pt2.y());
    if (CGAL::bounded_side_2(polygon.begin(), polygon.end(), pt) == CGAL::ON_BOUNDED_SIDE) {
        return true;
    }
    return false;
}
//- Explicit template instantiation
template bool geomutils::point_in_poly<Point_2>(const Point_2& pt2, const Polygon_2& polygon);
template bool geomutils::point_in_poly<Point_3>(const Point_3& pt2, const Polygon_2& polygon);

template <typename T>
bool geomutils::point_in_poly(const T& pt2, const Polygon_with_holes_2& polygon) {
    Point_2 pt(pt2.x(), pt2.y());

    //-- Check if the point falls within the outer surface
    if (point_in_poly(pt, polygon.outer_boundary())) {
        // Check if the point falls within one of the holes
        for (auto it_hole = polygon.holes_begin(); it_hole != polygon.holes_end(); ++it_hole) {
            if (point_in_poly(pt, *it_hole)) {
                return false;
            }
        }
        return true;
    }
    return false;
}
//- Explicit template instantiation
template bool geomutils::point_in_poly<Point_2>(const Point_2& pt2, const Polygon_with_holes_2& polygon);
template bool geomutils::point_in_poly<Point_3>(const Point_3& pt2, const Polygon_with_holes_2& polygon);

template <typename T>
void geomutils::make_round_poly(const Point_2& centre, double radius, T& poly) {
    const int nPts = 360; // Hardcoded
    const double angInt = 2 * M_PI / (double) nPts;
    double ang = 0;
    for (auto i = 0; i < nPts; ++i) {
        double xPt = centre.x() + radius * cos(ang + angInt);
        double yPt = centre.y() + radius * sin(ang + angInt);
        ang = ang + angInt;
        poly.push_back(Point_2(xPt, yPt));
    }
}
//- Explicit template instantiation
template void geomutils::make_round_poly<Polygon_2>(const Point_2& centre, double radius, Polygon_2& poly);
template void geomutils::make_round_poly<Polygon_with_holes_2>(const Point_2& centre, double radius, Polygon_with_holes_2& poly);

template <typename T>
void geomutils::make_round_poly(const Point_2& centre, double radius1, double radius2,
                                int nPts, double angInt, double ang, T& poly) {
    for (auto i = 0; i < nPts; ++i) {
        double xPt = centre.x() + radius1 * cos(ang + angInt);
        double yPt = centre.y() + radius2 * sin(ang + angInt);
        ang = ang + angInt;
        poly.push_back(Point_2(xPt, yPt));
    }
}
//- Explicit template instantiation
template void geomutils::make_round_poly<Polygon_2>(const Point_2& centre, double radius1, double radius2, int nPts, double angInt, double ang, Polygon_2& poly);
template void geomutils::make_round_poly<Polygon_with_holes_2>(const Point_2& centre, double radius1, double radius2, int nPts, double angInt, double ang, Polygon_with_holes_2& poly);

template <typename T, typename U>
void geomutils::smooth_dt(const Point_set_3& pointCloud, T& dt) {
    // Smooth heights with 5 successive Gaussian filters
#ifdef CGAL_LINKED_WITH_TBB
    using Concurrency_tag = CGAL::Parallel_tag;
#else
    using Concurrency_tag = CGAL::Sequential_tag;
#endif
    double spacing = CGAL::compute_average_spacing<Concurrency_tag>(pointCloud, 6);
    spacing *= 2;
    double gaussian_variance = 4 * spacing * spacing;
    for (typename T::Vertex_handle vh : dt.finite_vertex_handles())
    {
        double z = CGAL::to_double(vh->point().z());
        double total_weight = 1;
        typename T::Vertex_circulator circ = dt.incident_vertices (vh),
                start = circ;
        do
        {
            if (!dt.is_infinite(circ))
            {
                double sq_dist = CGAL::to_double(CGAL::squared_distance(vh->point(), circ->point()));
                double weight = std::exp(-sq_dist / gaussian_variance);
                z += weight * CGAL::to_double(circ->point().z());
                total_weight += weight;
            }
        }
        while (++ circ != start);
        z /= total_weight;
        vh->point() = CGAL::Point_3<U>(vh->point().x(), vh->point().y(), z);
    }
}
//-- Explicit template instantiation
template void geomutils::smooth_dt<DT, EPICK>(const Point_set_3& pointCloud, DT& dt);
template void geomutils::smooth_dt<CDT, EPECK>(const Point_set_3& pointCloud, CDT& dt);

template <typename T>
Polygon_2 geomutils::calc_bbox_poly(const T& inputPts) {
    double bxmin, bymin, bxmax, bymax;
    bxmin = global::largnum; bymin = global::largnum; bxmax = -global::largnum; bymax = -global::largnum;
    for (auto& pt : inputPts) {
        if (pt.x() > bxmax) bxmax = pt.x();
        if (pt.x() < bxmin) bxmin = pt.x();
        if (pt.y() > bymax) bymax = pt.y();
        if (pt.y() < bymin) bymin = pt.y();
    }
    Polygon_2 bboxPoly;
    bboxPoly.push_back(Point_2(bxmin, bymin));
    bboxPoly.push_back(Point_2(bxmax, bymin));
    bboxPoly.push_back(Point_2(bxmax, bymax));
    bboxPoly.push_back(Point_2(bxmin, bymax));

    return bboxPoly;
}
//- Explicit template instantiation
template Polygon_2 geomutils::calc_bbox_poly<std::vector<Point_2>>(const std::vector<Point_2>& inputPts);
template Polygon_2 geomutils::calc_bbox_poly<Polygon_2>(const Polygon_2& inputPts);

template <typename T>
void geomutils::pop_back_if_equal_to_front(CGAL::Polygon_2<T>& poly)
{
    auto it = poly.end();
    --it;
    if((*poly.begin()) == *it)
        poly.erase(it);
}
//- Explicit template instantiation
template void geomutils::pop_back_if_equal_to_front<EPICK>(CGAL::Polygon_2<EPICK>& poly);
template void geomutils::pop_back_if_equal_to_front<EPECK>(CGAL::Polygon_2<EPECK>& poly);
