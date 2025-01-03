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

#ifndef CITY4CFD_GEOMUTILS_H
#define CITY4CFD_GEOMUTILS_H

#include "Config.h"

namespace geomutils {
    double  avg(const std::vector<double>& values);
    double  percentile(std::vector<double> values, const double percentile);
//    bool    point_in_circle(const Point_3& pt, const Point_2& center, const double& radius);
    void    cdt_to_mesh(CDT& cdt, Mesh& mesh, const int surfaceLayerID = -9999);
    void    dt_to_mesh(DT& dt, Mesh& mesh);
    void    mark_domains(CDT& cdt);
    void    mark_domains(CDT& ct, const Face_handle& start, int index,
                         std::list<CDT::Edge>& border);
    void    shorten_long_poly_edges(Polygon_2& poly, double maxLen = Config::get().edgeMaxLen);
    Point_2 rotate_pt(const Point_2& pt, const double angle, Point_2 centerPt = Point_2(0, 0));
    Point_3 rotate_pt_xy(const Point_3& pt, const double angle, Point_2 centerPt = Point_2(0, 0));
    void    interpolate_poly_from_pc(const Polygon_2& poly, std::vector<double>& elevations, const Point_set_3& pointCloud);
    bool    polygons_in_contact(const Polygon_with_holes_2& firstPoly, const Polygon_with_holes_2& secondPoly);
    Polygon_2 offset_polygon(const Polygon_2& poly, double offset);
    Polygon_with_holes_2 offset_polygon_with_holes(const Polygon_with_holes_2& poly, double offset);
    Polygon_with_holes_2 exact_poly_to_poly(const CGAL::Polygon_with_holes_2<EPECK>& exactPoly);
    void    remove_self_intersections(Mesh& mesh);

    //-- Templated functions
    template <typename T> bool point_in_poly(const T& pt2, const Polygon_with_holes_2& polygon);
    template <typename T> bool point_in_poly(const T& pt2, const Polygon_2& polygon);
    template <typename T> bool point_in_poly_and_boundary(const T& pt2, const Polygon_with_holes_2& polygon);
    template <typename T> bool point_in_poly_and_boundary(const T& pt2, const Polygon_2& polygon);
    template <typename T> void make_round_poly(const Point_2& centre, double radius, T& poly);
    template <typename T> void make_round_poly(const Point_2& centre, double radius1, double radius2,
                                               int nPts, double angInt, double ang, T& poly);
//    template <typename T, typename U> void smooth_dt (const Point_set_3& pointCloud, T& dt);
    template <typename T> Polygon_2 calc_bbox_poly(const T& inputPts);
    template <typename T> T offset_polygon_geos(T poly, double offset);
    template <typename T> void pop_back_if_equal_to_front(CGAL::Polygon_2<T>& poly);

    struct Array_traits {
        typedef std::array<EPICK::FT, 3>  Custom_point;
        struct Equal_3 {
            bool operator()(const Custom_point& p, const Custom_point& q) const {
                return (p == q);
            }
        };
        struct Less_xyz_3 {
            bool operator()(const Custom_point& p, const Custom_point& q) const {
                return std::lexicographical_compare(p.begin(), p.end(), q.begin(), q.end());
            }
        };
        Equal_3 equal_3_object() const { return Equal_3(); }
        Less_xyz_3 less_xyz_3_object() const { return Less_xyz_3(); }
    };
}

#endif //CITY4CFD_GEOMUTILS_H