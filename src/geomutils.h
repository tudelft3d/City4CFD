#ifndef CITYCFD_GEOMUTILS_H
#define CITYCFD_GEOMUTILS_H

#include "config.h"

namespace geomutils {
    double  avg(const std::vector<double>& values);
    double  percentile(std::vector<double> values, const double percentile);
    bool    point_in_circle(const Point_3& pt, const Point_2& center, const double& radius);
    void    cdt_to_mesh(CDT& cdt, Mesh& mesh, const int surfaceLayerID = -9999);
    void    mark_domains(CDT& cdt, PolyFeatures features = {});
    void    mark_domains(CDT& ct, const Face_handle& start, int index,
                         std::list<CDT::Edge>& border, PolyFeatures& features);
    bool    check_layer_normal(const Face_handle& fh, int surfaceLayer);
    void    shorten_long_poly_edges(Polygon_2& poly, double maxLen = config::edgeMaxLen);
    Point_2 rotate_pt(Point_2& pt, const double angle, Point_2 centerPt = Point_2(0, 0));
    void    interpolate_poly_from_pc(const Polygon_2& poly, std::vector<double>& heights, const Point_set_3& pointCloud);
    bool    polygons_in_contact(const Polygon_with_holes_2& firstPoly, const Polygon_with_holes_2& secondPoly);

    //-- Templated functions
    template <typename T> bool point_in_poly(const T& pt2, const Polygon_with_holes_2& polygon);
    template <typename T> bool point_in_poly(const T& pt2, const Polygon_2& polygon);
    template <typename T> void make_round_poly(Point_2& centre, double radius, T& poly);
    template <typename T> void make_round_poly(Point_2& centre, double radius1, double radius2,
                                               int nPts, double angInt, double ang, T& poly);
    template <typename T, typename U> void smooth_dt (const Point_set_3& pointCloud, T& dt);
    template <typename T> Polygon_2 calc_bbox_poly(const T& inputPts);

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

#endif //CITYCFD_GEOMUTILS_H