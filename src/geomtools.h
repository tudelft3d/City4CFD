#ifndef CITYCFD_GEOMTOOLS_H
#define CITYCFD_GEOMTOOLS_H

#include "definitions.h"
#include "TopoFeature.h"

class PolyFeature;

namespace geomtools {
    double avg(const std::vector<double>& values);
    double percentile(std::vector<double> values, const double percentile);
    bool   point_in_circle(const Point_3& pt, const Point_2& center, const double& radius);
    void   cdt_to_mesh(CDT& cdt, Mesh& mesh, const int surfaceLayerID = -9999);
    void   mark_domains(CDT& cdt, PolyFeatures features = {});
    void   mark_domains(CDT& ct, const Face_handle& start, int index,
                        std::list<CDT::Edge>& border, PolyFeatures& features);
    void   shorten_long_poly_edges(Polygon_2& poly);

    //-- Templated functions
    template <typename T> bool point_in_poly(const T& pt2, const Polygon_with_holes_2& polygon);
    template <typename T> void make_round_poly(Point_2& centre, double radius, T& poly);
    template <typename T, typename U> void smooth_dt (const Point_set_3& pointCloud, T& dt);
}

#endif //CITYCFD_GEOMTOOLS_H