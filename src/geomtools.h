#ifndef CITYCFD_GEOMTOOLS_H
#define CITYCFD_GEOMTOOLS_H

#include "definitions.h"
#include "TopoFeature.h"

class PolyFeature;

namespace geomtools {
    double avg(const std::vector<double>& values);
    double median(std::vector<double> values);
    double percentile(std::vector<double> values, const double percentile);
    bool   check_inside(const Point_3& pt2, const Polygon_with_holes_2& polygon);
    bool   point_in_circle(const Point_3& pt, const Point_2& center, const double& radius);
    void   cdt_to_mesh(CDT& cdt, Mesh& mesh, const int surfaceLayerID = -9999);
    void   mark_domains(CDT& ct, Face_handle start, int index, std::list<CDT::Edge>& border);
    void   mark_domains(CDT& cdt);
    void   mark_surface_layer(CDT& cdt, const int surfaceLayerIdx);
    void   mark_surface_layer(CDT& cdt, std::vector<PolyFeature*> features);
    void   mark_surface_layer(CDT& ct, Face_handle start, int index, std::list<CDT::Edge>& border,
                            std::vector<PolyFeature*>& features);
}

#endif //CITYCFD_GEOMTOOLS_H