#ifndef CITYCFD_GEOMTOOLS_H
#define CITYCFD_GEOMTOOLS_H

#include "definitions.h"
#include "TopoFeature.h"

double avg(const std::vector<double>& values);
double median(std::vector<double> values);
double percentile(std::vector<double> values, const double percentile);
bool   check_inside(const Point_3& pt2, const Polygon_with_holes_2& polygon);
bool   point_in_circle(const Point_3& pt, const Point_2& center, const double& radius);
void   cdt_to_mesh(const CDT& cdt, Mesh& mesh);
void   mark_domains(CDT& ct, Face_handle start, int index, std::list<CDT::Edge>& border);
void   mark_domains(CDT& cdt);

#endif //CITYCFD_GEOMTOOLS_H