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

#include "Boundary.h"

#include "geomutils.h"

#include <CGAL/centroid.h>

Boundary::Boundary()
        : TopoFeature(), m_sideOutputPts() {}

Boundary::Boundary(const int outputLayerID)
        : TopoFeature(outputLayerID), m_sideOutputPts() {}

//-- Static member definition
std::vector<Point_3> Boundary::s_outerPts;
double               Boundary::s_outerBndHeight;

void Boundary::set_bnd_poly(Polygon_2& bndPoly, Polygon_2& pcBndPoly, Polygon_2& startBufferPoly) {
    Polygon_2 bufferPoly;
    Point_2 center = CGAL::centroid(bndPoly.begin(), bndPoly.end(),
                                    CGAL::Dimension_tag<0>());
    double bufferLen = 0.;

    //-- Set bndPoly and pcBndPoly depending on the buffer setup
    if (Config::get().domainBuffer > -global::largnum) {
        bufferLen = Config::get().domainBuffer / 100.;
        for (auto& pt: bndPoly)
            bufferPoly.push_back(pt + (pt - center) * bufferLen);

        if (Config::get().domainBuffer < 0) {
            startBufferPoly = bufferPoly;
        } else {
            startBufferPoly = bndPoly;
            bndPoly = bufferPoly;
        }
    } else {
        startBufferPoly = bndPoly;
        bufferLen = 1.; // arbitrary value for pcBndPoly
    }
    for (auto& pt : startBufferPoly)
        pcBndPoly.push_back(pt - (pt - center) * 0.001 * std::abs(bufferLen));
}

//-- Deactivate point cloud points that are out of bounds
void Boundary::set_bounds_to_buildings_pc(Point_set_3& pointCloud, const Polygon_2& pcBndPoly) {
    //-- Remove points out of the boundary region
    auto it = pointCloud.points().begin();
    int count = 0;
    while (it != pointCloud.points().end()) {
        if (!geomutils::point_in_poly(*it, pcBndPoly)) {
            pointCloud.remove(pointCloud.begin() + count);
        } else {
            ++it;
            ++count;
        }
    }
    pointCloud.collect_garbage(); // Free removed points from the memory
}

void Boundary::set_bounds_to_terrain_pc(Point_set_3& pointCloud, const Polygon_2& bndPoly,
                                        const Polygon_2& pcBndPoly, const Polygon_2& startBufferPoly) {
    //-- Remove points out of the boundary region
    Boundary::set_bounds_to_buildings_pc(pointCloud, pcBndPoly);

    std::vector<double> bndHeights;
    geomutils::interpolate_poly_from_pc(bndPoly, bndHeights, pointCloud);
    s_outerBndHeight = geomutils::avg(bndHeights); // Height for buffer (for now) - average of outer pts
    Config::get().logSummary << "Domain edge elevation: " << s_outerBndHeight << std::endl;

    if (Config::get().domainBuffer > -global::largnum + global::smallnum) {
        /*
        for (auto& pt : startBufferPoly) {
            pointCloud.insert(Point_3(pt.x(), pt.y(), s_outerBndHeight));
        }
         */
        for (auto& pt : bndPoly) {
            s_outerPts.emplace_back(Point_3(pt.x(), pt.y(), s_outerBndHeight));
            pointCloud.insert(s_outerPts.back());
        }
    } else {
        int i = 0;
        for (auto& pt : bndPoly) {
            s_outerPts.emplace_back(Point_3(pt.x(), pt.y(), bndHeights[i++]));
            pointCloud.insert(s_outerPts.back());
        }
    }
    s_outerPts.push_back(s_outerPts.front());
}

void Boundary::prep_output() {
    m_sideOutputPts = s_outerPts;
}

//-- Find all outerPts along the defined edge
void Boundary::prep_output(Vector_2 edge) {
    edge /= sqrt(edge.squared_length());

    //-- Search the outerPts for the same vector
    for (int i = 0; i < s_outerPts.size() - 1; ++i) {
        Vector_2 checkEdge(s_outerPts[i + 1].x() - s_outerPts[i].x(),
                           s_outerPts[i + 1].y() - s_outerPts[i].y());
        checkEdge /= sqrt(checkEdge.squared_length());

        if (edge * checkEdge > 1 - global::smallnum && edge * checkEdge < 1 + global::smallnum) {
            m_sideOutputPts.push_back(s_outerPts[i]);
            m_sideOutputPts.push_back(s_outerPts[i + 1]);

            bool collinear = true;
            while (collinear) {
                int j = i + 2;
                if (s_outerPts.begin() + j == s_outerPts.end()) break;

                Vector_2 nextEdge(s_outerPts[j].x() - s_outerPts[i + 1].x(),
                                  s_outerPts[j].y() - s_outerPts[i + 1].y());
                nextEdge /= sqrt(nextEdge.squared_length());

                if (nextEdge * edge > 1 - global::smallnum && nextEdge * edge < 1 + global::smallnum) {
                    m_sideOutputPts.push_back(s_outerPts[j]);
                } else {
                    collinear = false;
                }
                ++i;
            }
            return;
        }
    }
    throw city4cfd_error("Cannot find side for output!");
}

std::vector<double> Boundary::get_outer_bnd_bbox() {
    double maxx(-global::largnum), maxy(-global::largnum), maxz(-global::largnum);
    double minx(global::largnum),  miny(global::largnum),  minz(global::largnum);

    for (auto& pt : s_outerPts) {
        if (pt.x() > maxx) maxx = pt.x();
        else if (pt.x() < minx) minx = pt.x();

        if (pt.y() > maxy) maxy = pt.y();
        else if (pt.y() < miny) miny = pt.y();

        if (pt.z() > maxz) maxz = pt.z();
        else if (pt.z() < minz) minz = pt.z();
    }

    return std::vector<double> {minx, miny, minz, maxx, maxy, maxz};
}