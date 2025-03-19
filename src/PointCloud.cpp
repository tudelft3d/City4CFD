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

#include "PointCloud.h"

#include "io.h"
#include "geomutils.h"
#include "Config.h"
#include "Building.h"
#include "Quadtree/Quadtree.h"

#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/smooth_shape.h>
#include <CGAL/wlop_simplify_and_regularize_point_set.h>
#include <CGAL/create_offset_polygons_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/Straight_skeleton_builder_2.h>
#include <CGAL/compute_outer_frame_margin.h>

void PointCloud::random_thin_pts() {
    if (Config::get().terrainThinning > 0 + global::smallnum) {
        std::cout <<"\nRandomly thinning terrain points" << std::endl;
        m_pointCloudTerrain.remove(CGAL::random_simplify_point_set(m_pointCloudTerrain,
                                                                   Config::get().terrainThinning),
                                   m_pointCloudTerrain.end());
        m_pointCloudTerrain.collect_garbage();
        std::cout << "    Terrain points after thinning: " << m_pointCloudTerrain.size() << std::endl;
    }
}

void PointCloud::smooth_terrain() {
    typedef CGAL::Parallel_if_available_tag Concurrency_tag;

    std::cout << "Smoothing terrain" << std::endl;

    //-- WLOP simplification and regularization
    double retainPercentage = 100;
    int& maxTerrainPts = Config::get().maxSmoothPts;
    if (maxTerrainPts > 0 && m_pointCloudTerrain.size() > maxTerrainPts) {
        retainPercentage = (double)maxTerrainPts / (double)m_pointCloudTerrain.size() * 100.;
        std::cout << "    Performing additional (optimized) terrain thinning to " << maxTerrainPts << " points" << std::endl;
    }

    std::cout << "    Smoothing terrain 1/3..." << std::flush;
    const double neighborRadius = 0.5;   // neighbors size.
    Point_set_3 simplPts;
    CGAL::wlop_simplify_and_regularize_point_set<Concurrency_tag>
            (m_pointCloudTerrain, simplPts.point_back_inserter(),
             CGAL::parameters::select_percentage(retainPercentage).
                     neighbor_radius (neighborRadius));
    m_pointCloudTerrain.clear();

    std::cout << "\r    Smoothing terrain 2/3..." << std::flush;

    //-- Create CDT of current terrain pts
    DT dt;
    dt.insert(simplPts.points().begin(), simplPts.points().end());
    simplPts.clear(); // end of scope for simplePts

    //-- Make a mesh out of DT
    Mesh mesh;
    geomutils::dt_to_mesh(dt, mesh);
    dt.clear(); // end of scope for the dt

    //-- Isotropic remeshing
    const double targetEdgeLength = 0;
    const unsigned int nbIter =  10;
    PMP::remove_degenerate_faces(mesh);
    PMP::isotropic_remeshing(faces(mesh), targetEdgeLength, mesh,
                             PMP::parameters::number_of_iterations(nbIter)
                                     );

    //-- Smoothing
    std::cout << "\r    Smoothing terrain 3/3..." << std::flush;
    const double time = 1;
    PMP::smooth_shape(mesh, time, CGAL::parameters::number_of_iterations(Config::get().nSmoothIterations));

    std::cout << "\r    Smoothing terrain...done" << std::endl;

    //-- Mesh back to points
//    m_pointCloudTerrain.clear();
    for (auto& pt : mesh.points()) {
        m_pointCloudTerrain.insert(pt);
    }
}

void PointCloud::terrain_points_in_polygon(BuildingsPtr& features) {
    typedef Quadtree_node<EPICK, Point_set_3> Point_index;
    Point_index pointCloudIndex;
    auto& pointCloud = m_pointCloudTerrain;

    pointCloudIndex.compute_extent(pointCloud);
    for (auto pointIndex = pointCloud.begin();
         pointIndex != pointCloud.end();
         ++pointIndex) {
        pointCloudIndex.insert_point(pointCloud, *pointIndex);
    }
    pointCloudIndex.optimise(pointCloud, 100, 10);

    //-- Find points belonging to individual buildings
    for (auto& f: features) {
        auto poly = f->get_poly().get_cgal_type();
//        const double offset = 1.; // offset hardcoded
//        auto offsetPoly = geomutils::offset_polygon_geos(poly, offset);
        auto& offsetPoly = poly; // temp

        std::vector<Point_index*> intersected_nodes;
        pointCloudIndex.find_intersections(intersected_nodes, offsetPoly.bbox().xmin(), offsetPoly.bbox().xmax(),
                                           offsetPoly.bbox().ymin(), offsetPoly.bbox().ymax());
        for (auto const& node: intersected_nodes) {
            for (auto const& pointIdx: node->points) {
                if (geomutils::point_in_poly_and_boundary(pointCloud.point(pointIdx), poly)) {
                    pointCloud.remove(pointIdx);
                }
                //todo temp add to building
//                f->insert_terrain_point(pointCloud.point(pointIdx)); //use bbox for roofer
            }
        }
    }
    pointCloud.collect_garbage();
}

void PointCloud::create_flat_terrain(const PolyFeaturesPtr& lsFeatures) {
    std::cout << "\nCreating flat terrain" << std::endl;
    for (auto& f : lsFeatures) {
        if (f->get_poly().rings().empty()) {
//            std::cout << "Empty polygon?" << std::endl; //todo investigate this
            continue;
        }
        for (auto& pt : f->get_poly().outer_boundary()) {
            m_pointCloudTerrain.insert(Point_3(pt.x(), pt.y(), 0.0));
        }
    }
}

void PointCloud::set_flat_terrain() {
    Point_set_3 flatPC;
    for (auto& pt : m_pointCloudTerrain.points()) {
        flatPC.insert(Point_3(pt.x(), pt.y(), 0.));
    }
    m_pointCloudTerrain = flatPC;
}

void PointCloud::flatten_polygon_pts(const PolyFeaturesPtr& lsFeatures,
                                     std::vector<EPECK::Segment_3>& constrainedEdges,
                                     std::vector<std::pair<Polygon_with_holes_2, int>>& newPolys) {
    std::cout << "\n    Flattening surfaces" << std::endl;
    std::map<int, Point_3> flattenedPts;

    //-- Construct a connectivity map and remove duplicates along the way
    std::unordered_map<Point_3, int> pointCloudConnectivity;
    auto it = m_pointCloudTerrain.points().begin();
    int count = 0;
    while (it != m_pointCloudTerrain.points().end()) {
        auto itPC = pointCloudConnectivity.find(*it);
        if (itPC != pointCloudConnectivity.end()) {
            m_pointCloudTerrain.remove(m_pointCloudTerrain.begin() + count);
        } else {
            pointCloudConnectivity[*it] = count;
            ++it;
            ++count;
        }
    }
    m_pointCloudTerrain.collect_garbage();

    //-- Construct search tree from ground points
    SearchTree searchTree(m_pointCloudTerrain.points().begin(),
                          m_pointCloudTerrain.points().end(),
                          Config::get().searchtree_bucket_size);

    //-- Perform flattening
    PolyFeaturesPtr vertBorders;
    for (auto& f : lsFeatures) {
        auto ita = Config::get().flattenSurfaces.find(f->get_output_layer_id());
        if (ita != Config::get().flattenSurfaces.end()) {
            // flatten points
            bool isNextToBuilding = false;
            if(f->flatten_polygon_inner_points(m_pointCloudTerrain, flattenedPts, searchTree, pointCloudConnectivity,
                                               constrainedEdges, newPolys, isNextToBuilding)) {
                // add to list if constructing vertical borders
                if (std::find(Config::get().flattenVertBorder.begin(), Config::get().flattenVertBorder.end(),
                              f->get_output_layer_id()) != Config::get().flattenVertBorder.end() && !isNextToBuilding) {
                    vertBorders.push_back(f);
                }
            }
        }
    }
    //-- Handle border for flattened polys
    if (!vertBorders.empty()) this->buffer_flat_edges(vertBorders, constrainedEdges);

    //-- Change points with flattened values
    int pcOrigSize = m_pointCloudTerrain.points().size();
    for (auto& itp : flattenedPts) {
        m_pointCloudTerrain.insert(itp.second);
    }
    for (int i = 0; i < pcOrigSize; ++i) {
        auto itp = flattenedPts.find(i);
        if (itp != flattenedPts.end()) {
            m_pointCloudTerrain.remove(i);
            flattenedPts.erase(i);
        }
    }
    m_pointCloudTerrain.collect_garbage();
}

void PointCloud::buffer_flat_edges(const PolyFeaturesPtr& avgFeatures,
                                   std::vector<EPECK::Segment_3>& constrainedEdges) {
//    std::cout << "\n    Buffering flattening surfaces" << std::endl;
    //-- Add buffer around flattened polygons
    typedef CGAL::Straight_skeleton_builder_traits_2<EPICK> SsBuilderTraits;
    typedef CGAL::Straight_skeleton_2<EPICK> Ss;
    typedef CGAL::Straight_skeleton_builder_2<SsBuilderTraits, Ss> SsBuilder;
    typedef CGAL::Polygon_offset_builder_traits_2<EPICK> OffsetBuilderTraits;
    typedef CGAL::Polygon_offset_builder_2<Ss, OffsetBuilderTraits, Polygon_2> OffsetBuilder;

    typedef std::shared_ptr<Polygon_2> ContourPtr;
    typedef std::vector<ContourPtr> ContourSequence;
    // get info using the original point cloud
//    std::vector<double> offsets{0.001, 0.2};
    std::vector<double> offsets{0.001};
    std::vector<Polygon_3> polyList;
    #pragma omp parallel for
    for (int i = 0; i < avgFeatures.size(); ++i) {
    //for (auto& f: avgFeatures) { // MSVC doesn't like range loop with OMP
        auto& poly = avgFeatures[i]->get_poly().outer_boundary();
        // set the frame
        auto margin = CGAL::compute_outer_frame_margin(poly.begin(), poly.end(), offsets.back());
        CGAL::Bbox_2 bbox = CGAL::bbox_2(poly.begin(), poly.end());
        // Compute the boundaries of the frame
        double fxmin = bbox.xmin() - *margin;
        double fxmax = bbox.xmax() + *margin;
        double fymin = bbox.ymin() - *margin;
        double fymax = bbox.ymax() + *margin;
        // Create the rectangular frame
        Point_2 frame[4] = {Point_2(fxmin, fymin), Point_2(fxmax, fymin), Point_2(fxmax, fymax), Point_2(fxmin, fymax)
        };

        // Instantiate the skeleton builder
        SsBuilder ssb;
        // Enter the frame
        ssb.enter_contour(frame, frame + 4);
        // Enter the polygon as a hole of the frame (NOTE: as it is a hole we insert it in the opposite orientation)
        poly.reverse_orientation();
        ssb.enter_contour(poly.begin(), poly.end());
        // Construct the skeleton
        auto ss = ssb.construct_skeleton();
        // Proceed only if the skeleton was correctly constructed.
        if (ss) {
            for (auto& offset: offsets) {
                // Instantiate the container of offsetSmall contours
                ContourSequence offset_contours;
                // Instantiate the offsetSmall builder with the skeleton
                OffsetBuilder ob(*ss);
                // Obtain the offsetSmall contours
                ob.construct_offset_contours(offset, std::back_inserter(offset_contours));
                // Locate the offsetSmall contour that corresponds to the frame
                // That must be the outmost offsetSmall contour, which in turn must be the one
                // with the largetst unsigned area.
                auto offsend = offset_contours.end();
                double lLargestArea = 0.0;
                for (auto j = offset_contours.begin(); j != offset_contours.end(); ++j) {
                    double lArea = CGAL_NTS abs((*j)->area()); //Take abs() as  Polygon_2::area() is signed.
                    if (lArea > lLargestArea) {
                        offsend = j;
                        lLargestArea = lArea;
                    }
                }
                // Remove the offsetSmall contour that corresponds to the frame.
                offset_contours.erase(offsend);

#ifndef NDEBUG
                if (offset_contours.size() != 1) {
                    std::cout << "DEBUG: SOMETHING WRONG WITH SKELETON, NUMBER OF POLYS: " << offset_contours.size()
                              << std::endl;
                }
#endif
                std::vector<double> height;
                Polygon_2& offsetPoly2 = *(offset_contours.back());
                geomutils::interpolate_poly_from_pc(offsetPoly2, height, m_pointCloudTerrain);
                Polygon_3 offsetPoly3;
                for (auto j = 0; j < offsetPoly2.size(); ++j) {
                    offsetPoly3.push_back(ePoint_3(offsetPoly2[j].x(), offsetPoly2[j].y(), height[j]));
                }
                #pragma omp critical
                polyList.push_back(offsetPoly3);
            }
        }
    }
    //-- Add only non-adjacent edges to constrain
    CDT cdt;
    for (const auto& ring : polyList) {
        cdt.insert_constraint(ring.begin(), ring.end());
    }
    geomutils::mark_domains(cdt);

    for (CDT::Edge e : cdt.constrained_edges()) {
        CDT::Face_handle f1 = e.first;
        CDT::Face_handle f2 = f1->neighbor(e.second);
        if (f1->info().nesting_level == 0 || f2->info().nesting_level == 0) {
            ePoint_3 p1 = e.first->vertex((e.second + 1) % 3)->point();
            ePoint_3 p2 = e.first->vertex((e.second + 2) % 3)->point();
            constrainedEdges.emplace_back(p1, p2);
        }
    }
}

void PointCloud::read_point_clouds() {
    //- Read ground points
    if (!Config::get().ground_xyz.empty()) {
        std::cout << "Reading ground points" << std::endl;
        IO::read_point_cloud(Config::get().ground_xyz, m_pointCloudTerrain);

        std::cout << "    Points read: " << m_pointCloudTerrain.size() << std::endl;
    } else {
        std::cout << "INFO: Did not find any ground points! Will calculate ground as a flat surface." << std::endl;
        std::cout << "WARNING: Ground elevation of buildings can only be approximated. "
                  << "If you are using point cloud to reconstruct buildings, building height estimation can be wrong.\n"
                  << std::endl;
    }

    //- Read building points
    if (!Config::get().buildings_xyz.empty()) {
        std::cout << "Reading building points" << std::endl;
        IO::read_point_cloud(Config::get().buildings_xyz, m_pointCloudBuildings);
        if (m_pointCloudBuildings.empty()) throw std::invalid_argument("Didn't find any building points!");

        std::cout << "    Points read: " << m_pointCloudBuildings.size() << std::endl;
    }
}

Point_set_3& PointCloud::get_terrain() {
    return m_pointCloudTerrain;
}
Point_set_3& PointCloud::get_buildings() {
    return m_pointCloudBuildings;
}

const Point_set_3& PointCloud::get_terrain() const {
    return m_pointCloudTerrain;
}

const Point_set_3& PointCloud::get_buildings() const {
    return m_pointCloudBuildings;
}
