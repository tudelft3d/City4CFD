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

#include "PointCloud.h"

#include "io.h"
#include "geomutils.h"
#include "Config.h"
#include "PolyFeature.h"

#include <boost/locale.hpp>

PointCloud::PointCloud()  = default;
PointCloud::~PointCloud() = default;

void PointCloud::random_thin_pts() {
    if (Config::get().terrainThinning > 0 + global::smallnum && Config::get().las_files.empty()) {
        std::cout <<"\nRandomly thinning terrain points" << std::endl;
        _pointCloudTerrain.remove(CGAL::random_simplify_point_set(_pointCloudTerrain,
                                                                  Config::get().terrainThinning),
                                  _pointCloudTerrain.end());
        _pointCloudTerrain.collect_garbage();
        std::cout << "    Terrain points after thinning: " << _pointCloudTerrain.size() << std::endl;
    }
}

void PointCloud::smooth_terrain() {
    std::cout << "\nSmoothing terrain" << std::endl;
    DT dt(_pointCloudTerrain.points().begin(), _pointCloudTerrain.points().end());
    geomutils::smooth_dt<DT, EPICK>(_pointCloudTerrain, dt);

    //-- Return new pts to the point cloud
    _pointCloudTerrain.clear();
    for (auto& pt : dt.points()) _pointCloudTerrain.insert(pt);
    _pointCloudTerrain.add_property_map<bool> ("is_building_point", false);
}

void PointCloud::create_flat_terrain(const PolyFeatures& lsFeatures) {
    std::cout << "\nCreating flat terrain" << std::endl;
    for (auto& f : lsFeatures) {
        if (f->get_poly().rings().empty()) {
//            std::cout << "Empty polygon?" << std::endl; //todo investigate this
            continue;
        }
        for (auto& pt : f->get_poly().outer_boundary()) {
            _pointCloudTerrain.insert(Point_3(pt.x(), pt.y(), 0.0));
        }
    }
    _pointCloudTerrain.add_property_map<bool> ("is_building_point", false);
}

void PointCloud::set_flat_terrain() {
    Point_set_3 flatPC;
    for (auto& pt : _pointCloudTerrain.points()) {
        flatPC.insert(Point_3(pt.x(), pt.y(), 0.));
    }
    _pointCloudTerrain = flatPC;
    _pointCloudTerrain.add_property_map<bool> ("is_building_point", false);
}

void PointCloud::flatten_polygon_pts(const PolyFeatures& lsFeatures) {
    std::cout << "\n    Flattening surfaces" << std::endl;
    std::map<int, Point_3> flattenedPts;

    //-- Construct a connectivity map and remove duplicates along the way
    auto is_building_pt = _pointCloudTerrain.property_map<bool>("is_building_point").first;
    std::unordered_map<Point_3, int> pointCloudConnectivity;
    auto it = _pointCloudTerrain.points().begin();
    int count = 0;
    while (it != _pointCloudTerrain.points().end()) {
        auto itPC = pointCloudConnectivity.find(*it);
        if (itPC != pointCloudConnectivity.end()) {
            _pointCloudTerrain.remove(_pointCloudTerrain.begin() + count);
        } else {
            pointCloudConnectivity[*it] = count;
            ++it;
            ++count;
        }
    }
    _pointCloudTerrain.collect_garbage();

    //-- Construct search tree from ground points
    SearchTree searchTree(_pointCloudTerrain.points().begin(), _pointCloudTerrain.points().end());

    //-- Perform averaging
    for (auto& f : lsFeatures) {
        auto it = Config::get().flattenSurfaces.find(f->get_output_layer_id());
        if (it != Config::get().flattenSurfaces.end()) {
            f->flatten_polygon_inner_points(_pointCloudTerrain, flattenedPts, searchTree, pointCloudConnectivity);
        }
    }

    //-- Change points with flattened values
    int pcOrigSize = _pointCloudTerrain.points().size();
    for (auto& it : flattenedPts) {
        _pointCloudTerrain.insert(it.second);
    }
    for (int i = 0; i < pcOrigSize; ++i) {
        auto it = flattenedPts.find(i);
        if (it != flattenedPts.end()) {
            _pointCloudTerrain.remove(i);
            flattenedPts.erase(i);
        }
    }
    _pointCloudTerrain.collect_garbage();
}

SearchTreePtr PointCloud::make_search_tree_buildings() {
    return std::make_shared<SearchTree>(_pointCloudBuildings.points().begin(),
                                        _pointCloudBuildings.points().end());
}

void PointCloud::read_point_clouds() {
    //-- Translation matrix in relation to the point of interest
    CGAL::Aff_transformation_3<EPICK> translate(CGAL::TRANSLATION,
                                                CGAL::Vector_3<EPICK>(-Config::get().pointOfInterest.x(),
                                                                      -Config::get().pointOfInterest.y(),
                                                                      0));
    //-- Automatic input using LAS files or manually defining ground and buildings
    if (!Config::get().las_files.empty()) {
        //-- Add all used LAS classes to one vector
        std::vector<int> usedClasses;
        usedClasses.reserve(Config::get().las_classes_building.size() + Config::get().las_classes_ground.size());
        for (auto classif : Config::get().las_classes_ground)   usedClasses.push_back(classif);
        for (auto classif : Config::get().las_classes_building) usedClasses.push_back(classif);

        int readPts = 0;
        for (auto& pointFile: Config::get().las_files) {
            std::cout << "Reading LAS/LAZ file: " << pointFile << "\n";

            LASreadOpener lasreadopener;
            lasreadopener.set_file_name(pointFile.c_str());
            //-- set to compute bounding box
            lasreadopener.set_populate_header(true);
            LASreader* lasreader = lasreadopener.open();

            try {
                //-- check if file is open
                if (lasreader == nullptr)
                    throw(std::runtime_error(std::string("Error reading LAS/LAZ file.")));

                LASheader header = lasreader->header;

                //-- read each point 1-by-1
                uint32_t pointCount = header.number_of_point_records;

                std::cout << "    " << boost::locale::as::number << pointCount << " points in the file\n";
                double percentLeft = 1 - (Config::get().terrainThinning / 100);
                if (Config::get().terrainThinning > 0) {
                    std::cout << "    Skipping " << Config::get().terrainThinning << "% points" << std::endl;
                }
                else {
                    std::cout << "    all points used, no skipping" << std::endl;
                }

                //-- Read defined classes or use CSF to determine ground from non-ground
                if (!Config::get().las_classes_ground.empty() &&
                    !Config::get().las_classes_building.empty()) {
                    std::cout << "    Reading LAS classes: ";
                    for (int i : Config::get().las_classes_ground) {
                        std::cout << i << " ";
                    }
                    std::cout << "(ground) and ";
                    for (int i : Config::get().las_classes_building) {
                        std::cout << i << " ";
                    }
                    std::cout << "(buildings)" << std::endl;

                    int i = 0;
                    double currPercent = 0.;
                    IO::print_progress_bar(0);
                    while (lasreader->read_point()) {
                        LASpoint const& p = lasreader->point;
                        //-- set the thinning filter
                        if (currPercent >= 1.) {
                            //-- set the classification filter
                            if (std::find(usedClasses.begin(), usedClasses.end(), (int) p.classification)
                                != usedClasses.end()) {
                                this->add_elevation_point(p, translate);
                                ++readPts;
                            }
                            currPercent -= 1.;
                        } else {
                            currPercent += percentLeft;
                        }
                        if (i % (pointCount / 200) == 0) {
                            IO::print_progress_bar(100 * (i / double(pointCount)));
                        }
                        i++;
                    }
                    IO::print_progress_bar(100);
                    std::clog << std::endl;

                    lasreader->close();
                } else { //todo implement CSF here
                    std::cout << "    No point cloud classification information given. "
                            << "Will filter ground using the Cloth Simulation Filter" << std::endl;
                    CSF csf;
                    //todo csf params setting
                    csf.params.bSloopSmooth = true;
                    csf.params.class_threshold = 0.5;
                    csf.params.cloth_resolution = 2;
                    csf.params.interations = 500;
                    csf.params.rigidness = 3;
                    csf.params.time_step = 0.65;
                    //-- Read the data
                    int i = 0;
                    double currPercent = 0.;
                    IO::print_progress_bar(0);
                    while (lasreader->read_point()) {
                        LASpoint const& p = lasreader->point;
                        //-- set the thinning filter
                        if (currPercent >= 1.) {
                            //-- set the classification filter
                            this->add_elevation_point_to_csf(p, translate, csf);

                            ++readPts;
                            currPercent -= 1.;
                        } else {
                            currPercent += percentLeft;
                        }
                        if (i % (pointCount / 200) == 0) {
                            IO::print_progress_bar(100 * (i / double(pointCount)));
                        }
                        i++;
                    }
                    IO::print_progress_bar(100);
                    std::clog << std::endl;

                    lasreader->close();

                    //-- Perform CS filtering
                    std::cout << "    Applying CSF to the point cloud" << std::endl;
                    std::vector<int> groundIndices, offGroundIndices;
                    csf.do_filtering(groundIndices, offGroundIndices, false);

                    //-- Add points to point clouds
                    for (auto idx : groundIndices) {
                        Point_3 pt(csf.getPointCloud().at(idx).x,
                                   csf.getPointCloud().at(idx).z,
                                   -csf.getPointCloud().at(idx).y);
                        _pointCloudTerrain.insert(pt);
                    }
                    for (auto idx : offGroundIndices) {
                        Point_3 pt(csf.getPointCloud().at(idx).x,
                                   csf.getPointCloud().at(idx).z,
                                   -csf.getPointCloud().at(idx).y);
                        _pointCloudBuildings.insert(pt);
                    }
                }
                std::clog << "Points read: " << readPts << "\n" << std::endl;
                _pointCloudTerrain.add_property_map<bool>("is_building_point", false);
            }
            catch (std::exception& e) {
                if (lasreader != nullptr) lasreader->close();
                throw;
            }
        }
    } else {
        //-- Second input option is to explicitly define ground and/or building points
        //- Read explicitly defined ground points
        std::cout << "Reading explicitly defined ground and/or building points" << std::endl;
        if (!Config::get().ground_xyz.empty()) {
            std::cout << "Reading ground points" << std::endl;
            IO::read_point_cloud(Config::get().ground_xyz, _pointCloudTerrain);
            _pointCloudTerrain.add_property_map<bool>("is_building_point", false);

            std::cout << "    Points read: " << _pointCloudTerrain.size() << std::endl;
        } else {
            std::cout << "INFO: Did not find any ground points! Will calculate ground as a flat surface." << std::endl;
            std::cout << "WARNING: Ground height of buildings can only be approximated. "
                      << "If you are using point cloud to reconstruct buildings, building height estimation can be wrong.\n"
                      << std::endl;
        }

        //- Read explicitly defined building points
        if (!Config::get().buildings_xyz.empty()) {
            std::cout << "Reading building points" << std::endl;
            IO::read_point_cloud(Config::get().buildings_xyz, _pointCloudBuildings);
            if (_pointCloudBuildings.empty()) throw std::invalid_argument("Didn't find any building points!");

            std::cout << "    Points read: " << _pointCloudBuildings.size() << std::endl;
        }
    }
}

void PointCloud::add_elevation_point(const LASpoint& laspt, const CGAL::Aff_transformation_3<EPICK>& translate) {
    int classID = laspt.classification;
    Point_3 pt(laspt.get_x(), laspt.get_y(), laspt.get_z());
    std::vector<int> classList = Config::get().las_classes_ground;
    if (std::find(classList.begin(), classList.end(), classID)
                  != classList.end()) {
        _pointCloudTerrain.insert(pt.transform(translate));
    }
    classList = Config::get().las_classes_building;
    if (std::find(classList.begin(), classList.end(), classID)
        != classList.end()) {
        _pointCloudBuildings.insert(pt.transform(translate));
    }
}

void PointCloud::add_elevation_point_to_csf(const LASpoint& laspt,
                                            const CGAL::Aff_transformation_3<EPICK>& translate,
                                            CSF& csf) {
    Point_3 pt(laspt.get_x(), laspt.get_y(), laspt.get_z());
    pt = pt.transform(translate);
    csf.addPoint(pt.x(), pt.y(), pt.z());
}

/*
bool PointCloud::check_bounds(const double xmin, const double xmax, const double ymin, const double ymax) {
    if ((xmin < _maxxradius || xmax > _minxradius) &&
        (ymin < _maxyradius || ymax > _minyradius)) {
        return true;
    }
    return false;
}
*/

Point_set_3& PointCloud::get_terrain() {
    return _pointCloudTerrain;
}
Point_set_3& PointCloud::get_buildings() {
    return _pointCloudBuildings;
}

const Point_set_3& PointCloud::get_terrain() const {
    return _pointCloudTerrain;
}

const Point_set_3& PointCloud::get_buildings() const {
    return _pointCloudBuildings;
}