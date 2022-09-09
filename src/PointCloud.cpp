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
    //-- Read ground points
    if (!Config::get().las_files.empty()) { // add check for las/laz from config file
        for (auto& pointFile: Config::get().las_files) {
            std::clog << "Reading LAS/LAZ file: " << pointFile.filename << std::endl;

            LASreadOpener lasreadopener;
            lasreadopener.set_file_name(pointFile.filename.c_str());
            //-- set to compute bounding box
            lasreadopener.set_populate_header(true);
            LASreader* lasreader = lasreadopener.open();

            try {
                //-- check if file is open
                if (lasreader == nullptr) {
                    std::cerr << "\tERROR: could not open file: " << pointFile.filename << std::endl;
                    throw(std::string("\tERROR: could not open file: " + pointFile.filename));
                }
                LASheader header = lasreader->header;

//                if (this->check_bounds(header.min_x, header.max_x, header.min_y, header.max_y)) {
                //-- LAS classes to omit
                std::vector<int> lasomits;
                for (int i : pointFile.lasomits) {
                    lasomits.push_back(i);
                }

                //-- read each point 1-by-1
                uint32_t pointCount = header.number_of_point_records;

                //== IP: these are just info messages here==//
                //todo adapt thinning for this implementation
                std::clog << "\t(" << boost::locale::as::number << pointCount << " points in the file)\n";
                if ((pointFile.thinning > 1)) {
                    std::clog << "\t(skipping every " << pointFile.thinning << "th points, thus ";
                    std::clog << boost::locale::as::number << (pointCount / pointFile.thinning) << " are used)\n";
                }
                else
                    std::clog << "\t(all points used, no skipping)\n";

                if (!pointFile.lasomits.empty()) {
                    std::clog << "\t(omitting LAS classes: ";
                    for (int i : pointFile.lasomits)
                        std::clog << i << " ";
                    std::clog << ")\n";
                }
                //== IP: info messages up to here

                IO::print_progress_bar(0);
                int i = 0;
                //-- need if statement to make sure whether points go terrain or building
                //todo
                if (true) { //todo gotta include CSF somewhere here
                    while (lasreader->read_point()) {
                        LASpoint const& p = lasreader->point;
                        //-- set the thinning filter
                        if (i % pointFile.thinning == 0) {
                            //-- set the classification filter
                            if (std::find(lasomits.begin(), lasomits.end(), (int) p.classification) == lasomits.end()) {
                                //-- set the bounds filter
//                            if (this->check_bounds(p.X, p.X, p.Y, p.Y)) {
                                this->add_elevation_point(p, translate);
//                            }
                            }
                        }
                        if (i % (pointCount / 500) == 0)
                            IO::print_progress_bar(100 * (i / double(pointCount)));
                        i++;
                    }
                } else {

                }
                IO::print_progress_bar(100);
                std::clog << std::endl;
//                }
//                else {
//                    std::clog << "\tskipping file, bounds do not intersect polygon extent\n";
//                }
                _pointCloudTerrain.add_property_map<bool>("is_building_point", false);
                lasreader->close();
            }
            catch (std::exception& e) {
                lasreader->close();
                std::cerr << std::endl << e.what() << std::endl;
                throw (e.what());
            }
        }
    } else {
        std::cout << "Explicitly reading ground and/or building points" << std::endl;
        if (!Config::get().points_xyz.empty()) {
            std::cout << "Reading ground points" << std::endl;
            IO::read_point_cloud(Config::get().points_xyz, _pointCloudTerrain);
            _pointCloudTerrain.add_property_map<bool>("is_building_point", false);

            std::cout << "    Points read: " << _pointCloudTerrain.size() << std::endl;
        } else {
            std::cout << "INFO: Did not find any ground points! Will calculate ground as a flat surface." << std::endl;
            std::cout << "WARNING: Ground height of buildings can only be approximated. "
                      << "If you are using point cloud to reconstruct buildings, building height estimation can be wrong.\n"
                      << std::endl;
        }

        //-- Read building points
        if (!Config::get().buildings_xyz.empty()) {
            std::cout << "Reading building points" << std::endl;
            IO::read_point_cloud(Config::get().buildings_xyz, _pointCloudBuildings);
            if (_pointCloudBuildings.empty()) throw std::invalid_argument("Didn't find any building points!");

            std::cout << "    Points read: " << _pointCloudBuildings.size() << std::endl;
        }
    }
}

void PointCloud::add_elevation_point(const LASpoint& laspt, const CGAL::Aff_transformation_3<EPICK>& translate) {
    Point_set_3 transformPC;
    Point_3 pt(laspt.get_x(), laspt.get_y(), laspt.get_z());
    if (true) {//todo check if class belongs to the terrain
        _pointCloudTerrain.insert(pt.transform(translate));
    }
    if (true) {// todo check if the class belong to building
        _pointCloudBuildings.insert(pt.transform(translate));
    }
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