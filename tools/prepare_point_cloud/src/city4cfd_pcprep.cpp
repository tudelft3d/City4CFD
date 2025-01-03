/*
  City4CFD
 
  Copyright (c) 2021-2025, 3D Geoinformation Research Group, TU Delft

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

/*
 * Point cloud preparation tool:
 * Takes LAS/LAZ files and extracts ground and buildings
 * You can manually specify classes in case of classified point cloud
 * If the point cloud is not classified, the tool will try to extract the two automatically
 */

//todo: LAS/LAZ output?
//      Outlier and vegetation filter

#include <chrono>
#include <boost/locale.hpp>
#include <iomanip>

#include "lasreader.hpp"
#include "CSF/src/CSF.h"

#include "PCConfig.h"

using namespace PCPrep;

typedef std::array<double, 3> Point;

void printWelcome() {
    auto logo{
            R"(
     #==============================================================#
     #                        __                                    #
     #                   __  |''|                                   #
     #                  |""| |''|  _   /|__                         #
     #                __|""| |''|_| | | |""|/\_                     #
     #               |''|""| |''|'| __| |""|'''|  _____             #
     #          _ _  |''|""|^|''|'||""| |""|'''| |"""""|            #
     #         |"|"| |''|""|||''|'||""| |""|'''| |"""""|            #
     #     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     #
     #    ____   _   _                ___________________________   #
     #   / ___| /_\ | |_   _   _      __  ____/___  ____/___  __ \  #
     #  | |     |"| | __| | | | |     _  /     __  /_    __  / / /  #
     #  | |___  |"| | |_  | |_| |  4  / /___   _  __/    _  /_/ /   #
     #   \____| |"|  \__|  \__, |     \____/   /_/       /_____/    #
     #                     |___/                                    #
     #                                                              #
     #==============================================================#
)"
    };

    auto info{
            R"(#===============================================================================================================#
# Point cloud preparation tool                                                                                  #
# This tool takes LAS/LAZ point cloud tiles and exports building and ground points                              #
# In case the point cloud is not classified, it uses Cloth Simulation Filter (CSF) to perform ground filtering  #
#============================================================================================================== #

)"
    };
    std::cout << logo;
    std::cout << "City4CFD Copyright (C) 2021-2025 3D Geoinformation Research Group, TU Delft\n" << std::endl;
    std::cout << info;
}

void printHelp() {
    auto helpMsg{
            R"(USAGE:
    city4cfd_pcprep config.json
)"
    };
    std::cout << helpMsg;
}

void print_progress_bar(int percent) {
    std::string bar;
    for (int i = 0; i < 50; i++) {
        if (i < (percent / 2)) {
            bar.replace(i, 1, "=");
        }
        else if (i == (percent / 2)) {
            bar.replace(i, 1, ">");
        }
        else {
            bar.replace(i, 1, " ");
        }
    }
    std::clog << "\r" "    [" << bar << "] ";
    std::clog.width(3);
    std::clog << percent << "%     " << std::flush;
}

bool inBbox(const LASpoint& pt) {
    double ptx = pt.get_x();
    double pty = pt.get_y();
    if (ptx > Config::get().xmin && ptx < Config::get().xmax &&
        pty > Config::get().ymin && pty < Config::get().ymax) {
        return true;
    }
    return false;
}

void outputPts(const std::vector<Point>& pts, const std::string& filename) {
    std::cout << "Outputting " << filename << std::endl;

    std::ofstream of;
    of.open(filename.c_str());
    for (auto& pt : pts) {
        of << pt[0] << " " << pt[1] << " " << pt[2] << std::setprecision(16) << std::endl;
    }
    of.close();
}


int main(int argc, char** argv) {
    try {
        printWelcome();

        auto startTime = std::chrono::steady_clock::now();

        std::string config_path;
        if (argc >= 2) {
            config_path = fs::current_path().append(argv[1]).string();
        } else {
            printHelp();
            return EXIT_SUCCESS;
        }

        Config::get().read_config(config_path);

        //-- Point structures for building and ground?
        std::vector<Point> groundPts;
        std::vector<Point> buildingPts;

        int readPts = 0;
        for (auto& pointFile: Config::get().las_files) {
            std::cout << "Reading LAS/LAZ file: " << pointFile << "\n";

            LASreadOpener lasreadopener;
            lasreadopener.set_file_name(pointFile.c_str());
            //-- set to compute bounding box
            lasreadopener.set_populate_header(true);
            LASreader* lasreader = lasreadopener.open();

            //-- check if file is open
            if (lasreader == nullptr)
                throw(std::runtime_error(std::string("Error reading LAS/LAZ file.")));

            LASheader header = lasreader->header;

            //-- read each point 1-by-1
            uint64_t pointCount = header.number_of_point_records;

            std::vector<Point> pts;

            std::cout << "    " << boost::locale::as::number << pointCount << " points in the file\n";
            double percentLeft = 1 - (Config::get().pointCloudThinning / 100);
            if (Config::get().pointCloudThinning > 0) {
                std::cout << "    Skipping " << Config::get().pointCloudThinning << "% points" << std::endl;
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

                int i = -1;
                double currPercent = 1.01;
                print_progress_bar(0);
                while (lasreader->read_point()) {
                    if (++i % (pointCount / 200) == 0) {
                        print_progress_bar(100 * (i / double(pointCount)));
                    }

                    LASpoint const& p = lasreader->point;
                    if (p.return_number != p.number_of_returns) continue;
                    if (Config::get().checkBbox)
                        if (!inBbox(p)) continue;
                    //-- set the thinning filter
                    if (currPercent >= 1.) {
                        //-- add the point to the point cloud
                        bool isRead = false;
                        int classID = p.classification;
                        std::vector<int> classList = Config::get().las_classes_ground;
                        if (std::find(classList.begin(), classList.end(), classID)
                            != classList.end()) {
                            groundPts.push_back({p.get_x(), p.get_y(), p.get_z()});
                            isRead = true;
                        }
                        classList = Config::get().las_classes_building;
                        if (std::find(classList.begin(), classList.end(), classID)
                            != classList.end()) {
                            buildingPts.push_back({p.get_x(), p.get_y(), p.get_z()});
                            isRead = true;
                        }
                        if (isRead) ++readPts;
                        currPercent -= 1.;
                    }
                    currPercent += percentLeft;
                }
                print_progress_bar(100);
                std::clog << std::endl;

                lasreader->close();
            } else { // Cloth Simulation Filter
                std::cout << "    No point cloud classification information given. "
                          << "Will filter ground using the Cloth Simulation Filter" << std::endl;
                CSF csf;
                //-- Hardcoded CSF parameters
                csf.params.bSloopSmooth = true;
                csf.params.class_threshold = 0.5; //1
                csf.params.cloth_resolution = 2;  //3
                csf.params.interations = 500;
                csf.params.rigidness = 3;
                csf.params.time_step = 0.65;


                Point translatePt = {(header.max_x + header.min_x / 2),
                                     (header.max_y + header.min_y / 2),
                                     (header.max_z + header.min_z / 2)};

                //-- Read the data
                int i = -1;
                double currPercent = 0.;
                print_progress_bar(0);
                while (lasreader->read_point()) {
                    if (++i % (pointCount / 200) == 0) {
                        print_progress_bar(100 * (i / double(pointCount)));
                    }

                    LASpoint const& p = lasreader->point;
                    if (p.return_number != p.number_of_returns) continue;
                    if (Config::get().checkBbox)
                        if (!inBbox(p)) continue;
                    //-- set the thinning filter
                    if (currPercent >= 1.) {
                        //-- Add directly to CSF's data structure
                        csf.addPoint(p.get_x() - translatePt[0],
                                     p.get_y() - translatePt[1],
                                     p.get_z() - translatePt[2]);

                        ++readPts;
                        currPercent -= 1.;
                    }
                    currPercent += percentLeft;
                }
                print_progress_bar(100);
                std::clog << std::endl;

                lasreader->close();

                //-- Perform CS filtering
                std::cout << "    Applying CSF to the point cloud" << std::endl;
                std::vector<int> groundIndices, offGroundIndices;
                csf.do_filtering(groundIndices, offGroundIndices, false);

                //-- Add points to point clouds
                for (auto idx : groundIndices) {
                    groundPts.push_back({csf.getPointCloud().at(idx).x + translatePt[0],
                                         csf.getPointCloud().at(idx).z + translatePt[1],
                                        -csf.getPointCloud().at(idx).y + translatePt[2]});
                }
                for (auto idx : offGroundIndices) {
                    buildingPts.push_back({csf.getPointCloud().at(idx).x + translatePt[0],
                                           csf.getPointCloud().at(idx).z + translatePt[1],
                                          -csf.getPointCloud().at(idx).y + translatePt[2]});
                }
            }
        }
        std::clog << "\nPoints read: " << readPts << "\n" << std::endl;

        //-- Output points to respective files
        //- Output ground
        outputPts(groundPts, "ground.txt");

        //- Output buildings
        outputPts(buildingPts, "buildings.txt");

        auto endTime = std::chrono::steady_clock::now();
        auto diffTime = endTime - startTime;
        std::cout << "\nProgram executed in " << std::chrono::duration<double>(diffTime).count() << " s" << std::endl;
        std:: cout << "End" << std::endl;

        return EXIT_SUCCESS;
    } catch (std::exception& e) {

        std::cerr << "\nProgram failed! Reason: " << e.what() << std::endl;
        std::cout << "End" << std::endl;
        return EXIT_FAILURE;
    }
}
