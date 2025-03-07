/*
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

#include "Config.h"
#include "io.h"
#include "Map3d.h"

#include  <boost/algorithm/string/predicate.hpp>

std::string CITY4CFD_VERSION = "0.6.3";

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

    std::cout << logo;
    std::cout << "City4CFD Copyright (C) 2021-2025 3D Geoinformation Research Group, TU Delft\n" << std::endl;
}

void printHelp() {
    auto helpMsg{
R"(USAGE:
    City4CFD config_file.json OPTIONS

AVAILABLE OPTIONS:
    --help            Prints out this help message
    --version         Displays City4CFD version information
    --output_dir      Sets the directory where output files are stored
    --output_file     Overrides output file(s) name from the configuration file
)"
    };

    std::cout << helpMsg;
}

void printVersion() {
    std::cout << "Version: " << CITY4CFD_VERSION << std::endl;
}

int main(int argc, char** argv) {
    try {
        printWelcome();

        auto startTime = std::chrono::steady_clock::now();

        std::string config_path;
        //-- Path to config.json file
        if (argc >= 2) {
            config_path = fs::absolute(argv[1]).string();
        } else {
            printHelp();
            return EXIT_SUCCESS;
        }

        //-- Input arguments
        for (auto i = 1; i < argc; ++i) {
            if (boost::iequals(argv[i], "--help")) {
                printHelp();
                return EXIT_SUCCESS;
            } else if (boost::iequals(argv[i], "--version")) {
                printVersion();
                return EXIT_SUCCESS;
            } else if (boost::iequals(argv[i], "--output_dir")) {
                if (i + 1 == argc) throw std::invalid_argument("Missing argument for --output_dir");

                Config::get().outputDir = fs::absolute(argv[++i]);
                if (!fs::exists(Config::get().outputDir)) throw std::invalid_argument(std::string("Output directory '" + Config::get().outputDir.string() + "' does not exist!"));
            } else if (boost::iequals(argv[i], "--output_file")) {
                if (i + 1 == argc) throw std::invalid_argument("Missing argument for --output_file");

                Config::get().outputFileName = argv[++i];
            } else {
                if (i > 1) throw std::invalid_argument(std::string("Unknown option " + std::string(argv[i])));
            }
        }

        //-- Read configuration file
        IO::read_config(config_path);

        //-- Create the main class
        Map3d map3d;

        //-- Read point cloud and polygons, and store them in Map3d
        map3d.read_data();

        //-- Calculate elevations and triangulate
        map3d.reconstruct();

        //-- Output data
        map3d.output();

        //-- Output log
        IO::output_log(map3d.get_failed_buildings());

        auto endTime = std::chrono::steady_clock::now();
        auto diffTime = endTime - startTime;
        std::cout << "\nProgram executed in " << std::chrono::duration<double>(diffTime).count() << " s" << std::endl;
        std:: cout << "End" << std::endl;

        return EXIT_SUCCESS;
    } catch (std::exception& e) {
        //-- Output log
        IO::output_log();

        std::cerr << "\nProgram failed! Reason: " << e.what() << std::endl;
        std::cout << "End" << std::endl;
        return EXIT_FAILURE;
    }
}
