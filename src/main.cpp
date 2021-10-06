#include <iostream>
#include "nlohmann/json.hpp"
#include "definitions.h"
#include "io.h"
#include "Map3d.h"
#include <chrono>

int main() {
    try {
        auto startTime = std::chrono::steady_clock::now();

        //-- Data input - this needs to be sorted
        const char* config = "";

        //-- Read configuration file TODO
        IO::read_config(config);

        //-- Create the main class
        Map3d map3d;

        //-- Read point cloud and polygons, and store them in Map3d
        map3d.read_data();

        //-- Calculate elevations and triangulate
        map3d.reconstruct();

        //-- Output data
        map3d.output();

        auto endTime = std::chrono::steady_clock::now();
        auto diffTime = endTime - startTime;
        std::cout << "-> Program executed in " << std::chrono::duration<double>(diffTime).count() << " s" << std::endl;

        return EXIT_SUCCESS;
    } catch (std::exception& e) {
        std::cerr << "--> Program failed! Reason: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
}