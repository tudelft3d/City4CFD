#include <iostream>
#include "nlohmann/json.hpp"
#include "definitions.h"
#include "io.h"
#include "Map3d.h"
#include <chrono>

int main() {
    auto startTime = std::chrono::steady_clock::now();

    //-- Data input - this needs to be sorted
    const char* points_xyz = "data/input/ground_simplified.ply";
    const char* gisdata    = "data/input/tudcampus.geojson";
    const char* buildings_xyz    = "data/input/building_simplified.xyz";
    const char* config = "";

    //-- Read configuration file TODO
    if(!IO::read_config(config)){
        return 1;
    };

    //-- Create the main class
    Map3d map3d;

    //-- Read point cloud and polygons, and store them in Map3d
    if(!map3d.read_data()){
        return 1;
    };

    //-- Calculate elevations and triangulate
    map3d.reconstruct();

    //-- Output data
    map3d.output();

    auto endTime = std::chrono::steady_clock::now();
    auto diffTime = endTime - startTime;
    std::cout << "-> Program executed in " << std::chrono::duration<double> (diffTime).count() << " s" << std::endl;

    return 0;
}