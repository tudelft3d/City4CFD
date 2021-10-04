#include "config.h"

//-- Config values are hardcoded until I add config file reading
namespace config {
    //-- Input info
//    const char*              points_xyz    = "data/input/sampled_ground_simpl.ply";
    const char*              points_xyz    = "data/input/ground_simplified.ply";
    const char*              gisdata       = "data/input/tudcampus.geojson";
    const char*              buildings_xyz = "data/input/building_simplified.xyz";
    std::vector<const char*> topoLayers    = {"data/input/Vegetation.geojson", "data/input/Water.geojson"};
//    std::vector<const char*> topoLayers    = {};

    //-- Output info
    std::string outputFileName = "Mesh";

    //-- Domain dimensions
    Point_2 pointOfInterest = Point_2(85420., 446221.);
    double radiusOfInfluRegion = 350.0;
    double dimOfDomain = 1000.0;
    double topHeight = 300.;

    //-- Reconstruction related
    double lod                = 1.2;
    double buildingPercentile = 0.9;

    //-- Polygons related
    double edgeMaxLen = 3;
    std::vector<bool> avgHeights{false, false, false, true};

    //-- Output flags
    OutputFormat outputFormat = OBJ;
    bool outputSeparately = true;
}

bool config::read_config_file() {return true;}