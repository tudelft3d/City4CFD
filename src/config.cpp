#include "config.h"

//-- Config values are hardcoded until I add config file reading
namespace config {
    //-- Input info
    const char* points_xyz = "data/input/ground_simplified.ply";
    const char* gisdata    = "data/input/tudcampus.geojson";
    const char* buildings_xyz    = "data/input/building_simplified.xyz";
    const char* topoSem = "data/input/Water.geojson";

    //-- Output info
    std::string outputFileName = "Mesh";

    //-- Domain dimensions
    Point_2 pointOfInterest = Point_2(85420, 446221);
    double radiusOfInfluRegion = 350.0;
    double dimOfDomain = 1000.0;
    double topHeight = 300.;

    //-- Reconstruction related
    double lod                = 1.2;
    double buildingPercentile = 0.9;

    //-- Output flags
    OutputFormat outputFormat = OBJ;
    bool outputSeparately = false;
}

bool config::read_config_file() {return true;}