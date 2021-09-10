#include "config.h"

//-- Config values are hardcoded until I add config file reading
namespace config {
    //-- Domain dimensions
    Point_2 pointOfInterest = Point_2(85420, 446221);
    double radiusOfInfluRegion = 350.0;
    double dimOfDomain = 1000.0;
    double topHeight = 300.;

    //-- Reconstruction related
    double buildingPercentile = 0.9;

    //-- Output flags
    OutputFormat outputFormat = OBJ;
    bool outputSeparately = false;
}

bool config::read_config_file() {return true;}