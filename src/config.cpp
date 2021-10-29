#include "config.h"

namespace config {
    //-- Input info
    std::string              points_xyz;      // Ground
    std::string              buildings_xyz;   // Buildings
    std::string              gisdata;         // Building Polygons
    std::vector<std::string> topoLayers = {}; // Other polygons

    //-- Domain dimensions
    Point_2     pointOfInterest;
    double      topHeight = 0;
    boost::variant<bool, double, std::string, Polygon_2> influRegionConfig;
    boost::variant<bool, double, std::string, Polygon_2> domainBndConfig;
    DomainType            bpgDomainType;
    Vector_2              flowDirection;
    std::vector<double>   bpgDomainSize;
    double                domainBuffer = -g_largnum;

    //-- Reconstruction related
    double lod;
    double buildingPercentile;

    //-- Polygons related
    double edgeMaxLen;

    //-- Output
    fs::path                  outputDir = fs::current_path();
    std::string               outputFileName;
    OutputFormat              outputFormat;
    bool                      outputSeparately = false;
    std::vector<std::string>  outputSurfaces = {"Terrain", "Buildings", "Top"};
    int                       numSides = 1;
    std::vector<int>          surfaceLayerIDs;
}

void config::set_config(nlohmann::json& j) {
    //todo make it with JSON schema validator
    //-- Path to point cloud(s)
    points_xyz    = j["point_clouds"]["ground"];
    buildings_xyz = j["point_clouds"]["buildings"];

    //-- Path to polygons
    bool foundBuildingPoly = false;
    int i = 0;
    for (auto& poly : j["polygons"]) {
        if (poly["type"] == "Building") {
            gisdata = poly["path"];
            foundBuildingPoly = true;
        }
        if (poly["type"] == "SurfaceLayer") {
            topoLayers.push_back(poly["path"]);
            if (poly.contains("layer_name"))
                outputSurfaces.push_back(poly["layer_name"]);
            else
                outputSurfaces.push_back("SurfaceLayer" + std::to_string(++i));
        }
    }
    if (!foundBuildingPoly) throw std::invalid_argument("Didn't find a path to building polygons in configuration file!");

    //-- Domain setup
    pointOfInterest = Point_2(j["point_of_interest"][0], j["point_of_interest"][1]);

    //- Influence region
    std::string influRegionCursor = "influence_region";
    config::set_region(influRegionConfig, influRegionCursor, j);

    //- Domain boundaries
    std::string domainBndCursor = "domain_bnd";
    config::set_region(domainBndConfig, domainBndCursor, j);
    // Define domain type if using BPG
    if (domainBndConfig.type() == typeid(bool)) {
        if (!j.contains("flow_direction")) throw std::invalid_argument("Missing information on flow direction!");
        if (j["flow_direction"].size() != 2) throw std::invalid_argument("Flow direction array size is not of size 3!");
        flowDirection = Vector_2(j["flow_direction"][0], j["flow_direction"][1]);

        std::string bpgDomainConfig = j["bnd_type_bpg"];
        if (boost::iequals(bpgDomainConfig, "round")) {
            bpgDomainType = ROUND;
            if (j.contains("bpg_domain_size")) {
                if (j["bpg_domain_size"].size() == 2)
                    bpgDomainSize = {j["bpg_domain_size"][0], j["bpg_domain_size"][1]};
                else
                    throw std::invalid_argument("Wrong input for overriding BPG domain size");
            } else bpgDomainSize = {15, 5};
        } else if (boost::iequals(bpgDomainConfig, "rectangle")) {
            bpgDomainType = RECTANGLE;
        } else if (boost::iequals(bpgDomainConfig, "oval")) {
            bpgDomainType = OVAL;
        } else throw std::invalid_argument(std::string("'" + bpgDomainConfig+ "'" + " is unknown domain type!"));
    }
    // Sort out few specifics for domain types
    if (bpgDomainType == RECTANGLE || bpgDomainType == OVAL) {
        if (j.contains("bpg_domain_size")) {
            if (j["bpg_domain_size"].size() == 4)
                bpgDomainSize = {j["bpg_domain_size"][0], j["bpg_domain_size"][1], j["bpg_domain_size"][2], j["bpg_domain_size"][3]};
            else
                throw std::invalid_argument("Wrong input for overriding BPG domain size");
        } else bpgDomainSize = {5, 5, 15, 5};
    }

    // Set domain side
        if (domainBndConfig.type() == typeid(double) || bpgDomainType != RECTANGLE) {
            outputSurfaces.insert(outputSurfaces.begin() + 2, "Sides");
        } else if (bpgDomainType == RECTANGLE) { // Expand output surfaces with front and back
            numSides = 4;
            outputSurfaces.insert(outputSurfaces.begin() + 2, "Front");
            outputSurfaces.insert(outputSurfaces.begin() + 2, "Side_2");
            outputSurfaces.insert(outputSurfaces.begin() + 2, "Back");
            outputSurfaces.insert(outputSurfaces.begin() + 2, "Side_1");
        } else if (domainBndConfig.type() == typeid(Polygon_2)) {
            Polygon_2 poly = boost::get<Polygon_2>(domainBndConfig);
            numSides = poly.size();
            for (int i = 0; i < poly.size(); ++i) {
                outputSurfaces.insert(outputSurfaces.begin() + 2,
                                      std::string("Side_" + std::to_string(i % poly.size())));
            }
        } // If it is a JSON poly it has to be deferred until the polygon is parsed

    // Buffer region
    if (j.contains("buffer_region"))
        if (j["buffer_region"].is_number() && j["buffer_region"] > -99)
            domainBuffer = j["buffer_region"];

    // Top height
    if (j.contains("top_height") && j["top_height"].is_number())
        topHeight = j["top_height"].front();
    else if (domainBndConfig.type() != typeid(bool))
        throw std::invalid_argument("Invalid 'top_height' argument!");

    //-- Reconstruction related
    lod = j["lod"].front();
    buildingPercentile = (double)j["building_percentile"].front() / 100.;

    //-- Polygons related
    if (j.contains("edge_max_len"))
        edgeMaxLen = j["edge_max_len"].front();
    else throw std::invalid_argument("Missing 'edge_max_len' argument!");

    //-- Output
    //- File name
    if (outputFileName.empty() && j["output_file_name"].is_string())
        outputFileName = j["output_file_name"];
    else
        throw std::invalid_argument("Invalid output file name!");

    //- Output format
    std::string outputFormatConfig = j["output_format"];
    if (boost::iequals(outputFormatConfig, "obj")) {
        outputFormat = OBJ;
    } else if (boost::iequals(outputFormatConfig, "stl")) {
        outputFormat = STL;
    } else if (boost::iequals(outputFormatConfig, "cityjson")) {
        outputFormat = CityJSON;
    } else throw std::invalid_argument(std::string("'" + outputFormatConfig + "'" + " is unsupported file format!"));

    if (j.contains("output_separately") && j["output_separately"].is_boolean())
        outputSeparately = j["output_separately"];
}

//-- influRegion and domainBndConfig flow control
void config::set_region(boost::variant<bool, double, std::string, Polygon_2>& regionType,
                        std::string& regionName,
                        nlohmann::json& j) {
    if (j[regionName].is_string()) { //- Search for GeoJSON polygon
        regionType = (std::string)j[regionName];
    } else if (j[regionName].size() > 2) { // Explicitly defined region polygon with points
        Polygon_2 tempPoly;
        for (auto& pt : j[regionName]) tempPoly.push_back(Point_2(pt[0], pt[1]));
        regionType = tempPoly;
    } else if (j[regionName].size() == 2) {
        throw std::invalid_argument("Unknown setup of the region: " + regionName + "!");
    } else if (j[regionName].is_number() || j[regionName].is_array() && j[regionName][0].is_number()) { // Influ region radius
        regionType = (double)j[regionName].front();
    } else regionType = true; // Leave it to BPG
}