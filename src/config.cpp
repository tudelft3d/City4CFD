#include "config.h"

#include "valijson/adapters/nlohmann_json_adapter.hpp"
#include "valijson/schema.hpp"
#include "valijson/schema_parser.hpp"
#include "valijson/validator.hpp"

#include "io.h"

#include "configSchema.inc"

namespace config {
    //-- Input info
    std::string              points_xyz;         // Ground
    std::string              buildings_xyz;      // Buildings
    std::string              gisdata;            // Building Polygons
    std::vector<std::string> topoLayers = {};    // Other polygons
    std::string              importedBuildings;  // Additional pre-reconstructed buildings

    //-- Domain setup
    Point_2     pointOfInterest;
    double      topHeight = 0;
    //- Influ region and domain bnd
    boost::variant<bool, double, Polygon_2> influRegionConfig;
    boost::variant<bool, double, Polygon_2> domainBndConfig;
    DomainType            bpgDomainType;
    bool                  bpgBlockageRatioFlag = false;
    double                bpgBlockageRatio = 0.03;
    Vector_2              flowDirection(1, 0);
    std::vector<double>   bpgDomainSize;
    double                domainBuffer = -g_largnum;

    //-- Reconstruction
    //- Terrain
    double    terrainSimplification = 0.;
    bool      smoothTerrain = false;
    //- Buildings
    std::string lod;
    double      buildingPercentile;
    //- Imported buildings
    bool        importAdvantage;
    bool        importTrueHeight;
    std::string importLoD;

    //-- Polygons related
    double                edgeMaxLen;
    std::map<int, double> averageSurfaces;

    //-- Output
    fs::path                  workDir;
    fs::path                  outputDir = fs::current_path();
    std::string               outputFileName;
    OutputFormat              outputFormat;
    bool                      outputSeparately = false;
    std::vector<std::string>  outputSurfaces = {"Terrain", "Buildings"};
    int                       numSides = 1;
    std::vector<int>          surfaceLayerIDs;

    //-- Data log
    bool               outputLog = false;
    std::string        logName("log");
    std::ostringstream log;
    std::ostringstream logSummary;
    std::vector<int>   failedBuildings;
}

void config::validate(nlohmann::json& j) {
    using valijson::Schema;
    using valijson::SchemaParser;
    using valijson::Validator;
    using valijson::adapters::NlohmannJsonAdapter;
    using valijson::ValidationResults;

    //- Load schema
    Schema configSchema;
    SchemaParser parser;
    NlohmannJsonAdapter configSchemaAdapter(jsonschema::schema);
    // debug
//    nlohmann::json schema = nlohmann::json::parse(std::ifstream("../../data/input/schema.json"), nullptr, true, true);
//    NlohmannJsonAdapter configSchemaAdapter(schema);

    parser.populateSchema(configSchemaAdapter, configSchema);

    //- Validate with schema
    Validator validator;
    NlohmannJsonAdapter configAdapter(j);
    ValidationResults results;
    if (!validator.validate(configSchema, configAdapter, &results)) {
        std::stringstream err_oss;
        err_oss << "Validation failed." << std::endl;
        ValidationResults::Error error;
        int error_num = 1;
        while (results.popError(error))
        {
            std::string context;
            std::vector<std::string>::iterator itr = error.context.begin();
            for (; itr != error.context.end(); itr++)
                context += *itr;

            err_oss << "Error #" << error_num << std::endl
                    << "  context: " << context << std::endl
                    << "  desc:    " << error.description << std::endl;
            ++error_num;
        }
        throw std::runtime_error(err_oss.str());
    }
}

void config::set_config(nlohmann::json& j) {
    //-- Schema validation
    //-- Path to point cloud(s)
    if (j.contains("point_clouds")) {
        if (j["point_clouds"].contains("ground")) points_xyz = j["point_clouds"]["ground"];
        if (j["point_clouds"].contains("buildings")) buildings_xyz = j["point_clouds"]["buildings"];
    }

    //-- Path to polygons
    int i = 0;
    int surfLayerIdx = 2; // 0 - terrain, 1 - buildings, surface layers star from 2
    for (auto& poly : j["polygons"]) {
        if (poly["type"] == "Building") {
            gisdata = poly["path"];
        }
        if (poly["type"] == "SurfaceLayer") {
            topoLayers.push_back(poly["path"]);
            if (poly.contains("layer_name")) {
                outputSurfaces.push_back(poly["layer_name"]);
            } else {
                outputSurfaces.push_back("SurfaceLayer" + std::to_string(++i));
            }
            if (poly.contains("average_surface")) {
                if (poly["average_surface"]) {
                    averageSurfaces[surfLayerIdx] = poly["surface_percentile"];
                }
            }
            ++surfLayerIdx;
        }
    }

    //-- Additional geometries
    if (j.contains("import_geometries"))
        importedBuildings = j["import_geometries"]["path"];

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
        if (j.contains("flow_direction"))
            flowDirection = Vector_2(j["flow_direction"][0], j["flow_direction"][1]);

        std::string bpgDomainConfig = j["bnd_type_bpg"];
        if (boost::iequals(bpgDomainConfig, "round")) {
            bpgDomainType = ROUND;
            if (j.contains("bpg_domain_size")) {
                bpgDomainSize = {j["bpg_domain_size"][0], j["bpg_domain_size"][1]};
            } else bpgDomainSize = {15, 5}; // BPG
        } else if (boost::iequals(bpgDomainConfig, "rectangle")) {
            bpgDomainType = RECTANGLE;
        } else if (boost::iequals(bpgDomainConfig, "oval")) {
            bpgDomainType = OVAL;
        }
    }
    // Sort out few specifics for domain types
    if (bpgDomainType == RECTANGLE || bpgDomainType == OVAL) {
        if (j.contains("bpg_domain_size")) {
            bpgDomainSize = {j["bpg_domain_size"][0], j["bpg_domain_size"][1], j["bpg_domain_size"][2], j["bpg_domain_size"][3]};
        } else bpgDomainSize = {5, 5, 15, 5}; // BPG
    }

    // Set domain side and top
    if (domainBndConfig.type() == typeid(Polygon_2)) {
        Polygon_2 poly = boost::get<Polygon_2>(domainBndConfig);
        numSides = poly.size();
        for (int i = 0; i < poly.size(); ++i) {
            outputSurfaces.emplace_back("Side_" + std::to_string(i % poly.size()));
        }
    } else if (domainBndConfig.type() == typeid(double) || bpgDomainType != RECTANGLE) {
        outputSurfaces.emplace_back("Sides");
    } else if (bpgDomainType == RECTANGLE) { // Expand output surfaces with front and back
        numSides = 4;
        outputSurfaces.emplace_back("Front");
        outputSurfaces.emplace_back("Side_2");
        outputSurfaces.emplace_back("Back");
        outputSurfaces.emplace_back("Side_1");
    }
    outputSurfaces.emplace_back("Top");

    // Blockage ratio
    if (j.contains("bpg_blockage_ratio")) {
        if (j["bpg_blockage_ratio"].is_boolean()) {
            bpgBlockageRatioFlag = j["bpg_blockage_ratio"];
        } else if (j["bpg_blockage_ratio"].is_number()) {
            bpgBlockageRatioFlag = true;
            bpgBlockageRatio = (double)j["bpg_blockage_ratio"] / 100;
        }
    }

    // Top height
    if (j.contains("top_height"))
        topHeight = j["top_height"];

    // Buffer region
    if (j.contains("buffer_region"))
        domainBuffer = j["buffer_region"];

    //-- Reconstruction
    // Terrain
    if (j.contains("terrain_simplification"))
        terrainSimplification = j["terrain_simplification"];
    if (j.contains("smooth_terrain"))
        smoothTerrain = j["smooth_terrain"];

    // Buildings
    lod = j["lod"].front();
    buildingPercentile = (double)j["building_percentile"].front() / 100.;

    // Imported buildings
    if (j.contains("import_geometries")) {
        importAdvantage  = j["import_geometries"]["advantage"];
        importTrueHeight = j["import_geometries"]["true_height"];
        importLoD = j["import_geometries"]["lod"];
    }

    //-- Polygons related
    edgeMaxLen = j["edge_max_len"];

    //-- Output
    //- File name
    outputFileName = j["output_file_name"];

    //- Output format
    std::string outputFormatConfig = j["output_format"];
    if (boost::iequals(outputFormatConfig, "obj")) {
        outputFormat = OBJ;
    } else if (boost::iequals(outputFormatConfig, "stl")) {
        outputFormat = STL;
    } else if (boost::iequals(outputFormatConfig, "cityjson")) {
        outputFormat = CityJSON;
    } else throw std::invalid_argument(std::string("'" + outputFormatConfig + "'" + " is unsupported file format!"));

    outputSeparately = j["output_separately"];

    //-- Data log
    if (j.contains("output_log")) {
        outputLog = true;
        if (j.contains("log_file")) logName = j["log_file"];
    }
    log << "// ========================================= CITYCFD LOG ========================================= //" << std::endl;
    logSummary <<"\n// =========================================== SUMMARY =========================================== //" << std::endl;
}

//-- influRegion and domainBndConfig flow control
void config::set_region(boost::variant<bool, double, Polygon_2>& regionType,
                        std::string& regionName,
                        nlohmann::json& j) {
    if (j[regionName].is_string()) { // Search for GeoJSON polygon
        std::string polyFilePath = (std::string)j[regionName];
        if(!fs::exists(polyFilePath)) {
            throw std::invalid_argument(std::string("Cannot find polygon file '" +
                                        polyFilePath + "' for " + regionName));
        }
        //-- Read poly
        Polygon_2 tempPoly;
        JsonVector influJsonPoly;
        IO::read_geojson_polygons(polyFilePath, influJsonPoly);
        for (auto& coords : influJsonPoly.front()->front()) { // I know it should be only 1 polygon with 1 ring
            tempPoly.push_back(Point_2(coords[0], coords[1]));
        }
        //-- Prepare poly
        CGAL::internal::pop_back_if_equal_to_front(tempPoly);
        if (tempPoly.is_clockwise_oriented()) tempPoly.reverse_orientation();

        regionType = tempPoly;
    } else if (j[regionName].size() > 2) { // Explicitly defined region polygon with points
        Polygon_2 tempPoly;
        for (auto& pt : j[regionName]) tempPoly.push_back(Point_2(pt[0], pt[1]));
        regionType = tempPoly;
    } else if (j[regionName].is_number() || j[regionName].is_array() && j[regionName][0].is_number()) { // Influ region radius
        regionType = (double)j[regionName].front();
    } else regionType = true; // Leave it to BPG
}