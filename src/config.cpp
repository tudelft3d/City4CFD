#include "config.h"
#include "configSchema.inc"

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
    bool                  bpgBlockageRatioFlag = false;
    double                bpgBlockageRatio = 0.03;
    Vector_2              flowDirection;
    std::vector<double>   bpgDomainSize;
    double                domainBuffer = -g_largnum;

    //-- Reconstruction related
    std::string lod;
    double      buildingPercentile;

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
    points_xyz    = j["point_clouds"]["ground"];
    buildings_xyz = j["point_clouds"]["buildings"];

    //-- Path to polygons
    int i = 0;
    for (auto& poly : j["polygons"]) {
        if (poly["type"] == "Building") {
            gisdata = poly["path"];
        }
        if (poly["type"] == "SurfaceLayer") {
            topoLayers.push_back(poly["path"]);
            if (poly.contains("layer_name"))
                outputSurfaces.push_back(poly["layer_name"]);
            else
                outputSurfaces.push_back("SurfaceLayer" + std::to_string(++i));
        }
    }

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

    // Blockage ratio
    if (j.contains("bpg_blockage_ratio")) {
        if (j["bpg_blockage_ratio"].is_boolean()) {
            bpgBlockageRatioFlag = j["bpg_blockage_ratio"];
        } else if (j["bpg_blockage_ratio"].is_number()) {
            bpgBlockageRatioFlag = true;
            bpgBlockageRatio = j["bpg_blockage_ratio"];
        }
    }

    // Buffer region
    if (j.contains("buffer_region"))
        domainBuffer = j["buffer_region"];

    // Top height
    if (j.contains("top_height"))
        topHeight = j["top_height"];

    //-- Reconstruction related
    lod = j["lod"].front();
    buildingPercentile = (double)j["building_percentile"].front() / 100.;

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
    } else if (j[regionName].is_number() || j[regionName].is_array() && j[regionName][0].is_number()) { // Influ region radius
        regionType = (double)j[regionName].front();
    } else regionType = true; // Leave it to BPG
}