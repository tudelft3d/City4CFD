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

#ifndef CITY4CFD_PCCONFIG_H
#define CITY4CFD_PCCONFIG_H

#include <iostream>
#include <fstream>
#include <string>

#include "boost/filesystem.hpp"

#include "valijson/adapters/nlohmann_json_adapter.hpp"
#include "valijson/schema.hpp"
#include "valijson/schema_parser.hpp"
#include "valijson/validator.hpp"

#include "PCconfigSchema.inc"

namespace fs = boost::filesystem;

namespace PCPrep {

class Config {
public:
    static Config& get() {
        static Config inst;
        return inst;
    }

protected:
    Config() = default;
    Config(const Config&) = delete;
    Config& operator=(const Config&) = delete;
    ~Config() = default;

public:
    //-- Input 
    std::vector<std::string> las_files;                  // LAS point cloud
    std::vector<int>         las_classes_ground;         // LAS classes used by terrain
    std::vector<int>         las_classes_building;       // LAS classes used by buildings
    double                   xmin                 = 0.;
    double                   xmax                 = 0.;
    double                   ymin                 = 0.;
    double                   ymax                 = 0.;
    double                   pointCloudThinning   = 0.;

    //-- Output
    fs::path                  workDir;
    fs::path                  outputDir           = fs::current_path();

    //-- Other settings
    bool checkBbox                                = false;

    //== Public functions
    void read_config(std::string& config_path) {
        std::ifstream json_file(config_path);
        if (!json_file)
            throw std::invalid_argument(std::string("Configuration file " + config_path + " not found."));

        //-- Filepaths in the json file are relative to the location of the json file
        this->workDir = fs::path(config_path).parent_path();
        fs::current_path(this->workDir);
        std::cout << "Working directory: " << fs::canonical(this->workDir) << std::endl;

        nlohmann::json j;
        try {
            j = nlohmann::json::parse(json_file, nullptr, true, true);
        } catch (std::exception& e) {
            std::cerr << "ERROR: Configuration file is not a valid JSON file. "
                      << "Check if you're missing a comma or a brace." << std::endl;
            throw;
        }

        std::cout << "\nValidating JSON configuration file...";
        this->validate(j);
        std::cout <<"Configuration file is valid! \n" << std::endl;

        this->set_config(j);
    }

    void validate(nlohmann::json& j) {
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
    };

    void set_config(nlohmann::json& j) {
        if (j.contains("las_files")) {
            for (std::string filename: j["las_files"]) {
                las_files.push_back(filename);
            }
        }
        if (j.contains("ground_classes")) {
            for (int classID: j["ground_classes"]) {
                las_classes_ground.push_back(classID);
            }
        }
        if (j.contains("building_classes")) {
            for (int classID: j["building_classes"]) {
                las_classes_building.push_back(classID);
            }
        }
        if (j.contains("bbox")) {
            xmin = j["bbox"][0];
            ymin = j["bbox"][1];
            xmax = j["bbox"][2];
            ymax = j["bbox"][3];
            checkBbox = true;
        }
        if (j.contains("thinning")) {
            pointCloudThinning = j["thinning"];
        }
    };
};

}

#endif //CITY4CFD_PCCONFIG_H
