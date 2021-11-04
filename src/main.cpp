#include "config.h"
#include "io.h"
#include "Map3d.h"

void printWelcome() {
    std::cout << "cityCFD wooho" << std::endl;
}

void printHelp() {
    std::cout << "Help goes here" << std::endl;
}

int main(int argc, char** argv) {
    try {
        printWelcome();

        auto startTime = std::chrono::steady_clock::now();

        //-- TEMP
        const char* CONFIG_PATH = "data/input/config.json";
        std::string config_path = CONFIG_PATH;

        //-- Path to config.json file
        if (argc >= 2) {
            config_path = fs::current_path().append(argv[1]).string();
            std::cout << config_path << std::endl;
        } else {
//            throw std::invalid_argument("Missing path to config file!");
        }

        //-- TODO optional arguments like output dir and log file name, I don't know yet
        for (auto i = 1; i < argc; ++i) {
            if (boost::iequals(argv[i], "--help")) {
                printHelp();
                return EXIT_SUCCESS;
            } else if (boost::iequals(argv[i], "--output_dir") && (i + 1) != argc) {
                config::outputDir = fs::absolute(fs::current_path().append(argv[i+1]));
            } else if (boost::iequals(argv[i], "--output_file") && (i + 1) != argc) {
                config::outputFileName = argv[i+1];
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

        auto endTime = std::chrono::steady_clock::now();
        auto diffTime = endTime - startTime;
        std::cout << "-> Program executed in " << std::chrono::duration<double>(diffTime).count() << " s" << std::endl;

        return EXIT_SUCCESS;
    } catch (std::exception& e) {
        std::cerr << "--> Program failed! Reason: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
}