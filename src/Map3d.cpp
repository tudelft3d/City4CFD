#include "Map3d.h"

Map3d::~Map3d() {
    this->clear_features();
}

void Map3d::reconstruct() {
    //-- Prepare features
    this->set_features();

    //-- Define influence region, domain limits and boundaries
    this->set_boundaries();

    //-- Find polygon footprint elevation from point cloud
    this->set_footprint_elevation();

    //-- Reconstruct 3D features with respective algorithms
    this->threeDfy();

    //-- Add semantics

    //-- Add boundary
    this->add_boundary();
}

void Map3d::set_features() {
    //-- First feature is the terrain
    _terrain = new Terrain();

    //-- Add polygons as features -- ONLY BUILDING and BAG for now
    int pid = 1; // for now
    for (auto poly : _polygons["features"]) {
        Building* building = new Building(poly, pid++);
        _lsFeatures.push_back(building);
    }

    //-- Boundary
    _boundary = new Boundary();
}

void Map3d::set_boundaries() {
    //-- Set the influence region --//
    //- Define radius of interest
    double radiusOfInfluRegion;
    if (_configData->radiusOfInterest == -infty) {
        std::cout << "--> Radius of interest not defined in config, calculating automatically" << std::endl;
        //-- Find building where the point of interest lies in and define radius of interest with BPG
        for (auto& f : _lsFeatures) {
            if (f->get_class() != BUILDING) continue;
            //TODO function that searches for the building where point lies
            // no building found - throw exception
        }
    } else radiusOfInfluRegion = _configData->radiusOfInterest;

    //-- Deactivate buildings that are out of the influence region
    for (auto& f : _lsFeatures) {
        if (f->get_class() != BUILDING) continue;
        f->check_influ_region(_configData->pointOfInterest, radiusOfInfluRegion);
    }

    //-- Set the domain size --//
    // TODO: make it relative depending on round or rectangular domain
    double dimOfDomain;
    if (_configData->dimOfDomain == -infty) {
        std::cout << "--> Domain size not defined in config, calculating automatically" << std::endl;
        //- Domain boundaries deferred until all building heights in the influ region are determined
    } else {
        //-- Deactivate point cloud points that are out of bounds
        _boundary->set_bounds_to_pc(_pointCloud);
        _boundary->set_bounds_to_pc(_pointCloudBuildings);
//        _pointCloud.collect_garbage();
        _boundary->add_buffer(_pointCloud);
    }
}

void Map3d::set_footprint_elevation() {
    //-- Construct a searchTree to search for elevation
    SearchTree searchTree;
    searchTree.insert(_pointCloud.points().begin(), _pointCloud.points().end());

    for (auto f : _lsFeatures) {
        if (f->is_active() && f->get_class() != BUILDING) continue; // For now only building footprints
        try {
            f->calc_footprint_elevation(searchTree);
        } catch (std::exception& e) {
            std::cerr << std::endl << "Footprint elevation calculation failed for object \'" << f->get_id() << "\' (class " << f->get_class() << ") with error: " << e.what() << std::endl;
        }
    }
}

void Map3d::threeDfy() {
    //-- Measure execution time
    auto startTime = std::chrono::steady_clock::now();

    //-- CDT the terrain
    _terrain->threeDfy(_pointCloud, _lsFeatures);

    //-- Reconstruct buildings
    SearchTree searchTree;
    searchTree.insert(_pointCloudBuildings.points().begin(), _pointCloudBuildings.points().end());

    for (auto& f : _lsFeatures) {
        if (f->is_active() && f->get_class() != BUILDING) continue;
        f->threeDfy(searchTree);
    }

    //-- Measure execution time
    auto endTime = std::chrono::steady_clock::now();
    auto diffTime = endTime - startTime;
    std::cout << "-> Calculations executed in " << std::chrono::duration<double> (diffTime).count() << " s" << std::endl;
}

void Map3d::add_boundary() {
    //-- Add sides

    //-- Add top

}

//-- Input functions are temporary, will bemoved to IO
bool Map3d::read_config(const char *points_xyz) { //TODO still needs implementing
    _configData = new ConfigData();
    return true;
}

bool Map3d::read_point_cloud(const char* points_xyz) {
    std::ifstream ifile(points_xyz, std::ios_base::binary);
    ifile >> _pointCloud;
    std::cerr << "POINT CLOUD: "<< _pointCloud.size() << " point read" << std::endl;
    return true;
}

bool Map3d::read_point_cloud_buildings(const char* points_xyz) {
    std::ifstream ifile(points_xyz, std::ios_base::binary);
    ifile >> _pointCloudBuildings;
    std::cerr << "POINT CLOUD BUILDINGS: "<< _pointCloudBuildings.size() << " point read" << std::endl;
    return true;
}

bool Map3d::read_polygons(const char* gisdata) {
    std::ifstream ifs(gisdata);
    _polygons = json::parse(ifs);
    return true;
}

void Map3d::output() {
    bool output_separately = false;
    // TODO: config file to output as single model or separately

    //-- Create output file
    std::ofstream of;
    of.open("mesh.obj");
    std::string fs, bs;

    //-- Output terrain
    _terrain->output_feature(fs, bs, _dPts);

    //-- Output buildings
    bs += "\ng Building"; // temp
    for (auto &f : _lsFeatures) {
        if (!f->is_active()) continue;
        f->output_feature(fs, bs, _dPts);
    }

    //-- Output side and top boundary
    _boundary->get_mesh();
    _boundary->output_feature(fs, bs, _dPts);

    of << fs << bs;
    of.close();
}

void Map3d::clear_features() {
    delete _configData;
    delete _terrain;
    delete _boundary;
    _configData = nullptr;
    _terrain    = nullptr;
    _boundary   = nullptr;
    for (auto f : _lsFeatures) {
        delete f;
        f = nullptr;
    }
}
