#include "Map3d.h"

Map3d::Map3d() = default;

Map3d::~Map3d() {
    this->clear_features();
}

void Map3d::reconstruct() {
    //-- Prepare features
    this->set_features();

    std::cout << "Features done" << std::endl;
    std::cout << "Num of features: " << _lsFeatures.size() << std::endl;

    //-- Define influence region, domain limits and boundaries
    this->set_boundaries();

    std::cout << "Bnds done" << std::endl;
    this->collect_garbage();
    std::cout << "Num of features: " << _lsFeatures.size() << std::endl;

    //-- Find polygon footprint elevation from point cloud
    this->set_footprint_elevation();

    std::cout << "Elevation done" << std::endl;

    //-- Reconstruct 3D features with respective algorithms
    this->threeDfy();

    std::cout << "3dfy done" << std::endl;

    //-- Add semantics

}

void Map3d::set_features() {
    //-- First feature is the terrain
    _terrain = new Terrain();

    //-- Add polygons as features
    //- Buildings
    for (auto& poly : _polygonsBuildings["features"]) {
        Building* building = new Building(poly);
        _lsFeatures.push_back(building);
    }
    //- Other polygons
    int count = 0;
    for (auto& semanticLayer : _polygonsSemanticLayers) {
        for (auto& poly : semanticLayer["features"]) {
            if (poly["geometry"]["type"] != "Polygon") continue; // Make sure only polygons are added

            SemanticPoly* semanticPoly = new SemanticPoly(poly, count);
            _lsFeatures.push_back(semanticPoly);
        }
        ++count;
    }

    //-- Boundary
    Sides* sides = new Sides(); Top* top = new Top();
    _boundaries.push_back(sides); _boundaries.push_back(top);

    //-- Group all features in one data structure
    _allFeatures.push_back(_terrain);
    for (auto& f : _lsFeatures) {
        _allFeatures.push_back(f);
    }
    for (auto& b : _boundaries) {
        _allFeatures.push_back(b);
    }
}

void Map3d::set_boundaries() {
    //-- Set the influence region --//
    //- Define radius of interest
    if (config::radiusOfInfluRegion == -infty) {
        std::cout << "--> Radius of interest not defined in config, calculating automatically" << std::endl;
        //-- Find building where the point of interest lies in and define radius of interest with BPG
        for (auto& f : _lsFeatures) {
            if (f->get_class() != BUILDING) continue;
            //TODO function that searches for the building where point lies
            // no building found - throw exception
        }
    }

    //-- Deactivate features that are out of the influence region
    for (auto& f : _lsFeatures) {
//        if (f->get_class() != BUILDING) continue;
        f->check_influ_region();
    }

    //-- Set the domain size --//
    // TODO: make it relative depending on round or rectangular domain
    if (config::dimOfDomain == -infty) {
        std::cout << "--> Domain size not defined in config, calculating automatically" << std::endl;
        //- Domain boundaries deferred until all building heights in the influ region are determined
    } else {
        //- Deactivate point cloud points that are out of bounds - static function of Boundary
        Boundary::set_bounds_to_pc(_pointCloud);
        Boundary::set_bounds_to_pc(_pointCloudBuildings);
//        _pointCloud.collect_garbage();

        //- Add flat buffer zone between the terrain and boundary
        Boundary::add_buffer(_pointCloud);
    }
}

void Map3d::set_footprint_elevation() {
    //-- Construct a searchTree to search for elevation
    SearchTree searchTree;
    searchTree.insert(_pointCloud.points().begin(), _pointCloud.points().end());

    for (auto& f : _lsFeatures) {
//        if (!f->is_active() || f->get_class() != BUILDING) continue; // For now only building footprints
        if (!f->is_active()) continue; // Checking for surface layers
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
        if (!f->is_active() || f->get_class() != BUILDING) continue;
        f->threeDfy(searchTree);
    }

    //-- Reconstruct boundaries
    for (auto& b : _boundaries) {
        b->threeDfy();
    }

    //-- Measure execution time
    auto endTime = std::chrono::steady_clock::now();
    auto diffTime = endTime - startTime;
    std::cout << "-> Calculations executed in " << std::chrono::duration<double> (diffTime).count() << " s" << std::endl;
}

bool Map3d::read_data() { // This will change with time
    //-- Read ground points
    IO::read_point_cloud(config::points_xyz, _pointCloud);
    //-- Read building points
    IO::read_point_cloud(config::buildings_xyz, _pointCloudBuildings);

    //-- Read building polygons
    IO::read_polygons(config::gisdata, _polygonsBuildings);
    //-- Read semantic polygons - will add it to a vector later
    _polygonsSemanticLayers.emplace_back();
    IO::read_polygons(config::topoSem, _polygonsSemanticLayers.back());

    return true;
}

void Map3d::output() {
    switch (config::outputFormat) {
        case OBJ:
            IO::output_obj(_allFeatures);
            break;
        case STL: // Only ASCII stl for now
            IO::output_stl(_allFeatures);
            break;
        case CityJSON:
            //-- Remove inactives and add ID's to features - obj and stl don't need id
            // just temp for now
            this->prep_feature_output();
            IO::output_cityjson(_allFeatures);
            break;
    }
}

void Map3d::prep_feature_output() { // Temp impl, might change
    for (unsigned long i = 0; i < _allFeatures.size();) {
        if (_allFeatures[i]->is_active()) {
            _allFeatures[i]->set_id(i++);
        }
        else {
            _allFeatures[i] = nullptr;
            _allFeatures.erase(_allFeatures.begin() + i);
        }
    }
};

void Map3d::collect_garbage() { // Just a test for now, could be helpful when other polygons get implemented
    for (unsigned long i = 0; i < _allFeatures.size();) {
        if (_allFeatures[i]->is_active()) ++i;
        else {
            _allFeatures[i] = nullptr;
            _allFeatures.erase(_allFeatures.begin() + i);
        }
    }
    for (unsigned long i = 0; i < _lsFeatures.size();) {
        if (_lsFeatures[i]->is_active()) ++i;
        else {
            delete _lsFeatures[i]; _lsFeatures[i] = nullptr;
            _lsFeatures.erase(_lsFeatures.begin() + i);
        }
    }
}

void Map3d::clear_features() {
    for (auto f : _allFeatures) {
        f = nullptr;
    }
    delete _terrain;
    _terrain = nullptr;
    for (auto& f : _lsFeatures) {
        delete f; f = nullptr;
    }
    for (auto& b : _boundaries) {
        delete b; b = nullptr;
    }
}