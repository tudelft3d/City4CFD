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

    //-- Remove inactive features
    this->collect_garbage();
    std::cout << "Num of features: " << _lsFeatures.size() << std::endl;

    //-- Add PC points to DT
    this->triangulate_terrain();
    std::cout << "CDT terrain done" << std::endl;

    //-- Avoid having too long polygons
    this->polygon_processing();
    std::cout << "Checking edge length done" << std::endl;

    //-- Find polygon footprint elevation from point cloud
    this->set_footprint_elevation();
    std::cout << "Elevation done" << std::endl;

    //-- Reconstruct 3D features with respective algorithms
    this->threeDfy();
    std::cout << "3dfy done" << std::endl;
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
    int count = 4; //- Surface layer ID is 4 onwards
    for (auto& surfaceLayer : _polygonsSurfaceLayers) {
        for (auto& poly : surfaceLayer["features"]) {
            if (poly["geometry"]["type"] != "Polygon") continue; // Make sure only polygons are added

            SurfaceLayer* semanticPoly = new SurfaceLayer(poly, count);
            _lsFeatures.push_back(semanticPoly);
        }
        ++count;
    }

    //-- Boundary
    Sides* sides = new Sides(); Top* top = new Top();
    _boundaries.push_back(sides); _boundaries.push_back(top);
}

void Map3d::set_boundaries() {
    //-- Set the influence region --//
    //- Define radius of interest
    if (config::influenceRegionRadius == -infty) { // temp, will change the condition
        std::cout << "--> Radius of interest not defined in config, calculating automatically" << std::endl;
        //-- Find building where the point of interest lies in and define radius of interest with BPG
        for (auto& f : _lsFeatures) {
            if (f->get_class() != BUILDING) continue;
            //TODO function that searches for the building where point lies
            // no building found - throw exception
        }
    }

    //-- Deactivate features that are out of their scope
    //todo write it for CDT and PC respectively and go with the faster one
    for (auto& f : _lsFeatures) {
        f->check_feature_scope();
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

        //- Add flat buffer zone between the terrain and boundary
        Boundary::add_buffer(_pointCloud);
    }
}

void Map3d::triangulate_terrain() {
    _terrain->set_cdt(_pointCloud);
}

void Map3d::polygon_processing() {
    for (auto& f : _lsFeatures) {
        if (!f->is_active()) continue;
        std::cout << "New fieature!" << std::endl;
        for (auto& ring : f->get_poly().rings()) {
            geomtools::shorten_long_poly_edges(ring);
        }
    }
}

void Map3d::set_footprint_elevation() {
    //-- Construct a searchTree to search for elevation
//    SearchTree searchTree;
//    searchTree.insert(_pointCloud.points().begin(), _pointCloud.points().end());

    //-- Make a DT with inexact constructions for fast interpolation
    DT dt;
    dt.insert(_pointCloud.points().begin(), _pointCloud.points().end());

    for (auto& f : _lsFeatures) {
        if (!f->is_active()) continue;
        try {
#ifdef NDEBUG
            f->calc_footprint_elevation_nni(dt);
#else
            f->calc_footprint_elevation_linear(dt);  // NNI is quite slow in debug mode, better to use linear in that case
#endif
//            f->calc_footprint_elevation_from_pc(searchTree);
        } catch (std::exception& e) {
            std::cerr << std::endl << "Footprint elevation calculation failed for object \'" << f->get_id() << "\' (class " << f->get_class() << ") with error: " << e.what() << std::endl;
        }
    }
}

void Map3d::threeDfy() {
    //-- Measure execution time
    auto startTime = std::chrono::steady_clock::now();

    //-- Construct the terrain with surface layers
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
    //-- Read surface layer polygons
    for (auto& topoLayer : config::topoLayers) {
        _polygonsSurfaceLayers.emplace_back();
        IO::read_polygons(topoLayer, _polygonsSurfaceLayers.back());
    }

    return true;
}

void Map3d::output() {
    //-- Group all features for output
    this->prep_feature_output();

    switch (config::outputFormat) {
        case OBJ:
            IO::output_obj(_outputFeatures);
            break;
        case STL: // Only ASCII stl for now
            IO::output_stl(_outputFeatures);
            break;
        case CityJSON:
            //-- Remove inactives and add ID's to features - obj and stl don't need id
            // just temp for now
            this->prep_cityjson_output();
            IO::output_cityjson(_outputFeatures);
            break;
    }
}

void Map3d::prep_feature_output() {
    _outputFeatures.push_back(_terrain);
    for (auto& f : _lsFeatures) {
        if (!f->is_active() || f->get_class() != BUILDING) continue; // Surface layers are grouped in terrain
        _outputFeatures.push_back(f);
    }
    for (auto& b : _boundaries) {
        _outputFeatures.push_back(b);
    }
    for (auto& l : _terrain->get_surface_layers()) {
        _outputFeatures.push_back(l);
    }
}

void Map3d::prep_cityjson_output() { // Temp impl, might change
    for (unsigned long i = 0; i < _outputFeatures.size();) {
        if (_outputFeatures[i]->is_active()) {
            _outputFeatures[i]->set_id(i++);
        }
        else {
            _outputFeatures[i] = nullptr;
            _outputFeatures.erase(_outputFeatures.begin() + i);
        }
    }
};

void Map3d::collect_garbage() {
    for (unsigned long i = 0; i < _lsFeatures.size();) {
        if (_lsFeatures[i]->is_active()) ++i;
        else {
            delete _lsFeatures[i]; _lsFeatures[i] = nullptr;
            _lsFeatures.erase(_lsFeatures.begin() + i);
        }
    }
}

void Map3d::clear_features() {
    for (auto f : _outputFeatures) {
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