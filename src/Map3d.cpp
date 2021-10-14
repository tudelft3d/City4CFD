#include "Map3d.h"

Map3d::Map3d() = default;
Map3d::~Map3d() = default;

void Map3d::reconstruct() {
    //-- Prepare features
    this->set_features();
    std::cout << "Features done" << std::endl;
    std::cout << "Num of features: " << _lsFeatures.size() << std::endl;

    //-- Define influence region
    this->set_influ_region();

    //-- Store boundary polygon or prepare for BPG calculation
    this->set_bnd_calc();
    std::cout << "Bnds done" << std::endl;

    //-- Different flow if explicitly defining domain boundary or leaving it to BPG
    if (!_bndBPG) {
        //-- Avoid having too long polygons
        this->shorten_polygons(_lsFeatures);
        std::cout << "Checking edge length done" << std::endl;

        //-- Find polygon footprint elevation from point cloud
        this->set_footprint_elevation(_lsFeatures);
        std::cout << "Elevation done" << std::endl;

        //-- Constrain features, generate terrain mesh from it
        //-- todo combine terrain handling in one function
        this->triangulate_terrain();
        this->constrain_features(_lsFeatures);
        this->generate_terrain_mesh();
        std::cout << "Terrain mesh done" << std::endl;

        //-- Reconstruct buildings and boundary
        this->reconstruct_buildings();
        this->reconstruct_boundaries();
        std::cout << "3dfy done" << std::endl;
    } else {
        // todo tbd
        //-- Add PC points to DT
        this->triangulate_terrain();

        //-- Avoid having too long polygons
        this->shorten_polygons(_buildings);
        std::cout << "Checking edge length done" << std::endl;

        //-- Find polygon footprint elevation from point cloud
        this->set_footprint_elevation(_buildings);
        std::cout << "Elevation done" << std::endl;

        //-- Reconstruct 3D features with respective algorithms
        std::cout << "3dfy done" << std::endl;
    }
}

void Map3d::set_features() {
    //-- First feature is the terrain
    _terrain = std::make_shared<Terrain>();

    //-- Add polygons as features
    //- Buildings
    for (auto& poly : _polygonsBuildings) {
        auto building = std::make_shared<Building>(*poly);
        _lsFeatures.push_back(building);
        _buildings.push_back(building);
    }
    //- Other polygons
    int count = 4; //- Surface layer ID is 4 onwards
    for (auto& surfaceLayer : _polygonsSurfaceLayers) {
        for (auto& poly : surfaceLayer) {
            auto surfacePoly = std::make_shared<SurfaceLayer>(*poly, count);
            _lsFeatures.push_back(surfacePoly);
            _surfaceLayers.push_back(surfacePoly);
        }
        ++count;
    }

    //-- Boundary
    auto sides = std::make_shared<Sides>(); auto top = std::make_shared<Top>();
    _boundaries.push_back(sides); _boundaries.push_back(top);

    //-- Boundary calculation with BPG flag
    if (config::domainBndConfig.type() == typeid(bool)) _bndBPG = true;

    //-- Make a DT with inexact constructions for fast interpolation
    _dt.insert(_pointCloud.points().begin(), _pointCloud.points().end());
#ifdef SMOOTH
    geomtools::smooth_dt<DT, EPICK>(_pointCloud, _dt);
#endif
}

void Map3d::triangulate_terrain() {
    _terrain->set_cdt(_pointCloud);
}

void Map3d::set_influ_region() {
    //-- Set the influence region --//
    if (config::influRegionConfig.type() == typeid(bool)) { // Automatically calculate influ region with BPG
        std::cout << "--- INFO: Influence region not defined in config. "
                  << "Calculating with BPG. ---" << std::endl;
        _influRegion(_pointCloud, _pointCloudBuildings, _buildings);
    } else { // Define influ region either with radius or predefined polygon
        boost::apply_visitor(_influRegion, config::influRegionConfig);
    }

    //-- Deactivate buildings that are out of influ region
    for (auto& f : _buildings) {
        f->check_feature_scope(_influRegion.get_bounded_region());
    }
    this->collect_garbage();
}

void Map3d::set_bnd_calc() {
    if (_bndBPG) { // Automatically calculate boundary with BPG
        //-- BND calc deferred until buildings are reconstructed
        std::cout << "--- INFO: Domain boundaries not defined in config. "
                  << "Calculating with BPG. ---" << std::endl;
    } else { // Define boundary region either with radius or predefined polygon

        boost::apply_visitor(_domainBnd, config::domainBndConfig);
        this->bnd_sanity_check();

        //- Deactivate point cloud points that are out of bounds - static function of Boundary
        Boundary::set_bounds_to_pc(_pointCloud, _domainBnd.get_bounded_region());
        Boundary::set_bounds_to_pc(_pointCloudBuildings, _domainBnd.get_bounded_region());;

        _boundaries[0]->set_bnd_poly(_domainBnd.get_bounded_region(), _pointCloud);

        //-- Check feature scope for surface layers now that the full domain is known
        for (auto& f: _surfaceLayers) {
            f->check_feature_scope(_domainBnd.get_bounded_region());
        }
        this->collect_garbage();
    }
}

void Map3d::set_outer_bnd_bpg() {
    //todo
}

void Map3d::bnd_sanity_check() {
    auto& domainBndPoly = _domainBnd.get_bounded_region();
    for (auto& pt : _influRegion.get_bounded_region()) {
        if (!geomtools::point_in_poly(pt, domainBndPoly))
            throw std::domain_error("The influence region is larger than the domain boundary!");
    }
}

void Map3d::generate_terrain_mesh() {
    _terrain->create_mesh(_lsFeatures);
}

void Map3d::reconstruct_buildings() {
    SearchTree searchTree;
    searchTree.insert(_pointCloudBuildings.points().begin(), _pointCloudBuildings.points().end());

//    double failed = 0;
    for (auto& f : _buildings) {
        if (!f->is_active()) continue;
        try {
            f->reconstruct(searchTree);
        } catch (std::exception& e) {
//            std::cout << "Failed so far: " << ++failed << "; " << e.what() << std::endl;
            // Add to warning log when individual buildings don't reconstruct
        }
    }
}

void Map3d::reconstruct_boundaries() {
    for (auto& b : _boundaries) {
        b->reconstruct();
    }
}

void Map3d::read_data() { // This will change with time
    //-- Read ground points
    IO::read_point_cloud(config::points_xyz, _pointCloud);
    if (_pointCloud.size() == 0) {
        std::cout << "Didn't find any ground points! Calculating ground as flat surface" << std::endl;
    }

    //-- Read building points
    IO::read_point_cloud(config::buildings_xyz, _pointCloudBuildings);
    if (_pointCloudBuildings.empty()) throw std::invalid_argument("Didn't find any building points!");

    //-- Read building polygons
    IO::read_geojson_polygons(config::gisdata, _polygonsBuildings);
    if (_polygonsBuildings.empty()) throw std::invalid_argument("Didn't find any building polygons!");

    //-- Read surface layer polygons
    for (auto& topoLayer: config::topoLayers) {
        _polygonsSurfaceLayers.emplace_back();
        IO::read_geojson_polygons(topoLayer, _polygonsSurfaceLayers.back());
    }
}

void Map3d::output() {
    fs::current_path(config::outputDir);

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
    for (auto& f : _buildings) {
        if (!f->is_active()) continue;
        _outputFeatures.push_back(f);
    }
    for (auto& b : _boundaries) {
        _outputFeatures.push_back(b);
    }
    for (auto& l : _terrain->get_surface_layers()) { // Surface layers are grouped in terrain
        _outputFeatures.push_back(l);
    }
}

void Map3d::prep_cityjson_output() { // Temp impl, might change
    for (unsigned long i = 0; i < _outputFeatures.size();) {
        if (_outputFeatures[i]->is_active()) {
            _outputFeatures[i]->set_id(i++);
        }
        else {
            _outputFeatures.erase(_outputFeatures.begin() + i);
        }
    }
};

void Map3d::collect_garbage() {
    for (unsigned long i = 0; i < _buildings.size();) {
        if (_buildings[i]->is_active()) ++i;
        else {
            _buildings.erase(_buildings.begin() + i);
        }
    }
    for (unsigned long i = 0; i < _surfaceLayers.size();) {
        if (_surfaceLayers[i]->is_active()) ++i;
        else {
            _surfaceLayers.erase(_surfaceLayers.begin() + i);
        }
    }
    for (unsigned long i = 0; i < _lsFeatures.size();) {
        if (_lsFeatures[i]->is_active()) ++i;
        else {
            _lsFeatures.erase(_lsFeatures.begin() + i);
        }
    }
}

//-- Templated functions
template<typename T>
void Map3d::shorten_polygons(T& feature) {
    for (auto& f : feature) {
        if (!f->is_active()) continue;
        for (auto& ring : f->get_poly().rings()) {
            geomtools::shorten_long_poly_edges(ring);
        }
    }
}
//- Explicit template instantiation
template void Map3d::shorten_polygons<Buildings>    (Buildings& feature);
template void Map3d::shorten_polygons<SurfaceLayers>(SurfaceLayers& feature);
template void Map3d::shorten_polygons<PolyFeatures> (PolyFeatures& feature);

template<typename T>
void Map3d::set_footprint_elevation(T& features) {
    for (auto& f : features) {
        if (!f->is_active()) continue;
#ifdef NDEBUG
        f->calc_footprint_elevation_nni(_dt);
#else
        f->calc_footprint_elevation_linear(_dt);  // NNI is quite slow in debug mode, better to use linear in that case
#endif
    }
}
//- Explicit template instantiation
template void Map3d::set_footprint_elevation<Buildings>    (Buildings& feature);
template void Map3d::set_footprint_elevation<SurfaceLayers>(SurfaceLayers& feature);
template void Map3d::set_footprint_elevation<PolyFeatures> (PolyFeatures& feature);

template<typename T>
void Map3d::constrain_features(const T& features) {
    _terrain->constrain_features(features);
}
//-- Explicit template instantiation
template void Map3d::constrain_features<Buildings>    (const Buildings& feature);
template void Map3d::constrain_features<SurfaceLayers>(const SurfaceLayers& feature);
template void Map3d::constrain_features<PolyFeatures> (const PolyFeatures& feature);