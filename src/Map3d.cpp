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

    //-- Different flow if explicitly defining domain boundary or leaving it to BPG
    if (!_bndBPG) {
        //-- Set outer boundary
        this->set_bnd();
        std::cout << "Bnds done" << std::endl;

        //-- Avoid having too long polygons
        this->shorten_polygons(_lsFeatures);

        //-- Find footprint elevation of all polygons using smoothed DT
        this->set_footprint_elevation(_lsFeatures);

        //-- Reconstruct buildings in the influ region
        this->reconstruct_buildings();
    } else {
        //-- First the buildings are reconstructed
        //- Prepare polygons for buildings
        this->shorten_polygons(_buildings);

        //- Find building footprint elevation using smoothed DT
        this->set_footprint_elevation(_buildings);

        //- Reconstruct buildings in the influ region
        this->reconstruct_buildings();

        //-- Second, find the highest building and set domain according to BPG
        this->set_bnd();

        //-- Add surface layers now that the domain size is known
        this->shorten_polygons(_surfaceLayers);
        this->set_footprint_elevation(_surfaceLayers);
    }

    //-- Constrain features, generate terrain mesh from CDT
    std::cout << "---Reconstructing terrain ---" << std::endl;
    this->reconstruct_terrain();

    //-- Generate side and top boundaries
    this->reconstruct_boundaries();
    std::cout << "--- Reconstruction done ---" << std::endl;
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

    //-- Boundary
    for (int i = 0; i < config::numSides; ++i)
        _boundaries.push_back(std::make_shared<Sides>(TopoFeature::get_num_output_layers()));
    _boundaries.push_back(std::make_shared<Top>(TopoFeature::get_num_output_layers()));

    //- Other polygons
    for (auto& surfaceLayer : _polygonsSurfaceLayers) {
        int outputLayerID = TopoFeature::get_num_output_layers();
        config::surfaceLayerIDs.push_back(outputLayerID); // Need it for later
        for (auto& poly : surfaceLayer) {
            auto surfacePoly = std::make_shared<SurfaceLayer>(*poly, outputLayerID);
            _lsFeatures.push_back(surfacePoly);
            _surfaceLayers.push_back(surfacePoly);
        }
    }

    //-- BPG flags for influ region and domain boundary
    if (config::influRegionConfig.type() == typeid(bool)) _influRegionBPG = true;
    if (config::domainBndConfig.type() == typeid(bool))   _bndBPG = true;

    //-- Make a DT with inexact constructions for fast interpolation
    _dt.insert(_pointCloud.points().begin(), _pointCloud.points().end());
#ifdef SMOOTH
    geomtools::smooth_dt<DT, EPICK>(_pointCloud, _dt);
#endif
}

void Map3d::set_influ_region() {
    //-- Set the influence region --//
    if (_influRegionBPG) { // Automatically calculate influ region with BPG
        std::cout << "--- INFO: Influence region not defined in config. "
                  << "Calculating with BPG. ---" << std::endl;
        _influRegion.calc_influ_region_bpg(_dt, _pointCloudBuildings, _buildings);
    } else { // Define influ region either with radius or predefined polygon
        boost::apply_visitor(_influRegion, config::influRegionConfig);
    }

    //-- Deactivate buildings that are out of influ region
    for (auto& f : _buildings) {
        f->check_feature_scope(_influRegion.get_bounding_region());
    }
    this->clear_inactives();
}

void Map3d::set_bnd() {
    double hMax = 0;
    if (_bndBPG) { // Automatically calculate boundary with BPG
        std::cout << "--- INFO: Domain boundaries not defined in config. "
                  << "Calculating with BPG. ---" << std::endl;

        //-- Find the highest building in the influ region
        for (auto& b: _buildings) {
            if (b->get_height() > hMax) hMax = b->get_height();
        }
        //-- Calculate the boundary polygon according to BPG and defined boundary type
        _domainBnd.calc_bnd_bpg(hMax, _influRegion.get_bounding_region(), _buildings);
    } else {
        //-- Define boundary region with values set in config
        boost::apply_visitor(_domainBnd, config::domainBndConfig);
    }
    this->bnd_sanity_check(); // Check if outer bnd is larger than the influ region

    //-- Prepare the outer boundary polygon for sides and top, and polygon for feature scope
    Polygon_2 bndPoly, pcBndPoly, startBufferPoly; // Depends on the buffer region
    bndPoly = _domainBnd.get_bounding_region();
    if (_boundaries.size() > 2) {
        geomtools::shorten_long_poly_edges(bndPoly, hMax); // Outer poly edge size is hMax ATM
        Boundary::set_bnd_poly(bndPoly, pcBndPoly, startBufferPoly);
    } else
        Boundary::set_bnd_poly(bndPoly, pcBndPoly, startBufferPoly);

    //-- Deactivate point cloud points that are out of bounds
    Boundary::set_bounds_to_terrain(_pointCloud, bndPoly,
                                    pcBndPoly, startBufferPoly);
    Boundary::set_bounds_to_pc(_pointCloudBuildings, startBufferPoly);

    //-- Check feature scope for surface layers now that the full domain is known
    for (auto& f: _surfaceLayers) {
        f->check_feature_scope(pcBndPoly);
    }
    this->clear_inactives();

}

void Map3d::bnd_sanity_check() {
    auto& domainBndPoly = _domainBnd.get_bounding_region();
    for (auto& pt : _influRegion.get_bounding_region()) {
        if (!geomtools::point_in_poly(pt, domainBndPoly))
            throw std::domain_error("The influence region is larger than the domain boundary!");
    }
}

void Map3d::reconstruct_terrain() {
    _terrain->set_cdt(_pointCloud);
    _terrain->constrain_features(_lsFeatures);
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
    if (_boundaries.size() > 2) { // Means more than one side
        for (auto i = 0; i < _boundaries.size() - 1; ++i) {
            //-- Each boundary object is one side of the boundary
            _boundaries[i]->prep_output(_domainBnd.get_bounding_region().edge(i).to_vector());
        }
    } else {
        _boundaries.front()->prep_output();
    }

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
#ifndef NDEBUG
    assert(config::outputSurfaces.size() == TopoFeature::get_num_output_layers());
#endif
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
    for (unsigned long i = 0; i < _outputFeatures.size(); ++i) {
        if (_outputFeatures[i]->is_active()) {
            _outputFeatures[i]->set_id(i++);
        }
        else {
            _outputFeatures.erase(_outputFeatures.begin() + i);
        }
    }
};

void Map3d::clear_inactives() {
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