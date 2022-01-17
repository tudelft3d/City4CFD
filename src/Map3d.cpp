#include "Map3d.h"

#include "geomutils.h"
#include "io.h"
#include "Terrain.h"
#include "ReconstructedBuilding.h"
#include "ImportedBuilding.h"
#include "SurfaceLayer.h"
#include "Sides.h"
#include "Top.h"

Map3d::Map3d() = default;
Map3d::~Map3d() = default;

void Map3d::reconstruct() {
    //-- Prepare features
    this->set_features();

    //-- Define influence region
    this->set_influ_region();

    //-- Different flow if explicitly defining domain boundary or leaving it to BPG
    if (!_bndBPG) {
        //-- Set outer boundary
        this->set_bnd();

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
    this->reconstruct_terrain();

    //-- Generate side and top boundaries
    this->reconstruct_boundaries();
}

void Map3d::set_features() {
    //-- First feature is the terrain
    _terrain = std::make_shared<Terrain>();

    //-- Add features - order in _lsFeatures defines the advantage in marking terrain polygons
    //- Buildings
    int internalID = 0;
    for (auto& poly : _polygonsBuildings) {
        auto building = std::make_shared<ReconstructedBuilding>(*poly, internalID++);
        _reconstructedBuildings.push_back(building);
        _buildings.push_back(building);
        _lsFeatures.push_back(building);
    }
    //- Imported buildings
    if (!_importedBuildingsJson.empty())  std::cout << "Importing CityJSON geometries" << std::endl;
    std::vector<std::shared_ptr<ImportedBuilding>> appendingBuildings;
    internalID = 0;
    for (auto& importedBuilding : _importedBuildingsJson) {
        auto explicitGeom = std::make_shared<ImportedBuilding>(*importedBuilding, _importedBuildingsPts, internalID++);
        if (!explicitGeom->is_appending()) {
            _importedBuildings.push_back(explicitGeom);
            _buildings.push_back(explicitGeom);
            _lsFeatures.push_back(explicitGeom);
        } else {
            appendingBuildings.push_back(explicitGeom);
        }
    }
    //- Check for building parts that do not have footprint and append to another instance of the same building
    for (auto& b : appendingBuildings) {
        for (auto& importedBuilding :  _importedBuildings) {
            if (b->get_parent_building_id() == importedBuilding->get_parent_building_id()) {
                importedBuilding->append_nonground_part(b);
                break;
            }
        }
    }
    if (!_importedBuildings.empty()) {
        this->clear_inactives();
        std::cout << "    Geometries imported: " << _importedBuildings.size() << std::endl;
    }

    //- Other polygons
    for (auto& surfaceLayer : _polygonsSurfaceLayers) {
        int outputLayerID = TopoFeature::get_num_output_layers();
        config::surfaceLayerIDs.push_back(outputLayerID); // Need it for later
        for (auto& poly : surfaceLayer) {
            auto surfacePoly = std::make_shared<SurfaceLayer>(*poly, outputLayerID);
            _surfaceLayers.push_back(surfacePoly);
            _lsFeatures.push_back(surfacePoly);
        }
    }
    std::cout << "    Polygons read: " << _lsFeatures.size() << std::endl;

    //-- Boundaries
    for (int i = 0; i < config::numSides; ++i)
        _boundaries.push_back(std::make_shared<Sides>(TopoFeature::get_num_output_layers()));
    _boundaries.push_back(std::make_shared<Top>(TopoFeature::get_num_output_layers()));

    //-- Simplify terrain points
    if (config::terrainSimplification > 0 + g_smallnum) {
        std::cout <<"\nSimplyfing terrain points" << std::endl;
        _pointCloud.remove(CGAL::random_simplify_point_set(_pointCloud, config::terrainSimplification), _pointCloud.end());
        _pointCloud.collect_garbage();
        std::cout << "    Terrain points after simplification: " << _pointCloud.size() << std::endl;
    }

    //-- BPG flags for influ region and domain boundary
    if (config::influRegionConfig.type() == typeid(bool)) _influRegionBPG = true;
    if (config::domainBndConfig.type() == typeid(bool))   _bndBPG = true;

    //-- Make a DT with inexact constructions for fast interpolation
    _dt.insert(_pointCloud.points().begin(), _pointCloud.points().end());
    if (config::smoothTerrain) {
        geomutils::smooth_dt<DT, EPICK>(_pointCloud, _dt);
    }
}

void Map3d::set_influ_region() {
    std::cout << "\nDefining influence region" << std::endl;
    //-- Set the influence region --//
    if (_influRegionBPG) { // Automatically calculate influ region with BPG
        std::cout << "\nINFO: Influence region not defined in config. "
                  << "Calculating with BPG." << std::endl;

        //-- Check if imported and reconstructed buildings are overlapping
        if (!_importedBuildings.empty()) this->solve_building_conflicts();

        //-- Prepare search tree in case of reconstruction
        std::shared_ptr<SearchTree> searchTree;
        searchTree = std::make_shared<SearchTree>(_pointCloudBuildings.points().begin(), _pointCloudBuildings.points().end());
        for (auto& building : _reconstructedBuildings) building->set_search_tree(searchTree);

        //-- Calculate influ region
        _influRegion.calc_influ_region_bpg(_dt, _pointCloudBuildings, _buildings);
    } else { // Define influ region either with radius or predefined polygon
        boost::apply_visitor(_influRegion, config::influRegionConfig);
    }

    //-- Deactivate buildings that are out of influ region
    for (auto& f : _buildings) {
        f->check_feature_scope(_influRegion.get_bounding_region());
    }
    this->clear_inactives();

    //-- Check if imported and reconstructed buildings are overlapping
    if (!_importedBuildings.empty()) this->solve_building_conflicts();

    std::cout << "    Number of building geometries in the influence region: " << _buildings.size() << std::endl;
}

void Map3d::set_bnd() {
    if (_bndBPG) { // Automatically calculate boundary with BPG
        std::cout << "\nINFO: Domain boundaries not defined in config. "
                  << "Calculating with BPG." << std::endl;

        //-- Calculate the boundary polygon according to BPG and defined boundary type
        _domainBnd.calc_bnd_bpg(_influRegion.get_bounding_region(), _buildings);
    } else {
        //-- Define boundary region with values set in config
        boost::apply_visitor(_domainBnd, config::domainBndConfig);
    }
    this->bnd_sanity_check(); // Check if outer bnd is larger than the influ region

    //-- Prepare the outer boundary polygon for sides and top, and polygon for feature scope
    Polygon_2 bndPoly, pcBndPoly, startBufferPoly; // Depends on the buffer region
    bndPoly = _domainBnd.get_bounding_region();
    if (_boundaries.size() > 2) { 
        geomutils::shorten_long_poly_edges(bndPoly, 20 * config::edgeMaxLen); // Outer poly edge size is hardcoded atm
        Boundary::set_bnd_poly(bndPoly, pcBndPoly, startBufferPoly);
    } else
        // If it's only one side bnd, edge length is already okay
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
        if (!geomutils::point_in_poly(pt, domainBndPoly))
            throw std::domain_error("The influence region is larger than the domain boundary!");
    }
}

void Map3d::reconstruct_terrain() {
    std::cout << "\nReconstructing terrain" << std::endl;

    _terrain->prep_constraints(_lsFeatures, _pointCloud);
    if (!config::averageSurfaces.empty()) this->average_polygon_points();
    _terrain->set_cdt(_pointCloud);
    _terrain->constrain_features();

    std::cout << "\n    Creating terrain mesh" << std::endl;
    _terrain->create_mesh(_lsFeatures);
}

void Map3d::reconstruct_buildings() {
    std::cout << "\nReconstructing buildings" << std::endl;
    if (!_importedBuildings.empty()) {
        std::cout << "    Will try to reconstruct imported buildings in LoD: " << config::importLoD
                  << ". If I cannot find a geometry with that LoD, I will reconstruct in the highest LoD available"
                  << std::endl;
    }

    std::shared_ptr<SearchTree> searchTree;
    searchTree = std::make_shared<SearchTree>(_pointCloudBuildings.points().begin(), _pointCloudBuildings.points().end());
    for (auto& building : _reconstructedBuildings) building->set_search_tree(searchTree);

    int failed = 0;
    for (auto& f : _buildings) {
        if (!f->is_active()) continue;
        try {
            f->reconstruct();
        } catch (std::exception& e) {
            ++failed;
            //-- Add information to log file
            config::log << "Failed to reconstruct building ID: " << f->get_internal_id()
                        << " Reason: " << e.what() << std::endl;
            //-- Get JSON file ID for failed reconstructions output
            config::failedBuildings.push_back(f->get_internal_id());
        }
    }
    config::logSummary << "BUILDING RECONSTRUCTION SUMMARY: TOTAL FAILED RECONSTRUCTIONS: "
                       << failed << std::endl;

    this->clear_inactives();
}

void Map3d::reconstruct_boundaries() {
    std::cout << "\nReconstructing boundaries" << std::endl;
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

void Map3d::average_polygon_points() {
    std::cout << "\n    Averaging surfaces" << std::endl;
    std::map<int, Point_3> averagedPts;

    //-- Construct a connectivity map and remove duplicates along the way
    auto is_building_pt = _pointCloud.property_map<bool>("is_building_point").first;
    std::unordered_map<Point_3, int> pointCloudConnectivity;
    auto it = _pointCloud.points().begin();
    int count = 0;
    while (it != _pointCloud.points().end()) {
        auto itPC = pointCloudConnectivity.find(*it);
        if (itPC != pointCloudConnectivity.end()) {
            _pointCloud.remove(_pointCloud.begin() + count);
        } else {
            pointCloudConnectivity[*it] = count;
            ++it;
            ++count;
        }
    }
    _pointCloud.collect_garbage();

    //-- Construct search tree from ground points
    SearchTree searchTree(_pointCloud.points().begin(), _pointCloud.points().end());

    //-- Perform averaging
    for (auto& f : _lsFeatures) {
        auto it = config::averageSurfaces.find(f->get_output_layer_id());
        if (it != config::averageSurfaces.end()) {
            f->average_polygon_inner_points(_pointCloud, averagedPts, searchTree, pointCloudConnectivity);
        }
    }

    //-- Change points with averaged values
    int pcOrigSize = _pointCloud.points().size();
    for (auto& it : averagedPts) {
        _pointCloud.insert(it.second);
    }
    for (int i = 0; i < pcOrigSize; ++i) {
        auto it = averagedPts.find(i);
        if (it != averagedPts.end()) {
            _pointCloud.remove(i);
            averagedPts.erase(i);
        }
    }
    _pointCloud.collect_garbage();
}

void Map3d::solve_building_conflicts() {
    for (auto& importedBuilding : _importedBuildings) {
        for (auto& reconstructedBuilding : _reconstructedBuildings) {
            if (geomutils::polygons_in_contact(importedBuilding->get_poly(), reconstructedBuilding->get_poly())) {
                if (config::importAdvantage) {
                    reconstructedBuilding->deactivate();
                } else {
                    importedBuilding->deactivate();
                }
            }
        }
    }
    this->clear_inactives();

   // to check if conflicts are solved
//    for (auto& b : _importedBuildings) b->deactivate();
//    this->clear_inactives();
}

void Map3d::read_data() { // This will change with time
    //-- Read ground points
    if (!config::points_xyz.empty()) {
    std::cout << "Reading ground points" << std::endl;
    IO::read_point_cloud(config::points_xyz, _pointCloud);
    _pointCloud.add_property_map<bool> ("is_building_point", false);
    std::cout << "    Points read: " << _pointCloud.size() << std::endl;
    } else {
        std::cout << "INFO: Did not find any ground points! Calculating ground as flat surface\n" << std::endl;
        //todo needs to be implemented - handled with the schema for now
    }

    //-- Read building points
    if (!config::buildings_xyz.empty()) {
        std::cout << "Reading building points" << std::endl;
        IO::read_point_cloud(config::buildings_xyz, _pointCloudBuildings);
        if (_pointCloudBuildings.empty()) throw std::invalid_argument("Didn't find any building points!");

        std::cout << "    Points read: " << _pointCloudBuildings.size() << std::endl;
    }

    //-- Read building polygons
    if (!config::gisdata.empty()) {
        std::cout << "Reading polygons" << std::endl;
        IO::read_geojson_polygons(config::gisdata, _polygonsBuildings);
        if (_polygonsBuildings.empty()) throw std::invalid_argument("Didn't find any building polygons!");
    }

    //-- Read surface layer polygons
    for (auto& topoLayer: config::topoLayers) {
        _polygonsSurfaceLayers.emplace_back();
        IO::read_geojson_polygons(topoLayer, _polygonsSurfaceLayers.back());
    }

    //-- Read imported buildings
    if (!config::importedBuildings.empty()) {
//        std::cout << "Importing CityJSON geometries" << std::endl;
        IO::read_explicit_geometries(config::importedBuildings, _importedBuildingsJson, _importedBuildingsPts);
    }
}

void Map3d::output() {
#ifndef NDEBUG
    assert(config::outputSurfaces.size() == TopoFeature::get_num_output_layers());
#endif
    fs::current_path(config::outputDir);
    std::cout << "\nOutputting surface meshes "      << std::endl;
    std::cout << "    Folder: " << fs::canonical(fs::current_path()) << std::endl;
//    std::cout << "    Format: " << config::outputFormat << std::endl; //todo

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
    for (unsigned long i = 0; i < _reconstructedBuildings.size();) {
        if (_reconstructedBuildings[i]->is_active()) ++i;
        else {
            _reconstructedBuildings.erase(_reconstructedBuildings.begin() + i);
        }
    }
    std::vector<std::string> inactiveBuildingIdxs;
    for (auto& importedBuilding : _importedBuildings) {
        if (!importedBuilding->is_active()) {
            auto it = std::find(inactiveBuildingIdxs.begin(), inactiveBuildingIdxs.end(),
                                importedBuilding->get_parent_building_id());
            if (it == inactiveBuildingIdxs.end())
                inactiveBuildingIdxs.push_back(importedBuilding->get_parent_building_id());
        }
    }
    for (unsigned long i = 0; i < _importedBuildings.size();) {
        auto it = std::find(inactiveBuildingIdxs.begin(), inactiveBuildingIdxs.end(),
                            _importedBuildings[i]->get_parent_building_id());
        if (it == inactiveBuildingIdxs.end()) ++i;
        else {
            _importedBuildings[i]->deactivate();
            _importedBuildings.erase(_importedBuildings.begin() + i);
        }
    }
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
        if (!f->is_active() || f->is_imported()) continue;
        for (auto& ring : f->get_poly().rings()) {
            geomutils::shorten_long_poly_edges(ring);
        }
    }
}
//- Explicit template instantiation
template void Map3d::shorten_polygons<Buildings>              (Buildings& feature);
template void Map3d::shorten_polygons<SurfaceLayers>          (SurfaceLayers& feature);
template void Map3d::shorten_polygons<PolyFeatures>           (PolyFeatures& feature);

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