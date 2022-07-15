/*
  City4CFD
 
  Copyright (c) 2021-2022, 3D Geoinformation Research Group, TU Delft  

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

    //-- Clip building bottoms
    if (Config::get().clip) this->clip_buildings();

    //-- Flatten terrain with flag
    if (Config::get().flatTerrain) this->reconstruct_with_flat_terrain();

    //-- Constrain features, generate terrain mesh from CDT
    this->reconstruct_terrain();

    //-- Generate side and top boundaries
    if (Config::get().reconstructBoundaries) this->reconstruct_boundaries();
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
    if (Config::get().avoidBadPolys) this->clear_inactives(); // Remove buildings that potentially couldn't be imported
    //- Imported buildings
    if (!_importedBuildingsJSON.empty()) {
        std::cout << "Importing CityJSON geometries" << std::endl;

        std::vector<std::shared_ptr<ImportedBuilding>> appendingBuildings;
        internalID = 0;
        for (auto& importedBuilding: _importedBuildingsJSON) {
            auto explicitCityJSONGeom = std::make_shared<ImportedBuilding>(importedBuilding, _importedBuildingsPts,
                                                                           internalID++);
            if (!explicitCityJSONGeom->is_appending()) {
                _importedBuildings.push_back(explicitCityJSONGeom);
                _buildings.push_back(explicitCityJSONGeom);
                _lsFeatures.push_back(explicitCityJSONGeom);
            } else {
                appendingBuildings.push_back(explicitCityJSONGeom);
            }
        }
        //- Check for building parts that do not have footprint and append to another instance of the same building
        for (auto& b: appendingBuildings) {
            for (auto& importedBuilding: _importedBuildings) {
                if (b->get_parent_building_id() == importedBuilding->get_parent_building_id()) {
                    importedBuilding->append_nonground_part(b);
                    break;
                }
            }
        }
        _cityjsonInput = true;
        _importedBuildingsJSON.clear();
    } else if (!_importedBuildingsOther.empty()) {
        std::cout << "Importing geometries" << std::endl;
        for (auto& mesh : _importedBuildingsOther) {
            auto explicitOBJGeom = std::make_shared<ImportedBuilding>(mesh, internalID++);
            _importedBuildings.push_back(explicitOBJGeom);
            _buildings.push_back(explicitOBJGeom);
            _lsFeatures.push_back(explicitOBJGeom);
        }
        Config::get().logSummary << "Number of buildings not imported due to bad surface connectivity: "
                                 << ImportedBuilding::noBottom << std::endl;
        _importedBuildingsOther.clear();
    }
    if (!_importedBuildings.empty()) {
        this->clear_inactives();
        std::cout << "    Geometries imported: " << _importedBuildings.size() << std::endl;
    }
    //-- Boundaries
    for (int i = 0; i < Config::get().numSides; ++i)
        _boundaries.push_back(std::make_shared<Sides>(TopoFeature::get_num_output_layers()));
    _boundaries.push_back(std::make_shared<Top>(TopoFeature::get_num_output_layers()));

    //- Other polygons
    for (auto& surfaceLayer : _polygonsSurfaceLayers) {
        int outputLayerID = TopoFeature::get_num_output_layers();
        Config::get().surfaceLayerIDs.push_back(outputLayerID); // Need it for later
        for (auto& poly : surfaceLayer) {
            auto surfacePoly = std::make_shared<SurfaceLayer>(*poly, outputLayerID);
            _surfaceLayers.push_back(surfacePoly);
            _lsFeatures.push_back(surfacePoly);
        }
    }
    std::cout << "    Polygons read: " << _lsFeatures.size() << std::endl;

    //-- Set flat terrain or random thin terrain points
    if (_pointCloud.get_terrain().empty()) {
        _pointCloud.create_flat_terrain(_lsFeatures);
        Config::get().flatTerrain = false;
    } else {
        _pointCloud.random_thin_pts();
    }

    //-- BPG flags for influ region and domain boundary
    if (Config::get().influRegionConfig.type() == typeid(bool)) _influRegionBPG = true;
    if (Config::get().domainBndConfig.type() == typeid(bool))   _bndBPG = true;

    //-- Smooth terrain
    if (Config::get().smoothTerrain) {
//        geomutils::smooth_dt<DT, EPICK>(_pointCloud.get_terrain(), _dt);
        _pointCloud.smooth_terrain();
    }
    //-- Make a DT with inexact constructions for fast interpolation
    _dt.insert(_pointCloud.get_terrain().points().begin(),
               _pointCloud.get_terrain().points().end());
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
        std::shared_ptr<SearchTree> searchTree = _pointCloud.make_search_tree_buildings();
        for (auto& building : _reconstructedBuildings) building->set_search_tree(searchTree);

        //-- Calculate influ region
        _influRegion.calc_influ_region_bpg(_dt, _pointCloud.get_buildings(), _buildings);
    } else { // Define influ region either with radius or predefined polygon
        boost::apply_visitor(_influRegion, Config::get().influRegionConfig);
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
        boost::apply_visitor(_domainBnd, Config::get().domainBndConfig);
    }
    this->bnd_sanity_check(); // Check if outer bnd is larger than the influ region

    //-- Prepare the outer boundary polygon for sides and top, and polygon for feature scope
    Polygon_2 bndPoly, pcBndPoly, startBufferPoly; // Depends on the buffer region
    bndPoly = _domainBnd.get_bounding_region();
    if (_boundaries.size() > 2) { 
        geomutils::shorten_long_poly_edges(bndPoly, 20 * Config::get().edgeMaxLen); // Outer poly edge size is hardcoded atm
        Boundary::set_bnd_poly(bndPoly, pcBndPoly, startBufferPoly);
    } else
        // If it's only one side bnd, edge length is already okay
        Boundary::set_bnd_poly(bndPoly, pcBndPoly, startBufferPoly);

    //-- Deactivate point cloud points that are out of bounds
    Boundary::set_bounds_to_terrain(_pointCloud.get_terrain(),
                                    bndPoly, pcBndPoly, startBufferPoly);
    Boundary::set_bounds_to_pc(_pointCloud.get_buildings(), startBufferPoly);

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
    if (_terrain->get_cdt().number_of_vertices() == 0) {
        std::cout << "\nReconstructing terrain" << std::endl;
        _terrain->prep_constraints(_lsFeatures, _pointCloud.get_terrain());
        if (!Config::get().flattenSurfaces.empty()) _pointCloud.flatten_polygon_pts(_lsFeatures);
        _terrain->set_cdt(_pointCloud.get_terrain());
        _terrain->constrain_features();
    }

    std::cout << "\n    Creating terrain mesh" << std::endl;
    _terrain->create_mesh(_lsFeatures);
}

void Map3d::reconstruct_buildings() {
    std::cout << "\nReconstructing buildings" << std::endl;
    if (!_importedBuildings.empty() && _cityjsonInput) {
        std::cout << "    Will try to reconstruct imported buildings in LoD: " << Config::get().importLoD
                  << ". If I cannot find a geometry with that LoD, I will reconstruct in the highest LoD available"
                  << std::endl;
    }

    std::shared_ptr<SearchTree> searchTree = _pointCloud.make_search_tree_buildings();
    for (auto& building : _reconstructedBuildings) building->set_search_tree(searchTree);

    int failed = 0;
    for (auto& f : _buildings) {
        if (!f->is_active()) continue;
        try {
            f->reconstruct();
            //-- In case of hybrid boolean/constraining reconstruction
            if (Config::get().clip && !Config::get().handleSelfIntersect && f->has_self_intersections()) {
                f->set_clip_flag(false);
                f->reconstruct();
            }
            if (Config::get().refineBuildings) f->refine();
        } catch (std::exception& e) {
            ++failed;
            //-- Add information to log file
            Config::get().log << "Failed to reconstruct building ID: " << f->get_id()
                        << " Reason: " << e.what() << std::endl;
            //-- Get JSON file ID for failed reconstructions output
            //   For now only polygons (reconstructed buildings) are stored to GeoJSON
            if (!f->is_imported())
                Config::get().failedBuildings.push_back(f->get_internal_id());
        }
    }
    Config::get().logSummary << "Building reconstruction summary: total failed reconstructions : "
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

void Map3d::reconstruct_with_flat_terrain() {
    //-- Account for zero terrain height of surface layers
    for (auto& sl : _surfaceLayers) {
        sl->set_zero_borders();
    }
    //-- Account for zero terrain height of buildings
    for (auto& b : _buildings) {
        b->set_to_zero_terrain();
    }
    //-- Set terrain point cloud to zero height
    _pointCloud.set_flat_terrain();
}

void Map3d::solve_building_conflicts() {
    for (auto& importedBuilding : _importedBuildings) {
        for (auto& reconstructedBuilding : _reconstructedBuildings) {
            if (geomutils::polygons_in_contact(importedBuilding->get_poly(), reconstructedBuilding->get_poly())) {
                if (Config::get().importAdvantage) {
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

void Map3d::clip_buildings() {
    //-- Prepare terrain with subset
    std::cout << "\nReconstructing terrain" << std::endl;
    _terrain->prep_constraints(_lsFeatures, _pointCloud.get_terrain());
    if (!Config::get().flattenSurfaces.empty()) _pointCloud.flatten_polygon_pts(_lsFeatures);
    _terrain->set_cdt(_pointCloud.get_terrain());
    _terrain->constrain_features();
    _terrain->prepare_subset();

    //-- Do the clipping
    std::cout << "\n    Clipping buildings to terrain" << std::endl;
    int count = 0;
    for (auto& b : _buildings) {
        b->clip_bottom(_terrain);

        if ((count % 50) == 0) IO::print_progress_bar(100 * count / _buildings.size());
        ++count;
    }
    IO::print_progress_bar(100); std::cout << std::endl;
    _terrain->clear_subset();
}

void Map3d::read_data() {
    //-- Read point clouds
    _pointCloud.read_point_clouds();

    //-- Read building polygons
    if (!Config::get().gisdata.empty()) {
        std::cout << "Reading polygons" << std::endl;
        IO::read_geojson_polygons(Config::get().gisdata, _polygonsBuildings);
        if (_polygonsBuildings.empty()) throw std::invalid_argument("Didn't find any building polygons!");
    }
    //-- Read surface layer polygons
    for (auto& topoLayer: Config::get().topoLayers) {
        _polygonsSurfaceLayers.emplace_back();
        IO::read_geojson_polygons(topoLayer, _polygonsSurfaceLayers.back());
    }
    //-- Read imported buildings
    if (!Config::get().importedBuildingsPath.empty()) {
//        std::cout << "Importing CityJSON geometries" << std::endl;
        auto& inputfile = Config::get().importedBuildingsPath;
        if (IO::has_substr(inputfile, ".json")) {
            _importedBuildingsPts = std::make_shared<std::vector<Point_3>>();
            IO::read_cityjson_geometries(inputfile, _importedBuildingsJSON, _importedBuildingsPts);
        } else if (IO::has_substr(inputfile, ".obj") ||
                   IO::has_substr(inputfile, ".stl") ||
                   IO::has_substr(inputfile, ".vtp") ||
                   IO::has_substr(inputfile, ".ply") ||
                   IO::has_substr(inputfile, ".off")) {
            IO::read_other_geometries(inputfile, _importedBuildingsOther);
        } else {
            throw std::runtime_error(std::string("File " + inputfile + "contains unknown import format."
                                                                  " Available inputs: .obj, .stl, .vtp, "
                                                                  ".ply. .off, or .json (CityJSON)"));
        }
    }
}

void Map3d::output() {
#ifndef NDEBUG
    assert(Config::get().outputSurfaces.size() == TopoFeature::get_num_output_layers());
#endif
    fs::current_path(Config::get().outputDir);
    std::cout << "\nOutputting surface meshes "      << std::endl;
    std::cout << "    Folder: " << fs::canonical(fs::current_path()) << std::endl;
//    std::cout << "    Format: " << Config::get().outputFormat << std::endl; //todo

    //-- Group all features for output
    this->prep_feature_output();

    switch (Config::get().outputFormat) {
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
            _outputFeatures[i]->set_id(i);
            ++i;
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
    if (_cityjsonInput) {
        std::vector<std::string> inactiveBuildingIdxs;
        for (auto& importedBuilding: _importedBuildings) {
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
    } else {
        for (unsigned long i = 0; i < _importedBuildings.size();) {
            if (_importedBuildings[i]->is_active()) ++i;
            else {
                _importedBuildings.erase(_importedBuildings.begin() + i);
            }
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