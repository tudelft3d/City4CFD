/*
  City4CFD
 
  Copyright (c) 2021-2023, 3D Geoinformation Research Group, TU Delft

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

    //-- Attach points to buildings
    this->add_building_pts();

    //-- Define influence region
    this->set_influ_region();

    //-- Remove points belonging to buildings from terrain
    this->remove_extra_terrain_pts();

    //-- Different flow if explicitly defining domain boundary or leaving it to BPG
    if (!_bndBPG) {
        //-- Set outer boundary
        this->set_bnd();

        //-- Avoid having too long polygons
        this->shorten_polygons(_allFeaturesPtr);

        //-- Find footprint elevation of all polygons using smoothed DT
        this->set_footprint_elevation(_allFeaturesPtr);

        //-- Reconstruct buildings in the influ region
        this->reconstruct_buildings();
    } else {
        //-- First the buildings are reconstructed
        //- Prepare polygons for buildings
        this->shorten_polygons(_buildingsPtr);

        //- Find building footprint elevation using smoothed DT
        this->set_footprint_elevation(_buildingsPtr);

        //- Reconstruct buildings in the influ region
        this->reconstruct_buildings();

        //-- Second, find the highest building and set domain according to BPG
        this->set_bnd();

        //-- Add surface layers now that the domain size is known
        this->shorten_polygons(_surfaceLayersPtr);
        this->set_footprint_elevation(_surfaceLayersPtr);
    }

    //-- Clip building bottoms
    if (Config::get().clip) this->clip_buildings();

    //-- Flatten terrain with flag
    if (Config::get().flatTerrain) this->reconstruct_with_flat_terrain();

    //-- Constrain features, generate terrain mesh from CDT
    this->reconstruct_terrain();

    //-- Geometry wrap (experimental)
    if (Config::get().alphaWrap) this->wrap();

    //-- Generate side and top boundaries
    if (Config::get().reconstructBoundaries) this->reconstruct_boundaries();
}

void Map3d::set_features() {
    //-- First feature is the terrain
    _terrainPtr = std::make_shared<Terrain>();

    //-- Add features - order in _allFeaturesPtr defines the advantage in marking terrain polygons
    //- Buildings
    int internalID = 0;
    for (auto& poly : _polygonsBuildings) {
        auto building = std::make_shared<ReconstructedBuilding>(*poly, internalID++);
        _reconstructedBuildingsPtr.push_back(building);
        _buildingsPtr.push_back(building);
        _allFeaturesPtr.push_back(building);
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
                _importedBuildingsPtr.push_back(explicitCityJSONGeom);
                _buildingsPtr.push_back(explicitCityJSONGeom);
                _allFeaturesPtr.push_back(explicitCityJSONGeom);
            } else {
                appendingBuildings.push_back(explicitCityJSONGeom);
            }
        }
        //- Check for building parts that do not have footprint and append to another instance of the same building
        for (auto& b: appendingBuildings) {
            for (auto& importedBuilding: _importedBuildingsPtr) {
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
            _importedBuildingsPtr.push_back(explicitOBJGeom);
            _buildingsPtr.push_back(explicitOBJGeom);
            _allFeaturesPtr.push_back(explicitOBJGeom);
        }
        Config::get().logSummary << "Number of buildings not imported due to bad surface connectivity: "
                                 << ImportedBuilding::noBottom << std::endl;
        _importedBuildingsOther.clear();
    }
    if (!_importedBuildingsPtr.empty()) {
        this->clear_inactives();
        std::cout << "    Geometries imported: " << _importedBuildingsPtr.size() << std::endl;
    }
    //-- Boundaries
    for (int i = 0; i < Config::get().numSides; ++i)
        _boundariesPtr.push_back(std::make_shared<Sides>(TopoFeature::get_num_output_layers()));
    _boundariesPtr.push_back(std::make_shared<Top>(TopoFeature::get_num_output_layers()));

    //- Other polygons
    for (auto& surfaceLayer : _polygonsSurfaceLayers) {
        int outputLayerID = TopoFeature::get_num_output_layers();
        Config::get().surfaceLayerIDs.push_back(outputLayerID); // Need it for later
        for (auto& poly : surfaceLayer) {
            auto surfacePoly = std::make_shared<SurfaceLayer>(*poly, outputLayerID);
            _surfaceLayersPtr.push_back(surfacePoly);
            _allFeaturesPtr.push_back(surfacePoly);
        }
    }
    std::cout << "    Polygons read: " << _allFeaturesPtr.size() << std::endl;

    //-- Set flat terrain or random thin terrain points
    if (_pointCloud.get_terrain().empty()) {
        _pointCloud.create_flat_terrain(_allFeaturesPtr);
        Config::get().flatTerrain = false; // all points are already at 0 elevation
    } else {
        _pointCloud.random_thin_pts();
    }
    if (Config::get().flatTerrain) std::cout << "\nINFO: Reconstructing with flat terrain" << std::endl;

    //-- BPG flags for influ region and domain boundary
    if (Config::get().influRegionConfig.type() == typeid(bool)) _influRegionBPG = true;
    if (Config::get().domainBndConfig.type() == typeid(bool))   _bndBPG = true;

    //-- Smooth terrain
    if (Config::get().smoothTerrain) {
        _pointCloud.smooth_terrain();
    }
    //-- Make a DT with inexact constructions for fast interpolation
    _dt.insert(_pointCloud.get_terrain().points().begin(),
               _pointCloud.get_terrain().points().end());
}

void Map3d::add_building_pts() {
    if (_pointCloud.get_buildings().empty() || _reconstructedBuildingsPtr.empty()) return;
    std::cout << "\nAttaching point cloud points to buildings" << std::endl;

    //-- Construct a search tree from all building points
    SearchTree searchTree(_pointCloud.get_buildings().points().begin(),
                          _pointCloud.get_buildings().points().end(),
                          Config::get().searchtree_bucket_size);

    _pointCloud.get_buildings().clear(); // release the loaded building point cloud from memory

    //-- Find points belonging to individual buildings
    for (auto& b: _reconstructedBuildingsPtr) {
        auto& poly = b->get_poly();

        std::vector<Point_3> subsetPts;
        Point_2 bbox1(poly.bbox().xmin(), poly.bbox().ymin());
        Point_2 bbox2(poly.bbox().xmax(), poly.bbox().ymax());
        Fuzzy_iso_box pts_range(bbox1, bbox2);
        searchTree.search(std::back_inserter(subsetPts), pts_range);

        //-- Check if subset point lies inside the polygon
        for (auto& pt : subsetPts) {
            if (geomutils::point_in_poly(pt, poly)) {
                b->insert_point(pt);
            }
        }
    }
}

void Map3d::remove_extra_terrain_pts() {
    std::cout << "\nRemove extra terrain points" << std::endl;
    //-- Remove terrain points that lay in buildings
    _pointCloud.remove_points_in_polygon(_buildingsPtr);
    //-- Update DT for interpolation
    _dt.clear();
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
        if (!_importedBuildingsPtr.empty()) this->solve_building_conflicts(); // have to do it earlier if BPG

        //-- Calculate influ region
        _influRegion.calc_influ_region_bpg(_dt, _buildingsPtr);
    } else { // Define influ region either with radius or predefined polygon
        boost::apply_visitor(_influRegion, Config::get().influRegionConfig);
    }

    //-- Deactivate buildings that are out of influ region
    for (auto& f : _buildingsPtr) {
        f->check_feature_scope(_influRegion.get_bounding_region());
    }
    this->clear_inactives();

    //-- Check if imported and reconstructed buildings are overlapping
    if (!_importedBuildingsPtr.empty() && !_influRegionBPG) this->solve_building_conflicts();

    std::cout << "    Number of building geometries in the influence region: " << _buildingsPtr.size() << std::endl;
    if (_buildingsPtr.empty()) {
        throw std::runtime_error("No buildings were reconstructed in the influence region!"
                                 " If using polygons and point cloud, make sure they are aligned.");
    }
}

void Map3d::set_bnd() {
    if (_bndBPG) { // Automatically calculate boundary with BPG
        std::cout << "\nINFO: Domain boundaries not defined in config. "
                  << "Calculating with BPG." << std::endl;

        //-- Calculate the boundary polygon according to BPG and defined boundary type
        _domainBnd.calc_bnd_bpg(_influRegion.get_bounding_region(), _buildingsPtr);
    } else {
        //-- Define boundary region with values set in config
        boost::apply_visitor(_domainBnd, Config::get().domainBndConfig);
    }
    this->bnd_sanity_check(); // Check if outer bnd is larger than the influ region

    //-- Prepare the outer boundary polygon for sides and top, and polygon for feature scope
    Polygon_2 bndPoly, pcBndPoly, startBufferPoly; // Depends on the buffer region
    bndPoly = _domainBnd.get_bounding_region();
    if (_boundariesPtr.size() > 2) {
        geomutils::shorten_long_poly_edges(bndPoly, 20 * Config::get().edgeMaxLen); // Outer poly edge size is hardcoded atm
        Boundary::set_bnd_poly(bndPoly, pcBndPoly, startBufferPoly);
    } else
        // If it's only one side bnd, edge length is already okay
        Boundary::set_bnd_poly(bndPoly, pcBndPoly, startBufferPoly);

    //-- Deactivate point cloud points that are out of bounds
    Boundary::set_bounds_to_terrain_pc(_pointCloud.get_terrain(),
                                       bndPoly, pcBndPoly, startBufferPoly);

    //-- Check feature scope for surface layers now that the full domain is known
    for (auto& f: _surfaceLayersPtr) {
        f->check_feature_scope(bndPoly);
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
    this->clear_inactives();
    if (_terrainPtr->get_cdt().number_of_vertices() == 0) {
        std::cout << "\nReconstructing terrain" << std::endl;
        _terrainPtr->prep_constraints(_allFeaturesPtr, _pointCloud.get_terrain());
        if (!Config::get().flattenSurfaces.empty())
            _pointCloud.flatten_polygon_pts(_allFeaturesPtr, _terrainPtr->get_extra_constrained_edges());
        _terrainPtr->set_cdt(_pointCloud.get_terrain());
        _terrainPtr->constrain_features();
    }

    std::cout << "\n    Creating terrain mesh" << std::endl;
    _terrainPtr->create_mesh(_allFeaturesPtr);
}

void Map3d::reconstruct_buildings() {
    std::cout << "\nReconstructing buildings" << std::endl;
    if (!_importedBuildingsPtr.empty() && _cityjsonInput) {
        std::cout << "    Will try to reconstruct imported buildings in LoD: " << Config::get().importLoD
                  << ". If I cannot find a geometry with that LoD, I will reconstruct in the highest LoD available"
                  << std::endl;
    }
    int failed = 0;
    #pragma omp parallel for
    for (auto& f : _buildingsPtr) {
        if (!f->is_active()) continue;
        try {
            f->reconstruct();
            //-- In case of hybrid boolean/constraining reconstruction
            if (Config::get().clip && !Config::get().handleSelfIntersect && f->has_self_intersections()) {
                f->set_clip_flag(false);
                f->reconstruct();
            }
        } catch (std::exception& e) {
            #pragma omp atomic
            ++failed;
            // add information to log file
            Config::write_to_log("Building ID: " + f->get_id() + " Failed to reconstruct. Reason: " + e.what());
            //-- Get JSON file ID for failed reconstructions output
            //   For now only polygons (reconstructed buildings) are stored to GeoJSON
            if (!f->is_imported())
                #pragma omp critical
                Config::get().failedBuildings.push_back(f->get_internal_id());
        }
    }
    this->clear_inactives();
    std::cout << "    Number of successfully reconstructed buildings: " << _buildingsPtr.size() << std::endl;
    Config::get().logSummary << "Building reconstruction summary: successfully reconstructed buildings: "
                             << _buildingsPtr.size() << std::endl;
    Config::get().logSummary << "                                 num of failed reconstructions: "
                             << failed << std::endl;
}

void Map3d::reconstruct_boundaries() {
    std::cout << "\nReconstructing boundaries" << std::endl;
    if (_boundariesPtr.size() > 2) { // Means more than one side
        for (auto i = 0; i < _boundariesPtr.size() - 1; ++i) {
            //-- Each boundary object is one side of the boundary
            _boundariesPtr[i]->prep_output(_domainBnd.get_bounding_region().edge(i).to_vector());
        }
    } else {
        _boundariesPtr.front()->prep_output();
    }
    for (auto& b : _boundariesPtr) {
        b->reconstruct();
    }
}

void Map3d::reconstruct_with_flat_terrain() {
    //-- Account for zero terrain height of surface layers
    for (auto& sl : _surfaceLayersPtr) {
        sl->set_zero_borders();
    }
    //-- Account for zero terrain height of buildings
    for (auto& b : _buildingsPtr) {
        b->set_to_zero_terrain();
    }
    //-- Set terrain point cloud to zero height
    _pointCloud.set_flat_terrain();
}

void Map3d::solve_building_conflicts() {
    for (auto& importedBuilding : _importedBuildingsPtr) {
        for (auto& reconstructedBuilding : _reconstructedBuildingsPtr) {
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
//    for (auto& b : _importedBuildingsPtr) b->deactivate();
//    this->clear_inactives();
}

void Map3d::clip_buildings() {
    //-- Prepare terrain with subset
    std::cout << "\nReconstructing terrain" << std::endl;
    _terrainPtr->prep_constraints(_allFeaturesPtr, _pointCloud.get_terrain());
    if (!Config::get().flattenSurfaces.empty())
        _pointCloud.flatten_polygon_pts(_allFeaturesPtr, _terrainPtr->get_extra_constrained_edges());
    _terrainPtr->set_cdt(_pointCloud.get_terrain());
    _terrainPtr->constrain_features();
    _terrainPtr->prepare_subset();

    //-- Do the clipping
    std::cout << "\n    Clipping buildings to terrain" << std::endl;
    int count = 0;
    for (auto& b : _buildingsPtr) {
        b->clip_bottom(_terrainPtr);

        if ((count % 50) == 0) IO::print_progress_bar(100 * count / _buildingsPtr.size());
        ++count;
    }
    IO::print_progress_bar(100); std::cout << std::endl;
    _terrainPtr->clear_subset();
}

void Map3d::wrap() {
    std::cout << "\nAlpha wrapping all buildings..." << std::flush;

    //-- New mesh that will be output of wrapping
    Mesh newMesh;

    //-- Perform alpha wrapping
    Building::alpha_wrap(_buildingsPtr, newMesh);

    //-- Deactivate all individual buildings and add the new mesh
    for (auto& b : _buildingsPtr) b->deactivate();
    this->clear_inactives();
    _buildingsPtr.push_back(std::make_shared<ReconstructedBuilding>(newMesh));
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
            _importedBuildingsPts = std::make_shared<Point_set_3>();
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
    assert(Config::get().outputSurfaces.size() == TopoFeature::get_num_output_layers());
    fs::current_path(Config::get().outputDir);
    std::cout << "\nOutputting surface meshes "      << std::endl;
    std::cout << "    Folder: " << fs::canonical(fs::current_path()) << std::endl;
//    std::cout << "    Format: " << Config::get().outputFormat << std::endl; //todo

    //-- Group all features for output
    this->prep_feature_output();

    switch (Config::get().outputFormat) {
        case OBJ:
            IO::output_obj(_outputFeaturesPtr);
            break;
        case STL: // Only ASCII stl for now
            IO::output_stl(_outputFeaturesPtr);
            break;
        case CityJSON:
            //-- Remove inactives and add ID's to features - obj and stl don't need id
            // just temp for now
            this->prep_cityjson_output();
            IO::output_cityjson(_outputFeaturesPtr);
            break;
    }
}

void Map3d::prep_feature_output() {
    _outputFeaturesPtr.push_back(_terrainPtr);
    for (auto& f : _buildingsPtr) {
        if (!f->is_active()) continue;
        _outputFeaturesPtr.push_back(f);
    }
    for (auto& b : _boundariesPtr) {
        _outputFeaturesPtr.push_back(b);
    }
    for (auto& l : _terrainPtr->get_surface_layers()) { // Surface layers are grouped in terrain
        _outputFeaturesPtr.push_back(l);
    }
}

void Map3d::prep_cityjson_output() { // Temp impl, might change
    for (unsigned long i = 0; i < _outputFeaturesPtr.size(); ++i) {
        if (_outputFeaturesPtr[i]->is_active()) {
            _outputFeaturesPtr[i]->set_id(i);
            ++i;
        }
        else {
            _outputFeaturesPtr.erase(_outputFeaturesPtr.begin() + i);
        }
    }
};

void Map3d::clear_inactives() {
    for (unsigned long i = 0; i < _reconstructedBuildingsPtr.size();) {
        if (_reconstructedBuildingsPtr[i]->is_active()) ++i;
        else {
            _reconstructedBuildingsPtr.erase(_reconstructedBuildingsPtr.begin() + i);
        }
    }
    if (_cityjsonInput) {
        std::vector<std::string> inactiveBuildingIdxs;
        for (auto& importedBuilding: _importedBuildingsPtr) {
            if (!importedBuilding->is_active()) {
                auto it = std::find(inactiveBuildingIdxs.begin(), inactiveBuildingIdxs.end(),
                                    importedBuilding->get_parent_building_id());
                if (it == inactiveBuildingIdxs.end())
                    inactiveBuildingIdxs.push_back(importedBuilding->get_parent_building_id());
            }
        }
        for (unsigned long i = 0; i < _importedBuildingsPtr.size();) {
            auto it = std::find(inactiveBuildingIdxs.begin(), inactiveBuildingIdxs.end(),
                                _importedBuildingsPtr[i]->get_parent_building_id());
            if (it == inactiveBuildingIdxs.end()) ++i;
            else {
                _importedBuildingsPtr[i]->deactivate();
                _importedBuildingsPtr.erase(_importedBuildingsPtr.begin() + i);
            }
        }
    } else {
        for (unsigned long i = 0; i < _importedBuildingsPtr.size();) {
            if (_importedBuildingsPtr[i]->is_active()) ++i;
            else {
                _importedBuildingsPtr.erase(_importedBuildingsPtr.begin() + i);
            }
        }
    }
    for (unsigned long i = 0; i < _buildingsPtr.size();) {
        if (_buildingsPtr[i]->is_active()) ++i;
        else {
            _buildingsPtr.erase(_buildingsPtr.begin() + i);
        }
    }
    for (unsigned long i = 0; i < _surfaceLayersPtr.size();) {
        if (_surfaceLayersPtr[i]->is_active()) ++i;
        else {
            _surfaceLayersPtr.erase(_surfaceLayersPtr.begin() + i);
        }
    }
    for (unsigned long i = 0; i < _allFeaturesPtr.size();) {
        if (_allFeaturesPtr[i]->is_active()) ++i;
        else {
            _allFeaturesPtr.erase(_allFeaturesPtr.begin() + i);
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
template void Map3d::shorten_polygons<BuildingsPtr>    (BuildingsPtr& feature);
template void Map3d::shorten_polygons<SurfaceLayersPtr>(SurfaceLayersPtr& feature);
template void Map3d::shorten_polygons<PolyFeaturesPtr> (PolyFeaturesPtr& feature);

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
template void Map3d::set_footprint_elevation<BuildingsPtr>    (BuildingsPtr& feature);
template void Map3d::set_footprint_elevation<SurfaceLayersPtr>(SurfaceLayersPtr& feature);
template void Map3d::set_footprint_elevation<PolyFeaturesPtr> (PolyFeaturesPtr& feature);