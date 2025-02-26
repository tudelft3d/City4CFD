/*
  City4CFD
 
  Copyright (c) 2021-2025, 3D Geoinformation Research Group, TU Delft

  This file is part of City4CFD.

  City4CFD is free software: you can redistribute it and/or modify
  it under the terms of the GNU Affero General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  City4CFD is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Affero General Public License for more details.

  You should have received a copy of the GNU Affero General Public License
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
    if (!m_bndBPG) {
        //-- Set outer boundary
        this->set_bnd();

        //-- Avoid having too long polygons
        this->shorten_polygons(m_allFeaturesPtr);

        //-- Find footprint elevation of all polygons using smoothed DT
        this->set_footprint_elevation(m_surfaceLayersPtr);
        this->set_footprint_elevation(m_buildingsPtr);

        //-- Reconstruct buildings in the influ region
        this->reconstruct_buildings();
    } else {
        //-- First the buildings are reconstructed
        //- Prepare polygons for buildings
        this->shorten_polygons(m_buildingsPtr);

        //- Find building footprint elevation using smoothed DT
        this->set_footprint_elevation(m_buildingsPtr);

        //- Reconstruct buildings in the influ region
        this->reconstruct_buildings();

        //-- Second, find the highest building and set domain according to BPG
        this->set_bnd();

        //-- Add surface layers now that the domain size is known
        this->shorten_polygons(m_surfaceLayersPtr);
        this->set_footprint_elevation(m_surfaceLayersPtr);
    }

    //-- Clip building bottoms
    if (Config::get().clip) this->clip_buildings();

    //-- Flatten terrain with flag
    if (Config::get().flatTerrain) this->reconstruct_with_flat_terrain();

    //-- Constrain features, generate terrain mesh from CDT
    this->reconstruct_terrain();

    //-- Geometry wrap (experimental)
    if (Config::get().alphaWrapAll) this->wrap();

    //-- Generate side and top boundaries
    if (Config::get().reconstructBoundaries) this->reconstruct_boundaries();
}

void Map3d::set_features() {
    //-- First feature is the terrain
    m_terrainPtr = std::make_shared<Terrain>();

    //-- Add features - order in m_allFeaturesPtr defines the advantage in marking terrain polygons
    //- Buildings
    // Handle the number of building output layers
    TopoFeature::add_recon_region_output_layers(Config::get().reconRegions.size());
    // Initalize buildings
    for (auto& poly : m_polygonsBuildings) {
        auto building = std::make_shared<ReconstructedBuilding>(*poly);
        m_reconstructedBuildingsPtr.push_back(building);
        m_buildingsPtr.push_back(building);
        m_allFeaturesPtr.push_back(building);
    }
    this->clear_inactives(); // Remove buildings that potentially couldn't be imported
    //- Imported buildings
    if (!m_importedBuildingsJSON.empty()) {
        std::cout << "Importing CityJSON geometries" << std::endl;

        std::vector<std::shared_ptr<ImportedBuilding>> appendingBuildings;
        for (auto& importedBuilding: m_importedBuildingsJSON) {
            auto explicitCityJSONGeom = std::make_shared<ImportedBuilding>(importedBuilding, m_importedBuildingsPts);
            if (!explicitCityJSONGeom->is_appending()) {
                m_importedBuildingsPtr.push_back(explicitCityJSONGeom);
                m_buildingsPtr.push_back(explicitCityJSONGeom);
                m_allFeaturesPtr.push_back(explicitCityJSONGeom);
            } else {
                appendingBuildings.push_back(explicitCityJSONGeom);
            }
        }
        //- Check for building parts that do not have footprint and append to another instance of the same building
        for (auto& b: appendingBuildings) {
            for (auto& importedBuilding: m_importedBuildingsPtr) {
                if (b->get_id() == importedBuilding->get_id()) {
                    importedBuilding->append_nonground_part(b);
                    break;
                }
            }
        }
        m_cityjsonInput = true;
        m_importedBuildingsJSON.clear();
    } else if (!m_importedBuildingsOther.empty()) {
        std::cout << "Importing geometries" << std::endl;
        for (auto& mesh : m_importedBuildingsOther) {
            auto explicitOBJGeom = std::make_shared<ImportedBuilding>(mesh);
            m_importedBuildingsPtr.push_back(explicitOBJGeom);
            m_buildingsPtr.push_back(explicitOBJGeom);
            m_allFeaturesPtr.push_back(explicitOBJGeom);
        }
        Config::get().logSummary << "Number of buildings not imported due to bad surface connectivity: "
                                 << ImportedBuilding::noBottom << std::endl;
        m_importedBuildingsOther.clear();
    }
    if (!m_importedBuildingsPtr.empty()) {
        this->clear_inactives();
        std::cout << "    Geometries imported: " << m_importedBuildingsPtr.size() << std::endl;
    }
    //-- Boundaries
    for (int i = 0; i < Config::get().numSides; ++i)
        m_boundariesPtr.push_back(std::make_shared<Sides>(TopoFeature::get_num_output_layers()));
    m_boundariesPtr.push_back(std::make_shared<Top>(TopoFeature::get_num_output_layers()));

    //- Other polygons
    for (auto& surfaceLayer : m_polygonsSurfaceLayers) {
        int outputLayerID = TopoFeature::get_num_output_layers();
        Config::get().surfaceLayerIDs.push_back(outputLayerID); // Need it for later
        for (auto& poly : surfaceLayer) {
            auto surfacePoly = std::make_shared<SurfaceLayer>(*poly, outputLayerID);
            m_surfaceLayersPtr.push_back(surfacePoly);
            m_allFeaturesPtr.push_back(surfacePoly);
        }
    }
    std::cout << "Polygons read: " << m_allFeaturesPtr.size() << std::endl;

    //-- Set flat terrain or random thin terrain points
    if (m_pointCloud.get_terrain().empty()) {
        m_pointCloud.create_flat_terrain(m_allFeaturesPtr);
        Config::get().flatTerrain = false; // all points are already at 0 elevation
    } else {
        m_pointCloud.random_thin_pts();
    }
    if (Config::get().flatTerrain) std::cout << "\nINFO: Reconstructing with flat terrain" << std::endl;

    //-- Smooth terrain
    if (Config::get().smoothTerrain) {
        m_pointCloud.smooth_terrain();
    }
    //-- Make a DT with inexact constructions for fast interpolation
    m_dt.insert(m_pointCloud.get_terrain().points().begin(),
                m_pointCloud.get_terrain().points().end());

    //-- Initialize bounding regions (reconstruction/influence and domain boundary)
    //BPG flags for influ region and domain boundary
    for (auto& reconRegion : Config::get().reconRegions) {
        m_reconRegions.emplace_back(reconRegion);
    }
    if (std::holds_alternative<bool>(Config::get().domainBndConfig)) m_bndBPG = true;

    //-- Apply area filtering
    if (Config::get().minArea > 0.) {
        for (auto& b: m_buildingsPtr) {
            if (b->get_poly().outer_boundary().area() < Config::get().minArea) b->deactivate();
            Config::write_to_log("Building ID: " + b->get_id() + "  Skipped reconstruction; footprint smaller than defined value.");
        }
        this->clear_inactives();
    }
}

void Map3d::add_building_pts() {
    if (m_pointCloud.get_buildings().empty() || m_reconstructedBuildingsPtr.empty()) return;
    std::cout << "\nAttaching point cloud points to buildings" << std::endl;

    //-- Construct a search tree from all building points
    SearchTree searchTree(m_pointCloud.get_buildings().points().begin(),
                          m_pointCloud.get_buildings().points().end(),
                          Config::get().searchtree_bucket_size);

    m_pointCloud.get_buildings().clear(); // release the loaded building point cloud from memory

    //-- Find points belonging to individual buildings
    for (auto& b: m_reconstructedBuildingsPtr) {
        auto cgalPoly = b->get_poly().get_cgal_type();
        auto poly = geomutils::offset_polygon_geos(cgalPoly, 2.);

        std::vector<Point_3> subsetPts;
        Point_2 bbox1(poly.bbox().xmin(), poly.bbox().ymin());
        Point_2 bbox2(poly.bbox().xmax(), poly.bbox().ymax());
        Fuzzy_iso_box pts_range(bbox1, bbox2);
        searchTree.search(std::back_inserter(subsetPts), pts_range);

        //-- Check if subset point lies inside the offset polygon
        for (auto& pt : subsetPts) {
            if (geomutils::point_in_poly_and_boundary(pt, poly)) {
                b->insert_point(pt);
            }
        }
    }
}

void Map3d::remove_extra_terrain_pts() {
    std::cout << "\nRemoving extra terrain points" << std::endl;
    //-- Handle terrain points that lay in buildings
    m_pointCloud.terrain_points_in_polygon(m_buildingsPtr);
    //-- Update DT for interpolation
    m_dt.clear();
    m_dt.insert(m_pointCloud.get_terrain().points().begin(),
                m_pointCloud.get_terrain().points().end());
}

void Map3d::set_influ_region() {
    std::cout << "\nDefining influence region" << std::endl;
    //-- Set the reconstruction (influence) regions --//
    double maxDim = -1.; // this works if there's one point of interest
    for (int i = 0; i < m_reconRegions.size(); ++i) {
        if (std::holds_alternative<bool>(m_reconRegions[i].m_reconSettings->influRegionConfig)) {// bool defines BPG request
            std::cout << "INFO: Reconstruction region "<< i << " not defined in config. "
                      << "Calculating with BPG." << std::endl;
            if (maxDim < 0.)
                maxDim = m_reconRegions[i].calc_influ_region_bpg(m_dt, m_buildingsPtr);
            else
                m_reconRegions[i].calc_influ_region_bpg(maxDim);
        } else
            std::visit(m_reconRegions[i], Config::get().reconRegions[i]->influRegionConfig);
    }
    // Check if regions get larger with increasing index
    for (int i = 0; i < m_reconRegions.size(); ++i) {
        if (i == 0) continue;
        if (!m_reconRegions[i - 1].is_subset_of(m_reconRegions[i]))
            std::cout << "WARNING: Reconstruction region "
                    << i - 1 << " is not a full subset of region " << i << std::endl;
    }
    //-- Set the reconstruction rules from reconstruction regions to individual buildings
    //   also filter out buildings that are not being reconstructed
    for (auto& b: m_buildingsPtr) {
        for (auto& reconRegion : m_reconRegions) {
            if (!b->has_reconstruction_region() && b->is_part_of(reconRegion.get_bounding_region())) // first come, first served with region setup
                b->set_reconstruction_rules(reconRegion);
        }
        if (!b->has_reconstruction_region()) b->deactivate();
    }
    this->clear_inactives();

    //-- Check if imported and reconstructed buildings are overlapping
    if (!m_importedBuildingsPtr.empty()) this->solve_building_conflicts();

    std::cout << "    Number of building geometries in the influence region: " << m_buildingsPtr.size() << std::endl;
    if (m_buildingsPtr.empty()) {
        throw city4cfd_error("No buildings were reconstructed in the influence region!"
                                 " If using polygons and point cloud, make sure they are aligned.");
    }
}

void Map3d::set_bnd() {
    if (m_bndBPG) { // Automatically calculate boundary with BPG
        std::cout << "\nINFO: Domain boundaries not defined in config. "
                  << "Calculating with BPG." << std::endl;

        //-- Calculate the boundary polygon according to BPG and defined boundary type
        m_domainBnd.calc_bnd_bpg(m_reconRegions.back().get_bounding_region(), m_buildingsPtr);
    } else {
        //-- Define boundary region with values set in config
        std::visit(m_domainBnd, Config::get().domainBndConfig);
    }
    this->bnd_sanity_check(); // Check if outer bnd is larger than the influ region

    //-- Prepare the outer boundary polygon for sides and top, and polygon for feature scope
    Polygon_2 bndPoly, pcBndPoly, startBufferPoly; // Depends on the buffer region
    bndPoly = m_domainBnd.get_bounding_region();
    if (m_boundariesPtr.size() > 2) {
        geomutils::shorten_long_poly_edges(bndPoly, 20 * Config::get().edgeMaxLen); // Outer poly edge size is hardcoded atm
        Boundary::set_bnd_poly(bndPoly, pcBndPoly, startBufferPoly);
    } else
        // If it's only one side bnd, edge length is already okay
        Boundary::set_bnd_poly(bndPoly, pcBndPoly, startBufferPoly);

    //-- Deactivate point cloud points that are out of bounds
    Boundary::set_bounds_to_terrain_pc(m_pointCloud.get_terrain(),
                                       bndPoly, pcBndPoly, startBufferPoly);

    // update the terrain DT for interpolation
    m_dt.clear();
    m_dt.insert(m_pointCloud.get_terrain().points().begin(),
                m_pointCloud.get_terrain().points().end());

    //-- Check feature scope for surface layers now that the full domain is known
    for (auto& f: m_surfaceLayersPtr) {
        f->check_feature_scope(bndPoly);
    }
    this->clear_inactives();
}

void Map3d::bnd_sanity_check() {
    auto& domainBndPoly = m_domainBnd.get_bounding_region();
    for (auto& pt : m_reconRegions.back().get_bounding_region()) {
        if (!geomutils::point_in_poly(pt, domainBndPoly))
            throw city4cfd_error("The influence region is larger than the domain boundary!");
    }
}

void Map3d::reconstruct_terrain() {
    this->clear_inactives();
    if (m_terrainPtr->get_cdt().number_of_vertices() == 0) {
        std::cout << "\nReconstructing terrain" << std::endl;
        m_terrainPtr->prep_constraints(m_allFeaturesPtr, m_pointCloud.get_terrain());
        // Handle flattening
        if (!Config::get().flattenSurfaces.empty()) {
            std::vector<std::pair<Polygon_with_holes_2, int>> additionalPolys;
            m_pointCloud.flatten_polygon_pts(m_allFeaturesPtr, m_terrainPtr->get_extra_constrained_edges(), additionalPolys);
            if (!additionalPolys.empty()) {
                for (auto& polyToAdd : additionalPolys) {
                    Polygon_with_attr newPolyToAdd;
                    newPolyToAdd.polygon = polyToAdd.first;
                    auto surfacePoly
                            = std::make_shared<SurfaceLayer>(newPolyToAdd, polyToAdd.second);
                    m_surfaceLayersPtr.push_back(surfacePoly);
                    m_allFeaturesPtr.push_back(surfacePoly);
                }
            }
        }
        m_terrainPtr->set_cdt(m_pointCloud.get_terrain());
        m_terrainPtr->constrain_features();
    }

    std::cout << "\n    Creating terrain mesh" << std::endl;
    m_terrainPtr->create_mesh(m_allFeaturesPtr);
}

void Map3d::reconstruct_buildings() {
    std::cout << "\nReconstructing buildings" << std::endl;
    if (!m_importedBuildingsPtr.empty() && m_cityjsonInput) {
        std::cout << "    Will try to reconstruct imported buildings in LoD: " << Config::get().importLoD
                  << ". If I cannot find a geometry with that LoD, I will reconstruct in the highest LoD available"
                  << std::endl;
    }
    int count = 0;
    int tenPercent = m_buildingsPtr.size() / 10;
    #pragma omp parallel for
    for (int i = 0; i < m_buildingsPtr.size(); ++i) {
    //for (auto& f : m_buildingsPtr) { // MSVC doesn't like range loops with OMP
        auto& b = m_buildingsPtr[i];
        if (b->is_active()) this->reconstruct_one_building(b);

        if ((count % tenPercent) == 0)
            #pragma omp critical
            IO::print_progress_bar(100 * count / m_buildingsPtr.size());

        #pragma omp atomic
        ++count;
    }
    IO::print_progress_bar(100); std::cout << std::endl;
    /* //todo for groundPts
    // handle cases when a building is a multipart from roofer
    # pragma omp parallel for
    for (int i = 0; i < m_buildingsPtr.size(); ++i) {
        auto& b = m_buildingsPtr[i];
        if (!b->is_active()) continue;
        auto reconBuilding = std::dynamic_pointer_cast<ReconstructedBuilding>(b);
        if (reconBuilding) {
            auto& reconBuildingMeshes = reconBuilding->get_roofer_meshes();
            for (int j = 1; j < reconBuildingMeshes.size(); ++j) {
                auto newBuilding = std::make_shared<ReconstructedBuilding>(reconBuildingMeshes[j], reconBuilding);
                m_reconstructedBuildingsPtr.push_back(newBuilding);
                m_buildingsPtr.push_back(newBuilding);
                m_allFeaturesPtr.push_back(newBuilding);
            }
        }
    }
    */
    this->clear_inactives(); // renumber failed and in case of imported-reconstructed fallback
    // Gather failed reconstructions
    std::cout << "    Number of successfully reconstructed buildings: " << m_buildingsPtr.size() << std::endl;
    Config::get().logSummary << "Building reconstruction summary: successfully reconstructed buildings: "
                             << m_buildingsPtr.size() << std::endl;
    Config::get().logSummary << "                                 num of failed reconstructions: "
                             << m_failedBuildingsPtr.size() << std::endl;
}

void Map3d::reconstruct_one_building(std::shared_ptr<Building>& building) {
    try {
        building->reconstruct();
        //-- In case of hybrid boolean/constraining reconstruction
        if (Config::get().clip && !Config::get().handleSelfIntersect && building->has_self_intersections()) {
            building->set_clip_flag(false);
            building->reconstruct();
        }
    } catch (std::exception& e) {
        // add information to log file
        Config::write_to_log("Building ID: " + building->get_id() + " Failed to reconstruct. Reason: " + e.what());
        // fallback for failed reconstruction of imported buildings
        if (building->is_imported()) {
            building->deactivate(); // deactivate this and use reconstructed instead
            // try to recover by reconstructing LoD1.2 from geometry pts
                auto importToReconstructBuild =
                        std::make_shared<ReconstructedBuilding>(std::static_pointer_cast<ImportedBuilding>(building));
            #pragma omp critical
            {
                m_reconstructedBuildingsPtr.push_back(importToReconstructBuild);
                m_allFeaturesPtr.push_back(importToReconstructBuild);
                m_buildingsPtr.push_back(importToReconstructBuild);
            }
            std::shared_ptr<Building> buildToReconstruct = importToReconstructBuild;
            this->reconstruct_one_building(buildToReconstruct);
        } else {
            // mark for geojson output
            building->mark_as_failed();
        }
    }
}

void Map3d::reconstruct_boundaries() {
    std::cout << "\nReconstructing boundaries" << std::endl;
    if (m_boundariesPtr.size() > 2) { // Means more than one side
        for (auto i = 0; i < m_boundariesPtr.size() - 1; ++i) {
            //-- Each boundary object is one side of the boundary
            m_boundariesPtr[i]->prep_output(m_domainBnd.get_bounding_region().edge(i).to_vector());
        }
    } else {
        m_boundariesPtr.front()->prep_output();
    }
    for (auto& b : m_boundariesPtr) {
        b->reconstruct();
    }
}

void Map3d::reconstruct_with_flat_terrain() {
    //-- Account for zero terrain height of surface layers
    for (auto& sl : m_surfaceLayersPtr) {
        sl->set_zero_borders();
    }
    //-- Account for zero terrain height of buildings
    for (auto& b : m_buildingsPtr) {
        if (!b->is_active()) continue; // skip failed reconstructions
        b->set_to_zero_terrain();
    }
    //-- Set terrain point cloud to zero height
    m_pointCloud.set_flat_terrain();
}

void Map3d::solve_building_conflicts() {
    for (auto& importedBuilding : m_importedBuildingsPtr) {
        for (auto& reconstructedBuilding : m_reconstructedBuildingsPtr) {
            if (geomutils::polygons_in_contact(importedBuilding->get_poly(), reconstructedBuilding->get_poly())) {
                if (reconstructedBuilding->get_reconstruction_settings()->importAdvantage) {
                    reconstructedBuilding->deactivate();
                } else {
                    importedBuilding->deactivate();
                }
            }
        }
    }
    this->clear_inactives();
   // to check if conflicts are solved
//    for (auto& b : m_importedBuildingsPtr) b->deactivate();
//    this->clear_inactives();
}

void Map3d::clip_buildings() {
    //-- Prepare terrain with subset
    std::cout << "\nReconstructing terrain" << std::endl;
    m_terrainPtr->prep_constraints(m_allFeaturesPtr, m_pointCloud.get_terrain());
    // Handle flattening
    if (!Config::get().flattenSurfaces.empty()) {
        std::vector<std::pair<Polygon_with_holes_2, int>> additionalPolys;
        m_pointCloud.flatten_polygon_pts(m_allFeaturesPtr, m_terrainPtr->get_extra_constrained_edges(), additionalPolys);
        if (!additionalPolys.empty()) {
            for (auto& polyToAdd : additionalPolys) {
                Polygon_with_attr newPolyToAdd;
                newPolyToAdd.polygon = polyToAdd.first;
                auto surfacePoly = std::make_shared<SurfaceLayer>(newPolyToAdd, polyToAdd.second);
                m_surfaceLayersPtr.push_back(surfacePoly);
                m_allFeaturesPtr.push_back(surfacePoly);
            }
        }
    }
    m_terrainPtr->set_cdt(m_pointCloud.get_terrain());
    m_terrainPtr->constrain_features();
    m_terrainPtr->prepare_subset();

    //-- Do the clipping
    std::cout << "\n    Clipping buildings to terrain" << std::endl;
    int count = 0;
    for (auto& b : m_buildingsPtr) {
        if (!b->is_active()) continue; // skip failed reconstructions
        b->clip_bottom(m_terrainPtr);

        if ((count % 50) == 0) IO::print_progress_bar(100 * count / m_buildingsPtr.size());
        ++count;
    }
    IO::print_progress_bar(100); std::cout << std::endl;
    m_terrainPtr->clear_subset();
}

void Map3d::wrap() {
    std::cout << "\nAlpha wrapping all buildings..." << std::flush;

    //-- New mesh that will be output of wrapping
    Mesh newMesh;

    //-- Perform alpha wrapping
    Building::alpha_wrap_all(m_buildingsPtr, newMesh);

    //-- Deactivate all individual buildings and add the new mesh
    for (auto& b : m_buildingsPtr) b->deactivate();
    this->clear_inactives();
    m_buildingsPtr.push_back(std::make_shared<ReconstructedBuilding>(newMesh));
    // add reconstruction settings from the first region
    m_buildingsPtr.back()->set_reconstruction_rules(m_reconRegions.front());
}

void Map3d::read_data() {
    //-- Read point clouds
    m_pointCloud.read_point_clouds();

    //-- Read building polygons
    if (!Config::get().gisdata.empty()) {
        std::cout << "Reading polygons" << std::endl;
        IO::read_polygons(Config::get().gisdata, m_polygonsBuildings, &Config::get().crsInfo);
        if (m_polygonsBuildings.empty()) throw std::invalid_argument("Didn't find any building polygons!");
    }
    //-- Read surface layer polygons
    for (auto& topoLayer: Config::get().topoLayers) {
        m_polygonsSurfaceLayers.emplace_back();
        IO::read_polygons(topoLayer, m_polygonsSurfaceLayers.back(), nullptr);
    }
    //-- Read imported buildings
    if (!Config::get().importedBuildingsPath.empty()) {
//        std::cout << "Importing CityJSON geometries" << std::endl;
        auto& inputfile = Config::get().importedBuildingsPath;
        if (IO::has_substr(inputfile, ".json")) {
            m_importedBuildingsPts = std::make_shared<Point_set_3>();
            IO::read_cityjson_geometries(inputfile, m_importedBuildingsJSON, m_importedBuildingsPts);
        } else if (IO::has_substr(inputfile, ".obj") ||
                   IO::has_substr(inputfile, ".stl") ||
                   IO::has_substr(inputfile, ".vtp") ||
                   IO::has_substr(inputfile, ".ply") ||
                   IO::has_substr(inputfile, ".off")) {
            IO::read_other_geometries(inputfile, m_importedBuildingsOther);
        } else {
            throw city4cfd_error(std::string("File " + inputfile + "contains unknown import format."
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
    m_outputFeaturesPtr.push_back(m_terrainPtr);
    for (auto& f : m_buildingsPtr) {
        if (!f->is_active()) continue; // skip failed reconstructions
        m_outputFeaturesPtr.push_back(f);
    }
    for (auto& b : m_boundariesPtr) {
        m_outputFeaturesPtr.push_back(b);
    }
    for (auto& l : m_terrainPtr->get_surface_layers()) { // Surface layers are grouped in terrain
        m_outputFeaturesPtr.push_back(l);
    }

    switch (Config::get().outputFormat) {
        case OBJ:
            IO::output_obj(m_outputFeaturesPtr);
            break;
        case STL: // Only ASCII stl for now
            IO::output_stl(m_outputFeaturesPtr);
            break;
        case CityJSON:
            //-- Remove inactives and add ID's to features - obj and stl don't need id
            // just temp for now
            this->prep_cityjson_output();
            IO::output_cityjson(m_outputFeaturesPtr);
            break;
    }
}

void Map3d::prep_cityjson_output() { // Temp impl, might change
    for (unsigned long i = 0; i < m_outputFeaturesPtr.size(); ++i) {
        if (m_outputFeaturesPtr[i]->is_active()) {
            m_outputFeaturesPtr[i]->set_id(i);
            ++i;
        }
        else {
            m_outputFeaturesPtr.erase(m_outputFeaturesPtr.begin() + i);
        }
    }
}

void Map3d::clear_inactives() {
    m_reconstructedBuildingsPtr.erase(
            std::remove_if(
                    m_reconstructedBuildingsPtr.begin(),
                    m_reconstructedBuildingsPtr.end(),
                    [](const ReconstructedBuildingPtr& b) { return !b->is_active(); }
            ),
            m_reconstructedBuildingsPtr.end()
    );
    if (m_cityjsonInput) {
        std::vector<std::string> inactiveBuildingIdxs;
        for (auto& importedBuilding: m_importedBuildingsPtr) {
            if (!importedBuilding->is_active()) {
                auto it = std::find(inactiveBuildingIdxs.begin(), inactiveBuildingIdxs.end(),
                                    importedBuilding->get_id());
                if (it == inactiveBuildingIdxs.end())
                    inactiveBuildingIdxs.push_back(importedBuilding->get_id());
            }
        }
        for (unsigned long i = 0; i < m_importedBuildingsPtr.size();) {
            auto it = std::find(inactiveBuildingIdxs.begin(), inactiveBuildingIdxs.end(),
                                m_importedBuildingsPtr[i]->get_id());
            if (it == inactiveBuildingIdxs.end()) ++i;
            else {
                m_importedBuildingsPtr[i]->deactivate();
                m_importedBuildingsPtr.erase(m_importedBuildingsPtr.begin() + i);
            }
        }
    } else {
        m_importedBuildingsPtr.erase(
                std::remove_if(
                        m_importedBuildingsPtr.begin(),
                        m_importedBuildingsPtr.end(),
                        [](const ImportedBuildingPtr& b) { return !b->is_active(); }
                ),
                m_importedBuildingsPtr.end()
        );
    }
    for (auto& b : m_buildingsPtr)
        if (b->has_failed_to_reconstruct()) m_failedBuildingsPtr.push_back(b);
    m_buildingsPtr.erase(
            std::remove_if(
                    m_buildingsPtr.begin(),
                    m_buildingsPtr.end(),
                    [](const BuildingPtr & b) { return !b->is_active(); }
            ),
            m_buildingsPtr.end()
    );
    m_surfaceLayersPtr.erase(
            std::remove_if(
                    m_surfaceLayersPtr.begin(),
                    m_surfaceLayersPtr.end(),
                    [](const std::shared_ptr<SurfaceLayer>& s) { return !s->is_active(); }
            ),
            m_surfaceLayersPtr.end()
    );
    m_allFeaturesPtr.erase(
            std::remove_if(
                    m_allFeaturesPtr.begin(),
                    m_allFeaturesPtr.end(),
                    [](const std::shared_ptr<PolyFeature>& p) { return !p->is_active(); }
            ),
            m_allFeaturesPtr.end()
    );
}

const BuildingsPtr& Map3d::get_failed_buildings() const {
    return m_failedBuildingsPtr;
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
        f->calc_footprint_elevation_nni(m_dt);
#else
        f->calc_footprint_elevation_linear(m_dt);  // NNI is quite slow in debug mode, better to use linear in that case
#endif
    }
}
//- Explicit template instantiation
template void Map3d::set_footprint_elevation<BuildingsPtr>    (BuildingsPtr& feature);
template void Map3d::set_footprint_elevation<SurfaceLayersPtr>(SurfaceLayersPtr& feature);
template void Map3d::set_footprint_elevation<PolyFeaturesPtr> (PolyFeaturesPtr& feature);
