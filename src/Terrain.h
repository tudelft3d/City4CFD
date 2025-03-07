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

#ifndef CITY4CFD_TERRAIN_H
#define CITY4CFD_TERRAIN_H

#include "TopoFeature.h"

typedef std::unordered_map<Point_3, std::vector<face_descriptor>> vertex_face_map;

class Terrain : public TopoFeature {
public:
    using TopoFeature::TopoFeature;
    Terrain();
    Terrain(int pid);
    ~Terrain() = default;

    void set_cdt(const Point_set_3 &pointCloud);
    void prep_constraints(const PolyFeaturesPtr& features, Point_set_3& pointCloud);
    void constrain_features();
    void create_mesh(const PolyFeaturesPtr& features);
    void prepare_subset();
    Mesh mesh_subset(const Polygon_with_holes_2& poly) const;
    void clear_subset();
    void  tag_layers(const PolyFeaturesPtr& features);
    void  tag_layers(const Face_handle& start, int index,
                     std::list<CDT::Edge>& border, const PolyFeaturesPtr& features);
//    void  check_layer(const Face_handle& fh, int surfaceLayer);

    CDT&                     get_cdt();
    const CDT&               get_cdt() const;
    std::vector<Polygon_3>&  get_constrained_polys();
    const vertex_face_map&   get_vertex_face_map() const;
    const SearchTree&        get_mesh_search_tree() const;
    std::vector<EPECK::Segment_3>& get_extra_constrained_edges();

    void         get_cityjson_geomobj_info(nlohmann::json& g) const override;
    void         get_cityjson_cityobj_info(nlohmann::json& f) const override;
    TopoClass    get_class() const override;
    std::string  get_class_name() const override;

    const SurfaceLayersPtr& get_surface_layers() const;

protected:
    CDT                    m_cdt;
    SurfaceLayersPtr       m_surfaceLayersTerrain;
    vertex_face_map        m_vertexFaceMap;
    SearchTree             m_searchTree;
    std::vector<Polygon_3> m_constrainedPolys;
    std::vector<EPECK::Segment_3> m_extraConstrainedEdges;
};

#endif //CITY4CFD_TERRAIN_H