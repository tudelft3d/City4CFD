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

#ifndef CITY4CFD_BUILDING_H
#define CITY4CFD_BUILDING_H

#include "PolyFeature.h"
#include "BoundingRegion.h"
#include "Config.h"

class Building : public PolyFeature {
public:
    Building();
    Building(const nlohmann::json& poly);
    Building(const Polygon_with_attr& poly);
    virtual ~Building() = default;

    static void alpha_wrap_all(const BuildingsPtr& buildings, Mesh& newMesh);

    virtual void   calc_elevation() = 0;
    virtual void   reconstruct() = 0;
    virtual void   reconstruct_flat_terrain() = 0;
    virtual void   insert_terrain_point(const Point_3& pt) = 0;

    double get_elevation();
    double get_height();
    void   insert_point(const Point_3& pt);
    void   clip_bottom(const TerrainPtr& terrain);
    void   refine();
    void   alpha_wrap(double relativeAlpha, double relativeOffset);
    void   translate_footprint(const double h);
    bool   is_part_of(const Polygon_2& otherPoly) const;
    void   set_reconstruction_rules(const BoundingRegion& reconRegion);
    void   remove_reconstruction_rules();
    bool   has_reconstruction_region() const;
    std::shared_ptr<const Config::ReconRegion> get_reconstruction_settings() const;
    void   set_clip_flag (const bool flag);
    void   mark_as_failed();
    bool   has_failed_to_reconstruct() const;
    bool   has_self_intersections() const;
    void   set_to_zero_terrain();
    double sq_max_dim();
    PointSet3Ptr get_points() const;
    std::string get_lod() const;

    virtual void        get_cityjson_cityobj_info(nlohmann::json& f) const override;
    virtual void        get_cityjson_geomobj_info(nlohmann::json& g) const override;
    virtual void        get_cityjson_semantics(nlohmann::json& g) const override;
    virtual TopoClass   get_class() const override;
    virtual std::string get_class_name() const override;

protected:
    PointSet3Ptr         m_ptsPtr;
    double               m_elevation;
    double               m_height;
    bool                 m_hasFailed;
    bool                 m_clipBottom = Config::get().clip;
    std::shared_ptr<const Config::ReconRegion> m_reconSettings;
};

//-- Struct for clipping
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Exact_predicates_exact_constructions_kernel EK;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;
typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
typedef Mesh::Property_map<vertex_descriptor,EK::Point_3> Exact_point_map;

namespace params = PMP::parameters;

struct Exact_vertex_point_map {
    // typedef for the property map
    typedef boost::property_traits<Exact_point_map>::value_type value_type;
    typedef boost::property_traits<Exact_point_map>::reference reference;
    typedef boost::property_traits<Exact_point_map>::key_type key_type;
    typedef boost::read_write_property_map_tag category;
    // exterior references
    Exact_point_map exact_point_map;
    Mesh* tm_ptr;
    // Converters
    CGAL::Cartesian_converter<K, EK> to_exact;
    CGAL::Cartesian_converter<EK, K> to_input;

    Exact_vertex_point_map()
            : tm_ptr(nullptr) {}

    Exact_vertex_point_map(const Exact_point_map& ep, Mesh& tm)
            : exact_point_map(ep), tm_ptr(&tm) {
        for (Mesh::Vertex_index v: vertices(tm))
            exact_point_map[v] = to_exact(tm.point(v));
    }

    friend reference get(const Exact_vertex_point_map& map, key_type k) {
        CGAL_precondition(map.tm_ptr != nullptr);
        return map.exact_point_map[k];
    }

    friend void put(const Exact_vertex_point_map& map, key_type k, const EK::Point_3& p) {
        CGAL_precondition(map.tm_ptr != nullptr);
        map.exact_point_map[k] = p;
        // create the input point from the exact one
        map.tm_ptr->point(k) = map.to_input(p);
    }
};

#endif //CITY4CFD_BUILDING_H