/*
  Copyright (c) 2021-2022,
  Ivan Pađen <i.paden@tudelft.nl>
  3D Geoinformation,
  Delft University of Technology

  This file is part of City4CFD.

  City4CFD is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  City4CFD is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

  You should have received a copy of the GNU General Public License
  along with this program. If not, see <http://www.gnu.org/licenses/>
*/

#ifndef CITY4CFD_BUILDING_H
#define CITY4CFD_BUILDING_H

#include "PolyFeature.h"

class Building : public PolyFeature {
public:
    Building();
    Building(const int internalID);
    Building(const nlohmann::json& poly);
    Building(const nlohmann::json& poly, const int internalID);
    ~Building();

    virtual void reconstruct() = 0;

    void   clip_bottom(const Terrainptr& terrain);
    void   translate_footprint(const double h);
    void   check_feature_scope(const Polygon_2& influRegion);
    double max_dim();

    double get_height() const;

    virtual void        get_cityjson_info(nlohmann::json& b) const override;
    virtual void        get_cityjson_semantics(nlohmann::json& g) const override;
    virtual std::string get_cityjson_primitive() const override;
    virtual TopoClass   get_class() const override;
    virtual std::string get_class_name() const override;

protected:
    double _height;
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