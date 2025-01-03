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

#ifndef CITY4CFD_CGALTYPES_H
#define CITY4CFD_CGALTYPES_H

/*! CGAL Basics !*/
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Projection_traits_xy_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel  EPICK;
typedef CGAL::Exact_predicates_exact_constructions_kernel    EPECK;
typedef CGAL::Projection_traits_xy_3<EPICK>                  iProjection_traits;
typedef CGAL::Projection_traits_xy_3<EPECK>                  eProjection_traits;

//-- Kernel Converter
template<typename T, typename U>
using Converter = CGAL::Cartesian_converter<T, U>;

/*! CGAL Primitives !*/
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>

//-- CGAL Point
typedef EPICK::Point_2  Point_2;
typedef EPICK::Point_3  Point_3;
typedef EPECK::Point_2  ePoint_2;
typedef EPECK::Point_3  ePoint_3;

//-- CGAL Normal
typedef   EPICK::Vector_3 Vector_3;
typedef   EPICK::Vector_2 Vector_2;

//-- CGAL Polygon
typedef CGAL::Polygon_2<EPICK>              Polygon_2;
typedef CGAL::Polygon_2<eProjection_traits> Polygon_3;

//- CGAL's Polygon_with_holes container expanded
struct Polygon_with_holes_2 {
    Polygon_with_holes_2() = default;

    Polygon_with_holes_2(const CGAL::Polygon_with_holes_2<EPICK>& cgalPoly) {
        m_rings.push_back(cgalPoly.outer_boundary());
        for (auto& hole : cgalPoly.holes()) m_rings.push_back(hole);
    }

    std::vector<Polygon_2> m_rings;

    std::vector<Polygon_2>& rings() {return m_rings;}
    const std::vector<Polygon_2>& rings() const {return m_rings;}

    void push_back(Point_2 point) {
        m_rings.front().push_back(point);
    }
    void push_back(Point_2& point) {
        m_rings.front().push_back(point);
    }
    void insert_ring(const Polygon_2& ring) {
        m_rings.push_back(ring);
    }
    bool has_holes() const {
        if (m_rings.size() > 1) return true;
        return false;
    }
    Polygon_2& outer_boundary() {
        return m_rings.front();
    }
    const Polygon_2& outer_boundary() const {
        return m_rings.front();
    }

    CGAL::Polygon_with_holes_2<EPICK> get_cgal_type() const {
       CGAL::Polygon_with_holes_2<EPICK> cgalPoly;
       cgalPoly.outer_boundary() = m_rings.front();
       for (int i = 1; i < m_rings.size(); ++i) {
           cgalPoly.add_hole(m_rings[i]);
       }
       return cgalPoly;
    }

    CGAL::Polygon_2<EPECK> get_exact_outer_boundary() const {
        Converter<EPICK, EPECK> to_exact;
        CGAL::Polygon_2<EPECK> cgalOuterPoly;
        for (auto& pt : m_rings.front()) {
            cgalOuterPoly.push_back(to_exact(pt));
        }
        return cgalOuterPoly;
    }

    CGAL::Polygon_with_holes_2<EPECK> get_exact() const {
        Converter<EPICK, EPECK> to_exact;
        CGAL::Polygon_with_holes_2<EPECK> cgalPoly;
        cgalPoly.outer_boundary() = get_exact_outer_boundary();
        for (auto hole = holes_begin(); hole != holes_end(); ++hole) {
            CGAL::Polygon_2<EPECK> holePoly;
            for (auto& pt : *hole) {
                holePoly.push_back(to_exact(pt));
            }
            cgalPoly.add_hole(holePoly);
        }
        return cgalPoly;
    }

    std::vector<Polygon_2>::const_iterator holes_begin() const {
        if (has_holes()) return m_rings.begin() + 1; else return m_rings.end();
    }
    std::vector<Polygon_2>::const_iterator holes_end() const {
        return m_rings.end();
    }

    CGAL::Bbox_2 bbox() const {return m_rings.front().bbox();}
};

//- Polygon_with_holes expanded with attributes from GDAL
struct Polygon_with_attr {
    Polygon_with_holes_2 polygon;
    std::unordered_map<std::string, std::string> attributes;
//    std::string semanticClass;
//    std::string id;
};

//-- Minimum bbox
struct MinBbox {
    CGAL::Vector_2<EPICK> vec1;
    CGAL::Vector_2<EPICK> vec2;

    bool empty() const {return is_empty;};
    void calc() {is_empty = false;};
private:
    bool is_empty = true;
};

/*! CGAL Polygon Processing !*/
#include <CGAL/Polygon_2_algorithms.h>

/*! CGAL Point Set !*/
#include <CGAL/Point_set_3.h>
#include <CGAL/random_simplify_point_set.h>

typedef CGAL::Point_set_3<Point_3> Point_set_3;

/*! CGAL Mesh !*/
#include <CGAL/Surface_mesh.h>

typedef CGAL::Surface_mesh<Point_3>                      Mesh;
typedef Mesh::Vertex_index                               vertex_descriptor;
typedef Mesh::Face_index                                 face_descriptor;
typedef Mesh::Property_map<face_descriptor, std::string> Face_property;

namespace PMP = CGAL::Polygon_mesh_processing;

/*! CGAL Triangulation !*/
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_vertex_base_with_id_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>

struct FaceInfo2
{
    FaceInfo2() {}
    int nesting_level;
    bool in_domain() {
        return nesting_level%2 == 1;
    }
    bool in_domain_noholes() {
        return nesting_level > 0;
    }
    int surfaceLayer = -9999; // Face handle to output mesh for specific surface layer
};
typedef CGAL::Triangulation_vertex_base_with_id_2<eProjection_traits>              Vb;
typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo2, eProjection_traits>   Fbb;
typedef CGAL::Constrained_triangulation_face_base_2<eProjection_traits, Fbb>       Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>                                TDS;
typedef CGAL::Exact_intersections_tag                                              Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<eProjection_traits, TDS, Itag>  CDTt;
typedef CGAL::Constrained_triangulation_plus_2<CDTt>                               CDT;
typedef CDT::Point                                                                 Point;
typedef CDT::Face_handle                                                           Face_handle;
typedef CDT::Vertex_handle                                                         Vertex_handle;
//- Different versions of DT
typedef CGAL::Delaunay_triangulation_2<CGAL::Projection_traits_xy_3<EPICK>>        DT;

/*! CGAL Search !*/
#include <boost/graph/adjacency_list.hpp>
#include <CGAL/boost/graph/split_graph_into_polylines.h>

#include <CGAL/Kd_tree.h>
#include <CGAL/algorithm.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Fuzzy_iso_box.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>

// Point_3 to Point_2 projection on the fly
template<class K>
struct Projection_xy_property_map
{
    typedef typename K::Point_3 key_type;
    typedef typename K::Point_2 value_type;
    typedef value_type reference;
    typedef boost::readable_property_map_tag category;

    friend value_type get(Projection_xy_property_map<K>, const key_type& k)
    {
        return value_type(k.x(), k.y());
    }
};

typedef CGAL::Search_traits_2<EPICK>                                   Traits_base;
typedef CGAL::Search_traits_adapter<Point_3,
                                    Projection_xy_property_map<EPICK>,
                                    Traits_base>                       Traits;
typedef CGAL::Orthogonal_k_neighbor_search<Traits>                     Neighbor_search;
typedef Neighbor_search::Tree                                          SearchTree;
typedef CGAL::Fuzzy_iso_box<Traits>                                    Fuzzy_iso_box;

//-- Smart pointers with CGAL types
typedef std::shared_ptr<Point_set_3>                                  PointSet3Ptr;
typedef std::vector<std::unique_ptr<Polygon_with_attr>>               PolyVecPtr;

//-- Global constants
namespace global {
    const Point_2 nullPt(0, 0);
}

#endif //CITY4CFD_CGALTYPES_H