#ifndef CITYCFD_CGALTYPES_H
#define CITYCFD_CGALTYPES_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/IO/WKT.h>
#include <CGAL/centroid.h>

//-- CGAL Basics
typedef CGAL::Exact_predicates_inexact_constructions_kernel  EPICK;
typedef CGAL::Exact_predicates_exact_constructions_kernel    EPECK;
typedef CGAL::Projection_traits_xy_3<EPICK>                  iProjection_traits;
typedef CGAL::Projection_traits_xy_3<EPECK>                  Projection_traits;

//-- Kernel Converter
template<typename T, typename U>
using Converter = CGAL::Cartesian_converter<T, U>;

//-- CGAL Primitives
#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/convex_hull_2.h>

//-- CGAL Point
typedef EPICK::Point_2             Point_2;
typedef EPICK::Point_3             Point_3;
typedef EPECK::Point_3             ePoint_3;
typedef CGAL::Point_set_3<Point_3> Point_set_3;

//-- CGAL Normal
namespace PMP = CGAL::Polygon_mesh_processing;
typedef   EPICK::Vector_3 Vector_3;
typedef   EPICK::Vector_2 Vector_2;

//-- CGAL Polygon
typedef CGAL::Polygon_2<EPICK>             Polygon_2;
typedef CGAL::Polygon_2<Projection_traits> Polygon_3;

//-- CGAL's Polygon_with_holes container expanded with iterator over all rings
struct Polygon_with_holes_2 {
    std::vector<Polygon_2> _rings;

    std::vector<Polygon_2>& rings() {return _rings;}
    const std::vector<Polygon_2>& rings() const {return _rings;}

    void push_back(Point_2 point) {
        _rings.front().push_back(point);
    }
    void push_back(Point_2& point) {
        _rings.front().push_back(point);
    }
    bool has_holes() const {
        if (_rings.size() > 1) return true;
        return false;
    }
    Polygon_2& outer_boundary() {
        return _rings.front();
    }
    const Polygon_2& outer_boundary() const {
        return _rings.front();
    }

    std::vector<Polygon_2>::const_iterator holes_begin() const {
        if (has_holes()) return _rings.begin() + 1; else return _rings.end();
    }
    std::vector<Polygon_2>::const_iterator holes_end() const {
        return _rings.end();
    }

    CGAL::Bbox_2 bbox() const {return _rings.front().bbox();}
};

//-- CGAL Mesh
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>

typedef CGAL::Surface_mesh<Point_3>                      Mesh;
typedef Mesh::Vertex_index                               vertex_descriptor;
typedef Mesh::Face_index                                 face_descriptor;
typedef Mesh::Property_map<face_descriptor, std::string> Face_property;

namespace PMP = CGAL::Polygon_mesh_processing;

//-- CGAL Triangulation
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
typedef CGAL::Triangulation_vertex_base_with_id_2<Projection_traits>               Vb;
typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo2, Projection_traits>    Fbb;
typedef CGAL::Constrained_triangulation_face_base_2<Projection_traits, Fbb>        Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>                                TDS;
typedef CGAL::Exact_intersections_tag                                              Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<Projection_traits, TDS, Itag>   CDTt;
typedef CGAL::Constrained_triangulation_plus_2<CDTt>                               CDT;
typedef CDT::Point                                                                 Point;
typedef CDT::Face_handle                                                           Face_handle;
typedef CDT::Vertex_handle                                                         Vertex_handle;
//- Different versions of DT
typedef CGAL::Delaunay_triangulation_2<CGAL::Projection_traits_xy_3<EPICK>>        DT;

//-- CGAL Search
#include <boost/graph/adjacency_list.hpp>
#include <CGAL/boost/graph/split_graph_into_polylines.h>
#include <CGAL/compute_average_spacing.h>

#include <CGAL/Kd_tree.h>
#include <CGAL/algorithm.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Fuzzy_iso_box.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/natural_neighbor_coordinates_2.h>
#include <CGAL/Barycentric_coordinates_2/Triangle_coordinates_2.h>

typedef CGAL::Search_traits_3<EPICK>                 Traits;
//typedef CGAL::Kd_tree<Traits>                         SearchTree;
typedef CGAL::Orthogonal_k_neighbor_search<Traits>    Neighbor_search;
typedef Neighbor_search::Tree                         SearchTree;
typedef CGAL::Fuzzy_iso_box<Traits>                   Fuzzy_iso_box;
typedef CGAL::Fuzzy_sphere<Traits>                    Fuzzy_sphere;

#endif //CITYCFD_CGALTYPES_H