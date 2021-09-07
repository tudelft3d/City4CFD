#ifndef CITYCFD_DEFINITIONS_H
#define CITYCFD_DEFINITIONS_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_vertex_base_with_id_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/boost/graph/graph_traits_Delaunay_triangulation_2.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>
#include <CGAL/compute_average_spacing.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/locate.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <boost/graph/adjacency_list.hpp>
#include <CGAL/boost/graph/split_graph_into_polylines.h>
#include <CGAL/IO/WKT.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>

#include <CGAL/Kd_tree.h>
#include <CGAL/algorithm.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Fuzzy_iso_box.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>

#include <CGAL/Point_set_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/draw_polygon_with_holes_2.h> // Test to draw
#include <CGAL/Qt/Basic_viewer_qt.h> // Test to draw

#include <boost/algorithm/string.hpp>

#include "nlohmann/json.hpp"

//-- CGAL basics
typedef CGAL::Exact_predicates_inexact_constructions_kernel  Kernel;
typedef CGAL::Projection_traits_xy_3<Kernel>                 Projection_traits;

//-- CGAL point
typedef Kernel::Point_2            Point_2;
typedef Kernel::Point_3            Point_3;
typedef Kernel::Segment_3          Segment_3;
typedef CGAL::Point_set_3<Point_3> Point_set_3;

//-- CGAL tree search
typedef CGAL::Search_traits_3<Kernel>                 Traits;
//typedef CGAL::Kd_tree<Traits>                         SearchTree;
typedef CGAL::Orthogonal_k_neighbor_search<Traits>    Neighbor_search;
typedef Neighbor_search::Tree                         SearchTree;
typedef CGAL::Fuzzy_iso_box<Traits> Fuzzy_iso_box;

//-- CGAL CDT
struct FaceInfo2
{
    FaceInfo2(){}
    int nesting_level;
    bool in_domain(){
        return nesting_level%2 == 1;
    }
};

typedef CGAL::Triangulation_vertex_base_with_id_2<Projection_traits>                Vb;
typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo2, Projection_traits>     Fbb;
typedef CGAL::Constrained_triangulation_face_base_2<Projection_traits, Fbb>         Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>                                 TDS;
typedef CGAL::Exact_predicates_tag                                                  Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<Projection_traits, TDS, Itag>    CDT;
typedef CDT::Point                                                                  Point;
typedef CGAL::Polygon_2<Kernel>                                                     Polygon_2;
typedef CGAL::Polygon_with_holes_2<Kernel>                                          Polygon_with_holes_2;
typedef CDT::Face_handle                                                            Face_handle;
typedef CDT::Vertex_handle                                                          Vertex_handle;

//-- tesh
using Mesh = CGAL::Surface_mesh<Point_3>;

typedef nlohmann::json  json;

struct vertex {
    double x;
    double y;
    double z;
};

struct Triangle {
    int v0;
    int v1;
    int v2;
};

typedef enum {
    BUILDING    = 0,
    WATER       = 1,
    BRIDGE      = 2,
    ROAD        = 3,
    TERRAIN     = 4,
    FOREST      = 5,
    BOUNDARY    = 6
} TopoClass;

const double infty     = 1e6;
const double smallnum  = 1e-6;

#endif //CITYCFD_DEFINITIONS_H
