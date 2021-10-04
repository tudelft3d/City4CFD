/*
 Copyright (c) 2014,
 Ken Arroyo Ohori    g.a.k.arroyoohori@tudelft.nl
 All rights reserved.
 
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 SOFTWARE.
 */

#ifndef ENHANCED_TRIANGULATION_2_H
#define ENHANCED_TRIANGULATION_2_H

template <class T>
class Enhanced_constrained_triangulation_2;

template <class T>
class Is_Delaunay {
public:
  typedef T Triangulation;
  typedef typename Triangulation::List_edges List_edges;
  
  static void if_Delaunay_make_Delaunay(Enhanced_constrained_triangulation_2<T> &et, List_edges &e) {
    std::cout << "Unknown triangulation: cannot make Delaunay" << std::endl;
  }
};

template <class K, class TDS, class Itag>
class Is_Delaunay<CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag> > {
public:
  typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag> Triangulation;
  typedef typename Triangulation::List_edges List_edges;
  
  static void if_Delaunay_make_Delaunay(Enhanced_constrained_triangulation_2<Triangulation> &et, List_edges &e) {
    et.propagating_flip(e);
  }
};

template <class K, class TDS, class Itag>
class Is_Delaunay<CGAL::Constrained_triangulation_2<K, TDS, Itag> > {
public:
  typedef CGAL::Constrained_triangulation_2<K, TDS, Itag> Triangulation;
  typedef typename Triangulation::List_edges List_edges;
  
  static void if_Delaunay_make_Delaunay(Enhanced_constrained_triangulation_2<Triangulation> &et, List_edges &e) {
    // Do nothing
  }
};

template <class T>
class Enhanced_constrained_triangulation_2 : public T {
public:
  typedef typename T::Point Point;
  typedef typename T::Vertex Vertex;
  typedef typename T::Edge Edge;
  typedef typename T::Vertex_handle Vertex_handle;
  typedef typename T::Face_handle Face_handle;
  typedef typename T::Locate_type Locate_type;
  typedef typename T::List_faces List_faces;
  typedef typename T::List_edges List_edges;
  
  Vertex_handle insert(const Point &p, Face_handle f = Face_handle()) {
    // std::cout << "Enhanced_triangulation_2::insert(const Point &, Face_handle)" << std::endl;
    Locate_type location_type;
    int location_vertex;
    Face_handle location = T::locate(p, location_type, location_vertex, f);
    return insert(p, location_type, location, location_vertex);
  }
  
  Vertex_handle insert(const Point &p, Locate_type &lt, Face_handle loc, int li) {
    // std::cout << "Enhanced_triangulation_2::insert(const Point &, Locate_type &, Face_handle, int)" << std::endl;
    // TODO: Read and reset info on faces & edges that are split because of the insertion
    return T::insert(p, lt, loc, li);
  }
  
  void odd_even_insert_constraint(const Point& a, const Point& b) {
    Vertex_handle va = insert(a);
    Vertex_handle vb = insert(b);
    if (va != vb) odd_even_insert_constraint(va, vb);
  }
  
  void odd_even_insert_constraint(Vertex_handle va, Vertex_handle vb) {
    CGAL_triangulation_precondition(va != vb);
    
    // If [va, vb] lies on an existing edge
    Vertex_handle vertex_on_other_end;
    Face_handle incident_face;
    int vertex_opposite_to_incident_edge;
    if (T::includes_edge(va, vb, vertex_on_other_end, incident_face, vertex_opposite_to_incident_edge)) {
      if (T::is_constrained(Edge(incident_face, vertex_opposite_to_incident_edge))) {
        T::remove_constrained_edge(incident_face, vertex_opposite_to_incident_edge);
        List_edges possibly_non_Delaunay_edges;
        possibly_non_Delaunay_edges.push_back(Edge(incident_face, vertex_opposite_to_incident_edge));
        Is_Delaunay<T>::if_Delaunay_make_Delaunay(*this, possibly_non_Delaunay_edges);
      } else T::mark_constraint(incident_face, vertex_opposite_to_incident_edge);
      if (vertex_on_other_end != vb) odd_even_insert_constraint(vertex_on_other_end, vb);
      return;
    }
    
    // If [va, vb] intersects a constrained edge or an existing vertex
    List_faces intersected_faces;
    List_edges conflict_boundary_ab, conflict_boundary_ba;
    Vertex_handle intersection;
    if (T::find_intersected_faces(va, vb, intersected_faces, conflict_boundary_ab, conflict_boundary_ba, intersection)) {
      if (intersection != va && intersection != vb) {
        odd_even_insert_constraint(va, intersection);
        odd_even_insert_constraint(intersection, vb);
      } else odd_even_insert_constraint(va, vb);
      return;
    }
    
    // Otherwise
    T::triangulate_hole(intersected_faces, conflict_boundary_ab, conflict_boundary_ba);
    if (intersection != vb) {
      odd_even_insert_constraint(intersection, vb);
    }
  }
};

#endif
