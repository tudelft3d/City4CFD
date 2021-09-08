#include "geomtools.h"

double avg(const std::vector<double>& values) {
    if (values.size() == 0)  std::cerr <<"Can't calculate average of a zero-sized vector!" << std::endl;
    double average = 0;
    for (auto &value : values) {
        average += value;
    }
    return (average / (double)values.size());
}

//-- Check if the point is inside a polygon on a 2D projection
bool check_inside(const Point_3& pt2, const Polygon_with_holes_2& polygon) {
    Point_2 pt(pt2.x(), pt2.y());

    //-- Check if the point falls within the outer surface
    if (CGAL::bounded_side_2(polygon.outer_boundary().begin(), polygon.outer_boundary().end(), pt) == CGAL::ON_BOUNDED_SIDE) {
        // Check if the point falls within one of the holes - TODO: do I even need this?
        for (auto it_hole = polygon.holes_begin(); it_hole != polygon.holes_end(); ++it_hole) {
            if (!CGAL::bounded_side_2(it_hole->begin(), it_hole->end(), pt ) == CGAL::ON_UNBOUNDED_SIDE) {
                return false;
            }
        }
        return true;
    }
    return false;
}

bool point_in_circle(const Point_3& pt, const Point_2& center, const double& radius) {
    if (pow(pt.x() - center.x(), 2)
      + pow(pt.y() - center.y(), 2)
      < pow(radius, 2)) {
        return true;
    }
    return false;
}

void cdt_to_mesh(const CDT& cdt, Mesh& mesh) {
    std::map<CDT::Vertex_handle, int> indices;
    std::vector<Mesh::vertex_index> mesh_vertex;
    std::vector<Mesh::face_index> face_index;
    mesh_vertex.reserve(cdt.dimension());

    int counter = 0;
    for (auto it : cdt.finite_vertex_handles()) {
        mesh_vertex.emplace_back(mesh.add_vertex(it->point()));
        //        outstream << it->point() << std::endl;
        indices.insert(std::pair<CDT::Vertex_handle, int>(it, counter++));
    }

    for (auto it : cdt.finite_face_handles()) {
        int v1 = indices[it->vertex(0)];
        int v2 = indices[it->vertex(1)];
        int v3 = indices[it->vertex(2)];
        mesh.add_face(mesh_vertex[v1], mesh_vertex[v2], mesh_vertex[v3]);
    }
}

void mark_domains(CDT& ct,
             Face_handle start,
             int index,
             std::list<CDT::Edge>& border )
{
    if(start->info().nesting_level != -1){
        return;
    }
    std::list<Face_handle> queue;
    queue.push_back(start);
    while(! queue.empty()){
        Face_handle fh = queue.front();
        queue.pop_front();
        if(fh->info().nesting_level == -1){
            fh->info().nesting_level = index;
            for(int i = 0; i < 3; i++){
                CDT::Edge e(fh,i);
                Face_handle n = fh->neighbor(i);
                if(n->info().nesting_level == -1){
                    if(ct.is_constrained(e)) border.push_back(e);
                    else queue.push_back(n);
                }
            }
        }
    }
}

void mark_domains(CDT& cdt) {
    for(CDT::Face_handle f : cdt.all_face_handles()){
        f->info().nesting_level = -1;
    }
    std::list<CDT::Edge> border;
    mark_domains(cdt, cdt.infinite_face(), 0, border);
    while(! border.empty()){
        CDT::Edge e = border.front();
        border.pop_front();
        Face_handle n = e.first->neighbor(e.second);
        if(n->info().nesting_level == -1){
            mark_domains(cdt, n, e.first->info().nesting_level+1, border);
        }
    }
}