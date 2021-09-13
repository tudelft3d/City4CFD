#include "Boundary.h"

// TODO: put more thought into this constructor once the time comes
/*
Boundary::Boundary()
    :TopoFeature(), _outerPts() {
    if (configData.dimOfDomain != -infty) {
        _dimOfDomain = configData.dimOfDomain;
    }

    if (configData.topHeight != -infty) {
        _topHeight = configData.topHeight;
    }
}
*/

//-- Static member definition
std::vector<Point_3> Boundary::_outerPts;

//-- Deactivate point cloud points that are out of bounds
void Boundary::set_bounds_to_pc(Point_set_3& pointCloud) { // Will try to template it to include CDT
    //- 80% of the total domain size. The rest is left for the buffer zone
    auto it = pointCloud.points().begin();
    int count = 0;
    while (it != pointCloud.points().end()) {
        if (!geomtools::point_in_circle(*it, config::pointOfInterest, 0.8 * config::dimOfDomain)) {
            pointCloud.remove(pointCloud.begin() + count);
        } else {
            ++it;
            ++count;
        }
    }
    pointCloud.collect_garbage(); // Free removed points from the memory
}

//-- Add buffer between the terrain and domain edge
void Boundary::add_buffer(Point_set_3& pointCloud) {
    const int nPts      = 360; // Hardcoded
    const double angInt = 2 * M_PI / (double)nPts;
    double ang = 0;
    for (auto i = 0; i < nPts; ++i) {
        double xPt = config::pointOfInterest.x() + config::dimOfDomain * cos(ang + angInt);
        double yPt = config::pointOfInterest.y() + config::dimOfDomain * sin(ang + angInt);
        ang = ang + angInt;
        Point_3 pt(xPt, yPt, 0.0); // Height set at 0 for now. Could be average of the edge or average of the whole terrain
        _outerPts.push_back(pt); // for top and sides
        pointCloud.insert(pt);
    }
    _outerPts.push_back(_outerPts[0]); // Put first point at the end to close the loop
}

//-- Sides class
void Sides::threeDfy() {
    std::vector<Mesh::vertex_index> mesh_vertex_side;

    //-- Add mesh vertices and store them in a vector
    for (auto it = _outerPts.begin(); it != _outerPts.end(); ++it) {
        mesh_vertex_side.emplace_back(_mesh.add_vertex(*it));
        mesh_vertex_side.emplace_back(_mesh.add_vertex(Point_3(it->x(), it->y(), config::topHeight)));
    }

    //-- Add middle top point to mesh
//    Mesh::vertex_index topMiddlePoint = _meshTop.add_vertex(Point_3(bndInfo.xcent, bndInfo.ycent, bndInfo.height));

    //-- Add mesh faces for side
    for (auto i = 0; i < mesh_vertex_side.size() - 3; i= i + 2) {
        // -- i + 1 is i lifted up
        int v1 = i;
        int v2 = i + 2;

        _mesh.add_face(mesh_vertex_side[v2], mesh_vertex_side[v1], mesh_vertex_side[v1 + 1]);
        _mesh.add_face(mesh_vertex_side[v2 + 1], mesh_vertex_side[v2], mesh_vertex_side[v1 + 1]);
    }
}

TopoClass Sides::get_class() const {
    return SIDES;
}

std::string Sides::get_class_name() const {
    return "Sides";
}

//-- Top class
void Top::threeDfy() {
    std::vector<Mesh::vertex_index> mesh_vertex_top;

    //-- Top is done by making a CDT of outerPts
    CDT cdt_top;
    for (auto &pt : _outerPts) {
        cdt_top.insert(Point_3(pt.x(), pt.y(), config::topHeight));
    }

    //-- Add mesh faces for top
    geomtools::cdt_to_mesh(cdt_top, _mesh);
}

TopoClass Top::get_class() const {
    return TOP;
}

std::string Top::get_class_name() const {
    return "Top";
}
