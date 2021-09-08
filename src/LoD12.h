#ifndef CITYCFD_LOD12_H
#define CITYCFD_LOD12_H

#include "definitions.h"
#include "geomtools.h"

class LoD12 {
public:
    LoD12() = default;
    LoD12(const Polygon_with_holes_2& poly, const std::vector<double>& base_heights, const std::vector<double>& building_pts);
    ~LoD12() = default;

    void   lod12reconstruct(Mesh &mesh);
    double get_height();

private:
    double _height;
    const Polygon_with_holes_2& _poly;
    const std::vector<double>& _base_heights;
    const std::vector<double>& _building_pts;

    void create_mesh(Mesh& mesh);
};


#endif //CITYCFD_LOD12_H
