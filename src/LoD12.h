#ifndef CITY4CFD_LOD12_H
#define CITY4CFD_LOD12_H

#include "CGALTypes.h"

class LoD12 {
public:
    LoD12() = default;
    LoD12(const Polygon_with_holes_2& poly, const std::vector<std::vector<double>>& base_heights, const std::vector<double>& building_pts);
    ~LoD12() = default;

    void   lod12reconstruct(Mesh &mesh, double& height);
    double get_height() const;

private:
    double _height;
    const Polygon_with_holes_2&              _poly;
    const std::vector<std::vector<double>>&  _baseHeights;
    const std::vector<double>&               _buildingPts;

    void create_mesh(Mesh& mesh);
};

#endif //CITY4CFD_LOD12_H