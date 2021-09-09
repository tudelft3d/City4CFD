#ifndef CITYCFD_BOUNDARY_H
#define CITYCFD_BOUNDARY_H

#include "definitions.h"
#include "geomtools.h"
#include "config.h"
#include "io.h"
#include "TopoFeature.h"

class Boundary : public TopoFeature {
public:
    Boundary() = default;
    //TODO: get domain info from the config file
    ~Boundary() = default;

    static void set_bounds_to_pc(Point_set_3& pointCloud);
    void set_bounds_to_cdt(CDT& cdt) const;
    void add_buffer(Point_set_3& pointCloud);
    void threeDfy();

    const Mesh& get_top_mesh() const;

    TopoClass get_class() const override;

private:
    std::vector<Point_3> _outerPts;
    Mesh                 _meshTop; // _mesh from TopoFeature is for the sides
};

#endif //CITYCFD_BOUNDARY_H