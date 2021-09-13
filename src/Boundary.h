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
    ~Boundary() = default;

    virtual void threeDfy() = 0;

    static void set_bounds_to_pc(Point_set_3& pointCloud);
    static void add_buffer(Point_set_3& pointCloud);
    void set_bounds_to_cdt(CDT& cdt) const;

    TopoClass   get_class() const = 0;
    std::string get_class_name() const = 0;

protected:
    static std::vector<Point_3> _outerPts;
};

class Sides : public Boundary {
public:
    void threeDfy();

    TopoClass   get_class() const override;
    std::string get_class_name() const override;

};

class Top : public Boundary {
    void threeDfy();

    TopoClass   get_class() const override;
    std::string get_class_name() const override;
};

#endif //CITYCFD_BOUNDARY_H