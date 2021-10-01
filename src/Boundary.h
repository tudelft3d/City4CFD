#ifndef CITYCFD_BOUNDARY_H
#define CITYCFD_BOUNDARY_H

#include "config.h"
#include "geomtools.h"
#include "TopoFeature.h"

class Boundary : public TopoFeature {
public:
    Boundary();
    Boundary(const int outputLayerID);
    virtual ~Boundary();

    virtual void threeDfy() = 0;

    static void set_bounds_to_pc(Point_set_3& pointCloud);
    static void add_buffer(Point_set_3& pointCloud);
    static std::vector<double> get_domain_bbox ();

    void set_bounds_to_cdt(CDT& cdt) const;
    TopoClass   get_class() const = 0;
    std::string get_class_name() const = 0;

protected:
    static std::vector<Point_3> _outerPts;
};

class Sides : public Boundary {
public:
    Sides();
    ~Sides();

    void threeDfy();

    TopoClass   get_class() const override;
    std::string get_class_name() const override;

};

class Top : public Boundary {
public:
    Top();
    ~Top();

    void threeDfy();

    TopoClass   get_class() const override;
    std::string get_class_name() const override;
};

#endif //CITYCFD_BOUNDARY_H