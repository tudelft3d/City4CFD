#ifndef CITYCFD_BOUNDARY_H
#define CITYCFD_BOUNDARY_H

#include "geomtools.h"
#include "io.h"
#include "TopoFeature.h"
#include "Building.h"

class Boundary : public PolyFeature {
public:
    Boundary();
    Boundary(const int outputLayerID);
    virtual ~Boundary();

    static void set_bounds_to_pc(Point_set_3& pointCloud);
    void  set_bounds_to_cdt(CDT& cdt) const; // Plans to implement it
    static void add_buffer(Point_set_3& pointCloud);
    static std::vector<double> get_domain_bbox();

    virtual void threeDfy() = 0;

    virtual TopoClass   get_class() const = 0;
    virtual std::string get_class_name() const = 0;
    virtual void        get_cityjson_info(nlohmann::json& b) const;
    virtual void        get_cityjson_semantics(nlohmann::json& g) const;
    virtual std::string get_cityjson_primitive() const;

protected:
    static std::vector<Point_3> _outerPts;
    static double               _outerBndHeight;
    static Polygon_2            _bndPoly;
};

class Sides : public Boundary {
public:
    Sides();
    ~Sides();

    virtual void threeDfy() override;

    void set_bnd_poly(double radius, SearchTree& searchTree);

    virtual TopoClass   get_class() const override;
    virtual std::string get_class_name() const override;
};

class Top : public Boundary {
public:
    Top();
    ~Top();

    virtual void threeDfy() override;

    virtual TopoClass   get_class() const override;
    virtual std::string get_class_name() const override;
};

class InfluRegion {
public:
    InfluRegion();
    ~InfluRegion();

    void operator()(double radius);
    void operator()(Polygon_2& poly);
    void operator()(std::string& polyPath);
    void operator()(Point_set_3& pointCloud, Point_set_3& pointCloudBuildings, Buildings& buildings);

    const Polygon_2& get_influ_region() const;

protected:
    Polygon_2 _influRegion;

    double calc_influ_region_radius_bpg(Point_set_3& pointCloud, Point_set_3& pointCloudBuildings, Buildings& buildings);
};

#endif //CITYCFD_BOUNDARY_H