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

    static void set_bounds_to_pc(Point_set_3& pointCloud, const Polygon_2& bndPoly);
    static std::vector<double> get_domain_bbox();

    virtual void reconstruct() = 0;

    void set_bnd_poly(const Polygon_2& bndPoly, Point_set_3& pointCloud);

    virtual TopoClass   get_class() const = 0;
    virtual std::string get_class_name() const = 0;
    virtual void        get_cityjson_info(nlohmann::json& b) const;
    virtual void        get_cityjson_semantics(nlohmann::json& g) const;
    virtual std::string get_cityjson_primitive() const;

protected:
    static std::vector<Point_3> _outerPts;
    static double               _outerBndHeight;
};

class Sides : public Boundary {
public:
    Sides();
    ~Sides();

    virtual void reconstruct() override;

    virtual TopoClass   get_class() const override;
    virtual std::string get_class_name() const override;
};

class Top : public Boundary {
public:
    Top();
    ~Top();

    virtual void reconstruct() override;

    virtual TopoClass   get_class() const override;
    virtual std::string get_class_name() const override;
};

class BoundedRegion {
public:
    BoundedRegion();
    ~BoundedRegion();

    void operator()(double radius);
    void operator()(Polygon_2& poly);
    void operator()(std::string& polyPath);
    void operator()(Point_set_3& pointCloud, Point_set_3& pointCloudBuildings, Buildings& buildings);

    const Polygon_2& get_bounded_region() const;

protected:
    Polygon_2 _boundedRegion;

    double calc_influ_region_radius_bpg(Point_set_3& pointCloud, Point_set_3& pointCloudBuildings, Buildings& buildings);
};

#endif //CITYCFD_BOUNDARY_H