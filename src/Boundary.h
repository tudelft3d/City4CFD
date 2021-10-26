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

    static void set_bnd_poly(Polygon_2& bndPoly, Polygon_2& pcBndPoly, Polygon_2& startBufferPoly);
    static void set_bounds_to_pc(Point_set_3& pointCloud, const Polygon_2& pcBndPoly);
    static void set_bounds_to_terrain(Point_set_3& pointCloud, const Polygon_2& bndPoly,
                                      const Polygon_2& pcBndPoly, const Polygon_2& startBufferPoly);
    static std::vector<double> get_domain_bbox();

    virtual void reconstruct() = 0;

    void prep_output();
    void prep_output(Vector_2 edge);

    virtual TopoClass   get_class() const = 0;
    virtual std::string get_class_name() const = 0;
    virtual void        get_cityjson_info(nlohmann::json& b) const;
    virtual void        get_cityjson_semantics(nlohmann::json& g) const;
    virtual std::string get_cityjson_primitive() const;

protected:
    static std::vector<Point_3> _outerPts;
    static double               _outerBndHeight;
    std::vector<Point_3>        _sideOutputPts;
};

class Sides : public Boundary {
public:
    Sides(const int outputLayerID);
    ~Sides();

    virtual void reconstruct() override;

    virtual TopoClass   get_class() const override;
    virtual std::string get_class_name() const override;
};

class Top : public Boundary {
public:
    Top(const int outputLayerID);
    ~Top();

    virtual void reconstruct() override;

    virtual TopoClass   get_class() const override;
    virtual std::string get_class_name() const override;
};

class BoundingRegion {
public:
    BoundingRegion();
    ~BoundingRegion();

    void operator()(double radius);
    void operator()(Polygon_2& poly);
    void operator()(std::string& polyPath);

    void calc_influ_region_bpg(const DT& dt, const Point_set_3& pointCloudBuildings, Buildings& buildings);
    void calc_bnd_bpg(double hMax, const Polygon_2& influRegionPoly, const Buildings& buildings);

    Polygon_2& get_bounding_region();
    const Polygon_2& get_bounding_region() const;

protected:
    Polygon_2    _boundingRegion;

};

#endif //CITYCFD_BOUNDARY_H