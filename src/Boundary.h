#ifndef CITY4CFD_BOUNDARY_H
#define CITY4CFD_BOUNDARY_H

#include "TopoFeature.h"

class Boundary : public TopoFeature {
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
#endif //CITY4CFD_BOUNDARY_H