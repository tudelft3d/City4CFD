#ifndef CITYCFD_POLYFEATURE_H
#define CITYCFD_POLYFEATURE_H

#include "TopoFeature.h"

class PolyFeature : public TopoFeature {
public:
    PolyFeature();
    PolyFeature(const int outputLayerID);
    PolyFeature(const nlohmann::json& poly);
    PolyFeature(const nlohmann::json& poly, const int outputLayerID);
    virtual ~PolyFeature();

    void  calc_footprint_elevation_nni(const DT& dt);
    void  calc_footprint_elevation_linear(const DT& dt);
    void  clear_feature();

    virtual void        get_cityjson_info(nlohmann::json& b) const = 0;
    virtual void        get_cityjson_semantics(nlohmann::json& g) const = 0;
    virtual std::string get_cityjson_primitive() const = 0;
    virtual TopoClass   get_class() const = 0;
    virtual std::string get_class_name() const = 0;

    Polygon_with_holes_2&                    get_poly();
    const Polygon_with_holes_2&              get_poly() const;
    const std::vector<std::vector<double>>&  get_base_heights() const;

protected:
    Polygon_with_holes_2              _poly;
    std::vector<std::vector<double>>  _base_heights;
};

#endif //CITYCFD_POLYFEATURE_H