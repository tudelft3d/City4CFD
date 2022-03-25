#ifndef CITY4CFD_BUILDING_H
#define CITY4CFD_BUILDING_H

#include "PolyFeature.h"

class Building : public PolyFeature {
public:
    Building();
    Building(const int internalID);
    Building(const nlohmann::json& poly);
    Building(const nlohmann::json& poly, const int internalID);
    ~Building();

    virtual void reconstruct() = 0;

    void   check_feature_scope(const Polygon_2& influRegion);
    double max_dim();

    double get_height() const;

    virtual void        get_cityjson_info(nlohmann::json& b) const override;
    virtual void        get_cityjson_semantics(nlohmann::json& g) const override;
    virtual std::string get_cityjson_primitive() const override;
    virtual TopoClass   get_class() const override;
    virtual std::string get_class_name() const override;

protected:
    double _height;
};

#endif //CITY4CFD_BUILDING_H