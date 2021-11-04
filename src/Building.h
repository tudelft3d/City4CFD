#ifndef CITYCFD_BUILDING_H
#define CITYCFD_BUILDING_H

#include "PolyFeature.h"

class Building : public PolyFeature {
public:
    Building();
    Building(const nlohmann::json& poly);
    ~Building();

    void   check_feature_scope(const Polygon_2& influRegion);
    double max_dim();
    void   reconstruct(const SearchTree& searchTree);

    double get_height() const;

    virtual void        get_cityjson_info(nlohmann::json& b) const override;
    virtual void        get_cityjson_semantics(nlohmann::json& g) const override;
    virtual std::string get_cityjson_primitive() const override;
    virtual TopoClass   get_class() const override;
    virtual std::string get_class_name() const override;

private:
    double _height;
};

#endif //CITYCFD_BUILDING_H