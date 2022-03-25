#ifndef CITY4CFD_RECONSTRUCTEDBUILDING_H
#define CITY4CFD_RECONSTRUCTEDBUILDING_H

#include "Building.h"

class ReconstructedBuilding : public Building {
public:
    ReconstructedBuilding();
    ReconstructedBuilding(const int internalID);
    ReconstructedBuilding(const nlohmann::json& poly);
    ReconstructedBuilding(const nlohmann::json& poly, const int internalID);
    ~ReconstructedBuilding();

    void   set_search_tree(const std::shared_ptr<SearchTree>& searchTree);

    virtual void  reconstruct() override;
    virtual void  get_cityjson_info(nlohmann::json& b) const override;
    virtual void  get_cityjson_semantics(nlohmann::json& g) const override;

protected:
    std::shared_ptr<SearchTree> _searchTree;
};


#endif //CITY4CFD_RECONSTRUCTEDBUILDING_H