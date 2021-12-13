#ifndef CITYCFD_EXPLICITBUILDING_H
#define CITYCFD_EXPLICITBUILDING_H

#include "Building.h"

class ImportedBuilding : public Building {
public:
    ImportedBuilding() = delete;
    ImportedBuilding(nlohmann::json  poly, std::vector<Point_3>& importedBuildingPts, const int internalID);
    ~ImportedBuilding();

    virtual void  reconstruct() override;

    void append_nonground_part(const std::shared_ptr<ImportedBuilding>& other);

    const nlohmann::json& get_building_json() const;
    const std::string&    get_parent_building_id() const;
    const int             get_lod_idx() const;
    const bool            is_appending() const;

//    virtual void  get_cityjson_info(nlohmann::json& b) const override;
//    virtual void  get_cityjson_semantics(nlohmann::json& g) const override;

protected:
    nlohmann::json                   _buildingJson;
    double                           _avgFootprintHeight;
    std::vector<int>                 _footprintIdxList;
    std::vector<std::vector<int>>    _footprintPtsIdxList;
    std::string                      _parentBuildingID;
    bool                             _appendToBuilding;
    int                              _lodIdx;
    std::vector<Point_3>&            _dPts;
};

#endif //CITYCFD_EXPLICITBUILDING_H