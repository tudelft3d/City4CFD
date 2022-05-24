/*
  Copyright (c) 2021-2022,
  Ivan Pađen <i.paden@tudelft.nl>
  3D Geoinformation,
  Delft University of Technology

  This file is part of City4CFD.

  City4CFD is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  City4CFD is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

  You should have received a copy of the GNU General Public License
  along with this program. If not, see <http://www.gnu.org/licenses/>
*/

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

    void  set_search_tree(const std::shared_ptr<SearchTree>& searchTree);

    virtual void  reconstruct() override;
    virtual void  reconstruct_flat_terrain() override;
    virtual void  get_cityjson_info(nlohmann::json& b) const override;
    virtual void  get_cityjson_semantics(nlohmann::json& g) const override;

protected:
    std::shared_ptr<SearchTree> _searchTree;
    double _attributeHeight;
    bool   _attributeHeightAdvantage;

    //-- Hardcoded low height info
    double _lowHeight = 2.0;

    void reconstruct_from_attribute();
    bool reconstruct_again_from_attribute(const std::string& reason);
};


#endif //CITY4CFD_RECONSTRUCTEDBUILDING_H