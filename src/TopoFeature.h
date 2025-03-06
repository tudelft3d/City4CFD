/*
  City4CFD
 
  Copyright (c) 2021-2025, 3D Geoinformation Research Group, TU Delft

  This file is part of City4CFD.

  City4CFD is free software: you can redistribute it and/or modify
  it under the terms of the GNU Affero General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  City4CFD is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Affero General Public License for more details.

  You should have received a copy of the GNU Affero General Public License
  along with City4CFD.  If not, see <http://www.gnu.org/licenses/>.

  For any information or further details about the use of City4CFD, contact
  Ivan Pađen
  <i.paden@tudelft.nl>
  3D Geoinformation Research Group
  Delft University of Technology
*/

#ifndef CITY4CFD_TOPOFEATURE_H
#define CITY4CFD_TOPOFEATURE_H

#include "types.h"
#include "CGALTypes.h"

class TopoFeature {
public:
    TopoFeature();
    TopoFeature(std::string pid);
    TopoFeature(int outputLayerID);
    virtual ~TopoFeature() = default;

    static void          add_recon_region_output_layers(const int numLayers);
    static int           get_num_output_layers();

    virtual TopoClass    get_class() const = 0;
    virtual std::string  get_class_name() const = 0;

    virtual void         get_cityjson_cityobj_info(nlohmann::json& f) const;
    virtual void         get_cityjson_geomobj_info(nlohmann::json& g) const;
    virtual void         get_cityjson_semantics(nlohmann::json& g) const;

    Mesh&       get_mesh();
    const Mesh& get_mesh() const;
    void        set_id(unsigned long id);
    std::string get_id() const;
    int         get_output_layer_id() const;
    bool        is_active() const;
    bool        is_imported() const;
    void        deactivate();

protected:
    static int s_numOfOutputLayers;

    Mesh           m_mesh;
    std::string    m_id;
    bool           m_f_active;
    bool           m_f_imported;
    int            m_outputLayerID; // 0  Terrain
                                    //    Building regions
                                    //    Sides
                                    //    Top
                                    //    Surface Layers
};

#endif //CITY4CFD_TOPOFEATURE_H