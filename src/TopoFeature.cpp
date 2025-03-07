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

#include "TopoFeature.h"

//-- TopoFeature class
TopoFeature::TopoFeature()
        : m_mesh(), m_id(), m_f_active(true), m_f_imported(false), m_outputLayerID(-1) {}

TopoFeature::TopoFeature(std::string pid)
        : m_mesh(), m_id(std::move(pid)), m_f_active(true), m_f_imported(false), m_outputLayerID(-1) {}

TopoFeature::TopoFeature(int outputLayerID)
        : m_mesh(), m_id(), m_f_active(true), m_f_imported(false), m_outputLayerID(outputLayerID) {
    if (m_outputLayerID >= s_numOfOutputLayers) s_numOfOutputLayers = m_outputLayerID + 1;
}

int TopoFeature::s_numOfOutputLayers = 0;

void TopoFeature::add_recon_region_output_layers(const int numLayers) {
    s_numOfOutputLayers += numLayers;
}

int TopoFeature::get_num_output_layers() {
    return s_numOfOutputLayers;
}

Mesh& TopoFeature::get_mesh() {
    return m_mesh;
}

const Mesh& TopoFeature::get_mesh() const {
    return m_mesh;
}

void TopoFeature::set_id(unsigned long id) {
    m_id = std::to_string(id);
}

std::string TopoFeature::get_id() const {
    return m_id;
}

int TopoFeature::get_output_layer_id() const {
    return m_outputLayerID;
}

bool TopoFeature::is_active() const {
    return m_f_active;
}

bool TopoFeature::is_imported() const {
    return m_f_imported;
}

void TopoFeature::deactivate() {
    m_f_active = false;
}

void TopoFeature::get_cityjson_cityobj_info(nlohmann::json& /* f */) const { }

void TopoFeature::get_cityjson_geomobj_info(nlohmann::json& /* g */) const { }

void TopoFeature::get_cityjson_semantics(nlohmann::json& /* g */) const { }