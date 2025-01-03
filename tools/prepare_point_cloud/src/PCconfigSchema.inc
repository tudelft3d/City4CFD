/*
  City4CFD
 
  Copyright (c) 2021-2025, 3D Geoinformation Research Group, TU Delft

  This file is part of City4CFD.

  City4CFD is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  City4CFD is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with City4CFD.  If not, see <http://www.gnu.org/licenses/>.

  For any information or further details about the use of City4CFD, contact
  Ivan Pađen
  <i.paden@tudelft.nl>
  3D Geoinformation Research Group
  Delft University of Technology
*/

#ifndef CITY4CFD_PCCONFIGSCHEMA_INC
#define CITY4CFD_PCCONFIGSCHEMA_INC

namespace jsonschema {
    //-- JSON config schema defined here as a string literal
    const nlohmann::json schema = R"(
    {
      "$schema": "https://json-schema.org/draft-07/schema#",
      "$id": "https://json-schema.org/draft-07/schema",
      "title": "City4CFD Point Cloud Preparation Config Schema",
      "type": "object",
      "properties": {
        "las_files": {
          "description": "Schema that deals with import of LAS/LAZ tiles",
          "type": "array",
          "items": {
            "type": "string"
          }
        },
        "ground_classes": {
          "type": "array",
          "items": {
            "type": "number"
          }
        },
        "building_classes": {
          "type": "array",
          "items": {
            "type": "number"
          }
        },
        "bbox": {
          "description": "Optional: set the bounding box",
          "type": "array",
          "minItems": 4,
          "maxItems": 4,
          "items": {
            "type": "number"
          }
        },
        "thinning": {
          "description": "Optional: Random thinning",
          "type": "number",
          "minimum": 0,
          "exclusiveMaximum": 100
        }
      },
      "required": [ "las_files" ]
    }
        )"_json;
}

#endif //CITY4CFD_PCCONFIGSCHEMA_INC