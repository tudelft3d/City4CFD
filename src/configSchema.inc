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

#ifndef CITY4CFD_CONFIGSCHEMA_INC
#define CITY4CFD_CONFIGSCHEMA_INC

#include "types.h"

namespace jsonschema {
    //-- JSON config schema defined here as a string literal
    const nlohmann::json schema = R"(
    {
      "$schema": "https://json-schema.org/draft-07/schema#",
      "$id": "https://json-schema.org/draft-07/schema",
      "title": "City4CFD Configuration File Schema",
      "type": "object",
      "properties": {
        "point_clouds": {
          "type": "object",
          "properties": {
            "las_datasets": {
              "description": "Schema that deals with import of LAS/LAZ tiles",
              "description": "Deprecated, but keeping the schema for now",
              "type": "object",
              "properties": {
                "files": {
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
                }
              }
            },
            "ground": {
              "type": "string"
            },
            "buildings": {
              "type": "string"
            }
          },
          "additionalProperties": false
        },
        "polygons": {
          "description": "Input polygons in GeoJSON format",
          "type": "array",
          "items": {
            "type": "object",
            "properties": {
              "type": {
                "enum": ["Building", "SurfaceLayer"]
              },
              "path": {
                "type": "string"
              }
            },
            "allOf": [
              {
                "if": {
                  "properties": {"type": {"const": "Building"}}
                },
                "then": {
                  "properties": {
                    "height_attribute": {"type": "string"},
                    "floor_attribute": {"type": "string"},
                    "floor_height": {
                      "type": "number",
                      "exclusiveMinimum": 0
                    },
                    "height_attribute_advantage": {"type": "boolean"},
                    "avoid_bad_polys": {"type": "boolean"}
                  },
                  "dependencies": {
                    "floor_attribute": ["floor_height"]
                  },
                  "if": {
                    "not": {
                      "properties": {"height_attribute_advantage": {"const": false}}
                    }
                  },
                  "then": {
                    "anyOf": [
                      {"required": ["height_attribute"]},
                      {"required": ["floor_attribute"]}
                    ]
                  }
                }
              },
              {
                "if": {
                  "properties": {"type": {"const": "SurfaceLayer"}}
                },
                "then": {
                  "properties": {
                    "layer_name": {"type": "string"},
                    "flatten_surface": {"type": "boolean"},
                    "flatten_vertical_border": {"type": "boolean"},
                    "surface_percentile": {
                      "type": "number",
                      "exclusiveMinimum": 0,
                      "maximum": 100
                    }
                  },
                  "dependencies": {
                    "flatten_surface": ["surface_percentile"]
                  }
                }
              }
            ],
            "required": ["type", "path"]
          }
        },
        "reconstruction_regions": {
          "description": "Reconstruction regions/regions of interest/influence zones",
          "type": "array",
          "minItems": 1,
          "items": {
            "type": "object",
            "properties": {
              "influence_region": {
                "$ref": "#/$defs/influBndStruct"
              },
              "lod": {
                "type": "string",
                "enum": ["1.2", "1.3", "2.2"]
              },
              "complexity_factor": {
                "type": "number",
                "exclusiveMinimum": 0,
                "exclusiveMaximum": 1
              },
              "lod13_step_height": {
                "type": "number",
                "minimum": 0
              },
              "validate": {
                "type": "boolean"
              },
              "enforce_validity": {
                "type": "string",
                "enum": ["lod1.2", "surface_wrap"]
              },
              "relative_alpha": {
                "type": "number",
                "exclusiveMinimum": 0
              },
              "relative_offset": {
                "type": "number",
                "exclusiveMinimum": 0
              },
              "skip_gap_closing": {
                "type": "boolean"
              },
              "bpg_influence_region_extra": {
                "type": "number",
                "minimum": 0
              },
              "import_advantage": {
                "type": "boolean"
              },
              "refine": {
                "type": "boolean"
              }
            },
            "required": ["influence_region", "lod"],
            "oneOf": [
              {
                "properties": {
                  "enforce_validity": {
                    "not": {
                      "const": "surface_wrap"
                    }
                  }
                }
              },
              {
                "properties": {
                  "enforce_validity": {
                    "const": "surface_wrap"
                  },
                  "relative_alpha": {
                    "type": "number"
                  },
                  "relative_offset": {
                    "type": "number"
                  }
                },
                "required": ["relative_alpha", "relative_offset"]
              }
            ]
          }
        },
        "import_geometries": {
          "description": "Handling imported geometries",
          "type": "object",
          "properties": {
            "path": {
              "type": "string"
            },
            "true_height": {
              "type": "boolean"
            },
            "lod": {
              "type": "string"
            },
            "refine": {
              "type": "boolean"
            }
          },
          "required": ["path", "true_height"]
        },
        "point_of_interest": {
          "$ref": "#/$defs/point"
        },
        "domain_bnd":{
          "$ref": "#/$defs/influBndStruct"
        },
        "top_height": {
          "description": "Not required when domain_bnd is null",
          "type": "number",
          "exclusiveMinimum": 0
        },
        "bnd_type_bpg": {
          "enum": ["Round", "Rectangle", "Oval"]
        },
        "bpg_blockage_ratio": {
          "oneOf": [
            {
              "type": "boolean"
            },
            {
              "type": "number",
              "exclusiveMinimum": 0,
              "exclusiveMaximum": 100
            }
          ]
        },
        "flow_direction": {
          "description": "Required for rectangular and oval domain, and only for round when using bpg_blockage_ratio",
          "$ref": "#/$defs/point"
        },
        "bpg_domain_size": {
          "description": "Optional: BPG multiplication factor, defined down below as a conditional statement"
        },
        "buffer_region": {
          "description": "Optional: Buffer region that can extrude both inwards and outwards",
          "type": "number",
          "exclusiveMinimum": -99
        },
        "terrain_thinning": {
          "description": "Optional: Random thinning of the terrain point cloud",
          "type": "number",
          "minimum": 0,
          "exclusiveMaximum": 100
        },
        "smooth_terrain": {
          "description": "Optional: Terrain smoothing flag",
          "type": "object",
          "properties": {
            "iterations": {
              "type": "integer",
              "minimum": 0
            },
            "max_pts": {
              "type": "integer",
              "exclusiveMinimum": 0
            }
          }
        },
        "flat_terrain": {
          "description": "Optional: Force flat terrain",
          "type": "boolean"
        },
        "building_percentile": {
          "type": "number",
          "exclusiveMinimum": 0,
          "exclusiveMaximum": 100
        },
        "min_height": {
          "description": "Default: 2.; minimum allowed height to be reconstructed",
          "type": "number",
          "exclusiveMinimum": 0
        },
        "min_area": {
          "description": "Optional; minimum allowed footprint area",
          "type": "number",
          "exclusiveMinimum": 0
        },
        "reconstruct_failed": {
          "description": "Optional: Reconstruct failed buildings to minimum height",
          "type": "boolean"
        },
        "intersect_buildings_terrain": {
          "description": "Optional: Force buildings to protrude through terrain",
          "type": "boolean"
        },
        "reconstruct_boundaries": {
          "type": "boolean"
        },
        "edge_max_len": {
          "type": "number",
          "exclusiveMinimum": 0
        },
        "output_file_name": {
          "type": "string"
        },
        "output_format": {
          "type": "string"
        },
        "output_separately": {
          "type": "boolean"
        },
        "output_log": {
          "description": "Flag to output log file and GeoJSON of failed reconstructions, optional",
          "type": "boolean"
        },
        "log_name": {
          "description": "Optional, default to 'log'",
          "type": "string"
        }
      },
      "$defs": {
        "point": {
          "type": "array",
          "items": {
            "type": "number"
          },
          "minItems": 2,
          "maxItems": 2
        },
        "influBndStruct": {
          "description": "Radius, Poly points, GeoJSON polygon or null for BPG",
          "oneOf": [
            {
              "type": "number"
            },
            {
              "type": "array",
              "items": {
                "type": "number"
              },
              "maxItems": 1
            },
            {
              "type": "array",
              "items": {
                "$ref": "#/$defs/point"
                },
              "minItems": 3
            },
            {
              "type": ["string", "null"]
            }
          ]
        }
      },
      "allOf": [
        {
          "anyOf": [
            {
              "description": "Input for buildings can be point cloud, reconstructed buildings, attribute and their combination",
              "required": ["point_clouds"],
              "properties": {
                "point_clouds": {"required": ["buildings"]},
                "polygons": {
                  "contains": {
                    "description": "Need at least building polygons",
                    "type": "object",
                    "properties": {
                      "type": {
                        "const": "Building"
                      }
                    }
                  }
                }
              }
            },
            {
              "properties": {
                "polygons": {
                  "contains": {
                    "type": "object",
                    "properties": {
                      "type": {"const": "Building"}
                    },
                    "anyOf": [
                      {"required": ["height_attribute"]},
                      {"required": ["floor_attribute"]}
                    ]
                  }
                }
              }
            },
            {
              "required": ["import_geometries"]
            }
          ],
          "if": {
            "properties": {"domain_bnd": {"type": "null"}}
          },
          "then": {
            "required": ["bnd_type_bpg"],
            "anyOf": [
              {
                "properties": {"bnd_type_bpg": {"const": "Round"}}
              },
              {
                "required": ["flow_direction"]
              }
            ]
          },
          "else": {
            "required": ["top_height"]
          }
        },
        {
          "if": {
            "not": {
              "properties": {"bpg_blockage_ratio": {"const": false}}
            }
          },
          "then": {
            "required": ["flow_direction"]
          }
        },
        {
          "if": {
            "properties": {"bnd_type_bpg": {"const": "Round"}}
          },
          "then": {
            "properties": {"bpg_domain_size": {"$ref": "#/$defs/point"}}
          },
          "else": {
            "properties": {
              "bpg_domain_size": {
                "type": "array",
                "items": {
                  "type": "number"
                },
                "minItems": 4,
                "maxItems": 4
              }
            }
          }
        }
      ],
      "required": [ "point_of_interest", "domain_bnd", "reconstruction_regions",
        "edge_max_len", "output_file_name", "output_format", "building_percentile",
        "output_separately"]
    }
        )"_json;
}

#endif //CITY4CFD_CONFIGSCHEMA_INC