#ifndef CITY4CFD_CONFIGSCHEMA_INC
#define CITY4CFD_CONFIGSCHEMA_INC

#include "types.h"

namespace jsonschema {
    //- JSON config schema defined here as a string literal
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
            "ground": {
              "type": "string"
            },
            "buildings": {
              "type": "string"
            }
          },
          "required": ["ground"],
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
            "if": {
              "properties": {"type": {"const": "SurfaceLayer"}}
            },
            "then": {
              "properties": {
                "layer_name": {"type": "string"},
                "average_surface": {"type": "boolean"},
                "surface_percentile": {"type": "number"}
              },
              "dependencies": {
                "average_surface": ["surface_percentile"]
              }
            },
            "required": ["type", "path"]
          }
        },
        "import_geometries": {
          "description": "Handling imoprted CityJSON geometries",
          "type": "object",
          "properties": {
            "path": {
              "type": "string"
            },
            "advantage": {
              "type": "boolean"
            },
            "true_height": {
              "type": "boolean"
            },
            "lod": {
              "type": "string"
            }
          },
          "required": ["path", "advantage", "true_height", "lod"]
        },
        "point_of_interest": {
          "$ref": "#/$defs/point"
        },
        "influence_region": {
          "$ref": "#/$defs/influBndStruct"
        },
        "domain_bnd":{
          "$ref": "#/$defs/influBndStruct"
        },
        "top_height": {
          "description": "Not required when domain_bnd is null",
          "type": "number"
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
        "terrain_simplification": {
          "description": "Optional: Random simplification of the terrain point cloud",
          "type": "number",
          "minimum": 0,
          "exclusiveMaximum": 100
        },
        "smooth_terrain": {
          "description": "Optional: Terrain smoothing flag",
          "type": "boolean"
        },
        "lod": {
          "type": "string"
        },
        "building_percentile": {
          "type": "number",
          "exclusiveMinimum": 0,
          "exclusiveMaximum": 100
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
        "required": [ "point_of_interest", "influence_region", "domain_bnd",
                      "building_percentile", "edge_max_len", "output_file_name", "output_format",
                      "output_separately" ]
    }
        )"_json;
}

#endif //CITY4CFD_CONFIGSCHEMA_INC