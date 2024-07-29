/*
  val3dity 

  Copyright (c) 2011-2024, 3D geoinformation research group, TU Delft

  This file is part of val3dity.

  val3dity is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  val3dity is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with val3dity.  If not, see <http://www.gnu.org/licenses/>.

  For any information or further details about the use of val3dity, contact
  Hugo Ledoux
  <h.ledoux@tudelft.nl>
  Faculty of the Built Environment & Architecture
  Delft University of Technology
  Julianalaan 134, Delft 2628BL, the Netherlands
*/

#include "val3dity.h"
#include "Feature.h"
#include "CityObject.h"
#include "GenericObject.h"
#include "Solid.h"
#include "Surface.h"
#include "input.h"

#include <iostream>
#include <exception>      // std::exception


namespace val3dity
{

std::string VAL3DITY_VERSION = "2.5.0b1";

struct verror : std::exception {
  std::string whattext;
  verror(std::string s) : whattext(s){};
  const char* what() const noexcept {return whattext.c_str();}
};

json
validate_onegeom(json& j,
                 Parameters params)
{
  std::vector<Feature*> lsFeatures;
  parse_tu3djson_onegeom(j, lsFeatures, params._tol_snap);
  //-- validate
  for (auto& f : lsFeatures)
      f->validate(params._planarity_d2p_tol, params._planarity_n_tol, params._overlap_tol);
  //-- get report in json
  IOErrors ioerrs;
  ioerrs.set_input_file_type("tu3djson_geom");
  json jr = get_report_json("JSON object",
                            lsFeatures,
                            VAL3DITY_VERSION,
                            params._tol_snap,
                            params._overlap_tol,
                            params._planarity_d2p_tol,
                            params._planarity_n_tol,
                            ioerrs);
  return jr;
}

json validate_jsonfg(json& j,
                     Parameters params)
{
  IOErrors ioerrs;
  ioerrs.set_input_file_type("JSON-FG");
  std::vector<Feature*> lsFeatures;
  parse_jsonfg(j, lsFeatures, params._tol_snap, ioerrs);
  //-- validate
  for (auto& f : lsFeatures)
      f->validate(params._planarity_d2p_tol, params._planarity_n_tol, params._overlap_tol);
  //-- get report in json
  json jr = get_report_json("JSON object",
                            lsFeatures,
                            VAL3DITY_VERSION,
                            params._tol_snap,
                            params._overlap_tol,
                            params._planarity_d2p_tol,
                            params._planarity_n_tol,
                            ioerrs);
  return jr;
}

json validate_tu3djson(json& j,
                       Parameters params)
{
  std::vector<Feature*> lsFeatures;
  parse_tu3djson(j, lsFeatures, params._tol_snap);
  //-- validate
  for (auto& f : lsFeatures)
      f->validate(params._planarity_d2p_tol, params._planarity_n_tol, params._overlap_tol);
  //-- get report in json
  IOErrors ioerrs;
  ioerrs.set_input_file_type("tu3djson");
  json jr = get_report_json("JSON object",
                            lsFeatures,
                            VAL3DITY_VERSION,
                            params._tol_snap,
                            params._overlap_tol,
                            params._planarity_d2p_tol,
                            params._planarity_n_tol,
                            ioerrs);
  return jr;
}

json
validate_cityjson(json& j,
                  Parameters params)
{
  std::vector<Feature*> lsFeatures;
  //-- parse the cityjson object
  //-- compute (_minx, _miny)
  compute_min_xy(j);
  //-- read and store the GeometryTemplates
  std::vector<GeometryTemplate*> lsGTs;
  if (j.count("geometry-templates") == 1)
  {
      process_cityjson_geometrytemplates(j["geometry-templates"], lsGTs, params._tol_snap);
  }
  //-- process each CO
  for (json::iterator it = j["CityObjects"].begin(); it != j["CityObjects"].end(); ++it)
  {
      //-- BuildingParts geometries are put with those of a Building
      if (it.value()["type"] == "BuildingPart")
          continue;
      CityObject* co = new CityObject(it.key(), it.value()["type"]);
      process_json_geometries_of_co(it.value(), co, lsGTs, j, params._tol_snap);
      //-- if Building has Parts, put them here in _lsPrimitives
      if ( (it.value()["type"] == "Building") && (it.value().count("children") != 0) )
      {
          for (std::string bpid : it.value()["children"])
          {
              process_json_geometries_of_co(j["CityObjects"][bpid], co, lsGTs, j, params._tol_snap);
          }
      }
      lsFeatures.push_back(co);
  }
  //-- validate
  for (auto& f : lsFeatures)
      f->validate(params._planarity_d2p_tol, params._planarity_n_tol, params._overlap_tol);
  //-- compile errors
  std::set<int> errors;
  for (auto& f : lsFeatures)
      for (auto& p : f->get_primitives())
          for (auto& code : p->get_unique_error_codes())
              errors.insert(code);
  //-- get report in json
  IOErrors ioerrs;
  ioerrs.set_input_file_type("CityJSON");
  json jr = get_report_json("JSON object",
                            lsFeatures,
                            VAL3DITY_VERSION,
                            params._tol_snap,
                            params._overlap_tol,
                            params._planarity_d2p_tol,
                            params._planarity_n_tol,
                            ioerrs);
  return jr;
}

json
validate_cityjsonfeature(json& j,
                         Parameters params)
{
  // j["transform"] = jtransform;
  std::vector<Feature*> lsFeatures;
  //-- compute (_minx, _miny)
  compute_min_xy(j);
  //-- list empty GeometryTemplate TODO: populate this?
  std::vector<GeometryTemplate*> lsGTs;
  //-- process each CO
  for (json::iterator it = j["CityObjects"].begin(); it != j["CityObjects"].end(); ++it)
  {
      //-- BuildingParts geometries are put with those of a Building
      if (it.value()["type"] == "BuildingPart")
          continue;
      CityObject* co = new CityObject(it.key(), it.value()["type"]);
      process_json_geometries_of_co(it.value(), co, lsGTs, j, params._tol_snap);
      //-- if Building has Parts, put them here in _lsPrimitives
      if ( (it.value()["type"] == "Building") && (it.value().count("children") != 0) )
      {
          for (std::string bpid : it.value()["children"])
          {
              process_json_geometries_of_co(j["CityObjects"][bpid], co, lsGTs, j, params._tol_snap);
          }
      }
      lsFeatures.push_back(co);
  }
  //-- validate
  for (auto& f : lsFeatures)
      f->validate(params._planarity_d2p_tol, params._planarity_n_tol, params._overlap_tol);
  //-- compile errors
  std::set<int> errors;
  for (auto& f : lsFeatures)
      for (auto& p : f->get_primitives())
          for (auto& code : p->get_unique_error_codes())
              errors.insert(code);
  //-- get report in json
  IOErrors ioerrs;
  ioerrs.set_input_file_type("CityJSONFeature");
  json jr = get_report_json("JSON object",
                            lsFeatures,
                            VAL3DITY_VERSION,
                            params._tol_snap,
                            params._overlap_tol,
                            params._planarity_d2p_tol,
                            params._planarity_n_tol,
                            ioerrs);
  return jr;
}

json
validate_indoorgml(std::string& input,
                   Parameters params)
{
  IOErrors ioerrs;
  ioerrs.set_input_file_type("IndoorGML");
  pugi::xml_document doc;
  pugi::xml_parse_result result = doc.load_string(input.c_str());
  if (!result) {
      ioerrs.add_error(901, "Input value not valid XML");
  }
  std::vector<Feature*> lsFeatures;
  if (ioerrs.has_errors() == false) {
      //-- parse namespace
      pugi::xml_node ncm = doc.first_child();
      std::map<std::string, std::string> thens = get_namespaces(ncm); //-- results in global variable NS in this unit
      if ( (thens.count("indoorgml") != 0) && (ncm.name() == (thens["indoorgml"] + "IndoorFeatures")) ) {
          //-- find (_minx, _miny)
          compute_min_xy(doc);
          //-- build dico of xlinks for <gml:Polygon>
          std::map<std::string, pugi::xpath_node> dallpoly;
          build_dico_xlinks(doc, dallpoly, ioerrs);
          ioerrs.set_input_file_type("IndoorGML");
          process_gml_file_indoorgml(doc, lsFeatures, dallpoly, ioerrs, params._tol_snap);
      }
      else
      {
          ioerrs.add_error(904, "GML files not supported (yes that includes CityGML files ==> upgrade to CityJSON)");
      }
  }
  //-- start the validation
  if (ioerrs.has_errors() == false) {
      //-- validate
      for (auto& f : lsFeatures)
          f->validate(params._planarity_d2p_tol, params._planarity_n_tol, params._overlap_tol);
  }
  //-- get report in json
  json jr = get_report_json("JSON object",
                            lsFeatures,
                            VAL3DITY_VERSION,
                            params._tol_snap,
                            params._overlap_tol,
                            params._planarity_d2p_tol,
                            params._planarity_n_tol,
                            ioerrs);
  return jr;
}

json
validate_obj(std::string& input,
             Parameters params)
{
  IOErrors ioerrs;
  ioerrs.set_input_file_type("OBJ");
  std::vector<Feature*> lsFeatures;
  std::istringstream iss(input);
  parse_obj(iss, lsFeatures, SOLID, ioerrs, params._tol_snap);
  //-- start the validation
  if (ioerrs.has_errors() == false) {
      //-- validate
      for (auto& f : lsFeatures)
          f->validate(params._planarity_d2p_tol, params._planarity_n_tol, params._overlap_tol);
  }
  //-- get report in json
  json jr = get_report_json("OBJ object",
                            lsFeatures,
                            VAL3DITY_VERSION,
                            params._tol_snap,
                            params._overlap_tol,
                            params._planarity_d2p_tol,
                            params._planarity_n_tol,
                            ioerrs);
  return jr;
}

json
validate_off(std::string& input,
             Parameters params)
{
  IOErrors ioerrs;
  ioerrs.set_input_file_type("OFF");
  std::vector<Feature*> lsFeatures;
  GenericObject* o = new GenericObject("none");
  std::istringstream iss(input);
  Surface* sh = parse_off(iss, 0, ioerrs, params._tol_snap);
  std::cout << "1" << std::endl;
  Solid* s = new Solid;
  s->set_oshell(sh);
  o->add_primitive(s);
  lsFeatures.push_back(o);
  //-- start the validation
  // std::cout << "errors: " << ioerrs.has_errors() << std::endl;
  for (auto& each : ioerrs.get_unique_error_codes())
      std::cout << each << std::endl;

  // if (ioerrs.has_errors() == false) {
  //-- validate
  for (auto& f : lsFeatures)
      f->validate(params._planarity_d2p_tol, params._planarity_n_tol, params._overlap_tol);
  // }
  std::cout << "2" << std::endl;
  //-- get report in json
  json jr = get_report_json("OFF object",
                            lsFeatures,
                            VAL3DITY_VERSION,
                            params._tol_snap,
                            params._overlap_tol,
                            params._planarity_d2p_tol,
                            params._planarity_n_tol,
                            ioerrs);
  std::cout << "3" << std::endl;
  return jr;
}

bool
is_valid(json& j,
         Parameters params)
{
  json re = validate(j, params);
  return re["validity"];
}

json 
validate(json& j,
         Parameters params)
{
  spdlog::set_level(spdlog::level::off);
  json re;
  //-- CityJSON
  if (j["type"] == "CityJSON") {
    json jr = validate_cityjson(j, params);
    re = jr;
  
  //-- CityJSONFeature
  } else if (j["type"] == "CityJSONFeature") {
    json jr = validate_cityjsonfeature(j, params);
    re = jr;
  
  //-- tu3djson
  } else if (j["type"] == "tu3djson") {
    json jr = validate_tu3djson(j, params);
    re = jr;
  
  //-- JSON-FG
  } else if (j["type"] == "Feature") { 
    json jr = validate_jsonfg(j, params);
    re = jr;
  } else if (j["type"] == "FeatureCollection") { 
    json jr = validate_jsonfg(j, params);
    re = jr;

  //-- tu3djson onegeom
  } else if ( (j["type"] == "MultiSurface") || 
              (j["type"] == "CompositeSurface") ||
              (j["type"] == "Solid") ||
              (j["type"] == "MultiSolid") ||
              (j["type"] == "CompositeSolid") ) { 
    json jr = validate_onegeom(j, params);
    re = jr;

  //-- then we don't support it   
  } else {
    throw verror("Flavour of JSON not supported");
  }  
  return re;
}

bool 
is_valid(const std::vector<std::array<double, 3>>& vertices,
         const std::vector<std::vector<int>>& faces,
         Parameters params)
{
  json re = validate(vertices, faces, params);
  return re["validity"];
}

json
validate(const std::vector<std::array<double, 3>>& vertices,
         const std::vector<std::vector<int>>& faces,
         Parameters params)
{
  spdlog::set_level(spdlog::level::off);
  double _minx = 9e15;
  double _miny = 9e15; 
  //-- find (minx, miny)
  for (auto& v: vertices) {
    if (v[0] < _minx)
      _minx = v[0];
    if (v[1] < _miny)
      _miny = v[1];
  }
  //-- create a Surface (a 2-manifold)
  Surface* sh = new Surface("0", params._tol_snap);
  std::vector<Point3*> allvertices;
  GenericObject* o = new GenericObject("none");
  //-- read all the vertices
  for (auto& v: vertices) {
    Point3 *p = new Point3(v[0] - _minx, v[1] - _miny, v[2]);
    allvertices.push_back(p);
  }

  //-- read all the faces (0-indexed!)
  for (auto& vids: faces) {
    std::vector<int> r;
    for (auto& vid: vids) {
      Point3* tp = allvertices[vid];
      r.push_back(sh->add_point(*tp));
    }
    std::vector< std::vector<int> > pgnids;
    pgnids.push_back(r);
    sh->add_face(pgnids);
  }
  //-- we assume it's a Solid (TODO: should this be a param?)
  Solid* sol = new Solid("");
  sol->set_oshell(sh);
  o->add_primitive(sol);
  o->validate(params._planarity_d2p_tol, params._planarity_n_tol, params._overlap_tol);
  IOErrors ioerrs;
  ioerrs.set_input_file_type("std::vectors");
  std::vector<Feature*> lsFeatures;
  lsFeatures.push_back(o); 
  json jr = get_report_json("std::vectors",
                            lsFeatures,
                            VAL3DITY_VERSION,
                            params._tol_snap,
                            params._overlap_tol,
                            params._planarity_d2p_tol,
                            params._planarity_n_tol,
                            ioerrs);
  return jr;
}

//-- for ASCII + XML formats
bool 
is_valid(std::string& input,
         std::string format,
         Parameters params)
{
  json re = validate(input, format, params);
  return re["validity"];
}

json
validate(std::string& input,
         std::string format,
         Parameters params)
{
  spdlog::set_level(spdlog::level::off);
  json re;
  if (format == "IndoorGML") {
    json j = validate_indoorgml(input, params);
    re = j;
  }
  else if (format == "OBJ") {
    json j = validate_obj(input, params);
    re = j;
  }
  else if (format == "OFF") {
    json j = validate_off(input, params);
    re = j;
  }
  else { 
    throw verror("File type not supported");
  }
  return re;
}

}
