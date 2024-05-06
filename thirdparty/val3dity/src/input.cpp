/*
  val3dity 

  Copyright (c) 2011-2023, 3D geoinformation research group, TU Delft

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
  Faculty of Architecture & the Built Environment
  Delft University of Technology
  Julianalaan 134, Delft 2628BL, the Netherlands
*/

#include "input.h"
#include "Primitive.h"
#include "Feature.h"
#include "CityObject.h"
#include "GenericObject.h"
#include "IndoorModel.h"
#include "IndoorGraph.h"
#include "Surface.h"
#include "MultiSurface.h"
#include "CompositeSurface.h"
#include "Solid.h"
#include "CompositeSolid.h"
#include "MultiSolid.h"
#include "GeometryTemplate.h"
#include "val3dity_ostream.h"


using namespace std;
using json = nlohmann::json;

namespace val3dity
{

double _minx = 9e15;
double _miny = 9e15;  

//-- XML namespaces map
std::map<std::string, std::string> NS; 

bool IOErrors::has_errors()
{
  if (_errors.size() == 0)
    return false;
  else
    return true;
}

bool IOErrors::has_specific_error(int i) 
{
  if (_errors.count(i) == 0)
    return false;
  else
    return true;
}

std::set<int> IOErrors::get_unique_error_codes()
{
  std::set<int> errs;
  for (auto& err : _errors)
    for (auto i : err.second)
      errs.insert(std::get<0>(err));
  return errs;
}


void IOErrors::add_error(int code, std::string info)
{
  _errors[code].push_back(info);
  *val3ditycout << "ERROR " << code << " : " << info << std::endl;
}


std::string IOErrors::get_report_text()
{
  std::stringstream ss;
  for (auto& err : _errors)
  {
    for (auto i : err.second)
    {
      ss << err.first << " -- " << ALL_ERRORS[err.first] << std::endl;
      ss << "\tInfo: " << i << std::endl;
    }
  }
  return ss.str();
}


json IOErrors::get_report_json()
{
  json j;
  for (auto& err : _errors)
  {
    for (auto i : err.second)
    {
      json jj;
      jj["type"] = "Error";
      jj["code"] = err.first;
      jj["description"] = ALL_ERRORS[err.first];
      jj["info"] = i;
      j.push_back(jj);
    }
  }
  return j;
}


std::string IOErrors::get_input_file_type() {
  return _inputfiletype;
}


void IOErrors::set_input_file_type(std::string s) {
  _inputfiletype = s;
}



//-- ignore XML namespace
std::string localise(std::string s)
{
  return "*[local-name(.) = '" + s + "']";
}

  
std::string remove_xml_namespace(const char* input)
{
  std::string s = input;
  return s.substr(s.find_first_of(":") + 1);
}


vector<int> process_gml_ring(const pugi::xml_node& n, Surface* sh, IOErrors& errs) {
  std::string s = "./" + NS["gml"] + "LinearRing" + "/" + NS["gml"] + "pos";
  pugi::xpath_node_set npos = n.select_nodes(s.c_str());
  std::vector<int> r;
  if (npos.size() > 0) //-- <gml:pos> used
  {
    for (pugi::xpath_node_set::const_iterator it = npos.begin(); it != npos.end(); ++it) {
      std::string buf;
      std::stringstream ss(it->node().child_value());
      std::vector<std::string> tokens;
      while (ss >> buf)
        tokens.push_back(buf);
      long double x = std::stold(tokens[0]);
      x -= _minx;
      long double y = std::stold(tokens[1]);
      y -= _miny;
      Point3 p(double(x), double(y), std::stod(tokens[2]));
      r.push_back(sh->add_point(p));
    }
  }
  else //-- <gml:posList> used
  {
    std::string s = "./" + NS["gml"] + "LinearRing" + "/" + NS["gml"] + "posList";
    pugi::xpath_node pl = n.select_node(s.c_str());
    if (pl == NULL)
    {
      throw 901;
    }
    std::string buf;
    std::stringstream ss(pl.node().child_value());
    std::vector<std::string> coords;
    while (ss >> buf)
      coords.push_back(buf);
    if (coords.size() % 3 != 0)
    {
      errs.add_error(901, "Error: <gml:posList> has bad coordinates.");
      return r;
    }
    for (int i = 0; i < coords.size(); i += 3)
    {
      long double x = std::stold(coords[i]);
      x -= _minx;
      long double y = std::stold(coords[i+1]);
      y -= _miny;
      Point3 p(double(x), double(y), std::stod(coords[i+2]));
      r.push_back(sh->add_point(p));
    }
  }
  return r;
}


Surface* process_gml_surface(const pugi::xml_node& n, int id, std::map<std::string, pugi::xpath_node>& dallpoly, double tol_snap, IOErrors& errs) 
{
  std::string s = ".//" + NS["gml"] + "surfaceMember";
  pugi::xpath_node_set nsm = n.select_nodes(s.c_str());
  Surface* sh = new Surface(id, tol_snap);
  int i = 0;
  for (pugi::xpath_node_set::const_iterator it = nsm.begin(); it != nsm.end(); ++it)
  {
    std::vector< std::vector<int> > oneface;
    bool bxlink = false;
    pugi::xml_node tmpnode = it->node();
    pugi::xpath_node p;
    bool fliporientation = false;
    for (pugi::xml_attribute attr = tmpnode.first_attribute(); attr; attr = attr.next_attribute())
    {
      if (strcmp(attr.value(), "xlink:href") != 0) {
        bxlink = true;
        break;
      }
    }
    if (bxlink == true) 
    {
      std::string k = it->node().attribute("xlink:href").value();
      if (k[0] == '#')
        k = k.substr(1);
      p = dallpoly[k];
    }
    else
    {
      for (pugi::xml_node child : it->node().children()) 
      {
        if (std::string(child.name()).find("Polygon") != std::string::npos) {
          p = child;
          break;
        }
        else if (std::string(child.name()).find("OrientableSurface") != std::string::npos) {
          if (std::strncmp(child.attribute("orientation").value(), "-", 1) == 0)
            fliporientation = true;
          for (pugi::xml_node child2 : child.children()) 
          {
            if (std::string(child2.name()).find("baseSurface") != std::string::npos) 
            {
              std::string k = child2.attribute("xlink:href").value();
              if (k != "")
              {
                if (k[0] == '#')
                  k = k.substr(1);
                p = dallpoly[k];
                break;
              }
              for (pugi::xml_node child3 : child2.children())
              {
                if (std::string(child3.name()).find("Polygon") != std::string::npos)
                {
                  p = child;
                  break;
                }
              }
            }
          }
          break;
        }
        else if (std::string(child.name()).find("CompositeSurface") != std::string::npos) 
          break;
        else {
          throw 901;
        }
      }
    }

    //-- this is to handle CompositeSurfaces part of MultiSurfaces
    if (p == NULL) 
      continue;
    
    if (std::strncmp(p.node().attribute("orientation").value(), "-", 1) == 0)
      fliporientation = true;
    //-- exterior ring (only 1)
    s = ".//" + NS["gml"] + "exterior";
    pugi::xpath_node ring = p.node().select_node(s.c_str());
    std::vector<int> r = process_gml_ring(ring.node(), sh, errs);
    if (fliporientation == true) 
      std::reverse(r.begin(), r.end());
    if (r.front() != r.back())
    {
      if (p.node().attribute("gml:id") == 0)
        sh->add_error(103, std::to_string(i));
      else
        sh->add_error(103, p.node().attribute("gml:id").value());
    }
    else
      r.pop_back(); 
    oneface.push_back(r);
    //-- interior rings
    s = ".//" + NS["gml"] + "interior";
    pugi::xpath_node_set nint = p.node().select_nodes(s.c_str());
    for (pugi::xpath_node_set::const_iterator it = nint.begin(); it != nint.end(); ++it) {
      std::vector<int> r = process_gml_ring(it->node(), sh, errs);
      if (fliporientation == true) 
        std::reverse(r.begin(), r.end());
      if (r.front() != r.back())
      {
        if (p.node().attribute("gml:id") == 0)
          sh->add_error(103, std::to_string(i));
        else
          sh->add_error(103, p.node().attribute("gml:id").value());
      }
      else
        r.pop_back(); 
      oneface.push_back(r);
    }
    sh->add_face(oneface, p.node().attribute("gml:id").value());
    i++;
  }
  return sh;
}


Solid* process_gml_solid(const pugi::xml_node& nsolid, std::map<std::string, pugi::xpath_node>& dallpoly, double tol_snap, IOErrors& errs)
{
  //-- exterior shell
  Solid* sol = new Solid;
  if (nsolid.attribute("gml:id") != 0)
    sol->set_id(std::string(nsolid.attribute("gml:id").value()));
  std::string s = "./" + NS["gml"] + "exterior";
  pugi::xpath_node next = nsolid.select_node(s.c_str());
  sol->set_oshell(process_gml_surface(next.node(), 0, dallpoly, tol_snap, errs));
  //-- interior shells
  s = "./" + NS["gml"] + "interior";
  pugi::xpath_node_set nint = nsolid.select_nodes(s.c_str());
  int id = 1;
  for (pugi::xpath_node_set::const_iterator it = nint.begin(); it != nint.end(); ++it)
  {
    sol->add_ishell(process_gml_surface(it->node(), id, dallpoly, tol_snap, errs));
    id++;
  }
  return sol;
}


MultiSolid* process_gml_multisolid(const pugi::xml_node& nms, std::map<std::string, pugi::xpath_node>& dallpoly, double tol_snap, IOErrors& errs)
{
  MultiSolid* ms = new MultiSolid;
  if (nms.attribute("gml:id") != 0)
    ms->set_id(std::string(nms.attribute("gml:id").value()));
  std::string s = ".//" + NS["gml"] + "Solid";
  pugi::xpath_node_set nn = nms.select_nodes(s.c_str());
  for (pugi::xpath_node_set::const_iterator it = nn.begin(); it != nn.end(); ++it)
  {
    Solid* s = process_gml_solid(it->node(), dallpoly, tol_snap, errs);
    if (s->get_id() == "")
      s->set_id(std::to_string(ms->number_of_solids()));
    ms->add_solid(s);
  }
  return ms;
}


CompositeSolid* process_gml_compositesolid(const pugi::xml_node& nms, std::map<std::string, pugi::xpath_node>& dallpoly, double tol_snap, IOErrors& errs)
{
  CompositeSolid* cs = new CompositeSolid;
  if (nms.attribute("gml:id") != 0)
    cs->set_id(std::string(nms.attribute("gml:id").value()));
  std::string s = ".//" + NS["gml"] + "Solid";
  pugi::xpath_node_set nn = nms.select_nodes(s.c_str());
  for (pugi::xpath_node_set::const_iterator it = nn.begin(); it != nn.end(); ++it)
  {
    Solid* s = process_gml_solid(it->node(), dallpoly, tol_snap, errs);
    if (s->get_id() == "")
      s->set_id(std::to_string(cs->number_of_solids()));
    cs->add_solid(s);
  }

  return cs;
}


MultiSurface* process_gml_multisurface(const pugi::xml_node& nms, std::map<std::string, pugi::xpath_node>& dallpoly, double tol_snap, IOErrors& errs)
{
  MultiSurface* ms = new MultiSurface;
  if (nms.attribute("gml:id") != 0)
    ms->set_id(std::string(nms.attribute("gml:id").value()));
  Surface* s = process_gml_surface(nms, 0, dallpoly, tol_snap, errs);
  ms->set_surface(s);
  return ms;
}


CompositeSurface* process_gml_compositesurface(const pugi::xml_node& nms, std::map<std::string, pugi::xpath_node>& dallpoly, double tol_snap, IOErrors& errs)
{
  CompositeSurface* cs = new CompositeSurface;
  if (nms.attribute("gml:id") != 0)
    cs->set_id(std::string(nms.attribute("gml:id").value()));
  Surface* s = process_gml_surface(nms, 0, dallpoly, tol_snap, errs);
  cs->set_surface(s);
  return cs;
}


void process_json_surface(std::vector< std::vector<int> >& pgn, json& j, Surface* sh)
{
  std::vector< std::vector<int> > pgnids;
  for (auto& r : pgn)
  {
    std::vector<int> newr;
    for (auto& i : r)
    {
      double x;
      double y;
      double z;
      if (j.count("transform") == 0)
      {
        x = double(j["vertices"][i][0]);
        y = double(j["vertices"][i][1]);
        z = double(j["vertices"][i][2]);
      }
      else
      {
        x = (double(j["vertices"][i][0]) * double(j["transform"]["scale"][0])) + double(j["transform"]["translate"][0]);
        y = (double(j["vertices"][i][1]) * double(j["transform"]["scale"][1])) + double(j["transform"]["translate"][1]);
        z = (double(j["vertices"][i][2]) * double(j["transform"]["scale"][2])) + double(j["transform"]["translate"][2]);
      }
      x -= _minx;
      y -= _miny;
      Point3 p3(x, y, z);
      newr.push_back(sh->add_point(p3));
    }
    pgnids.push_back(newr);
  }
  sh->add_face(pgnids);
}


void process_jsonfg_surface(std::vector<std::vector<std::vector<double>>>& pgn, Surface* sh, IOErrors& errs)
{
  std::vector<std::vector<int>> pgnids;
  for (auto& r : pgn)
  {
    std::vector<int> newr;
    for (auto& p : r) {
      double x = p[0];
      double y = p[1];
      double z = p[2];
      Point3 p3(x, y, z);
      newr.push_back(sh->add_point(p3));
    }
    if ( newr.front() != newr.back() ) {
      sh->add_error(103, "First-last must be the same for JSON-FG rings");
    }
    newr.pop_back(); //-- delete last one because jsonfg "closes" its rings
    pgnids.push_back(newr);
  }
  sh->add_face(pgnids);
}


void process_json_geometries_of_co(json& jco, CityObject* co, std::vector<GeometryTemplate*>& lsGTs, json& j, double tol_snap)
{
  int idgeom = co->number_of_primitives();
  for (auto& g : jco["geometry"]) {
    std::string theid = std::to_string(idgeom);
    if  (g["type"] == "Solid")
    {
      Solid* s = new Solid(theid);
      std::string thelod = "";
      if (g["lod"].is_number()) {
        thelod = std::to_string(g["lod"].get<double>());
      } else {
        thelod = g["lod"].get<std::string>();
      }
      s->set_lod(thelod);
      bool oshell = true;
      int c = 0;
      for (auto& shell : g["boundaries"]) 
      {
        Surface* sh = new Surface(c, tol_snap);
        c++;
        for (auto& polygon : shell) { 
          std::vector< std::vector<int> > pa = polygon;
          process_json_surface(pa, j, sh);
        }
        if (oshell == true)
        {
          oshell = false;
          s->set_oshell(sh);
        }
        else
          s->add_ishell(sh);
      }
      co->add_primitive(s);
    }
    else if ( (g["type"] == "MultiSurface") || (g["type"] == "CompositeSurface") ) 
    {
      Surface* sh = new Surface(-1, tol_snap);
      for (auto& p : g["boundaries"]) 
      { 
        std::vector< std::vector<int> > pa = p;
        process_json_surface(pa, j, sh);
      }
      std::string thelod = "";
      if (g["lod"].is_number()) {
        thelod = std::to_string(g["lod"].get<double>());
      } else {
        thelod = g["lod"].get<std::string>();
      }
      if (g["type"] == "MultiSurface")
      {
        MultiSurface* ms = new MultiSurface(theid);
        ms->set_lod(thelod);
        ms->set_surface(sh);
        co->add_primitive(ms);
      }
      else
      {
        CompositeSurface* cs = new CompositeSurface(theid);
        cs->set_lod(thelod);
        cs->set_surface(sh);
        co->add_primitive(cs);
      }
    }
    else if (g["type"] == "MultiSolid") 
    {
      MultiSolid* ms = new MultiSolid(theid);
      std::string thelod = "";
      if (g["lod"].is_number()) {
        thelod = std::to_string(g["lod"].get<double>());
      } else {
        thelod = g["lod"].get<std::string>();
      }
      ms->set_lod(thelod);
      for (auto& solid : g["boundaries"]) 
      {
        Solid* s = new Solid();
        bool oshell = true;
        for (auto& shell : solid) 
        {
          Surface* sh = new Surface(-1, tol_snap);
          for (auto& polygon : shell) { 
            std::vector< std::vector<int> > pa = polygon;
            process_json_surface(pa, j, sh);
          }
          if (oshell == true)
          {
            oshell = false;
            s->set_oshell(sh);
          }
          else
            s->add_ishell(sh);
        }
        ms->add_solid(s);
      }
      co->add_primitive(ms);
    }
    else if (g["type"] == "CompositeSolid") 
    {
      CompositeSolid* cs = new CompositeSolid(theid);
      std::string thelod = "";
      if (g["lod"].is_number()) {
        thelod = std::to_string(g["lod"].get<double>());
      } else {
        thelod = g["lod"].get<std::string>();
      }
      cs->set_lod(thelod);
      for (auto& solid : g["boundaries"]) 
      {
        Solid* s = new Solid();
        bool oshell = true;
        for (auto& shell : solid) 
        {
          Surface* sh = new Surface(-1, tol_snap);
          for (auto& polygon : shell) { 
            std::vector< std::vector<int> > pa = polygon;
            process_json_surface(pa, j, sh);
          }
          if (oshell == true)
          {
            oshell = false;
            s->set_oshell(sh);
          }
          else
            s->add_ishell(sh);
        }
        cs->add_solid(s);
      }
      co->add_primitive(cs);
    }
    else if (g["type"] == "GeometryInstance") 
    {
      int gti = g["template"];
      GeometryTemplate* g2 = lsGTs[gti];
      co->add_primitive(g2);
    }
    idgeom++;
  }
}

void read_file_json(std::string &ifile, std::vector<Feature*>& lsFeatures, IOErrors& errs, double tol_snap)
{
  std::ifstream input(ifile);
  json j;
  try 
  {
    input >> j;
  }
  catch (nlohmann::detail::parse_error e) 
  {
    errs.add_error(901, "Input file not a valid JSON file.");
    return;
  }
  // TODO: other validation for CityJSON or just let it crash?
  if (j["type"] == "CityJSON") {
    errs.set_input_file_type("CityJSON");
    parse_cityjson(j, lsFeatures, tol_snap);
  } 
  else if (j["type"] == "tu3djson") {
    errs.set_input_file_type("tu3djson");
    *val3ditycout << "tu3djson input file" << std::endl;
    *val3ditycout << "# Features found: " << j["features"].size() << std::endl;
    parse_tu3djson(j, lsFeatures, tol_snap);
  }
  else if ( (j["type"] == "Feature") || (j["type"] == "FeatureCollection") ) {
    errs.set_input_file_type("JSON-FG");
    *val3ditycout << "JSON-FG input file" << std::endl;
    if (j["type"] == "Feature")
      *val3ditycout << "# Features found: " << j["features"].size() << std::endl;
    else
      *val3ditycout << "# Features found: " << j["features"].size() << std::endl;
    parse_jsonfg(j, lsFeatures, tol_snap, errs);
  }
  else {
    errs.add_error(904, "Input file type not a supported JSON file (only: CityJSON, JSON-FG, and tu3djson ).");
    return;  
  }
}

void read_file_jsonl(std::string &ifile, std::vector<Feature*>& lsFeatures, IOErrors& errs, double tol_snap)
{
  *val3ditycout << "CityJSONL input file" << std::endl;
  std::ifstream infile(ifile.c_str(), std::ifstream::in);
  if (!infile)
  {
    errs.add_error(901, "Input file not found.");
    return;
  }
  //-- read and store the GeometryTemplates
  std::vector<GeometryTemplate*> lsGTs;
  //-- transform
  json jtransform;
  std::string l;
  int linecount = 0;
  while (std::getline(infile, l)) 
  {
    std::istringstream iss(l);
    json j;
    try 
    {
      iss >> j;
    }
    catch (nlohmann::detail::parse_error e) 
    {
      std::string s = "Input file has invalid JSON at line #" + std::to_string(linecount);
      errs.add_error(901, s);
      break;
    }
    if (j["type"] == "CityJSON") {
      if (j.count("geometry-templates") == 1)
        process_cityjson_geometrytemplates(j["geometry-templates"], lsGTs, tol_snap);
      if (j.count("transform") == 0) {
        std::string s = "Input file first line has no \"transform\" property";
        errs.add_error(901, s);
        break;
      } else {
        jtransform = j["transform"];
      }
    }
    if (j["type"] == "CityJSONFeature") {
      errs.set_input_file_type("CityJSONL");
      j["transform"] = jtransform; //-- add transform b/c BuildingPart overlap uses a tolerance
      parse_cityjsonl(j, lsFeatures, tol_snap, lsGTs);
    }
    linecount++;
  }
}

void parse_cityjson(json& j, std::vector<Feature*>& lsFeatures, double tol_snap)
{
  *val3ditycout << "CityJSON input file" << std::endl;
  *val3ditycout << "# City Objects found: " << j["CityObjects"].size() << std::endl;
  //-- compute (_minx, _miny)
  compute_min_xy(j);
  //-- read and store the GeometryTemplates
  std::vector<GeometryTemplate*> lsGTs;
  if (j.count("geometry-templates") == 1)
  {
    process_cityjson_geometrytemplates(j["geometry-templates"], lsGTs, tol_snap);
  }
  //-- process each CO
  for (json::iterator it = j["CityObjects"].begin(); it != j["CityObjects"].end(); ++it) 
  {
    //-- BuildingParts geometries are put with those of a Building
    if (it.value()["type"] == "BuildingPart")
      continue;
    CityObject* co = new CityObject(it.key(), it.value()["type"]);
    process_json_geometries_of_co(it.value(), co, lsGTs, j, tol_snap);
    //-- if Building has Parts, put them here in _lsPrimitives
    if ( (it.value()["type"] == "Building") && (it.value().count("children") != 0) ) 
    {
      for (std::string bpid : it.value()["children"])
      {
        process_json_geometries_of_co(j["CityObjects"][bpid], co, lsGTs, j, tol_snap);
      }
    }
    lsFeatures.push_back(co);
  }
}

void parse_cityjsonl(json& j, std::vector<Feature*>& lsFeatures, double tol_snap, std::vector<GeometryTemplate*>& lsGTs)
{
  //-- compute (_minx, _miny)
  compute_min_xy(j);
  //-- process each CO
  for (json::iterator it = j["CityObjects"].begin(); it != j["CityObjects"].end(); ++it) 
  {
    //-- BuildingParts geometries are put with those of a Building
    if (it.value()["type"] == "BuildingPart")
      continue;
    CityObject* co = new CityObject(it.key(), it.value()["type"]);
    process_json_geometries_of_co(it.value(), co, lsGTs, j, tol_snap);
    //-- if Building has Parts, put them here in _lsPrimitives
    if ( (it.value()["type"] == "Building") && (it.value().count("children") != 0) ) 
    {
      for (std::string bpid : it.value()["children"])
      {
        process_json_geometries_of_co(j["CityObjects"][bpid], co, lsGTs, j, tol_snap);
      }
    }
    lsFeatures.push_back(co);
  }
}


void process_cityjson_geometrytemplates(json& j, std::vector<GeometryTemplate*>& lsGTs, double tol_snap)
{
  int count = 0;
  for (auto& jt : j["templates"])
  {
    GeometryTemplate* gt = new GeometryTemplate(std::to_string(count));
    if  (jt["type"] == "Solid")
    {
      Solid* s = new Solid("0");
      bool oshell = true;
      int c = 0;
      for (auto& shell : jt["boundaries"]) 
      {
        Surface* sh = new Surface(c, tol_snap);
        c++;
        for (auto& polygon : shell) { 
          std::vector< std::vector<int> > pa = polygon;
          process_json_surface_geometrytemplate(pa, j, sh);
        }
        if (oshell == true)
        {
          oshell = false;
          s->set_oshell(sh);
        }
        else
          s->add_ishell(sh);
      }
      gt->add_primitive(s);
    }
    else if ( (jt["type"] == "MultiSurface") || (jt["type"] == "CompositeSurface") ) 
    {
      Surface* sh = new Surface(-1, tol_snap);
      for (auto& p : jt["boundaries"]) 
      { 
        std::vector< std::vector<int> > pa = p;
        process_json_surface_geometrytemplate(pa, j, sh);
      }
      if (jt["type"] == "MultiSurface")
      {
        MultiSurface* ms = new MultiSurface("0");
        ms->set_surface(sh);
        gt->add_primitive(ms);
      }
      else
      {
        CompositeSurface* cs = new CompositeSurface("0");
        cs->set_surface(sh);
        gt->add_primitive(cs);
      }
    }
    lsGTs.push_back(gt);
    count++;
  }
}


void process_json_surface_geometrytemplate(std::vector< std::vector<int> >& pgn, json& j, Surface* sh)
{
  std::vector< std::vector<int> > pgnids;
  for (auto& r : pgn)
  {
    std::vector<int> newr;
    for (auto& i : r)
    {
      double x;
      double y;
      double z;
      x = double(j["vertices-templates"][i][0]);
      y = double(j["vertices-templates"][i][1]);
      z = double(j["vertices-templates"][i][2]);
      Point3 p3(x, y, z);
      newr.push_back(sh->add_point(p3));
    }
    pgnids.push_back(newr);
  }
  sh->add_face(pgnids);
}


void process_gml_file_indoorgml(pugi::xml_document& doc, std::vector<Feature*>& lsFeatures, std::map<std::string, pugi::xpath_node>& dallpoly, IOErrors& errs, double tol_snap)
{
  //-- 0. read the header of the file and find its gml:name, if any
  std::string nameim = "";
  if (doc.first_child().attribute("gml:id") != 0) 
    nameim = doc.first_child().attribute("gml:id").value();
  else 
    nameim = "MyIndoorModel";
  IndoorModel* im = new IndoorModel(nameim);
  lsFeatures.push_back(im);
    
  //-- 1. read each cellSpaceMember in the file (the primal objects)
  //--    these can have different names, depending on the Extensions/ADEs used
  std::string s = ".//" + NS["indoorgml"] + "cellSpaceMember";
  pugi::xpath_node_set nn = doc.select_nodes(s.c_str());
  int pcounter = 0;
  for (pugi::xpath_node_set::const_iterator it = nn.begin(); it != nn.end(); ++it)
  {
    pugi::xml_node cs = it->node().first_child();
    std::string theid   = "";
    std::string duality = "";
    //-- get the ID
    if (cs.attribute("gml:id") != 0) {
      // *val3ditycout << "\t" << it->node().attribute("gml:id").value() << std::endl;
      theid = cs.attribute("gml:id").value();
    }
    else 
      theid = ("MISSING_ID_" + std::to_string(pcounter));
    //-- get the duality pointer (max one, sweet)
    s = NS["indoorgml"] + "duality";
    for (pugi::xml_node child : cs.children(s.c_str()))
    {
      if (child.attribute("xlink:href") != 0) {
        std::string s = child.attribute("xlink:href").value();
        if (s[0] == '#')
          s = s.substr(1);
        duality = s;
      }
    }
    // IndoorCell* cell = new IndoorCell(theid, duality);
    //-- get the geometry, either Solid or Surface
    s = NS["indoorgml"] + "cellSpaceGeometry";
    Solid* sol;
    for (pugi::xml_node child : cs.children(s.c_str()))
    {
      s = NS["indoorgml"] + "Geometry3D";
      for (pugi::xml_node child2 : child.children(s.c_str()))
      {
        s = NS["gml"] + "Solid";
        for (pugi::xml_node child3 : child2.children(s.c_str()))
        {
          // *val3ditycout << "Solid: " << child3.attribute("gml:id").value() << std::endl;
          sol = process_gml_solid(child3, dallpoly, tol_snap, errs);
          if (sol->get_id() == "")
            sol->set_id("MISSING_ID");
          // cell->add_primitive(sol);
        }
      }
    }
    pcounter++;
    im->add_cell(theid, sol, duality, cs.name());
  }

  //-- 2. read the dual graphs (yes there can be more than one) 
  s = ".//" + NS["indoorgml"] + "SpaceLayer";
  nn = doc.select_nodes(s.c_str());
  for (pugi::xpath_node_set::const_iterator it = nn.begin(); it != nn.end(); ++it)
  {
    std::string idg;
    if (it->node().attribute("gml:id") != 0) {
      idg = it->node().attribute("gml:id").value();
    }
    else 
      idg = "";
    // IndoorGraph ig(idg);
    IndoorGraph* ig = new IndoorGraph(idg);
    //-- fetch all the edges
    std::map<std::string, std::tuple<std::string,std::string>> edges;
    s = ".//" + NS["indoorgml"] + "Transition";
    pugi::xpath_node_set ntr = it->node().select_nodes(s.c_str());
    for (pugi::xpath_node_set::const_iterator it = ntr.begin(); it != ntr.end(); ++it)
    {
      std::string theid = it->node().attribute("gml:id").value();
      s = NS["indoorgml"] + "connects";
      std::vector<std::string> connects;
      for (pugi::xml_node child : it->node().children(s.c_str()))
      {
        if (child.attribute("xlink:href") != 0) {
          std::string s = child.attribute("xlink:href").value();
          if (s[0] == '#')
            s = s.substr(1);
          connects.push_back(s);
        }
      }
      edges[theid] = std::make_tuple(connects[0], connects[1]);
    }
    //-- fetch all the nodes
    s = ".//" + NS["indoorgml"] + "State";
    pugi::xpath_node_set nstate = it->node().select_nodes(s.c_str());
//    pugi::xpath_node_set nstate = doc.select_nodes(s.c_str());
    for (pugi::xpath_node_set::const_iterator it = nstate.begin(); it != nstate.end(); ++it)
    {
      // *val3ditycout << "---\n" << it->node().attribute("gml:id").value() << std::endl;
      std::string vid = it->node().attribute("gml:id").value();
      s = NS["indoorgml"] + "duality";
      std::string vdual;
      pugi::xml_node child = it->node().child(s.c_str());
      if (child.attribute("xlink:href") != 0) {
        vdual = child.attribute("xlink:href").value();
        if (vdual[0] == '#')
          vdual = vdual.substr(1);
        // *val3ditycout << "dual node: " << vdual << std::endl;
      }
      s = NS["indoorgml"] + "connects";
      std::vector<std::string> vadj;
      for (pugi::xml_node child : it->node().children(s.c_str()))
      {
        if (child.attribute("xlink:href") != 0) {
          std::string s = child.attribute("xlink:href").value();
          if (s[0] == '#')
            s = s.substr(1);
          if (std::get<1>(edges[s]) != vid)
            vadj.push_back(std::get<1>(edges[s]));
        }
      }
      s = ".//" + NS["gml"] + "pos";
      pugi::xpath_node n = it->node().select_node(s.c_str());
      // *val3ditycout << n.node().child_value() << std::endl;
      
      std::string buf;
      std::stringstream ss(n.node().child_value());
      std::vector<std::string> tokens;
      while (ss >> buf)
        tokens.push_back(buf);
      // long double x = std::stold(tokens[0]);
      // long double y = std::stold(tokens[1]);
      // long double z = std::stold(tokens[2]);
      // TODO: use long double?
      ig->add_vertex(vid, 
                     std::stold(tokens[0]), 
                     std::stold(tokens[1]), 
                     std::stold(tokens[2]),
                     vdual,
                     vadj);
    }
    im->add_graph(ig);
  }
  // if (im->get_number_graphs() == 0)
  //   *val3ditycout << "\t!!! No dual graph was found in the input file" << std::endl;
}


void set_min_xy(double minx, double miny)
{
  _minx = minx;
  _miny = miny;
  // *val3ditycout << "Translating all coordinates by (-" << _minx << ", -" << _miny << ")" << std::endl;
  Primitive::set_translation_min_values(_minx, _miny);
  Surface::set_translation_min_values(_minx, _miny);
}

void compute_min_xy(json& j)
{
  for (auto& v : j["vertices"])
  {
    if (v[0] < _minx)
      _minx = v[0];
    if (v[1] < _miny)
      _miny = v[1];
  }
  if (j.count("transform") != 0)
  {
    _minx = (_minx * double(j["transform"]["scale"][0])) + double(j["transform"]["translate"][0]);
    _miny = (_miny * double(j["transform"]["scale"][1])) + double(j["transform"]["translate"][1]);
  }
  // *val3ditycout << "Translating all coordinates by (-" << _minx << ", -" << _miny << ")" << std::endl;
  Primitive::set_translation_min_values(_minx, _miny);
  Surface::set_translation_min_values(_minx, _miny);
}


void compute_min_xy(pugi::xml_document& doc)
{
  std::string s = "//" + NS["gml"] + "posList";
  pugi::xpath_node_set nall = doc.select_nodes(s.c_str());
  for (auto& each : nall) 
  {
    std::string buf;
    std::stringstream ss(each.node().child_value());
    std::vector<std::string> coords;
    while (ss >> buf)
      coords.push_back(buf);
    for (int i = 0; i < coords.size(); i += 3)
    {
      if (std::stod(coords[0]) < _minx)
        _minx = std::stod(coords[0]);
      if (std::stod(coords[1]) < _miny)
        _miny = std::stod(coords[1]);
    }
  }
  s = "//" + NS["gml"] + "pos";
  nall = doc.select_nodes(s.c_str());
  for (auto& each : nall) 
  {
    std::string buf;
    std::stringstream ss(each.node().child_value());
    std::vector<std::string> tokens;
    while (ss >> buf)
      tokens.push_back(buf);
    if (std::stod(tokens[0]) < _minx)
      _minx = std::stod(tokens[0]);
    if (std::stod(tokens[1]) < _miny)
      _miny = std::stod(tokens[1]);
  }
  // *val3ditycout << "Translating all coordinates by (-" << _minx << ", -" << _miny << ")" << std::endl;
  Primitive::set_translation_min_values(_minx, _miny);
  Surface::set_translation_min_values(_minx, _miny);
}


void read_file_gml(std::string &ifile, std::vector<Feature*>& lsFeatures, IOErrors& errs, double tol_snap)
{
  *val3ditycout << "Reading file: " << ifile << std::endl;
  pugi::xml_document doc;
  pugi::xml_parse_result result = doc.load_file(ifile.c_str());
  if (!result) {
    if (result.status == pugi::status_file_not_found)
      errs.add_error(901, "Input file not found.");
    else
      errs.add_error(901, result.description());
    return;
  }
  //-- parse namespace
  pugi::xml_node ncm = doc.first_child();
  get_namespaces(ncm); //-- results in global variable NS in this unit
  if ( (NS.count("indoorgml") != 0) && (ncm.name() == (NS["indoorgml"] + "IndoorFeatures")) ) {
    *val3ditycout << "IndoorGML input file" << std::endl;
    //-- find (_minx, _miny)
    compute_min_xy(doc);
    //-- build dico of xlinks for <gml:Polygon>
    std::map<std::string, pugi::xpath_node> dallpoly;
    build_dico_xlinks(doc, dallpoly, errs);
    errs.set_input_file_type("IndoorGML");
    process_gml_file_indoorgml(doc, lsFeatures, dallpoly, errs, tol_snap);
  }
  else
  {
    errs.add_error(904, "GML files not supported (yes that includes CityGML files ==> upgrade to CityJSON)");
  }
}

std::map<std::string, std::string> get_namespaces(pugi::xml_node& root) {
  for (pugi::xml_attribute attr = root.first_attribute(); attr; attr = attr.next_attribute()) {
    std::string name = attr.name();
    if (name.find("xmlns") != std::string::npos) {
      // *val3ditycout << attr.name() << "=" << attr.value() << std::endl;
      std::string value = attr.value();
      std::string sns;
      if (value.find("http://www.opengis.net/citygml/") != std::string::npos) {
        sns = "citygml";
      }
      else if ( (value.find("http://www.opengis.net/gml") != std::string::npos) &&
                (value.find("http://www.opengis.net/gmlcov") == std::string::npos) ) {
        sns = "gml";
      }
      else if ( (value.find("http://www.opengis.net/indoorgml/1") != std::string::npos) &&
                (value.find("core") != std::string::npos) ) {
        sns = "indoorgml";
      }
      else if (value.find("http://www.w3.org/1999/xlink") != std::string::npos)
        sns = "xlink";
      else
        sns = "";
      if (sns != "") {
        size_t pos = name.find(":");
        if (pos == std::string::npos) 
          NS[sns] = "";
        else 
          NS[sns] = name.substr(pos + 1) + ":";
      }    
    }
  }
  std::map<std::string, std::string> tmp = NS;
  return tmp;
}


void build_dico_xlinks(pugi::xml_document& doc, std::map<std::string, pugi::xpath_node>& dallpoly, IOErrors& errs)
{
  std::string s = "//" + NS["gml"] + "Polygon" + "[@" + NS["gml"] + "id]";
  pugi::xpath_node_set nallpoly = doc.select_nodes(s.c_str());
  if (nallpoly.size() > 0)
   *val3ditycout << "XLinks found, resolving them..." << std::flush;
  for (pugi::xpath_node_set::const_iterator it = nallpoly.begin(); it != nallpoly.end(); ++it)
    dallpoly[it->node().attribute("gml:id").value()] = *it;
  //-- for <gml:OrientableSurface>
  s = "//" + NS["gml"] + "OrientableSurface" + "[@" + NS["gml"] + "id" + "]";
  pugi::xpath_node_set nallosurf = doc.select_nodes(s.c_str());
  for (pugi::xpath_node_set::const_iterator it = nallosurf.begin(); it != nallosurf.end(); ++it)
    dallpoly[it->node().attribute("gml:id").value()] = *it;
  //-- checking xlinks validity now, not to be bitten later
  s = "//" + NS["gml"] + "surfaceMember" + "[@" + NS["xlink"] + "href" + "]";
  pugi::xpath_node_set nsmxlink = doc.select_nodes(s.c_str());
  for (pugi::xpath_node_set::const_iterator it = nsmxlink.begin(); it != nsmxlink.end(); ++it) 
  {
    std::string k = it->node().attribute("xlink:href").value();
    if (k[0] == '#')
      k = k.substr(1);
    if (dallpoly.count(k) == 0) 
    {
      std::string r = "One XLink couldn't be resolved (";
      r += it->node().attribute("xlink:href").value();
      r += ")";
      errs.add_error(901, r);
      return;
    }
  }
  if (nallpoly.size() > 0)
    *val3ditycout << "done." << std::endl;
}


Surface* parse_poly(std::istream &input, int shellid, IOErrors& errs)
{
  //-- read the points
  int num, tmpint;
  float tmpfloat;
  double x, y, z;
  input >> num >> tmpint >> tmpint >> tmpint;
  //-- compute (_minx, _miny)
  for (int i = 0; i < num; i++)
  {
    input >> tmpint >> x >> y >> z;
    if (x < _minx)
      _minx = x;
    if (y < _miny)
      _miny = y;
  }
  *val3ditycout << "Translating all coordinates by (-" << _minx << ", -" << _miny << ")" << std::endl;
  Primitive::set_translation_min_values(_minx, _miny);
  Surface::set_translation_min_values(_minx, _miny);
  input.clear();
  input.seekg(0);
  input >> num >> tmpint >> tmpint >> tmpint;
  //-- read verticess
  Surface* sh = new Surface(shellid);  
  for (int i = 0; i < num; i++)
  {
    input >> tmpint >> x >> y >> z;
    x -= _minx;
    y -= _miny;
    Point3 p(x, y, z);
    sh->add_point(p);
  }
  //-- read the facets
  input >> num >> tmpint;
  int numf, numpt, numholes;
  string s;
  for (int i = 0; i < num; i++)
  {
    numholes = 0;
    input >> numf;
    while(true) {
      if (input.peek() == '\n')
        break;
      else if (input.peek() == ' ')
        input.ignore();
      else
        input >> numholes;
    }
    //-- read oring (there's always one and only one)
    input >> numpt;
    if (numpt == -1) {
      sh->add_error(103, std::to_string(i));
      continue;
    }
    std::vector<int> ids(numpt);
    for (int k = 0; k < numpt; k++)
      input >> ids[k];
    std::vector< std::vector<int> > pgnids;
    pgnids.push_back(ids);
    //-- check for irings
    for (int j = 1; j < numf; j++)
    {
      input >> numpt;
      if (numpt == -1) {
        sh->add_error(103, std::to_string(i));
        continue;
      }
      std::vector<int> ids(numpt);
      for (int l = 0; l < numpt; l++)
        input >> ids[l];
      pgnids.push_back(ids);
    }
    //-- skip the line about points defining holes (if present)
    for (int j = 0; j < numholes; j++)
      input >> tmpint >> tmpfloat >> tmpfloat >> tmpfloat;
    sh->add_face(pgnids);
  }
  return sh;
}

void printProgressBar(int percent) {
  std::string bar;
  for(int i = 0; i < 50; i++){
    if( i < (percent / 2)) {
      bar.replace(i, 1, "=");
    }
    else if( i == (percent / 2)) {
      bar.replace(i, 1, ">");
    }
    else{
      bar.replace(i, 1, " ");
    }
  }
  *val3ditycout << "\r" "[" << bar << "] ";
  (*val3ditycout).width(3);
  *val3ditycout << percent << "%     " << std::flush;
}


Surface* parse_off(std::istream &input, int shellid, IOErrors& errs, double tol_snap)
{
  //-- read the points
  int numpt, numf, tmpint;
  std::string s;
  input >> s;
  input >> numpt >> numf >> tmpint;
  if ( (s != "OFF") || (numpt <= 0) ) {
    errs.add_error(901, "Input file not a valid OFF file.");
    return NULL;
  }
  //-- compute (_minx, _miny)
  for (int i = 0; i < numpt; i++)
  {
    double x, y, z;
    input >> x >> y >> z;
    if (x < _minx)
      _minx = x;
    if (y < _miny)
      _miny = y;
  }
  *val3ditycout << "Translating all coordinates by (-" << _minx << ", -" << _miny << ")" << std::endl;
  Primitive::set_translation_min_values(_minx, _miny);
  Surface::set_translation_min_values(_minx, _miny);
  //-- reset the file
  input.clear();
  input.seekg(0);
  input >> s;
  input >> numpt >> numf >> tmpint;
  //-- read the points
  std::vector<int> newi;
  Surface* sh = new Surface(shellid, tol_snap);
  for (int i = 0; i < numpt; i++)
  {
    double x, y, z;
    input >> x >> y >> z;
    x -= _minx;
    y -= _miny;
    Point3 p(x, y, z);
    newi.push_back(sh->add_point(p));
  }
  //-- read the facets
  for (int i = 0; i < numf; i++)
  {
    input >> tmpint;
    std::vector<int> ids(tmpint);
    if (ids.empty() == true)
    {
      errs.add_error(901, "Some surfaces not defined correctly or are empty");
      return NULL;
    }
    for (int k = 0; k < tmpint; k++) {
      int t;
      input >> t;
      ids[k] = newi[t];
    }
    std::vector< std::vector<int> > pgnids;
    pgnids.push_back(ids);
    sh->add_face(pgnids);
  }
  return sh;
}


void parse_obj(std::istream &input, std::vector<Feature*>& lsFeatures, Primitive3D prim3d, IOErrors& errs, double tol_snap)
{
  //-- find (minx, miny)
  std::string l;
  while (std::getline(input, l)) 
  {
    std::istringstream iss(l);
    if (l.substr(0, 2) == "v ") {
      std::string tmp;
      double x, y, z;
      iss >> tmp >> x >> y >> z;
      if (x < _minx)
        _minx = x;
      if (y < _miny)
        _miny = y;
    }
  }
  *val3ditycout << "Translating all coordinates by (-" << _minx << ", -" << _miny << ")" << std::endl;
  Primitive::set_translation_min_values(_minx, _miny);
  Surface::set_translation_min_values(_minx, _miny);
  //-- read again file and parse everything
  input.clear();
  input.seekg(0);
  int primid = 0;
  std::string primids = "";
  Surface* sh = new Surface(0, tol_snap);
  std::vector<Point3*> allvertices;
  GenericObject* o = new GenericObject("none");
  while (std::getline(input, l)) {
    std::istringstream iss(l);
    if (l.substr(0, 2) == "v ") {
      std::string tmp;
      double x, y, z;
      iss >> tmp >> x >> y >> z;
      x -= _minx;
      y -= _miny;
      Point3 *p = new Point3(x, y, z);
      allvertices.push_back(p);
    }
    else if (l.substr(0, 2) == "o ") {
      if (sh->is_empty() == true) { //-- for the first "o" before any faces
        std::string tmp;
        iss >> tmp >> primids; 
      }
      else{
        std::string theid = "";
        if (primids != "")
          theid = primids;
        else
          theid = std::to_string(primid);
        if (prim3d == SOLID)
        {
          Solid* sol = new Solid(theid);
          sol->set_oshell(sh);
          o->add_primitive(sol);
        }
        else if ( prim3d == COMPOSITESURFACE)
        {
          CompositeSurface* cs = new CompositeSurface(theid);
          cs->set_surface(sh);
          o->add_primitive(cs);
        }
        else if (prim3d == MULTISURFACE)
        {
          MultiSurface* ms = new MultiSurface(theid);
          ms->set_surface(sh);
          o->add_primitive(ms);
        }
        primid++;
        std::string tmp;
        iss >> tmp >> primids; 
        sh = new Surface(0, tol_snap);
      }
    }
    else if (l.substr(0, 2) == "f ") {
      std::vector<int> r;
      std::string tmp;
      iss >> tmp;
      while (iss)
      {
        tmp.clear();
        iss >> tmp;
        if (tmp.compare("\\") == 0) {
          std::getline(input, l);
          iss.str(l);
          continue;
        }
        if (tmp.empty() == false) {
          std::size_t k = tmp.find("/");
          int index = std::stoi(tmp.substr(0, k));
          if (index == 0) {
            errs.add_error(901, "OBJ files are 1-indexed, vertex '0' found.");
            return;
          }
          if (index > allvertices.size()) {
            std::string r = "Vertex #";
            r += std::to_string(index);
            r += " doesn't exist in the input file.";
            errs.add_error(901, r);
            return;
          }
          Point3* tp = allvertices[index - 1];
          r.push_back(sh->add_point(*tp));
        }
      }
      std::vector< std::vector<int> > pgnids;
      pgnids.push_back(r);
      sh->add_face(pgnids);
    }
  }
  if (sh->is_empty() == true) {
    errs.add_error(902, "Some surfaces are not defined correctly or are empty");
    return;
  }
  std::string theid = "";
  if (primids != "")
    theid = primids;
  else
    theid = std::to_string(primid);
  if (prim3d == SOLID)
  {
    Solid* sol = new Solid(theid);
    sol->set_oshell(sh);
    o->add_primitive(sol);
  }
  else if ( prim3d == COMPOSITESURFACE)
  {
    CompositeSurface* cs = new CompositeSurface(theid);
    cs->set_surface(sh);
    o->add_primitive(cs);
  }
  else if (prim3d == MULTISURFACE)
  {
    MultiSurface* ms = new MultiSurface(theid);
    ms->set_surface(sh);
    o->add_primitive(ms);
  }
  for (auto& each : allvertices)
    delete each;
  allvertices.clear();
  lsFeatures.push_back(o);  
}


void parse_jsonfg(json& j, std::vector<Feature*>& lsFeatures, double tol_snap, IOErrors& errs)
{
  //-- TODO: not translation for json-fg, is that okay?
  set_min_xy(0.0, 0.0);
  if (j["type"] == "Feature") {
    parse_jsonfg_onefeature(j, lsFeatures, tol_snap, 0, errs);
    return;
  }  
  else { //-- then it's FeatureCollection
    int counter = 0;
    for (auto& f : j["features"]) {
      parse_jsonfg_onefeature(f, lsFeatures, tol_snap, counter, errs);
      counter++;
    }
  }
}


void parse_jsonfg_onefeature(json& j, std::vector<Feature*>& lsFeatures, double tol_snap, int counter, IOErrors& errs) 
{
  //-- find the id, counter is not present
  std::string sid = std::to_string(counter);
  if (j.contains("id")) {
    if (j["id"].is_number())
      sid = std::to_string(j["id"].get<double>());
    else
      sid = j["id"];
  }
  //-- create the object
  GenericObject* go = new GenericObject(sid);
  if (j["place"]["type"] == "Polyhedron") {
    Solid* s = new Solid();
    bool oshell = true;
    int c = 0;
    for (auto& shell : j["place"]["coordinates"]) 
    {
      Surface* sh = new Surface(c, tol_snap);
      c++;
      for (auto& polygon : shell) { 
        std::vector<std::vector<std::vector<double>>> pa = polygon;
        process_jsonfg_surface(pa, sh, errs);
      }
      if (oshell == true)
      {
        oshell = false;
        s->set_oshell(sh);
      }
      else
        s->add_ishell(sh);
    }
    go->add_primitive(s);
  }
  else if (j["place"]["type"] == "MultiPolyhedron") {
    MultiSolid* ms = new MultiSolid();
    for (auto& solid : j["place"]["coordinates"]) 
    {
      Solid* s = new Solid();
      bool oshell = true;
      for (auto& shell : solid) 
      {
        Surface* sh = new Surface(-1, tol_snap);
        for (auto& polygon : shell) { 
          std::vector<std::vector<std::vector<double>>> pa = polygon;
          process_jsonfg_surface(pa, sh, errs);
        }
        if (oshell == true)
        {
          oshell = false;
          s->set_oshell(sh);
        }
        else
          s->add_ishell(sh);
      }
      ms->add_solid(s);
    }
    go->add_primitive(ms);
  }
  lsFeatures.push_back(go);
}

void parse_tu3djson(json& j, std::vector<Feature*>& lsFeatures, double tol_snap)
{

  //-- TODO: not translation for tu3djson, is that okay?
  set_min_xy(0.0, 0.0);
  int counter = 0;
  for (auto& f : j["features"]) {
    GenericObject* go = new GenericObject(std::to_string(counter));
    counter++;
    if  (f["geometry"]["type"] == "Solid")
    {
      Solid* s = new Solid();
      bool oshell = true;
      int c = 0;
      for (auto& shell : f["geometry"]["boundaries"]) 
      {
        Surface* sh = new Surface(c, tol_snap);
        c++;
        for (auto& polygon : shell) { 
          std::vector< std::vector<int> > pa = polygon;
          process_json_surface(pa, f["geometry"], sh);
        }
        if (oshell == true)
        {
          oshell = false;
          s->set_oshell(sh);
        }
        else
          s->add_ishell(sh);
      }
      go->add_primitive(s);
    }
    else if ( (f["geometry"]["type"] == "MultiSurface") || (f["geometry"]["type"] == "CompositeSurface") ) 
    {
      Surface* sh = new Surface(-1, tol_snap);
      for (auto& p : f["geometry"]["boundaries"]) 
      { 
        std::vector< std::vector<int> > pa = p;
        process_json_surface(pa, f["geometry"], sh);
      }
      if (f["geometry"]["type"] == "MultiSurface")
      {
        MultiSurface* ms = new MultiSurface();
        ms->set_surface(sh);
        go->add_primitive(ms);
      }
      else
      {
        CompositeSurface* cs = new CompositeSurface();
        cs->set_surface(sh);
        go->add_primitive(cs);
      }
    }
    else if (f["geometry"]["type"] == "MultiSolid") 
    {
      MultiSolid* ms = new MultiSolid();
      for (auto& solid : f["geometry"]["boundaries"]) 
      {
        Solid* s = new Solid();
        bool oshell = true;
        for (auto& shell : solid) 
        {
          Surface* sh = new Surface(-1, tol_snap);
          for (auto& polygon : shell) { 
            std::vector< std::vector<int> > pa = polygon;
            process_json_surface(pa, f["geometry"], sh);
          }
          if (oshell == true)
          {
            oshell = false;
            s->set_oshell(sh);
          }
          else
            s->add_ishell(sh);
        }
        ms->add_solid(s);
      }
      go->add_primitive(ms);
    }
    else if (f["geometry"]["type"] == "CompositeSolid") 
    {
      CompositeSolid* cs = new CompositeSolid();
      for (auto& solid : f["geometry"]["boundaries"]) 
      {
        Solid* s = new Solid();
        bool oshell = true;
        for (auto& shell : solid) 
        {
          Surface* sh = new Surface(-1, tol_snap);
          for (auto& polygon : shell) { 
            std::vector< std::vector<int> > pa = polygon;
            process_json_surface(pa, f["geometry"], sh);
          }
          if (oshell == true)
          {
            oshell = false;
            s->set_oshell(sh);
          }
          else
            s->add_ishell(sh);
        }
        cs->add_solid(s);
      }
      go->add_primitive(cs);
    } 
    lsFeatures.push_back(go);
  }
}


void parse_tu3djson_onegeom(json& j, std::vector<Feature*>& lsFeatures, double tol_snap)
{
  //-- TODO: not translation for tu3djson, is that okay?
  set_min_xy(0.0, 0.0);
  GenericObject* go = new GenericObject("0");
  if  (j["type"] == "Solid")
  {
    Solid* s = new Solid();
    bool oshell = true;
    int c = 0;
    for (auto& shell : j["boundaries"]) 
    {
      Surface* sh = new Surface(c, tol_snap);
      c++;
      for (auto& polygon : shell) { 
        std::vector< std::vector<int> > pa = polygon;
        process_json_surface(pa, j, sh);
      }
      if (oshell == true)
      {
        oshell = false;
        s->set_oshell(sh);
      }
      else
        s->add_ishell(sh);
    }
    go->add_primitive(s);
  }
  else if ( (j["type"] == "MultiSurface") || (j["type"] == "CompositeSurface") ) 
  {
    Surface* sh = new Surface(-1, tol_snap);
    for (auto& p : j["boundaries"]) 
    { 
      std::vector< std::vector<int> > pa = p;
      process_json_surface(pa, j, sh);
    }
    if (j["type"] == "MultiSurface")
    {
      MultiSurface* ms = new MultiSurface();
      ms->set_surface(sh);
      go->add_primitive(ms);
    }
    else
    {
      CompositeSurface* cs = new CompositeSurface();
      cs->set_surface(sh);
      go->add_primitive(cs);
    }
  }
  else if (j["type"] == "MultiSolid") 
  {
    MultiSolid* ms = new MultiSolid();
    for (auto& solid : j["boundaries"]) 
    {
      Solid* s = new Solid();
      bool oshell = true;
      for (auto& shell : solid) 
      {
        Surface* sh = new Surface(-1, tol_snap);
        for (auto& polygon : shell) { 
          std::vector< std::vector<int> > pa = polygon;
          process_json_surface(pa, j, sh);
        }
        if (oshell == true)
        {
          oshell = false;
          s->set_oshell(sh);
        }
        else
          s->add_ishell(sh);
      }
      ms->add_solid(s);
    }
    go->add_primitive(ms);
  }
  else if (j["type"] == "CompositeSolid") 
  {
    CompositeSolid* cs = new CompositeSolid();
    for (auto& solid : j["boundaries"]) 
    {
      Solid* s = new Solid();
      bool oshell = true;
      for (auto& shell : solid) 
      {
        Surface* sh = new Surface(-1, tol_snap);
        for (auto& polygon : shell) { 
          std::vector< std::vector<int> > pa = polygon;
          process_json_surface(pa, j, sh);
        }
        if (oshell == true)
        {
          oshell = false;
          s->set_oshell(sh);
        }
        else
          s->add_ishell(sh);
      }
      cs->add_solid(s);
    }
    go->add_primitive(cs);
  } 
  lsFeatures.push_back(go);
}


json get_report_json(std::string ifile, 
                     std::vector<Feature*>& lsFeatures,
                     std::string val3dity_version,
                     double snap_tol,
                     double overlap_tol,
                     double planarity_d2p_tol,
                     double planarity_n_tol,
                     IOErrors ioerrs)
{
  json jr;
  jr["type"] = "val3dity_report";
  jr["val3dity_version"] = val3dity_version; 
  jr["input_file"] = ifile;
  jr["input_file_type"] = ioerrs.get_input_file_type();
  //-- time
  std::time_t currentTime;
  struct tm * localTime;
  std::time (&currentTime);
  localTime = std::localtime ( &currentTime );
  char buffer[80];
  std::strftime(buffer, 80, "%FT%T (%Z)", localTime);
  jr["time"] = buffer;
  //-- user-defined params
  jr["parameters"];
  jr["parameters"]["snap_tol"] = snap_tol;
  jr["parameters"]["overlap_tol"] = overlap_tol;
  jr["parameters"]["planarity_d2p_tol"] = planarity_d2p_tol;
  jr["parameters"]["planarity_n_tol"] = planarity_n_tol;

  //-- primitives overview
  std::map<int, std::tuple<int,int> > prim_o; //-- <primID, total, valid>
  std::set<int> theprimitives;
  for (auto& f : lsFeatures)
    for (auto& p : f->get_primitives())
      theprimitives.insert(p->get_type());
  for (auto& each : theprimitives)
    prim_o[each] = std::make_tuple(0, 0);
  for (auto& f : lsFeatures) {
    for (auto& p : f->get_primitives()) {
      std::get<0>(prim_o[p->get_type()]) += 1;
      if (p->is_valid() == true) {
        std::get<1>(prim_o[p->get_type()]) += 1;
      }
    }
  }
  jr["primitives_overview"] = json::array();
  for (auto& each : prim_o) {
    json j;
    switch(each.first)
    {
      case 0: j["type"] = "Solid"; break;
      case 1: j["type"] = "CompositeSolid"; break;
      case 2: j["type"] = "MultiSolid"; break;
      case 3: j["type"] = "CompositeSurface"; break;
      case 4: j["type"] = "MultiSurface"; break;
      case 5: j["type"] = "GeometryTemplate"; break;
      case 9: j["type"] = "ALL"; break;
    }
    j["total"] = std::get<0>(each.second);
    j["valid"] = std::get<1>(each.second);
    jr["primitives_overview"].push_back(j);
  }

  //-- features overview
  std::map<std::string, std::tuple<int,int> > feat_o; //-- <featureID, total, valid>
  std::set<std::string> thefeatures;
  for (auto& f : lsFeatures)
    thefeatures.insert(f->get_type());
  for (auto& each : thefeatures)
    feat_o[each] = std::make_tuple(0, 0);
  for (auto& f : lsFeatures) {
    std::get<0>(feat_o[f->get_type()]) += 1;
    if (f->is_valid() == true) {
      std::get<1>(feat_o[f->get_type()]) += 1;
    }
  }
  jr["features_overview"] = json::array();
  for (auto& each : feat_o) {
    json j;
    j["type"] = each.first; 
    j["total"] = std::get<0>(each.second);
    j["valid"] = std::get<1>(each.second);
    jr["features_overview"].push_back(j);
  }

  //-- each of the features with their primitives listed
  jr["features"] = json::array();
  for (auto& f : lsFeatures)
    jr["features"].push_back(f->get_report_json());
  
  //-- dataset errors (9xx)
  jr["dataset_errors"] = json::array();
  if (ioerrs.has_errors() == true)
    jr["dataset_errors"] = ioerrs.get_report_json();

  //-- overview of errors
  std::set<int> unique_errors;
  for (auto& f : lsFeatures)
    for (auto& code : f->get_unique_error_codes())
      unique_errors.insert(code);
  jr["all_errors"] = json::array();
  for (auto& e : unique_errors)
    jr["all_errors"].push_back(e);
  for (auto& e : ioerrs.get_unique_error_codes())
    jr["all_errors"].push_back(e);

  bool bValid = true;
  for (auto& f : lsFeatures) {
    if (f->is_valid() == false) {
      bValid = false;
      break;
    }
  }
  if (ioerrs.has_errors() == true)
    bValid = false;
  jr["validity"] = bValid;

  return jr;
}


} // namespace val3dity
