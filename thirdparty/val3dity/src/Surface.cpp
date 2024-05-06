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

#include "Surface.h"

#include "geomtools.h"
#include "input.h"
#include "validate_shell.h"
#include "val3dity_ostream.h"
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <geos_c.h>
#include <sstream>

using namespace std;



namespace val3dity
{

double Surface::_shiftx = 0.0;
double Surface::_shifty = 0.0;

Surface::Surface(int id, double tol_snap)
{
  _id = id;
  _is_valid_2d = -1;
  _vertices_added = 0;
  _tol_snap = tol_snap;
}

Surface::~Surface()
{
  // TODO: clear memory properly
  _lsPts.clear();
}

int Surface::get_id()
{
  return _id;
}

CgalPolyhedron* Surface::get_cgal_polyhedron()
{
  return _polyhedron;
}

bool Surface::has_errors()
{
return !(_errors.empty());
}

bool Surface::is_empty()
{
  return (_lsPts.empty() || _lsFaces.empty());
}


void Surface::add_error(int code, std::string faceid, std::string info)
{
  std::tuple<std::string, std::string> a(faceid, info);
  _errors[code].push_back(a);
  *val3dityclog << "\tERROR " << code << ": " << ALL_ERRORS[code];
  if (faceid.empty() == false)
    *val3dityclog << " (face " << faceid << ")";
  *val3dityclog << std::endl;
  if (info.empty() == false)
    *val3dityclog << "\t[" << info << "]" << std::endl;
}

std::set<int> Surface::get_unique_error_codes()
{
  std::set<int> errs;
  for (auto& err : _errors)
  {
    errs.insert(std::get<0>(err));
  }
  return errs;
}


json Surface::get_report_json(std::string preid)
{
  json j;
  for (auto& err : _errors)
  {
    for (auto& e : _errors[std::get<0>(err)])    
    {
      json jj;
      jj["type"] = "Error";
      jj["code"] = std::get<0>(err);
      jj["description"] = ALL_ERRORS[std::get<0>(err)];
      if (preid != "") {
        jj["id"] = preid + " | " + "face:" + std::get<0>(e);
      } else {
        jj["id"] = "face:" + std::get<0>(e);
      }
      jj["info"] = std::get<1>(e);
      j.push_back(jj);
    }
  }
  return j;
}


std::string Surface::get_poly_representation()
{
  std::ostringstream s;
  //-- points
  s << _lsPts.size() << " 3 0 0" << std::endl;
  int i = 0;
  for (auto& p : _lsPts)
  {
    s << i << " " << setprecision(15) << p.x() << " " << p.y() << " " << p.z() << std::endl;
    i++;
  }
  //-- faces
  s << _lsFaces.size() << " 0" << std::endl;
  for (auto& f : _lsFaces)
  {
    s << f.size() << " " << (f.size() - 1) << std::endl;
    for (auto& r : f)
    {
      s << r.size() << " ";
      for (auto p : r)
      {
        s << p << " ";
      }
      s << std::endl;
    }
  }
  s << "0" << std::endl;
  s << "0" << std::endl;
  return s.str();
}

std::string Surface::get_off_representation()
{
  std::stringstream ss;
  ss << "OFF" << std::endl;
  int notr = 0;
  for (auto& f: _lsTr) 
    for (auto& t: f) 
      notr++;
  ss << _lsPts.size() << " " << notr << " 0" << std::endl;
  //-- points
  for (auto& p : _lsPts)
    ss << setprecision(15) << p.x() << " " << p.y() << " " << p.z() << std::endl;
  //-- triangles
  for (auto& f: _lsTr)
  {
    for (auto& t: f) 
    {
      ss << "3 " << t[0] << " " << t[1] << " " << t[2] << std::endl;
    }
  }
  return ss.str();
}


std::string Surface::get_coords_key(Point3* p)
{
  std::string s;
  if (_tol_snap < 0.0)
  {
    s = std::to_string(p->x()) + std::to_string(p->y()) + std::to_string(p->z());
  }
  else
  {
    int tol = int(1 / _tol_snap);
    s = std::to_string(int64(p->x() * tol)) + " " + std::to_string(int64(p->y() * tol)) + " " + std::to_string(int64(p->z() * tol));
  }
  return s;
}

bool Surface::were_vertices_merged_during_parsing()
{
  if (this->number_vertices() == this->_vertices_added)
    return false;
  else 
    return true;
}


int Surface::get_number_parsed_vertices()
{
  return _vertices_added;
}


int Surface::add_point(Point3 pi)
{
  _vertices_added += 1;
  int i = 0;
  for (auto& p : _lsPts)
  {
    auto d = CGAL::squared_distance(pi, p);
    if (d < (_tol_snap*_tol_snap)) {
      return i;
    }
    i++;
  }
  _lsPts.push_back(pi);
  return (_lsPts.size() - 1);
}

void Surface::add_face(std::vector< std::vector<int> > f, std::string id)
{
  _lsFaces.push_back(f);
  if (id.empty() == true)
    _lsFacesID.push_back(std::to_string(_lsFaces.size() - 1));
  else
    _lsFacesID.push_back(id);
}


int Surface::number_vertices()
{
  return _lsPts.size();
}


int Surface::number_faces()
{
  return _lsFaces.size();
}


bool Surface::triangulate_shell()
{
  *val3dityclog << "-----Triangulation of each surface" << std::endl;
  //-- read the facets
  size_t num = _lsFaces.size();
  for (int i = 0; i < static_cast<int>(num); i++)
  {
    // These are the number of rings on this facet
    size_t numf = _lsFaces[i].size();
    std::vector<int> &idsob = _lsFaces[i][0]; // helpful alias for the outer boundary
    if ( (numf == 1) && (idsob.size() == 3)) 
    {
      std::vector<int*> oneface;
      int* tr = new int[3];
      tr[0] = idsob[0] ;
      tr[1] = idsob[1];
      tr[2] = idsob[2];
      oneface.push_back(tr);
      _lsTr.push_back(oneface);
      continue;
    }
    //-- get projected CT
    std::vector<int*> tris = construct_ct_one_face(_lsFaces[i]);
    if (tris.size() == 0)
    {
      this->add_error(999, _lsFacesID[i], "face does not have an outer boundary.");
      return false;
    }
    _lsTr.push_back(tris);
  }
  return true;
}


std::vector<int*> 
Surface::construct_ct_one_face(const std::vector<std::vector<int>>& pgnids)
{
  std::vector<int*> re;

  std::vector<Point3> planepts;
  for (auto& each : pgnids[0]) {
    planepts.push_back(_lsPts[each]);
  }
  Plane bestfitplane = get_best_fitted_plane(planepts);
  //-- check orientation (for good normals for the output, pointing outwards)
  bool reversed = false;
  Polygon pgn;
  for (auto& each : pgnids[0]) {
    Point3 p = _lsPts[each];
    pgn.push_back(bestfitplane.to_2d(p));
  }
  //-- check orientation, this works all must be ccw
  if (pgn.is_counterclockwise_oriented() == false) {
    reversed = true;
  }
  CT ct;
  for (auto& ring : pgnids) {
    //-- make another *closed* ring for simplicity
    std::vector<int> r2 = ring;
    r2.push_back(ring.front());
    std::vector<int>::const_iterator it;
    std::vector<int>::iterator it2 = std::prev(r2.end());
    for (it = r2.begin(); it != it2; it++) {
      Point2 p0 = bestfitplane.to_2d(_lsPts[*it]);
      CT::Vertex_handle v0 = ct.insert(p0);
      v0->id() = *it;
      auto it3 = it;
      it3++;
      Point2 p1 = bestfitplane.to_2d(_lsPts[*it3]);
      CT::Vertex_handle v1 = ct.insert(p1);
      v1->id() = *it3;
      if (v0 != v1) {
        ct.insert_constraint(v0, v1);
      }
    }
  }
  mark_domains(ct); 
  if (!ct.is_valid()) 
    return re;
  for (CT::Finite_faces_iterator fit = ct.finite_faces_begin();
       fit != ct.finite_faces_end(); 
       ++fit) 
  {
    if (fit->info().in_domain()) {
      int* tr = new int[3];
      tr[0] = fit->vertex(0)->id();
      if (reversed) {
        tr[1] = fit->vertex(2)->id();
        tr[2] = fit->vertex(1)->id();
      } else {
        tr[1] = fit->vertex(1)->id();
        tr[2] = fit->vertex(2)->id();
      }
      re.push_back(tr);
    }
  }
  return re;
}


void Surface::get_min_bbox(double& x, double& y)
{
  std::vector<Point3>::iterator it = _lsPts.begin();
  K::FT minx = 9e10;
  K::FT miny = 9e10;
  for ( ; it != _lsPts.end(); it++)
  {
    if (it->x() < minx)
      minx = it->x();
    if (it->y() < miny)
      miny = it->y();
  }
  x = minx;
  y = miny;
}


void Surface::translate_vertices()
{
  std::vector<Point3>::iterator it = _lsPts.begin();
  for (it = _lsPts.begin(); it != _lsPts.end(); it++)
  {
    Point3 tp(CGAL::to_double(it->x() - Surface::_shiftx), CGAL::to_double(it->y() - Surface::_shifty), CGAL::to_double(it->z()));
    *it = tp;
  }
}


void Surface::set_translation_min_values(double minx, double miny)
{
  Surface::_shiftx = minx;
  Surface::_shifty = miny;
}


bool Surface::validate_2d_primitives(double tol_planarity_d2p, double tol_planarity_normals)
{
  *val3dityclog << "-----2D validation of each surface" << std::endl;
  bool isValid = true;
  size_t num = _lsFaces.size();
  for (int i = 0; i < static_cast<int>(num); i++)
  {
    //-- test for too few points (<3 for a ring)
    if (has_face_rings_toofewpoints(_lsFaces[i]) == true)
    {
      this->add_error(101, _lsFacesID[i]);
      isValid = false;
      continue;
    }
    //-- test for 2 repeated consecutive points
    if (has_face_2_consecutive_repeated_pts(_lsFaces[i]) == true)
    {
      this->add_error(102, _lsFacesID[i]);
      isValid = false;
      continue;
    }
    size_t numf = _lsFaces[i].size();
    std::vector<int> &ids = _lsFaces[i][0]; // helpful alias for the outer boundary

    //-- if only 3 pts it's not valid, no need to process further
    if ( (numf == 1) && (ids.size() == 3)) 
    {
      if (CGAL::collinear(_lsPts[ids[0]], _lsPts[ids[1]], _lsPts[ids[2]]) == true) {
        this->add_error(104, _lsFacesID[i], " outer ring (a triangle) is collapsed to a line");
        isValid = false;
      }
      continue;
    }

    std::vector< Point3 > allpts;
    std::vector<int>::const_iterator itp = ids.begin();
    for ( ; itp != ids.end(); itp++)
    {
      allpts.push_back(_lsPts[*itp]);
    }
    //-- irings
    for (int j = 1; j < static_cast<int>(numf); j++)
    {
      std::vector<int> &ids2 = _lsFaces[i][j]; // helpful alias for the inner boundary
      std::vector<int>::const_iterator itp2 = ids2.begin();
      for ( ; itp2 != ids2.end(); itp2++)
      {
        allpts.push_back(_lsPts[*itp2]);
      }
    }
    double value;
    CgalPolyhedron::Plane_3 bestfitplane = get_best_fitted_plane(allpts);
    if (false == is_face_planar_distance2plane(allpts, bestfitplane, value, tol_planarity_d2p))
    {
      std::stringstream msg;
      msg << "distance to fitted plane: " << value << " (tolerance=" << tol_planarity_d2p << ")";
      this->add_error(203, _lsFacesID[i], msg.str());
      isValid = false;
      continue;
    }
    //-- get projected oring
    Polygon pgn;
    std::vector<Polygon> lsRings;
    create_cgal_polygon(_lsPts, ids, bestfitplane, pgn);
    if (validate_projected_ring(pgn, _lsFacesID[i]) == false)
    {
      isValid = false;
      continue;
    }
    lsRings.push_back(pgn);
    //-- check for irings
    for (int j = 1; j < static_cast<int>(numf); j++)
    {
      std::vector<int> &ids2 = _lsFaces[i][j]; // helpful alias for the inner boundary
      //-- get projected iring
      Polygon pgn;
      create_cgal_polygon(_lsPts, ids2, bestfitplane, pgn);
      if (validate_projected_ring(pgn, _lsFacesID[i]) == false)
      {
        isValid = false;
        continue;
      }
      lsRings.push_back(pgn);
    }
    //-- use GEOS to validate projected polygon
    if (!validate_polygon(lsRings, _lsFacesID[i]))
      isValid = false;
  }
  if (isValid)
  {
    //-- triangulate faces of the shell
    triangulate_shell();
    //-- check planarity by normal deviation method (of all triangle)
    *val3dityclog << "-----Planarity of surfaces (with normals deviation)" << std::endl;
    std::vector< std::vector<int*> >::iterator it = _lsTr.begin();
    int j = 0;
    double deviation;
    for ( ; it != _lsTr.end(); it++)
    { 
      if (is_face_planar_normals(*it, _lsPts, deviation, tol_planarity_normals) == false)
      {
        std::ostringstream msg;
        msg << "deviation normals: " << deviation << " (tolerance=" << tol_planarity_normals << ")";
        this->add_error(204, _lsFacesID[j], msg.str());
        isValid = false;
      }
      j++;
    }
  }
  _is_valid_2d = isValid;
  return isValid;
}


bool Surface::validate_as_multisurface(double tol_planarity_d2p, double tol_planarity_normals)
{
  *val3dityclog << "--- MultiSurface validation ---" << std::endl;
  if (_is_valid_2d == -1)
    return validate_2d_primitives(tol_planarity_d2p, tol_planarity_normals);
  else
  {
    if (_is_valid_2d == 1)
      return true;
    else
      return false;
  }
}


bool Surface::is_shell(double tol_planarity_d2p, double tol_planarity_normals) 
{
  return validate_as_shell(tol_planarity_d2p, tol_planarity_normals);
}


bool Surface::validate_as_compositesurface(double tol_planarity_d2p, double tol_planarity_normals)
{
  *val3dityclog << "--- CompositeSurface validation ---" << std::endl;
//-- 1. Each surface should individually be valid
  if (_is_valid_2d == -1)
    validate_2d_primitives(tol_planarity_d2p, tol_planarity_normals);
  if (_is_valid_2d == 0)
    return false;
//-- 2. Combinatorial consistency
  *val3dityclog << "--Combinatorial consistency" << std::endl;
  _polyhedron = construct_CgalPolyhedron_incremental(&(_lsTr), &(_lsPts), this);
  if (this->has_errors() == true)
    return false;
  if (_polyhedron != NULL)
  {
    if (_polyhedron->is_valid() == true)
    {
      if (_polyhedron->keep_largest_connected_components(1) > 0)
      {
        this->add_error(305);
        return false;
      }
    }
    else 
    {
      this->add_error(300, "", "Something went wrong during construction of the shell, reason is unknown.");
      return false;
    }
    //-- check for non-manifold vertices
    //-- errors will be added to the Surface
    this->contains_nonmanifold_vertices();
  }
  else
  {
    this->add_error(300, "", "Something went wrong during construction of the shell, reason is unknown.");
    return false;
  }
  if (this->has_errors() == true)
    return false;
//-- 2. Geometrical consistency (aka intersection tests between faces)
  *val3dityclog << "--Geometrical consistency" << std::endl;
  if (does_self_intersect() == false)
    return false;
  return true;
}


bool Surface::contains_nonmanifold_vertices()
{
  //-- check for non-manifold vertices (which oddly CGAL doesn't expose as a function)
  //-- the trick is to check if a given vertex has >2 incident edges (not half-edges, but edges)
  //-- that are border edges
  // count border edges 
  std::map<CgalPolyhedron::Vertex_handle, int> noincident;
  for (CgalPolyhedron::Halfedge_handle heh = _polyhedron->halfedges_begin(); heh != _polyhedron->halfedges_end(); heh++)
  {
    if (heh->is_border() == true)
    {
      auto it = noincident.find(heh->vertex());
      if (it == noincident.end())
        noincident[heh->vertex()] = 1;
      else
        noincident[heh->vertex()] += 1;
      it = noincident.find(heh->opposite()->vertex());
      if (it == noincident.end())
        noincident[heh->opposite()->vertex()] = 1;
      else
        noincident[heh->opposite()->vertex()] += 1;
    }
  }
  for (auto& v : noincident)
  {
    if (v.second > 2)
    {
      std::stringstream st;
      st << "Non-manifold vertex at (";;
      st << (v.first->point().x() + _shiftx);
      st << ", ";
      st << (v.first->point().y() + _shifty);
      st << ", ";
      st << v.first->point().z();
      st << ")";
      this->add_error(303, "", st.str());
    }
  }
  return (this->has_errors());
}

bool Surface::does_self_intersect()
{
  if (CGAL::Polygon_mesh_processing::does_self_intersect(*_polyhedron) == true)
  {
    std::vector<std::pair<CgalPolyhedron::Facet_const_handle, CgalPolyhedron::Facet_const_handle> > intersected_tris;
    CGAL::Polygon_mesh_processing::self_intersections(*_polyhedron, std::back_inserter(intersected_tris));
    std::set<CgalPolyhedron::Facet_const_handle> uniquetr;
    for (auto& each : intersected_tris)
    {
      uniquetr.insert(each.first);
      uniquetr.insert(each.second);
    }
    for (auto& each : uniquetr)
    {
      //-- report the triangle centroid, as an approximation of the intersection
      //-- faster and less error-prone than reporting exact location
      Point3 c = CGAL::centroid(each->halfedge()->vertex()->point(),
                                each->halfedge()->next()->vertex()->point(),
                                each->halfedge()->next()->next()->vertex()->point()); 
      std::stringstream st;
      st << "Location close to: (";
      st << c.x() + _shiftx;
      st << ", ";
      st << c.y() + _shifty;
      st << ", ";
      st << c.z();
      st << ")"; 
      this->add_error(306, "", st.str());
    }
    return false;
  }
  return true;
}


bool Surface::validate_as_shell(double tol_planarity_d2p, double tol_planarity_normals)
{
  *val3dityclog << "--- Shell validation (#" << _id << ") ---" << std::endl;
  if (_is_valid_2d == -1)
    validate_2d_primitives(tol_planarity_d2p, tol_planarity_normals);
  if (_is_valid_2d == 0)
    return false;
//-- 1. minimum number of faces = 4
  if (_lsTr.size() < 4) 
  {
    this->add_error(301);
    return false;
  }
//-- 2. Combinatorial consistency
  *val3dityclog << "-----Combinatorial consistency" << std::endl;
  _polyhedron = construct_CgalPolyhedron_incremental(&(_lsTr), &(_lsPts), this);
  if (this->has_errors() == true)
    return false;
  if (_polyhedron != NULL)
  {
    if (_polyhedron->is_valid() == true)
    {
      if (_polyhedron->is_closed() == false)
      {
        _polyhedron->normalize_border();
        //-- check for unconnected faces
        if (_polyhedron->keep_largest_connected_components(1) > 0)
        {
          //TODO: how to report what face is not connected? a bitch of a problem...
          this->add_error(305);
          return false;
        }
        else
        {
          //-- check if there are holes in the surface
          if (_polyhedron->is_closed() == false)
          {
            std::stringstream st;
            _polyhedron->normalize_border();
            while (_polyhedron->size_of_border_edges() > 0) {
              CgalPolyhedron::Halfedge_handle he = ++(_polyhedron->border_halfedges_begin());
              st << "Location hole: (";
              st << (he->vertex()->point().x() + _shiftx);
              st << ", ";
              st << (he->vertex()->point().y() + _shifty);
              st << ", ";
              st << he->vertex()->point().z();
              st << ")";
              this->add_error(302, "", st.str());
              st.str("");
              _polyhedron->fill_hole(he);
              _polyhedron->normalize_border();
            }
            return false;
          }
        }
      }
      else //-- check if >1 connected components exist (both valid)
      {
        if (_polyhedron->keep_largest_connected_components(1) > 0) 
        {
          this->add_error(305, "", "More than one connected components.");
          return false;
        }
      }
    }
    else 
    {
      this->add_error(300, "", "Something went wrong during construction of the shell, reason is unknown.");
      return false;
    }
  }
  else
  {
    this->add_error(300, "", "Something went wrong during construction of the shell, reason is unknown.");
    return false;
  }
//-- 3. Geometrical consistency (aka intersection tests between faces)
  *val3dityclog << "-----Geometrical consistency" << std::endl;
  if (does_self_intersect() == false)
    return false;
  return true;
}


bool Surface::validate_projected_ring(Polygon &pgn, std::string id)
{
  if ( (!pgn.is_simple()) || (pgn.orientation() == CGAL::COLLINEAR) )
  {
    this->add_error(104, id, "ring self-intersects or is collapsed to a line");
    return false;
  }
  return true;
}


bool Surface::validate_polygon(std::vector<Polygon> &lsRings, std::string polygonid)
{
  initGEOS(NULL, NULL);
  //-- check the orientation of the rings: oring != irings
  //-- we don't care about CCW or CW at this point, just opposite is important
  //-- GEOS doesn't do its job, so we have to do it here. Shame on you GEOS.
  bool isvalid = true;
  if (lsRings.size() > 1)
  {
    CGAL::Orientation ooring = lsRings[0].orientation();
    std::vector<Polygon>::iterator it = lsRings.begin();
    it++;
    for ( ; it != lsRings.end(); it++)
    {
      if (it->orientation() == ooring)
      {
        this->add_error(208, polygonid, "same orientation for outer and inner rings");
        isvalid = false;
        break;
      }
    }
  }
  if (isvalid == false)
    return isvalid;
  //-- check 2D validity of the surface by (1) projecting them; (2) use GEOS IsValid()
  stringstream wkt;
  wkt << setprecision(15);
  wkt << "POLYGON(";
  std::vector<Polygon>::iterator it = lsRings.begin();
  //-- oring
  wkt << "(";
  Polygon::Vertex_iterator vit = it->vertices_begin();
  for ( ; vit != it->vertices_end(); vit++)
  {
    wkt << vit->x() << " " << vit->y() << ",";
  }
  vit = it->vertices_begin();
  wkt << vit->x() << " " << vit->y() << ")";
  it++;
  //-- irings
  for ( ; it != lsRings.end(); it++)
  {
    wkt << ", (";
    it->reverse_orientation();
    Polygon::Vertex_iterator vit = it->vertices_begin();
    for ( ; vit != it->vertices_end(); vit++)
    {
      wkt << vit->x() << " " << vit->y() << ",";
    }
    vit = it->vertices_begin();
    wkt << vit->x() << " " << vit->y() << ")";
    it->reverse_orientation();
  }
  wkt << ")";
  GEOSWKTReader* r;
  r = GEOSWKTReader_create();
  GEOSGeometry* mygeom;
  mygeom = GEOSWKTReader_read(r, wkt.str().c_str());
  string reason = (string)GEOSisValidReason(mygeom);
  if (reason.find("Valid Geometry") == string::npos)
  {
    isvalid = false;
    if (reason.find("Self-intersection") != string::npos)
      this->add_error(201, polygonid, reason.c_str());
    else if (reason.find("Duplicate Rings") != string::npos)
      this->add_error(202, polygonid, reason.c_str());
    else if (reason.find("Interior is disconnected") != string::npos)
      this->add_error(205, polygonid, reason.c_str());
    else if (reason.find("Hole lies outside shell") != string::npos)
      this->add_error(206, polygonid, reason.c_str());
    else if (reason.find("Holes are nested") != string::npos)
      this->add_error(207, polygonid, reason.c_str());
    else
      this->add_error(999, polygonid, reason.c_str());
  }
  GEOSGeom_destroy( mygeom );
  finishGEOS();
  return isvalid;
}

bool Surface::has_face_2_consecutive_repeated_pts(const std::vector< std::vector<int> >& theface)
{
  bool bDuplicates = false;
  std::vector< std::vector<int> >::const_iterator itr = theface.begin();
  for ( ; itr != theface.end(); itr++) {
    size_t numv = itr->size();
    //-- first-last not the same (they are not in GML format anymore)
    if ((*itr)[0] == (*itr)[numv - 1]) {
      bDuplicates = true;
      break;
    }
    for (int i = 0; i < (static_cast<int>(numv) - 1); i++) {
      if ((*itr)[i] == (*itr)[i+1]) {
        bDuplicates = true;
        break;
      }
    }
  }
  return bDuplicates;
}

bool Surface::has_face_rings_toofewpoints(const std::vector< std::vector<int> >& theface)
{
  bool bErrors = false;
  std::vector< std::vector<int> >::const_iterator itr = theface.begin();
  for ( ; itr != theface.end(); itr++) {
    if (itr->size() < 3) {
      bErrors = true;
      break;
    }
  }
  return bErrors;
}

int Surface::side_of_triangle_surface(Point3& p)
  /*
   -2 = not valid polyhedron
   -1 = outside
   0 = directly on the boundary of polyhedron
   1 = inside
   */
{
  int re = -2;
  if ( (_polyhedron != NULL) && (CGAL::is_triangle_mesh(*_polyhedron) == true) )
  {
    CGAL::Side_of_triangle_mesh<CgalPolyhedron, K> inside(*_polyhedron);
    Point3 p_translated(p.x() - Surface::_shiftx, p.y() - Surface::_shifty, p.z());
    re = inside(p_translated);
  }
  return re;
}


} // namespace val3dity
