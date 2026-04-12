/*
  City4CFD

  Copyright (c) 2021-2026, 3D Geoinformation Research Group, TU Delft

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
*/


#include "BuildingFootprintFilter.h"

#include "Config.h"
#include "geomutils.h"

#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/geometry/strategies/index/cartesian.hpp>

#include <iostream>

namespace {
    namespace bg  = boost::geometry;
    namespace bgi = boost::geometry::index;
    typedef bg::model::point<double, 2, bg::cs::cartesian> BgPt;
    typedef bg::model::box<BgPt>             BgBox;
    typedef std::pair<BgBox, std::size_t>    BgVal;
    typedef bgi::rtree<BgVal, bgi::quadratic<16>> BgRtree;
}

struct BuildingFootprintFilter::Impl {
    std::vector<PolygonEntry> entries; // one per successfully indexed footprint
    BgRtree                   rtree;   // R-tree value is index into entries
};

BuildingFootprintFilter::BuildingFootprintFilter()
    : m_impl(std::make_unique<Impl>()) {}

BuildingFootprintFilter::~BuildingFootprintFilter() = default;

void BuildingFootprintFilter::build(const PolyVecPtr& polygons, double buffer) {
    m_impl->entries.reserve(polygons.size());
    std::size_t skipped = 0;
    for (std::size_t srcIdx = 0; srcIdx < polygons.size(); ++srcIdx) {
        const Polygon_2& outer = polygons[srcIdx]->polygon.outer_boundary();
        if (outer.size() < 3) {
            Config::write_to_log("BuildingFootprintFilter: skipped polygon at index "
                                 + std::to_string(srcIdx) + " (fewer than 3 vertices).");
            ++skipped;
            continue;
        }
        try {
            Polygon_2 buffered = geomutils::offset_polygon_geos(outer, buffer);
            if (buffered.is_empty()) {
                Config::write_to_log("BuildingFootprintFilter: skipped polygon at index "
                                     + std::to_string(srcIdx) + " (buffer produced empty polygon).");
                ++skipped;
                continue;
            }
            auto bb = buffered.bbox();
            BgBox box(BgPt(bb.xmin(), bb.ymin()), BgPt(bb.xmax(), bb.ymax()));
            m_impl->rtree.insert({box, m_impl->entries.size()});
            m_impl->entries.push_back({std::move(buffered), srcIdx, {}});
        } catch (const std::exception& e) {
            Config::write_to_log("BuildingFootprintFilter: skipped polygon at index "
                                 + std::to_string(srcIdx) + " (buffering failed: "
                                 + e.what() + ").");
            ++skipped;
        }
    }
    std::cout << "    Footprint filter: " << m_impl->entries.size()
              << " footprints indexed (buffer: " << buffer << " m)";
    if (skipped > 0) {
        std::cout << "  [WARNING: " << skipped << " degenerate footprint(s) skipped — see log]";
        Config::get().logSummary << "BuildingFootprintFilter: " << skipped
                                 << " footprint(s) skipped during indexing (degenerate geometry)."
                                 << std::endl;
    }
    std::cout << std::endl;
}

void BuildingFootprintFilter::build_from_polygons(const std::vector<Polygon_2>& polys, double buffer) {
    m_impl->entries.reserve(polys.size());
    std::size_t skipped = 0;
    for (std::size_t i = 0; i < polys.size(); ++i) {
        if (polys[i].size() < 3) { ++skipped; continue; }
        try {
            Polygon_2 buffered = (buffer > 0.)
                ? geomutils::offset_polygon_geos(polys[i], buffer)
                : polys[i];
            if (buffered.is_empty()) { ++skipped; continue; }
            auto bb = buffered.bbox();
            BgBox box(BgPt(bb.xmin(), bb.ymin()), BgPt(bb.xmax(), bb.ymax()));
            m_impl->rtree.insert({box, m_impl->entries.size()});
            m_impl->entries.push_back({std::move(buffered), i, {}});
        } catch (const std::exception&) {
            ++skipped;
        }
    }
    if (skipped > 0)
        std::cout << "  [WARNING: " << skipped << " degenerate footprint(s) skipped]" << std::endl;
}

bool BuildingFootprintFilter::contains(double x, double y) const {
    BgPt pt(x, y);
    std::vector<BgVal> candidates;
    m_impl->rtree.query(bgi::intersects(BgBox(pt, pt)), std::back_inserter(candidates));
    const Point_2 cgalPt(x, y);
    for (const auto& [box, idx] : candidates) {
        if (geomutils::point_in_poly_and_boundary(cgalPt, m_impl->entries[idx].polygon))
            return true;
    }
    return false;
}

bool BuildingFootprintFilter::collect_building_point(double x, double y, const Point_3& pt) {
    BgPt bgPt(x, y);
    std::vector<BgVal> candidates;
    m_impl->rtree.query(bgi::intersects(BgBox(bgPt, bgPt)), std::back_inserter(candidates));
    const Point_2 cgalPt(x, y);
    for (const auto& [box, idx] : candidates) {
        if (geomutils::point_in_poly_and_boundary(cgalPt, m_impl->entries[idx].polygon)) {
            m_impl->entries[idx].bucket.push_back(pt);
            return true;
        }
    }
    return false;
}

const std::vector<BuildingFootprintFilter::PolygonEntry>&
BuildingFootprintFilter::get_entries() const {
    return m_impl->entries;
}

bool BuildingFootprintFilter::empty() const {
    return m_impl->entries.empty();
}
