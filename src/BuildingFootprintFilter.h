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


#ifndef CITY4CFD_BUILDINGFOOTPRINTFILTER_H
#define CITY4CFD_BUILDINGFOOTPRINTFILTER_H

#include "types.h"
#include "CGALTypes.h"
#include "PolyFeature.h"

#include <memory>
#include <vector>

// Spatial index over building footprint polygons.
    // Used in two contexts:
    //   1. Read-time: routes building-class LAS points into per-footprint buckets and
    //      filters stray building points that fall outside every footprint.
    //   2. Post-influence-region: identifies terrain points that fall inside active
    //      building footprints so they can be removed before terrain reconstruction.
    class BuildingFootprintFilter {
    public:
        // One entry per successfully indexed footprint polygon.
        struct PolygonEntry {
            Polygon_2            polygon;        // buffered footprint (used for PIP)
            std::size_t          original_index; // position in the source collection passed to build()
            std::vector<Point_3> bucket;         // building points collected at read time
        };

        BuildingFootprintFilter();
        ~BuildingFootprintFilter();   // defined in .cpp where Impl is complete

        // Build from raw polygon data (outer boundary of each entry, dilated by buffer m).
        void build(const PolyVecPtr& polygons, double buffer = 0.);

        // Build from a vector of shared_ptr to any PolyFeature-derived type (Building,
        // ReconstructedBuilding, etc.).  buffer=0 skips the GEOS offset.
        template<typename T>
        void build(const std::vector<std::shared_ptr<T>>& features, double buffer = 0.) {
            static_assert(std::is_base_of<PolyFeature, T>::value,
                          "T must derive from PolyFeature");
            std::vector<Polygon_2> polys;
            polys.reserve(features.size());
            for (const auto& f : features)
                polys.push_back(f->get_poly().outer_boundary());
            build_from_polygons(polys, buffer);
        }

        // Returns true if (x, y) lies inside (or on the boundary of) any indexed footprint.
        bool contains(double x, double y) const;

        // Appends pt to the bucket of the first matching footprint.
        // Returns true if accepted, false if dropped as stray (outside every footprint).
        // Non-const — mutates internal buckets.
        bool collect_building_point(double x, double y, const Point_3& pt);

        // All entries (polygon + original_index + collected bucket) after a read.
        const std::vector<PolygonEntry>& get_entries() const;

        bool empty() const;

        // Release all indexed polygons and collected buckets.
        void clear();

    private:
        struct Impl;
        std::unique_ptr<Impl> m_impl;

        // Common implementation used by the template build() overload.
        void build_from_polygons(const std::vector<Polygon_2>& polys, double buffer);
    };

#endif //CITY4CFD_BUILDINGFOOTPRINTFILTER_H
