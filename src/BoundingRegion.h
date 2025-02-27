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

#ifndef CITY4CFD_BOUNDINGREGION_H
#define CITY4CFD_BOUNDINGREGION_H

#include "types.h"
#include "CGALTypes.h"
#include "Config.h"

class BoundingRegion {
public:
    CDT m_projCDT;
    std::shared_ptr<const Config::ReconRegion> m_reconSettings;

    BoundingRegion() = default;
    BoundingRegion(std::shared_ptr<Config::ReconRegion> reconRegion);
    ~BoundingRegion() = default;

    void operator()(double radius);
    void operator()(Polygon_2& poly);

    double calc_influ_region_bpg(const DT& dt, BuildingsPtr& buildings);
    void   calc_influ_region_bpg(const double maxDim);
    void calc_bnd_bpg(const Polygon_2& influRegionPoly, const BuildingsPtr& buildings);
    bool is_subset_of(const BoundingRegion& otherRegion) const;

    Polygon_2& get_bounding_region();
    const Polygon_2& get_bounding_region() const;

protected:
    Polygon_2  m_boundingRegion;

    Polygon_2  calc_bnd_poly(const std::vector<Point_2>& candidatePts, const double hMax,
                             const double angle, const double enlargeRatio = 1);

//    double     calc_blockage_ratio_from_chull (const BuildingsPtr& buildings, const double angle,
//                                               Polygon_2& localPoly) const;
//    double     calc_blockage_ratio_from_ashape(const BuildingsPtr& buildings, const double angle,
//                                               Polygon_2& localPoly) const;
    double     calc_blockage_ratio_comb(const BuildingsPtr& buildings, const double angle,
                                        Polygon_2& localPoly) const;
//    double     calc_blockage_ratio_from_ashape_alt(const BuildingsPtr& buildings, const double angle,
//                                                   Polygon_2& localPoly) const;
//    double     calc_blockage_ratio_from_edges(const BuildingsPtr& buildings, const double angle,
//                                              Polygon_2& localPoly) const;

    void       project_mesh_pts(const Mesh& mesh, const double angle, std::vector<Point_2>& buildingPts) const;
    void       chull_to_cdt(const std::vector<Point_2>& buildingPts, CDT& projCDT) const;
    void       ashape_to_cdt(const std::vector<Point_2>& buildingPts, CDT& projCDT, const double aVal) const;
    void       calc_cross_sec_areas(CDT& projCDT, Polygon_2& localPoly, const BuildingsPtr& buildings,
                                    double& blockArea, double& domainCrossArea) const;
};

#endif //CITY4CFD_BOUNDINGREGION_H