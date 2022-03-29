/*
  Copyright (c) 2021-2022,
  Ivan Pađen <i.paden@tudelft.nl>
  3D Geoinformation,
  Delft University of Technology

  This file is part of City4CFD.

  City4CFD is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  City4CFD is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

  You should have received a copy of the GNU General Public License
  along with this program. If not, see <http://www.gnu.org/licenses/>
*/

#ifndef CITY4CFD_TOP_H
#define CITY4CFD_TOP_H

#include "Boundary.h"

class Top : public Boundary {
public:
    Top(const int outputLayerID);
    ~Top();

    virtual void reconstruct() override;

    virtual TopoClass   get_class() const override;
    virtual std::string get_class_name() const override;
};

#endif //CITY4CFD_TOP_H