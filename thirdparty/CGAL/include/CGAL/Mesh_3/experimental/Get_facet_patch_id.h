// Copyright (c) 2010,2012  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL: https://github.com/CGAL/cgal/blob/v5.4/Mesh_3/include/CGAL/Mesh_3/experimental/Get_facet_patch_id.h $
// $Id: Get_facet_patch_id.h bf128b3 2021-01-25T16:19:23+01:00 Sébastien Loriot
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent Rineau
//

#ifndef CGAL_MESH_3_GET_FACET_PATCH_ID_H
#define CGAL_MESH_3_GET_FACET_PATCH_ID_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/Mesh_3/experimental/Facet_patch_id_map.h>

namespace CGAL { namespace Mesh_3 {

// backward compatibility with user code
template <typename MeshDomain>
using Get_facet_patch_id_sm = Facet_patch_id_map<MeshDomain, typename MeshDomain::AABB_tree>;

}} // end namespace CGAL::Mesh_3

#endif // CGAL_MESH_3_GET_FACET_PATCH_ID_H
