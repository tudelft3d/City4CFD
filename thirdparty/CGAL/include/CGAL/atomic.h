// Copyright (c) 2016 GeometryFactory Sarl (France)
//  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL: https://github.com/CGAL/cgal/blob/v5.4/Installation/include/CGAL/atomic.h $
// $Id: atomic.h 6481cb2 2021-08-13T16:44:53+02:00 Sébastien Loriot
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial

#ifndef CGAL_ATOMIC_H
#define CGAL_ATOMIC_H

#define CGAL_DEPRECATED_HEADER "<CGAL/atomic.h>"
#define CGAL_REPLACEMENT_HEADER "<CGAL/config.h>"

#include <CGAL/config.h>

#ifdef CGAL_HAS_THREADS

#    include <atomic>

namespace CGAL {
namespace cpp11 {
using std::atomic;

using std::memory_order_relaxed;
using std::memory_order_consume;
using std::memory_order_acquire;
using std::memory_order_release;
using std::memory_order_acq_rel;
using std::memory_order_seq_cst;

using std::atomic_thread_fence;
} }
#else
#  define CGAL_NO_ATOMIC "No atomic because CGAL_NO_THREADS is defined."
#endif // CGAL_HAS_THREADS

#endif // CGAL_ATOMIC_H
