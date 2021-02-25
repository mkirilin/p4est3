/*
  This file is part of p4est, version 3.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2019 individual authors
  Originally written by Carsten Burstedde, Lucas C. Wilcox, and Tobin Isaac

  p4est is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  p4est is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with p4est; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/** \file p4est3_quadrant_zyx.h
 *
 * 3D quadrant implementation using a 128 bit AVX accelerated data type.
 *
 * \ingroup p4est3
 */

#ifndef P4EST_QUADRANT_ZYX_H
#define P4EST_QUADRANT_ZYX_H

#include <p4est3_quadrant_vtable.h>

#ifdef __cplusplus
extern              "C"
{
#if 0
}
#endif
#endif

#define P4EST3_ZYX_MAXLEVEL 30
#define P4EST3_ZYX_QMAXLEVEL 30

/** Populate a 3D quadrant virtual table with an AVX implementation.
 * We use the level and x, y, z coordinates inside a 4x32 bit hardware type.
 * \param [out] qvt     Members populated with virtual functions.
 * \return              NULL on success, error object otherwise.
 *                      If AVX hardware support is not available,
 *                      return an error of kind SC3_ERROR_RUNTIME.
 */
sc3_error_t        *p4est3_quadrant_zyx_vtable (p4est3_quadrant_vtable_t *
                                                qvt);

#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif /* !P4EST_QUADRANT_ZYX_H */
