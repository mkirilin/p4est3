/*
  This file is part of p4est, version 3
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

/** \file p4est3_convert_p8est.h
 * Convert \ref p8est_t data type into \ref p4est3_t
 *
 * \ingroup p4est3
 */

#ifndef P4EST3_CONVERT_P8EST_H
#define P4EST3_CONVERT_P8EST_H

#include <p8est.h>
#include <p4est3.h>

#ifdef __cplusplus
extern              "C"
{
#if 0
}
#endif
#endif

/** Convert all the available data from \ref p8est_t object into preallocated
 *  \ref p4est3_t one. The \ref p4est3_t is return setup. Thus all preliminary
 * p4est3's manipulation parameters should be set. If they conflict with the
 * donor p4est_t object, then they will be overwritten.
 * \param [in] p            Ponter to \ref p8est_t object that must be setup.
 * \param [in,out] p3       Pointer to new allocated \ref p4est3_t object.
 *
 * \return                  NULL on success, error object otherwise.
 */
sc3_error_t        *p4est3_convert_p8est (p8est_t * p, p4est3_t * p3);

#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif /* !P4EST3_CONVERT_P8EST_H */
