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

/** \file p4est3_ghost_p4est.h
 * Fill p4est structure p4est_ghost_t by iteration over p4est3 structure.
 *
 * \ingroup p4est3
 */

#ifndef P4EST3_GHOST_P4EST_H
#define P4EST3_GHOST_P4EST_H

#include <p4est3.h>
#include <p4est_ghost.h>
#include <p4est3_iterate.h>

#ifdef __cplusplus
extern              "C"
{
#if 0
}
#endif
#endif

typedef struct p4est3_ghost_p4est
{
  p4est_ghost_t      *ghost;        /**< Pointer to the ghost structure */
  sc_hash_t          *gid_to_pos;   /**< Maps global IDs to a position in
                                         a ghost array */
}
p4est3_ghost_p4est_t;

/** Fill the \a p4est \a ghost layer structure by iterating over \a p4est3.
 * \param [in] p3       Forest passed for reference.
 * \param [out] ptr_ghost On output a pointer to a valid ghost structure.
 * \return              NULL on success, error object otherwise.
 */
sc3_error_t        *p4est3_ghost_fill_p4est (p4est3_t * p3,
                                             p4est3_ghost_p4est_t **
                                             ptr_ghost);

#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif /* !P4EST3_GHOST_P4EST_H */
