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

/** \file p4est3_iterate.h
 * Iterate through the leaves of a forest as well as the neighbor connections.
 *
 * \ingroup p4est3
 */

#ifndef P4EST3_ITERATE_H
#define P4EST3_ITERATE_H

#include <p4est3.h>

#ifdef __cplusplus
extern              "C"
{
#if 0
}
#endif
#endif

typedef struct p4est3_iterate_face_side
{
  p4est3_topidx       ntree;
  int                 nface;
  void               *quadrant;
}
p4est3_iterate_face_side_t;

typedef struct p4est3_iterate_face_info
{
  p4est3_t           *p3;
  void               *user_data;
  int                 orientation; /**< the orientation of the sides to each
                                        other, as in the definition of
                                        p4est_connectivity_t */
  int                 tree_boundary; /**< boolean */
  p4est3_iterate_face_side_t sides[2];
}
p4est3_iterate_face_info_t;

typedef sc3_error_t *(*p4est3_iterate_face_t) (p4est3_iterate_face_info_t *
                                               fi);

/** Iterate through the forest.
 * \return              NULL on success, error object otherwise.
 */
sc3_error_t        *p4est3_iterate (p4est3_t * p3,
                                    p4est3_iterate_face_t * face_callback,
                                    void *user_data);

#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif /* !P4EST3_ITERATE_H */
