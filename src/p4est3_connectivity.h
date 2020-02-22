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

#ifndef P4EST3_CONNECTIVITY_H
#define P4EST3_CONNECTIVITY_H

#include <p4est3_base.h>

#ifdef __cplusplus
extern              "C"
{
#if 0
}
#endif
#endif

typedef int         (*p4est3_connectivity_is_t) (void *slf, char *reason);
typedef sc3_error_t *(*p4est3_connectivity_inout_t) (void *slf);
typedef sc3_error_t *(*p4est3_connectivity_get_topidx_t)
                    (void *slf, p4est3_topidx * ptopidx);

typedef struct p4est3_connectivity_vtable
{
  p4est3_connectivity_get_topidx_t get_num_trees;
  p4est3_connectivity_inout_t destroy;
}
p4est3_connectivity_vtable_t;

typedef struct p4est3_connectivity p4est3_connectivity_t;

int                 p4est3_connectivity_is_valid (const p4est3_connectivity_t
                                                  * c, char *reason);
int                 p4est3_connectivity_is_new (const p4est3_connectivity_t *
                                                c, char *reason);
int                 p4est3_connectivity_is_setup (const p4est3_connectivity_t
                                                  * c, char *reason);

sc3_error_t        *p4est3_connectivity_new (sc3_allocator_t * alloc,
                                             p4est3_connectivity_t ** pc);

sc3_error_t        *p4est3_connectivity_set_vtable
  (p4est3_connectivity_t * c, p4est3_connectivity_vtable_t * cvt, void *slf);
sc3_error_t        *p4est3_connectivity_set_num_trees
  (p4est3_connectivity_t * c, p4est3_topidx num_trees);

sc3_error_t        *p4est3_connectivity_setup (p4est3_connectivity_t * c);
sc3_error_t        *p4est3_connectivity_ref (p4est3_connectivity_t * c);
sc3_error_t        *p4est3_connectivity_unref (p4est3_connectivity_t ** c);
sc3_error_t        *p4est3_connectivity_destroy (p4est3_connectivity_t ** c);

sc3_error_t        *p4est3_connectivity_get_num_trees
  (const p4est3_connectivity_t * c, p4est3_topidx * pnum_trees);

sc3_error_t        *p4est3_connectivity_new_num_trees
  (sc3_allocator_t * alloc, p4est3_topidx num_trees,
   p4est3_connectivity_t ** pc);
sc3_error_t        *p4est3_connectivity_new_unitcube
  (sc3_allocator_t * alloc, p4est3_connectivity_t ** pc);

#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif /* !P4EST3_CONNECTIVITY_H */
