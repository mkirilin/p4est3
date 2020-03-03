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

#ifndef P4EST3_QUADRANT_VTABLE
#define P4EST3_QUADRANT_VTABLE

#include <p4est3_base.h>

#ifdef __cplusplus
extern              "C"
{
#if 0
}
#endif
#endif

/* *INDENT-OFF* */
typedef int         (*p4est3_quadrant_int_t) (void);
typedef size_t      (*p4est3_quadrant_size_t) (void);
typedef int         (*p4est3_quadrant_is_t) (const void * q, char *reason);
typedef sc3_error_t *(*p4est3_quadrant_out_t) (void *r);
typedef sc3_error_t *(*p4est3_quadrant_in_j_t) (const void *q, int *j);
typedef sc3_error_t *(*p4est3_quadrant_in_i_j_t) (const void *q, int i, int *j);
typedef sc3_error_t *(*p4est3_quadrant_in_out_t) (const void *q, void *r);
typedef sc3_error_t *(*p4est3_quadrant_in_i_out_t) (const void *q, int i,
                                                    void *r);
typedef sc3_error_t *(*p4est3_quadrant_morton_t) (int l, p4est3_gloidx i,
                                                  void *r);

typedef p4est3_quadrant_in_j_t p4est3_quadrant_level_t;
typedef p4est3_quadrant_in_j_t p4est3_quadrant_child_id_t;
typedef p4est3_quadrant_in_i_j_t p4est3_quadrant_ancestor_id_t;
typedef p4est3_quadrant_out_t p4est3_quadrant_root_t;
typedef p4est3_quadrant_in_out_t p4est3_quadrant_parent_t;
typedef p4est3_quadrant_in_out_t p4est3_quadrant_predecessor_t;
typedef p4est3_quadrant_in_out_t p4est3_quadrant_successor_t;
typedef p4est3_quadrant_in_i_out_t p4est3_quadrant_ancestor_t;
typedef p4est3_quadrant_in_i_out_t p4est3_quadrant_child_t;
typedef p4est3_quadrant_in_i_out_t p4est3_quadrant_first_descendant_t;
typedef p4est3_quadrant_in_i_out_t p4est3_quadrant_last_descendant_t;
/* *INDENT-ON* */

typedef struct p4est3_quadrant_vtable
{
  int                 id;
  int                 dim;
  p4est3_quadrant_int_t max_level;
  p4est3_quadrant_int_t num_children;
  p4est3_quadrant_size_t quadrant_size;
  p4est3_quadrant_is_t quadrant_is_valid;
  p4est3_quadrant_level_t quadrant_level;
  p4est3_quadrant_child_id_t quadrant_child_id;
  p4est3_quadrant_ancestor_id_t quadrant_ancestor_id;
  p4est3_quadrant_root_t quadrant_root;
  p4est3_quadrant_parent_t quadrant_parent;
  p4est3_quadrant_predecessor_t quadrant_predecessor;
  p4est3_quadrant_successor_t quadrant_successor;
  p4est3_quadrant_child_t quadrant_child;
  p4est3_quadrant_ancestor_t quadrant_ancestor;
  p4est3_quadrant_first_descendant_t quadrant_first_descendant;
  p4est3_quadrant_last_descendant_t quadrant_last_descendant;
  p4est3_quadrant_morton_t quadrant_morton;
}
p4est3_quadrant_vtable_t;

int                 p4est3_max_level (p4est3_quadrant_vtable_t * qvt);
int                 p4est3_num_children (p4est3_quadrant_vtable_t * qvt);
size_t              p4est3_quadrant_size (p4est3_quadrant_vtable_t * qvt);
int                 p4est3_quadrant_is_valid (p4est3_quadrant_vtable_t * qvt,
                                              const void *q, char *reason);
sc3_error_t        *p4est3_quadrant_level (p4est3_quadrant_vtable_t * qvt,
                                           const void *q, int *l);
sc3_error_t        *p4est3_quadrant_child_id (p4est3_quadrant_vtable_t * qvt,
                                              const void *q, int *j);
sc3_error_t        *p4est3_quadrant_ancestor_id (p4est3_quadrant_vtable_t *
                                                 qvt, const void *q, int l,
                                                 int *j);
sc3_error_t        *p4est3_quadrant_root (p4est3_quadrant_vtable_t * qvt,
                                          void *r);
sc3_error_t        *p4est3_quadrant_parent (p4est3_quadrant_vtable_t * qvt,
                                            const void *q, void *r);
sc3_error_t        *p4est3_quadrant_predecessor (p4est3_quadrant_vtable_t * qvt,
                                                 const void *q, void *r);
sc3_error_t        *p4est3_quadrant_successor (p4est3_quadrant_vtable_t * qvt,
                                               const void *q, void *r);
sc3_error_t        *p4est3_quadrant_child (p4est3_quadrant_vtable_t * qvt,
                                           const void *q, int i, void *r);
sc3_error_t        *p4est3_quadrant_ancestor (p4est3_quadrant_vtable_t * qvt,
                                              const void *q, int l, void *r);
sc3_error_t        *p4est3_quadrant_first_descendant (p4est3_quadrant_vtable_t
                                                      * qvt, const void *q,
                                                      int l, void *r);
sc3_error_t        *p4est3_quadrant_last_descendant (p4est3_quadrant_vtable_t
                                                     * qvt, const void *q,
                                                     int l, void *r);
sc3_error_t        *p4est3_quadrant_morton (p4est3_quadrant_vtable_t * qvt,
                                            int level, p4est3_gloidx id,
                                            void *r);

#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif /* !P4EST3_QUADRANT_VTABLE */
