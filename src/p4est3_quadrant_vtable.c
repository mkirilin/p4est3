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

#include <p4est3_quadrant_vtable.h>

int
p4est3_quadrant_vtable_is_valid (p4est3_quadrant_vtable_t * qvt, char *reason)
{
  /* arguments */
  SC3E_TEST (qvt != NULL, reason);

  /*** test member variables ***/
  SC3E_TEST (0 < qvt->dim && qvt->dim <= 3, reason);
  SC3E_TEST (0 < qvt->max_level, reason);
  SC3E_TEST (0 < qvt->max_children, reason);

  /*** test member functions ***/
  SC3E_TEST (qvt->quadrant_tree_boundary != NULL, reason);
  SC3E_TEST (qvt->quadrant_num_uniform != NULL, reason);
  SC3E_TEST (qvt->quadrant_level != NULL, reason);
  SC3E_TEST (qvt->quadrant_child_id != NULL, reason);
  SC3E_TEST (qvt->quadrant_ancestor_id != NULL, reason);
  SC3E_TEST (qvt->quadrant_coordinates != NULL, reason);
  SC3E_TEST (qvt->quadrant_compare != NULL, reason);
  SC3E_TEST (qvt->quadrant_root != NULL, reason);
  SC3E_TEST (qvt->quadrant_copy != NULL, reason);
  SC3E_TEST (qvt->quadrant_parent != NULL || qvt->quadrant_ancestor != NULL,
             reason);
  SC3E_TEST (qvt->quadrant_face_neighbor != NULL, reason);
  SC3E_TEST (qvt->quadrant_face_neighbor_extra != NULL, reason);
  SC3E_TEST (qvt->quadrant_transform_face != NULL, reason);
  SC3E_TEST (qvt->quadrant_predecessor != NULL, reason);
  SC3E_TEST (qvt->quadrant_successor != NULL, reason);
  SC3E_TEST (qvt->quadrant_child != NULL, reason);
  SC3E_TEST (qvt->quadrant_first_descendant != NULL, reason);
  SC3E_TEST (qvt->quadrant_last_descendant != NULL, reason);
  SC3E_TEST (qvt->quadrant_morton != NULL, reason);

  /* this is it */
  SC3E_YES (reason);
}

int
p4est3_quadrant_id (p4est3_quadrant_vtable_t * qvt)
{
  if (qvt == NULL || qvt->id < 0) {
    return -1;
  }
  return qvt->id;
}

int
p4est3_quadrant_dim (p4est3_quadrant_vtable_t * qvt)
{
  if (qvt == NULL || qvt->dim <= 0 || qvt->dim > 3) {
    return -1;
  }
  return qvt->dim;
}

int
p4est3_quadrant_max_level (p4est3_quadrant_vtable_t * qvt)
{
  if (qvt == NULL || qvt->max_level < 0) {
    return -1;
  }
  return qvt->max_level;
}

int
p4est3_quadrant_max_children (p4est3_quadrant_vtable_t * qvt)
{
  if (qvt == NULL || qvt->max_children <= 0) {
    return -1;
  }
  return qvt->max_children;
}

size_t
p4est3_quadrant_size (p4est3_quadrant_vtable_t * qvt)
{
  if (qvt == NULL) {
    return 0;
  }
  return qvt->quadrant_size;
}

p4est3_gloidx
p4est3_quadrant_num_uniform (p4est3_quadrant_vtable_t * qvt, int level)
{
  if (qvt == NULL || qvt->quadrant_num_uniform == NULL) {
    return -1;
  }
  return qvt->quadrant_num_uniform (level);
}

int
p4est3_quadrant_is2_valid (p4est3_quadrant_vtable_t * qvt,
                           const void *q, char *reason)
{
  SC3E_TEST (qvt != NULL, reason);
  if (qvt->quadrant_is_valid != NULL) {
    SC3E_IS (qvt->quadrant_is_valid, q, reason);
  }
  SC3E_YES (reason);
}

int
p4est3_quadrant_is2_inside_root (p4est3_quadrant_vtable_t * qvt,
                                const void *q, char *reason)
{
  SC3E_TEST (qvt != NULL, reason);
  if (qvt->quadrant_is_inside_root != NULL) {
    SC3E_IS (qvt->quadrant_is_inside_root, q, reason);
  }
  SC3E_YES (reason);
}

int
p4est3_quadrant_is3_equal (p4est3_quadrant_vtable_t * qvt,
                           const void *q1, const void *q2, char *reason)
{
  SC3E_TEST (qvt != NULL, reason);
  if (qvt->quadrant_is_equal != NULL) {
    SC3E_IS2 (qvt->quadrant_is_equal, q1, q2, reason);
  }
  else {
    SC3E_TEST (!memcmp (q1, q2, qvt->quadrant_size), reason);
  }
  SC3E_YES (reason);
}

sc3_error_t        *
p4est3_quadrant_tree_boundary (p4est3_quadrant_vtable_t * qvt,
                              const void *q, sc3_array_t * nf)
{
  SC3A_CHECK (qvt != NULL);
  SC3E (qvt->quadrant_tree_boundary (q, nf));
  return NULL;
}

sc3_error_t        *
p4est3_quadrant_level (p4est3_quadrant_vtable_t * qvt, const void *q, int *l)
{
  SC3A_CHECK (qvt != NULL);
  if (qvt->quadrant_level != NULL) {
    SC3E (qvt->quadrant_level (q, l));
  }
  else {
    *l = qvt->max_level;
  }
  return NULL;
}

sc3_error_t        *
p4est3_quadrant_child_id (p4est3_quadrant_vtable_t * qvt,
                          const void *q, int *j)
{
  SC3A_CHECK (qvt != NULL && qvt->quadrant_child_id != NULL);
  SC3E (qvt->quadrant_child_id (q, j));
  return NULL;
}

sc3_error_t        *
p4est3_quadrant_ancestor_id (p4est3_quadrant_vtable_t * qvt,
                             const void *q, int l, int *j)
{
  SC3A_CHECK (qvt != NULL && qvt->quadrant_ancestor_id != NULL);
  SC3E (qvt->quadrant_ancestor_id (q, l, j));
  return NULL;
}

sc3_error_t        *
p4est3_quadrant_coordinates (p4est3_quadrant_vtable_t * qvt,
                             const void *q, int n, void *j)
{
  SC3A_CHECK (qvt != NULL && qvt->quadrant_coordinates != NULL);
  SC3E (qvt->quadrant_coordinates (q, n, j));
  return NULL;
}

sc3_error_t        *
p4est3_quadrant_num_children (p4est3_quadrant_vtable_t * qvt,
                              const void *q, int *n)
{
  SC3A_CHECK (qvt != NULL);
  if (qvt->quadrant_num_children != NULL) {
    SC3E (qvt->quadrant_num_children (q, n));
  }
  else {
    *n = qvt->max_children;
  }
  return NULL;
}

sc3_error_t        *
p4est3_quadrant_compare (p4est3_quadrant_vtable_t * qvt,
                         const void *q1, const void *q2, int *j)
{
  SC3A_CHECK (qvt != NULL && qvt->quadrant_compare != NULL);
  SC3E (qvt->quadrant_compare (q1, q2, j));
  return NULL;
}

sc3_error_t        *
p4est3_quadrant_root (p4est3_quadrant_vtable_t * qvt, void *r)
{
  SC3A_CHECK (qvt != NULL && qvt->quadrant_root != NULL);
  SC3E (qvt->quadrant_root (r));
  return NULL;
}

sc3_error_t        *
p4est3_quadrant_copy (p4est3_quadrant_vtable_t * qvt, const void *q, void *r)
{
  SC3A_CHECK (qvt != NULL);

  if (qvt->quadrant_copy != NULL) {
    SC3E (qvt->quadrant_copy (q, r));
  }
  else {
    memcpy (r, q, qvt->quadrant_size);
  }
  return NULL;
}

sc3_error_t        *
p4est3_quadrant_parent (p4est3_quadrant_vtable_t * qvt,
                        const void *q, void *r)
{
  SC3A_CHECK (qvt != NULL);

  if (qvt->quadrant_parent != NULL) {
    SC3E (qvt->quadrant_parent (q, r));
  }
  else {
    int                 level;

    SC3A_CHECK (qvt->quadrant_level != NULL);
    SC3E (qvt->quadrant_level (q, &level));
    SC3A_CHECK (level > 0);

    SC3A_CHECK (qvt->quadrant_ancestor != NULL);
    SC3E (qvt->quadrant_ancestor (q, level - 1, r));
  }
  return NULL;
}

sc3_error_t        *
p4est3_quadrant_face_neighbor (p4est3_quadrant_vtable_t * qvt,
                               const void *q, int i, void *r)
{
  SC3A_CHECK (qvt != NULL && qvt->quadrant_face_neighbor != NULL);
  SC3E (qvt->quadrant_face_neighbor (q, i, r));
  return NULL;
}

sc3_error_t        *
p4est3_quadrant_transform_face (p4est3_quadrant_vtable_t * qvt,
                                const void *q, sc3_array_t * transform,
                                void *r)
{
  SC3A_CHECK (qvt != NULL && qvt->quadrant_transform_face != NULL);
  SC3E (qvt->quadrant_transform_face (q, transform, r));
  return NULL;
}

sc3_error_t        *
p4est3_quadrant_face_neighbor_extra (p4est3_quadrant_vtable_t * qvt,
                                     void *c, const void *q, int i,
                                     int *j, int *k, void *r)
{
  SC3A_CHECK (qvt != NULL && qvt->quadrant_face_neighbor_extra != NULL);
  SC3E (qvt->quadrant_face_neighbor_extra (c, q, i, j, k, r));
  return NULL;
}

sc3_error_t        *
p4est3_quadrant_predecessor (p4est3_quadrant_vtable_t * qvt,
                             const void *q, void *r)
{
  SC3A_CHECK (qvt != NULL && qvt->quadrant_predecessor != NULL);
  SC3E (qvt->quadrant_predecessor (q, r));
  return NULL;
}

sc3_error_t        *
p4est3_quadrant_successor (p4est3_quadrant_vtable_t * qvt,
                           const void *q, void *r)
{
  SC3A_CHECK (qvt != NULL && qvt->quadrant_successor != NULL);
  SC3E (qvt->quadrant_successor (q, r));
  return NULL;
}

sc3_error_t        *
p4est3_quadrant_child (p4est3_quadrant_vtable_t * qvt,
                       const void *q, int i, void *r)
{
  SC3A_CHECK (qvt != NULL && qvt->quadrant_child != NULL);
  SC3E (qvt->quadrant_child (q, i, r));
  return NULL;
}

sc3_error_t        *
p4est3_quadrant_ancestor (p4est3_quadrant_vtable_t * qvt,
                          const void *q, int l, void *r)
{
  SC3A_CHECK (qvt != NULL);

  if (qvt->quadrant_ancestor != NULL) {
    SC3E (qvt->quadrant_ancestor (q, l, r));
  }
  else {
    int                 level;

    SC3A_CHECK (l >= 0);
    SC3E (qvt->quadrant_level (q, &level));
    SC3A_CHECK (level >= l);

    SC3E (p4est3_quadrant_copy (qvt, q, r));

    SC3A_CHECK (qvt->quadrant_parent != NULL);
    while (level > l) {
      SC3E (qvt->quadrant_parent (r, r));
      SC3E (qvt->quadrant_level (r, &level));
    }
  }
  return NULL;
}

sc3_error_t        *
p4est3_quadrant_first_descendant (p4est3_quadrant_vtable_t * qvt,
                                  const void *q, int l, void *r)
{
  SC3A_CHECK (qvt != NULL && qvt->quadrant_first_descendant != NULL);
  SC3E (qvt->quadrant_first_descendant (q, l, r));
  return NULL;
}

sc3_error_t        *
p4est3_quadrant_last_descendant (p4est3_quadrant_vtable_t * qvt,
                                 const void *q, int l, void *r)
{
  SC3A_CHECK (qvt != NULL && qvt->quadrant_last_descendant != NULL);
  SC3E (qvt->quadrant_last_descendant (q, l, r));
  return NULL;
}

sc3_error_t        *
p4est3_quadrant_morton (p4est3_quadrant_vtable_t * qvt,
                        int level, p4est3_gloidx id, void *r)
{
  SC3A_CHECK (qvt != NULL && qvt->quadrant_morton != NULL);
  SC3E (qvt->quadrant_morton (level, id, r));
  return NULL;
}

sc3_error_t        *
p4est3_nearest_common_ancestor (p4est3_quadrant_vtable_t * qvt,
                                const void *q1, const void *q2, void *r)
{
  SC3A_CHECK (qvt != NULL && qvt->nearest_common_ancestor != NULL);
  SC3E (qvt->nearest_common_ancestor (q1, q2, r));
  return NULL;
}

sc3_error_t        *
p4est3_quadrant_linear_id (p4est3_quadrant_vtable_t * qvt,
                           const void *q, int l, p4est3_gloidx * id)
{
  SC3A_CHECK (qvt != NULL && qvt->quadrant_linear_id != NULL);
  SC3E (qvt->quadrant_linear_id (q, l, id));
  return NULL;
}

sc3_error_t        *
p4est3_quadrant_is_ancestor (p4est3_quadrant_vtable_t * qvt,
                             const void *q1, const void *q2, int *j)
{
  SC3A_CHECK (qvt != NULL && qvt->quadrant_is_ancestor != NULL);
  SC3E (qvt->quadrant_is_ancestor (q1, q2, j));
  return NULL;
}

sc3_error_t        *
p4est3_quadrant_array_new (sc3_allocator_t * alloc,
                           p4est3_quadrant_vtable_t * qvt,
                           p4est3_locidx n, sc3_array_t ** arr)
{
  SC3E_RETVAL (arr, NULL);
  SC3A_IS (sc3_allocator_is_setup, alloc);
  SC3A_CHECK (qvt != NULL);
  SC3A_CHECK (n >= 0);

  SC3E (sc3_array_new (alloc, arr));
  SC3E (sc3_array_set_elem_size (*arr, qvt->quadrant_size));
  SC3E (sc3_array_set_elem_count (*arr, n));
  SC3E (sc3_array_setup (*arr));

  return NULL;
}
