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

#include <p4est3_internal.h>
#include <sc3_array.h>

#ifdef __cplusplus
extern              "C"
{
#if 0
}
#endif
#endif

static sc3_error_t *
p4est3_search_array_new (sc3_allocator_t * alloc, size_t esize,
                         int ealloc, int ecount, sc3_array_t ** arr)
{
  SC3E_RETVAL (arr, NULL);
  SC3A_IS (sc3_allocator_is_setup, alloc);
  SC3A_CHECK (ealloc >= 0);

  SC3E (sc3_array_new (alloc, arr));
  SC3E (sc3_array_set_elem_size (*arr, esize));
  SC3E (sc3_array_set_elem_alloc (*arr, ealloc));
  SC3E (sc3_array_set_elem_count (*arr, ecount));
  SC3E (sc3_array_set_initzero (*arr, 1));
  SC3E (sc3_array_setup (*arr));

  return NULL;
}

/** A callback function that describes the search window.
 *  The idea is to define the type of an array entry as type 1, if
 *  my_begin <= array[i], my_end > array[i] and as type 2, if the entry
 *  is not of type 1 and my_end <= array[i]. The remaining cases are
 *  indicated by type 0. This function is passed into sc_array_split by
 *  \ref p4est_find_partition. Note that this function, as well as
 *  p4est_find_partition, is dimension-independent; still we duplicate
 *  it in 3D in the usual way.
 */
static sc3_error_t *
type_fn_global_quad_index (sc3_array_t * array, int index,
                           void *data_array, int *type)
{
  p4est3_gloidx      *my_begin_end, *entry;

  SC3A_CHECK (data_array != NULL);
  my_begin_end = (p4est3_gloidx *) data_array;
  SC3E (sc3_array_index (array, index, &entry));

  if (*entry < my_begin_end[0]) {
    *type = 0;
  }
  else if ((my_begin_end[0] <= *entry) && (my_begin_end[1] > *entry)) {
    *type = 1;
  }
  else if (my_begin_end[1] <= *entry) {
    *type = 2;
  }
  else {
    *type = 0;
  }
  return NULL;
}

sc3_error_t        *
p4est3_find_partition (sc3_allocator_t * alloc, int num_entities,
                       p4est3_gloidx * search_in,
                       p4est3_gloidx my_begin, p4est3_gloidx my_end,
                       p4est3_gloidx * begin, p4est3_gloidx * end)
{
  int                *iptr;
  sc3_array_t        *view, *offsets;
  p4est3_gloidx       my_begin_end[2];

  SC3A_CHECK (my_begin <= my_end);
  SC3E (p4est3_search_array_new (alloc, sizeof (int), 0, 0, &offsets));
  SC3E (sc3_array_new_data (alloc, &view, search_in,
                            sizeof (p4est3_gloidx), 0, num_entities));

  my_begin_end[0] = my_begin;
  my_begin_end[1] = my_end;

  SC3E (sc3_array_split
        (view, offsets, 3, type_fn_global_quad_index, my_begin_end));

  SC3E (sc3_array_index (offsets, 1, &iptr));
  *begin = *iptr;
  SC3E (sc3_array_index (offsets, 2, &iptr));
  *end = *iptr;

  sc3_array_destroy (&offsets);
  sc3_array_destroy (&view);

  return NULL;
}

sc3_error_t        *
p4est3_search_lower_bound64 (int64_t target, const int64_t * array,
                             ssize_t nmemb, ssize_t * guess)
{
  ssize_t             k_low, k_high;
  int64_t             cur;

  if (nmemb == 0) {
    *guess = -1;
    return NULL;
  }

  k_low = 0;
  k_high = nmemb - 1;
  for (;;) {
    SC3A_CHECK (k_low <= k_high);
    SC3A_CHECK (k_low < nmemb && k_high < nmemb);
    SC3A_CHECK (k_low <= *guess && *guess <= k_high);

    /* compare two quadrants */
    cur = array[*guess];

    /* check if guess is higher or equal target and there's room below it */
    if (target <= cur && (*guess > 0 && target <= array[*guess - 1])) {
      k_high = *guess - 1;
      *guess = (k_low + k_high + 1) / 2;
      continue;
    }

    /* check if guess is lower than target */
    if (target > cur) {
      k_low = *guess + 1;
      if (k_low > k_high) {
        *guess = -1;
        return NULL;
      }
      *guess = (k_low + k_high) / 2;
      continue;
    }

    /* otherwise guess is the correct position */
    break;
  }

  SC3A_CHECK (*guess < nmemb);
  SC3A_CHECK (array[*guess] >= target);
  SC3A_CHECK (*guess == 0 || array[*guess - 1] < target);
  return NULL;
}

/* Wrap cont qvt since sc3's sc3_array_is_sorted takes non-const user data. */
typedef struct p4est3_qvt_non_const_wrapper
{
  const p4est3_quadrant_vtable_t *qvt;
} p4est3_qvt_non_const_wrapper_t;

#ifdef P4EST_ENABLE_DEBUG
/* Custom quadrant_compare function to align with sc3_array_is_sorted interface. */
static sc3_error_t *
p4est3_array_split_compare (const void * q1, const void * q2,
                            void * wqvt, int * j)
{
  SC3E (p4est3_quadrant_compare
          (((p4est3_qvt_non_const_wrapper_t *) wqvt)->qvt, q1, q2, j));
  return NULL;
}
#endif

/* Wrap user data to pass into custom p4est3_array_split_ancestor_id function. */
typedef struct p4est3_array_split_data
{
  p4est3_quadrant_ancestor_id_t quadrant_ancestor_id;
  int                *level;
}
p4est3_array_split_data_t;

/* Custom quadrant ancestor_id function to align with sc3_array_split interface. */
static sc3_error_t *
p4est3_array_split_ancestor_id (sc3_array_t * a, int index, void *data,
                                int *type)
{
  SC3A_CHECK (data != NULL);

  void               *q;
  p4est3_array_split_data_t *d = (p4est3_array_split_data_t *) data;
  SC3E (sc3_array_index (a, index, &q));

  SC3E (d->quadrant_ancestor_id (q, *(d->level), type));
  return NULL;
}

sc3_error_t         *
p4est3_quadrant_array_split (const p4est3_quadrant_vtable_t * qvt,
                             sc3_array_t * array, int level,
                             sc3_array_t * indices)
{
  p4est3_array_split_data_t data;
  p4est3_qvt_non_const_wrapper_t sqvtw, *qvtw = &sqvtw;
#ifdef P4EST_ENABLE_DEBUG
  void               *q1, *q2;
  int                 l, count;
#endif

  qvtw->qvt = qvt;
  SC3A_IS (sc3_array_is_setup, array);
  SC3A_IS (sc3_array_is_setup, indices);
  SC3A_CHECK (qvt != NULL);
  SC3A_CHECK (0 <= level && level < qvt->max_level);
  SC3A_IS3 (sc3_array_is_sorted, array, p4est3_array_split_compare, qvtw);

#ifdef P4EST_ENABLE_DEBUG
  SC3E (sc3_array_get_elem_count (array, &count));
  if (count > 0) {
    SC3E (sc3_array_index (array, 0, &q1));
    SC3E (p4est3_quadrant_level (qvt, q1, &l));
    SC3A_CHECK (l > level);
    SC3E (sc3_array_index (array, count - 1, &q2));
    SC3E (p4est3_quadrant_level (qvt, q2, &l));
    SC3A_CHECK (l > level);
  }
  /*TODO: check if l >= level, where l is a level of nearest
     common ancestor of q1 and q2.
   */
#endif

  level++;
  data.quadrant_ancestor_id = qvt->quadrant_ancestor_id;
  data.level = &level;
  SC3E (sc3_array_split (array, indices, p4est3_quadrant_num_children (qvt),
                         p4est3_array_split_ancestor_id, &data));
  return NULL;
}

#ifdef __cplusplus
#if 0
{
#endif
}
#endif
