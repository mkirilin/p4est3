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

#include <p4est3_iterate.h>
#include <p4est3_internal.h>

#ifdef __cplusplus
extern              "C"
{
#if 0
}
#endif
#endif

typedef struct p4est3_search_area
{
  p4est3_locidx      *begin;
  p4est3_locidx      *end;
  int                *level2nchildren;  /* Array specifing the number n
                                           of children processed on
                                           the particular level; n can't be
                                           greater than p3->qvt->max_children */
  int                 Level;
  int                 start_level;
  p4est3_tree_t      *tree;
  sc3_array_t        *idx_stack;
  sc3_array_t        *view;             /* array ptr to pass tree's quads
                                           into array_split */
  sc3_array_t        *split_offsets;    /* temp array of size max_children + 1
                                           for split_array output */
  p4est3_iterate_volume_info_t *vinfo;
}
p4est3_search_area_t;

sc3_error_t        *
p4est3_iterate_face (p4est3_t * p3,
                     p4est3_iterate_volume_t cvolume,
                     p4est3_iterate_face_t cface, void *user_data)
{
  SC3E (p4est3_iterate_codim (p3, 0x01, cvolume, cface, NULL, user_data));
  return NULL;
}

static sc3_error_t *
p4est3_internal_iterate_volume_simple (p4est3_t * p3,
                                       p4est3_iterate_volume_t cvolume,
                                       void *user_data)
{
  p4est3_topidx       ntree;
  p4est3_tree_t      *tree;
  p4est3_iterate_volume_info_t info;
  p4est3_locidx       si;

  info.p3 = p3;
  info.user_data = user_data;
  for (ntree = p3->fltree; ntree <= p3->lltree; ++ntree) {
    info.ntree = ntree;
    SC3E (p4est3_tree_index (p3, ntree, &tree));
    for (si = 0; si < tree->num_quads; ++si) {
      info.quadrant = tree->tquads + si * p3->qvt->quadrant_size;
      info.nquad = si + tree->quad_offset;
      SC3E (cvolume (&info));
    }
  }
  return NULL;
}

static sc3_error_t *
p4est3_internal_iterate_volume (p4est3_t * p3,
                                p4est3_iterate_volume_t cvolume,
                                p4est3_iterate_face_t cface,
                                p4est3_iterate_codim_t ccodim,
                                void *user_data,
                                p4est3_search_area_t * search_area)
{
  int                 i;
  int                 is_refine = 1;
  void               *first_quad;     /*first quadrant in this search area */
  int                 level;
  p4est3_locidx      *offsets_it, *stack_it;

  const int           max_children = p3->qvt->max_children;
  int                *l2nch = search_area->level2nchildren;
  int                *Level = &search_area->Level;
  const int           begin = *(search_area->begin);
  const int           end = *(search_area->end);
  p4est3_tree_t      *tree = search_area->tree;
  sc3_array_t        *view = search_area->view;
  sc3_array_t        *split_offsets = search_area->split_offsets;
  sc3_array_t        *idx_stack = search_area->idx_stack;
  p4est3_iterate_volume_info_t *vinfo = search_area->vinfo;

/*Are trees with no quads possible?*/
  first_quad = (void *) (tree->tquads + p3->qvt->quadrant_size * begin);
  SC3E (p4est3_quadrant_level (p3->qvt, first_quad, &level));
  if (level == *Level) {
    is_refine = 0;
    if (cvolume != NULL) {
      vinfo->quadrant = first_quad;
      vinfo->nquad = begin;
      SC3E (cvolume (vinfo));
    }
    l2nch[*Level]++;
    /*no return here because we want to iter_face a bit later */
  }

  if (is_refine) {
    SC3E (sc3_array_new_data (tree->tquads, p3->qvt->quadrant_size, begin,
                              end - begin, &view));
    SC3E (p4est3_quadrant_array_split (p3->qvt, view, *Level, split_offsets));
    l2nch[++(*Level)] = 0;
  }

  if (l2nch[search_area->start_level] > 0) {
    return NULL;
  }
  if (l2nch[*Level] == p3->qvt->max_children) {
    /*iterate faces */
    l2nch[--(*Level)]++;
  }
  else {
    SC3E (sc3_array_index (split_offsets, max_children, &offsets_it));
    for (i = 0; i <= max_children; ++i, --offsets_it) {
      SC3E (sc3_array_push (idx_stack, &stack_it));
      *(++stack_it) = *offsets_it;
    }
    for (i = 0; i < max_children; ++i) {
      search_area->begin = stack_it;
      search_area->end = --stack_it;
      SC3E (sc3_array_pop (idx_stack));
      if (*(search_area->begin) == *(search_area->end) - 1) {
        l2nch[*Level]++;
        continue;
      }
      SC3E (p4est3_internal_iterate_volume (p3, cvolume, cface, ccodim,
                                            user_data, search_area));
    }
    SC3E (sc3_array_pop (idx_stack));
  }
  return NULL;
}

sc3_error_t        *
p4est3_iterate_codim (p4est3_t * p3, int codims,
                      p4est3_iterate_volume_t cvolume,
                      p4est3_iterate_face_t cface,
                      p4est3_iterate_codim_t ccodim, void *user_data)
{
  p4est3_topidx       tree;

  /* This iteration is w/o ghost layer and for volumes only */
  if (codims < 0 || codims >= P4EST3_ITERATE_LAST) {
    return NULL;
  }
  SC3A_IS (p4est3_is_setup, p3);
  if (p3->fltree < 0 ||
     (cvolume == NULL && cface == NULL && ccodim == NULL)) {
    return NULL;
  }

  if (codims == P4EST3_ITERATE_VOLUME) {
    SC3E (p4est3_internal_iterate_volume_simple (p3, cvolume, user_data));
    return NULL;
  }

  /** ???Should not we loop over all trees and not just local trees because of the
   * ghost layer??? */
  for (tree = p3->fltree; tree <= p3->lltree; ++tree) {
    /*alloc memory that is necessary*/
    /*set possible fields in vinfo*/
  }
  return NULL;
}

#ifdef __cplusplus
#if 0
{
#endif
}
#endif
