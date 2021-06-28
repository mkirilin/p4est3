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

typedef struct p4est3_array_split_data
{
  p4est3_quadrant_ancestor_id_t quadrant_ancestor_id;
  int                *level;
}
p4est3_array_split_data_t;

typedef struct p4est3_search_area
{
  /* general section */
  p4est3_locidx      *begin;
  p4est3_locidx      *end;
  int                 start_level;
  p4est3_tree_t      *tree;
  sc3_array_t        *view_quads;     /* array ptr to pass tree's quads
                                   into array_split */
  sc3_array_t        *view_indices;
  sc3_array_t        *split_offsets;    /* temp array of size max_children + 1
                                           for split_array output */

  /* volume section */
  int                 Level;
  int                *level2nchildren;  /* Array specifing the number n
                                           of children processed on
                                           the particular level; n can't be
                                           greater than p3->qvt->max_children */
  sc3_array_t        *idx_vol_stack;
  p4est3_iterate_volume_info_t *vinfo;

  /* face section */
  int                 dir;
  int                *Level_sides;      /* array of 2, storing the level of
                                           current size */
  int                *is_refine;        /* array of 2, storing the information
                                           about neccesity of refenement. Must
                                           be initiated by 0 */
  p4est3_locidx     **begin_face;
  p4est3_locidx     **end_face;
  sc3_array_t        *inner_quadrants;  /* array of max_children, storing
                                           result split_array in face_init */
  sc3_array_t       **idx_face_stack;
  int                *dir_order;        /*Array of array of array (side -> dir -> order) */

  int                 orient;
  int                *faces;
  p4est3_tree_t     **trees_face;
  p4est3_iterate_face_info_t *finfo;
}
p4est3_search_area_t;

static sc3_error_t *
p4est3_array_is_sorted (const void * q1, const void * q2, void * qvt, int *j)
{
  p4est3_quadrant_vtable_t *qvtable;
  SC3A_IS (p4est3_quadrant_vtable_is_valid, qvt);

  qvtable = (p4est3_quadrant_vtable_t *) qvt;
  SC3A_IS (qvtable->quadrant_is_valid, q1);
  SC3A_IS (qvtable->quadrant_is_valid, q2);
  SC3A_CHECK (j != NULL);

  SC3E (qvtable->quadrant_compare (q1, q2, j));
  return NULL;
}

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

static sc3_error_t *
p4est3_quadrant_array_split (p4est3_quadrant_vtable_t * qvt,
                             sc3_array_t * array, int level,
                             sc3_array_t * indices)
{
  p4est3_array_split_data_t data;
#ifdef P4EST_ENABLE_DEBUG
  void               *q1, *q2;
  int                 l, count;
#endif

  SC3A_IS (sc3_array_is_setup, array);
  SC3A_IS (sc3_array_is_setup, indices);
  SC3A_CHECK (qvt != NULL);
  SC3A_CHECK (0 <= level && level < qvt->max_level);
  SC3A_IS3 (sc3_array_is_sorted, array, p4est3_array_is_sorted, qvt);

#ifdef P4EST_ENABLE_DEBUG
  SC3E (sc3_array_get_elem_count (array, &count));
  SC3E (sc3_array_index (array, 0, &q1));
  SC3E (p4est3_quadrant_level (qvt, q1, &l));
  SC3A_CHECK (l > level);
  SC3E (sc3_array_index (array, count - 1, &q2));
  SC3E (p4est3_quadrant_level (qvt, q2, &l));
  SC3A_CHECK (l > level);
  /*TODO: check if l >= level, where l is a level of nearest
     common ancestor of q1 and q2.
   */
#endif

  level++;
  data.quadrant_ancestor_id = qvt->quadrant_ancestor_id;
  data.level = &level;
  SC3E (sc3_array_split (array, indices, qvt->max_children,
                         p4est3_array_split_ancestor_id, &data));
  return NULL;
}

static inline int *
p4est3_direction_order (const int side, const int dir,
                        const int dim, int *dir_ord)
{                               /*TODO: optimize */
  return dir_ord + (1 << (dim - 1)) * (side * dim + dir);
}

static inline sc3_error_t *
p4est3_reverse_copy_vol (const int n, sc3_array_t * src, sc3_array_t * dst,
                         p4est3_locidx ** ptr_stack_it)
{
  int                 i;
  p4est3_locidx      *offsets_it, *stack_it;

  SC3E_RETVAL (ptr_stack_it, NULL);
  SC3E (sc3_array_index (src, n, &offsets_it));
  for (i = 0; i <= n; ++i, --offsets_it) {
    SC3E (sc3_array_push (dst, &stack_it));
    *(++stack_it) = *offsets_it;
  }
  *ptr_stack_it = stack_it;
  return NULL;
}

static inline sc3_error_t *
p4est3_reverse_copy_face (const int max_children, const int *order,
                          sc3_array_t * src, sc3_array_t * dst,
                          p4est3_locidx ** ptr_stack_it)
{
  int                 i;
  p4est3_locidx      *offsets_it, *stack_it;

  SC3E_RETVAL (ptr_stack_it, NULL);
  for (i = max_children / 2; i > 0; --i) {
    SC3E (sc3_array_index (src, order[i], &offsets_it));
    SC3E (sc3_array_push (dst, &stack_it));
    *(++stack_it) = *(offsets_it);
    SC3E (sc3_array_push (dst, &stack_it));
    *(++stack_it) = *(++offsets_it);
  }
  *ptr_stack_it = stack_it;
  return NULL;
}

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
p4est3_internal_iterate_face_rec (p4est3_t * p3, p4est3_iterate_face_t cface,
                                  p4est3_iterate_codim_t ccodim,
                                  p4est3_search_area_t * search_area)
{
  int                 i, side, level;
  int                *is_refine = search_area->is_refine;
  int                *dir_ord_arr = search_area->dir_order, *dir_ord_it;
  int                *Level = search_area->Level_sides;
  int                 max_children = p3->qvt->max_children;
  void               *first_quad;
  p4est3_locidx     **begin_face = search_area->begin_face,
    **end_face = search_area->end_face;
  p4est3_locidx      *stack_it;
  sc3_array_t        *sides = search_area->finfo->sides;
  sc3_array_t        *view = search_area->view;
  sc3_array_t        *split_offsets = search_area->split_offsets;
  sc3_array_t       **idx_face_stack = search_area->idx_face_stack;
  p4est3_tree_t      *tree = search_area->tree;
  p4est3_iterate_face_side_t *fside;

  for (side = 0; side < 2; ++side) {
    first_quad =
      (void *) (tree->tquads + p3->qvt->quadrant_size * (*begin_face[side]));
    SC3E (p4est3_quadrant_level (p3->qvt, first_quad, &level));
    if (level == Level[side]) {
      is_refine[side] = 0;
      SC3E (sc3_array_index (sides, side, &fside));
      fside->nquad = *begin_face[side];
      fside->quadrant = first_quad;
    }
    if (!is_refine[0] && !is_refine[1]) {
      if (cface != NULL) {
        SC3E (cface (search_area->finfo));
      }
      return NULL;
    }
    if (is_refine[side]) {
      SC3E (sc3_array_new_data (p3->alloc, &view, tree->tquads,
                                p3->qvt->quadrant_size, *(begin_face[side]),
                                *(end_face[side]) - *(begin_face[side])));
      SC3E (p4est3_quadrant_array_split
            (p3->qvt, view, Level[side], split_offsets));
      Level[side]++;
      dir_ord_it = p4est3_direction_order (side, search_area->dir,
                                           p3->qvt->dim, dir_ord_arr);
      SC3E (p4est3_reverse_copy_face (max_children, dir_ord_it,
                                      split_offsets, idx_face_stack[side],
                                      &stack_it));
      SC3A_CHECK (*stack_it == *(begin_face[side]) + 1);
      for (i = 0; i < max_children / 2; ++i) {
        begin_face[side] = stack_it;
        end_face[side] = --stack_it;
        SC3E (p4est3_internal_iterate_face_rec (p3, cface, ccodim,
                                                search_area));
        SC3E (sc3_array_pop (idx_face_stack[side]));
        SC3E (sc3_array_pop (idx_face_stack[side]));
      }
      Level[side]--;
    }
  }
  return NULL;
}

/* rearrange some data in the search ares so as to process faces */
static sc3_error_t *
p4est3_internal_iterate_face (p4est3_t * p3, p4est3_iterate_face_t cface,
                              p4est3_iterate_codim_t ccodim,
                              p4est3_search_area_t * search_area)
{
  const int           max_children = p3->qvt->max_children;
  int                 i, side, dir;
  int                 fn_shift, /* face neighbour shift */
                      dir_shift;        /* shift along directory */
  int                 idx = 0;
  int                *Level = &search_area->Level;
  int                *Level_sides = search_area->Level_sides;
  int                *is_refine = search_area->is_refine;

  p4est3_locidx      *begin_face, *end_face;
  p4est3_locidx      *idx_it;

  sc3_array_t        *sides = search_area->finfo->sides;
  sc3_array_t        *view = search_area->view;
  sc3_array_t        *inner_quadrants = search_area->inner_quadrants;
  p4est3_tree_t      *tree = search_area->tree;
  p4est3_iterate_face_side_t *fside;

  Level_sides[0] = Level_sides[1] = *Level;
  is_refine[0] = is_refine[1] = 1;

  end_face = search_area->end - 2;
  begin_face = search_area->begin - 2;

  if (*(begin_face) == *(end_face) - 1) {
    return NULL;
  }

  /* some preliminary definitions */
  search_area->finfo->orientation = 0 /*??? */ ;
  search_area->finfo->tree_boundary = 0;
  for (side = 0; side < 2; ++side) {
    SC3E (sc3_array_index (sides, side, &fside));
    fside->is_ghost = 0;
    fside->ntree = search_area->vinfo->ntree;
    /*ToDo: add #faces at the connection tracker */
    //fside->nface = 2; /*???*/ 
  }

  SC3E (sc3_array_new_data (p3->alloc, &search_area->view, tree->tquads,
                            p3->qvt->quadrant_size,
                            *begin_face, *end_face - *begin_face));
  SC3E (p4est3_quadrant_array_split (p3->qvt, view, *Level, inner_quadrants));

  /* try to procces quads of current level */
  for (dir = 0; dir < p3->qvt->dim; ++dir) {
    search_area->dir = dir;
    fn_shift = 1 << dir;
    dir_shift = fn_shift << 1;
    for (i = 0; i < max_children / 2; ++i) {
      SC3E (sc3_array_index (inner_quadrants, idx, &idx_it));
      idx = (idx + dir_shift) % (max_children - 1);
      for (side = 0; side < 2; ++side) {
        search_area->begin_face[side] = &(idx_it[side * fn_shift]);
        search_area->end_face[side] = &(idx_it[side * fn_shift + 1]);
        Level_sides[side]++;
      }
      SC3E (p4est3_internal_iterate_face_rec
            (p3, cface, ccodim, search_area));
    }
  }
  return NULL;
}

static sc3_error_t *
p4est3_iterate_face_frame (p4est3_t * p3,
                           p4est3_iterate_face_t cface,
                           p4est3_iterate_codim_t ccodim,
                           p4est3_search_area_t * search_area)
{
  const int max_children = p3->qvt->max_children;
  int *face = search_area->faces;
  p4est3_tree_t **trees = search_area->trees_face;
  p4est3_locidx **b_f = search_area->begin_face,
                **e_f = search_area->end_face;
  p4est3_locidx *stack_iter;
  sc3_array_t *sides = search_area->finfo->sides;
  sc3_array_t *view_q = search_area->view_quads;
  sc3_array_t *view_i = search_area->view_indices;
  sc3_array_t **idx_f_stack = search_area->idx_face_stack;
  p4est3_iterate_face_side_t *fside;
  int *Level_sides = search_area->Level_sides;
  int *is_refine = search_area->is_refine;
  int *Level = search_area->Level_sides;
  int i, side, level, acount, idx;
  int ori = search_area->orient;
  void *first_quad;

  /* first check if the whole quadrant passes */
  for (side = 0; side < 2; ++side) {
    if (!is_refine[side]) {
      continue;
    }
    first_quad =
      (void *) (trees[side]->tquads + p3->qvt->quadrant_size * (*b_f[side]));
    SC3E (p4est3_quadrant_level (p3->qvt, first_quad, &level));
    if (level == Level[side]) {
      is_refine[side] = 0;
      SC3E (sc3_array_index (sides, side, &fside));
      fside->nquad = *b_f[side];
      fside->quadrant = first_quad;
      fside->ntree = trees[side]->treeid;
      fside->nface = face[side];
    }
  }
  if (!is_refine[0] && !is_refine[1]) {
    if (cface != NULL) {
      SC3E (cface (search_area->finfo));
    }
    is_refine[0] = is_refine[1] = 1;
    return NULL;
  }
  for (i = 0; i < max_children / 2; ++i) {
    for (side = 0; side < 2; ++side) {
      if (!is_refine[side]) {
        continue;
      }
      SC3E (sc3_array_get_elem_count (idx_f_stack[side], &acount));
      SC3E (sc3_array_push_count (idx_f_stack[side], max_children + 1,
                                  &stack_iter));
      SC3E (sc3_array_new_data (p3->alloc, &view_q, trees[side]->tquads,
                                p3->qvt->quadrant_size, *(b_f[side]),
                                *(e_f[side]) - *(b_f[side])));
      SC3E (sc3_array_new_view (p3->alloc, &view_i, idx_f_stack[side], acount,
                                max_children + 1));
      SC3E (p4est3_quadrant_array_split (p3->qvt, view_q, Level[side],
                                         view_i));
      Level[side]++;
      idx = i;
      if (side == 1) {
        SC3E (p4est3_connectivity_face_neighbor_face_corner (p3->conn, &idx,
                                                             face[0],
                                                             face[1], ori));
      }
      SC3E (p4est3_connectivity_get_face_child_id (p3->conn, face[side], &idx));
      b_f[side] = stack_iter + idx;
      e_f[side] = stack_iter + idx + 1;
      SC3E (p4est3_iterate_face_frame (p3, cface, ccodim, search_area));
      Level[side]--;
    }
  }
  for (side = 0; side < 2; ++side) {
    if (!is_refine[side]) {
      continue;
    }
    for (i = 0; i < max_children + 1; ++i) {
      SC3E (sc3_array_pop (idx_f_stack[side]));
    }
  }

  return NULL;
}

static sc3_error_t *
p4est3_iterate_face_frame_init (p4est3_t * p3,
                                p4est3_iterate_face_t cface,
                                p4est3_iterate_codim_t ccodim,
                                p4est3_search_area_t * search_area)
{
  p4est3_locidx *idx;
  int *Level_sides = search_area->Level_sides;
  int *is_refine = search_area->is_refine;

  sc3_array_t ** idx_f_stack = search_area->idx_face_stack;
  p4est3_locidx **b_f = search_area->begin_face,
                **e_f = search_area->end_face;
  int side;

  Level_sides[0] = Level_sides[1] = 0;
  is_refine[0] = is_refine[1] = 1;

  for (side = 0; side < 2; ++side) {
    SC3E (sc3_array_index (idx_f_stack[side], 0, &b_f[side]));
    *(b_f[side]) = 0;
    SC3E (sc3_array_index (idx_f_stack[side], 1, &e_f[side]));
    *(e_f[side]) = search_area->trees_face[side]->num_quads;
  }
  return NULL;
}

static sc3_error_t *
p4est3_iterate_face_init (p4est3_t * p3,
                          p4est3_iterate_face_t cface,
                          p4est3_iterate_codim_t ccodim,
                          p4est3_search_area_t * search_area)
{
  p4est3_locidx *idx;
  int *Level_sides = search_area->Level_sides;
  int *is_refine = search_area->is_refine;

  sc3_array_t ** idx_f_stack = search_area->idx_face_stack;
  p4est3_locidx **b_f = search_area->begin_face,
                **e_f = search_area->end_face;
  p4est3_locidx *b_v = search_area->begin,
                *e_v = search_area->end;
  int side;

  Level_sides[0] = Level_sides[1] = 0;
  is_refine[0] = is_refine[1] = 1;

  for (side = 0; side < 2; ++side) {
    SC3E (sc3_array_index (idx_f_stack[side], 0, &b_f));
    *b_f[side] = *b_v;
    SC3E (sc3_array_index (idx_f_stack[side], 1, &e_f));
    *e_f[side] = *e_v;
  }

  return NULL;
}

static sc3_error_t *
p4est3_internal_iterate_volume (p4est3_t * p3,
                                p4est3_iterate_volume_t cvolume,
                                p4est3_iterate_face_t cface,
                                p4est3_iterate_codim_t ccodim,
                                p4est3_search_area_t * search_area)
{
  int                 i;
  int                 is_refine = 1;
  void               *first_quad;       /*first quadrant in this search area */
  int                 level;
  p4est3_locidx      *stack_it;

  const int           max_children = p3->qvt->max_children;
  int                *l2nch = search_area->level2nchildren;
  int                *Level = &search_area->Level;
  int                 icount;
  const int           begin = *(search_area->begin);
  const int           end = *(search_area->end);
  p4est3_tree_t      *tree = search_area->tree;
  sc3_array_t        *view_q = search_area->view_quads;
  sc3_array_t        *view_i = search_area->view_indices;
  sc3_array_t        *split_offsets = search_area->split_offsets;
  sc3_array_t        *idx_vol_stack = search_area->idx_vol_stack;
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
    SC3E (sc3_array_get_elem_count (idx_vol_stack, &icount));
    SC3E (sc3_array_push_count (idx_vol_stack, max_children + 1, &stack_it));
    SC3E (sc3_array_new_data (p3->alloc, &view_q, tree->tquads,
                              p3->qvt->quadrant_size, begin, end - begin));
    SC3E (sc3_array_new_view (p3->alloc, &view_i, idx_vol_stack, icount,
                              max_children + 1));
    SC3E (p4est3_quadrant_array_split (p3->qvt, view_q, *Level, view_i));
    l2nch[++(*Level)] = 0;
  }

  if (l2nch[search_area->start_level] > 0) {
    return NULL;
  }
  if (l2nch[*Level] == max_children) {
    SC3E (p4est3_iterate_face_frame_init (p3, cface, ccodim, search_area));
    SC3E (p4est3_iterate_face_frame (p3, cface, ccodim, search_area));
    l2nch[--(*Level)]++;
  }
  else {
    SC3E (p4est3_reverse_copy_vol (max_children, split_offsets,
                                   idx_vol_stack, &stack_it));
    for (i = 0; i < max_children; ++i) {
      search_area->begin = stack_it;
      search_area->end = ++stack_it;
      if (*(search_area->begin) == *(search_area->end) - 1) {
        l2nch[*Level]++;
        continue;
      }
      SC3E (p4est3_internal_iterate_volume (p3, cvolume, cface, ccodim,
                                            search_area));
    }
    for (i = 0; i < max_children + 1; ++i) {
      SC3E (sc3_array_pop (idx_vol_stack));
    }
  }
  return NULL;
}

sc3_error_t        *
p4est3_iterate_codim (p4est3_t * p3, int codims,
                      p4est3_iterate_volume_t cvolume,
                      p4est3_iterate_face_t cface,
                      p4est3_iterate_codim_t ccodim, void *user_data)
{
  p4est3_search_area_t ssa, *search_area = &ssa;
  p4est3_topidx       tree, tree_neighbor;
  int                 nfaces = p3->conn->num_faces;
  int                 i, face, orient;

  /* This iteration is w/o ghost layer and for volumes only */
  if (codims < 0 || codims >= P4EST3_ITERATE_LAST) {
    return NULL;
  }
  SC3A_IS (p4est3_is_setup, p3);
  if (p3->fltree < 0 || (cvolume == NULL && cface == NULL && ccodim == NULL)) {
    return NULL;
  }

  if (codims == P4EST3_ITERATE_VOLUME) {
    SC3E (p4est3_internal_iterate_volume_simple (p3, cvolume, user_data));
    return NULL;
  }

  /** ???Should not we loop over all trees and not just local trees because of the
   * ghost layer??? */
  for (tree = p3->fltree; tree <= p3->lltree; ++tree) {
    SC3E (p4est3_tree_index (p3, tree, &search_area->trees_face[0]));
    /*alloc memory that is necessary */
    /*set possible fields in vinfo */

    /* frame faces part */
    for (i = 0; i < nfaces; ++i) {
      face = i;
      tree_neighbor = tree;
      search_area->faces[0] = face;
      SC3E (p4est3_connectivity_get_face (p3->conn, &tree_neighbor, &face,
                                          &orient));
      SC3E (p4est3_tree_index (p3, tree_neighbor,
                               &search_area->trees_face[1]));
      search_area->faces[1] = face;
      search_area->orient = orient;
      if (tree_neighbor < tree) {
        continue;
      }
      else if (tree_neighbor == tree) {
        /* physical boundary */
      }
      else {
        SC3E (p4est3_iterate_face_frame_init (p3, cface, ccodim, search_area));
        SC3E (p4est3_iterate_face_frame (p3, cface, ccodim, search_area));
      }
    }
  }
  return NULL;
}

#ifdef __cplusplus
#if 0
{
#endif
}
#endif
