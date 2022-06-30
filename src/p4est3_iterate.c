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
  int                 max_children;
  int                 nfaces;   /*Global number of faces */
  int                 start_level;
  sc3_array_t        *view_quads;       /* array ptr to pass tree's quads
                                           into array_split */
  int                *children_face_neighbors;
  int                *face_dual;

  /* volume section */
  p4est3_tree_t      *tree;
  p4est3_locidx      *begin;    /*or p4est3_gloidx?? */
  p4est3_locidx      *end;
  int                 Level;
  int                *level2nchildren;  /* Array specifing the number n
                                           of children processed on
                                           the particular level; n can't be
                                           greater than p3->num_children */
  sc3_array_t        *idx_vol_stack;    /* 2D stack storing arrays of indices,
                                           that are output of split_array */
  p4est3_iterate_volume_info_t *vinfo;

  /* face section */
  int                 nsides;
  p4est3_tree_t      *tree_face[2];
  p4est3_locidx      *begin_face[2];
  p4est3_locidx      *end_face[2];
  int                 Level_face[2];    /* array of 2, storing the level of
                                           current size */
  int                 is_refine[2];     /* array of 2, storing the information
                                           about neccesity of refenement. Must
                                           be initiated by 0 */
  sc3_array_t        *idx_face_stack[2];        /* 2D stacks storing arrays of indices,
                                                   that are output of split_array */
  p4est3_iterate_face_info_t *finfo;
}
p4est3_search_area_t;

static sc3_error_t *
p4est3_array_new (sc3_allocator_t * alloc, size_t esize, int ealloc,
                  int ecount, int is_resizable, sc3_array_t ** arr)
{
  SC3E_RETVAL (arr, NULL);
  SC3A_IS (sc3_allocator_is_setup, alloc);
  SC3A_CHECK (ealloc >= 0);

  SC3E (sc3_array_new (alloc, arr));
  SC3E (sc3_array_set_elem_size (*arr, esize));
  SC3E (sc3_array_set_elem_alloc (*arr, ealloc));
  SC3E (sc3_array_set_elem_count (*arr, ecount));
  SC3E (sc3_array_set_initzero (*arr, 1));
  SC3E (sc3_array_set_resizable (*arr, is_resizable));
  SC3E (sc3_array_set_tighten (*arr, 0));
  SC3E (sc3_array_setup (*arr));

  return NULL;
}

#ifdef P4EST_ENABLE_DEBUG
static sc3_error_t *
p4est3_array_set_zero (sc3_array_t * arr)
{
  int                 ecount;
  size_t              esize;
  void               *idx;
  SC3E (sc3_array_get_elem_count (arr, &ecount));
  SC3E (sc3_array_get_elem_size (arr, &esize));
  SC3E (sc3_array_index (arr, 0, &idx));
  memset (idx, 0, ecount * esize);
  return NULL;
}
#endif

static sc3_error_t *
p4est3_set_children_face_neighbors (p4est3_t * p3, p4est3_search_area_t * sa)
{
  int                *cfn;
  SC3E (sc3_allocator_calloc (p3->alloc, sa->nfaces * (1 << p3->qvt->dim),
                              sizeof (int), &sa->children_face_neighbors));
  cfn = sa->children_face_neighbors;
  if (p3->qvt->dim == 2) {
    /* *INDENT-OFF* */
    cfn[0] = -1; cfn[1] = 1; cfn[2] = -1; cfn[3] = 2;
    cfn[4] = 0; cfn[5] = -1; cfn[6] = -1; cfn[7] = 3;
    cfn[8] = -1; cfn[9] = 3; cfn[10] = 0; cfn[11] = -1;
    cfn[12] = 2; cfn[13] = -1; cfn[14] = 1; cfn[15] = -1;
    /* *INDENT-ON* */
  }
  else if (p3->qvt->dim == 3) {
    /* *INDENT-OFF* */
    cfn[0] = -1; cfn[1] = 1; cfn[2] = -1; cfn[3] = 2; cfn[4] = -1; cfn[5] = 4;
    cfn[6] = 0; cfn[7] = -1; cfn[8] = -1; cfn[9] = 3; cfn[10] = -1; cfn[11] = 5;
    cfn[12] = -1; cfn[13] = 3; cfn[14] = 0; cfn[15] = -1; cfn[16] = -1; cfn[17] = 6;
    cfn[18] = 2; cfn[19] = -1; cfn[20] = 1; cfn[21] = -1; cfn[22] = -1; cfn[23] = 7;
    cfn[24] = -1; cfn[25] = 5; cfn[26] = -1; cfn[27] = 6; cfn[28] = 0; cfn[29] = -1;
    cfn[30] = 4; cfn[31] = -1; cfn[32] = -1; cfn[33] = 7; cfn[34] = 1; cfn[35] = -1;
    cfn[36] = -1; cfn[37] = 7; cfn[38] = 4; cfn[39] = -1; cfn[40] = 2; cfn[41] = -1;
    cfn[42] = 6; cfn[43] = -1; cfn[44] = 5; cfn[45] = -1; cfn[46] = 3; cfn[47] = -1;
    /* *INDENT-ON* */
  }
  return NULL;
}

static sc3_error_t *
p4est3_set_face_dual (p4est3_t * p3, p4est3_search_area_t * sa)
{
  int                *fd;
  SC3E (sc3_allocator_calloc
        (p3->alloc, 2 * p3->qvt->dim, sizeof (int), &sa->face_dual));
  fd = sa->face_dual;
  if (p3->qvt->dim == 2) {
    fd[0] = 1;
    fd[1] = 0;
    fd[2] = 3;
    fd[3] = 2;
  }
  else if (p3->qvt->dim == 3) {
    fd[0] = 1;
    fd[1] = 0;
    fd[2] = 3;
    fd[3] = 2;
    fd[4] = 5;
    fd[5] = 4;
  }
  return NULL;
}

static sc3_error_t *
p4est3_set_outer_data (p4est3_t * p3, p4est3_search_area_t * sa,
                       void *user_data)
{
  const int           ntypes = p3->num_children + 1;
  int                 i, side;
  void               *arr;
  /*set general section of sa */
  sa->max_children = p3->num_children;
  sa->nfaces = 2 * p3->qvt->dim;
  sa->start_level = 0;
  SC3E (p4est3_set_children_face_neighbors (p3, sa));
  SC3E (p4est3_set_face_dual (p3, sa));

  /*set volume section of sa */
  SC3E (p4est3_tree_index (p3, p3->fltree, &sa->tree));
  sa->begin = NULL;
  sa->end = NULL;
  sa->Level = 0;
  SC3E (sc3_allocator_calloc (p3->alloc, p3->qvt->max_level, sizeof (int),
                              &sa->level2nchildren));
  memset (sa->level2nchildren, 0, sizeof (int) * p3->qvt->max_level);

  /* set 2d stack (array of arrays) */
  SC3E (p4est3_array_new (p3->alloc, sizeof (sc3_array_t *),
                          p3->qvt->max_level, p3->qvt->max_level, 1,
                          &sa->idx_vol_stack));
  SC3E (sc3_array_index (sa->idx_vol_stack, 0, &arr));
  SC3E (p4est3_array_new (p3->alloc, sizeof (p4est3_locidx), 2, 2, 1,
                          (sc3_array_t **) arr));
  for (i = 1; i < p3->qvt->max_level; ++i) {
    SC3E (sc3_array_index (sa->idx_vol_stack, i, &arr));
    SC3E (p4est3_array_new (p3->alloc, sizeof (p4est3_locidx), ntypes, ntypes,
                            1, (sc3_array_t **) arr));
  }
  SC3E (sc3_array_resize (sa->idx_vol_stack, 1));

  SC3E (sc3_allocator_calloc_one (p3->alloc,
                                  sizeof (p4est3_iterate_volume_info_t),
                                  &sa->vinfo));
  sa->vinfo->p3 = p3;
  sa->vinfo->user_data = user_data;

  /*set face section of sa */
  sa->nsides = 2;
  SC3E (sc3_allocator_calloc_one (p3->alloc,
                                  sizeof (p4est3_iterate_face_info_t),
                                  &sa->finfo));
  SC3E (p4est3_array_new (p3->alloc, sizeof (p4est3_iterate_face_side_t), 2,
                          2, 1, &sa->finfo->sides));
  sa->finfo->p3 = p3;
  sa->finfo->user_data = user_data;
  for (side = 0; side < 2; ++side) {
    sa->tree_face[side] = NULL;
    sa->begin_face[side] = NULL;
    sa->end_face[side] = NULL;
    sa->Level_face[side] = 0;
    sa->is_refine[side] = 1;

    /* set 2d stack (array of arrays) */
    SC3E (p4est3_array_new (p3->alloc, sizeof (sc3_array_t *),
                            p3->qvt->max_level, p3->qvt->max_level, 1,
                            &sa->idx_face_stack[side]));
    SC3E (sc3_array_index (sa->idx_face_stack[side], 0, &arr));
    SC3E (p4est3_array_new
          (p3->alloc, sizeof (p4est3_locidx), 2, 2, 1, (sc3_array_t **) arr));
    for (i = 1; i < p3->qvt->max_level; ++i) {
      SC3E (sc3_array_index (sa->idx_face_stack[side], i, &arr));
      SC3E (p4est3_array_new
            (p3->alloc, sizeof (p4est3_locidx), ntypes, ntypes, 1,
             (sc3_array_t **) arr));
    }
    SC3E (sc3_array_resize (sa->idx_face_stack[side], 1));
  }

  /*temporarily set view on level2children array to be able to renew it later
     instead of making new / removing */
  SC3E (sc3_array_new_data (p3->alloc, &sa->view_quads, sa->tree->tquads,
                            p3->qvt->quadrant_size, 0, 0));
  return NULL;
}

static sc3_error_t *
p4est3_destroy_outer_data (p4est3_t * p3, p4est3_search_area_t * sa)
{
  int                 i, side;
  void               *arr;
  SC3E (sc3_allocator_free (p3->alloc, sa->children_face_neighbors));
  SC3E (sc3_allocator_free (p3->alloc, sa->face_dual));
  SC3E (sc3_allocator_free (p3->alloc, sa->level2nchildren));

  SC3E (sc3_array_resize (sa->idx_vol_stack, p3->qvt->max_level));
  for (i = 0; i < p3->qvt->max_level; ++i) {
    SC3E (sc3_array_index (sa->idx_vol_stack, i, &arr));
    SC3E (sc3_array_destroy ((sc3_array_t **) arr));
  }
  SC3E (sc3_array_destroy (&sa->idx_vol_stack));

  for (side = 0; side < 2; ++side) {
    SC3E (sc3_array_resize (sa->idx_face_stack[side], p3->qvt->max_level));
    for (i = 0; i < p3->qvt->max_level; ++i) {
      SC3E (sc3_array_index (sa->idx_face_stack[side], i, &arr));
      SC3E (sc3_array_destroy ((sc3_array_t **) arr));
    }
  }
  SC3E (sc3_array_destroy (&sa->idx_face_stack[0]));
  SC3E (sc3_array_destroy (&sa->idx_face_stack[1]));
  SC3E (sc3_array_destroy (&sa->view_quads));
  SC3E (sc3_array_destroy (&sa->finfo->sides));
  SC3E (sc3_allocator_free (p3->alloc, sa->vinfo));
  SC3E (sc3_allocator_free (p3->alloc, sa->finfo));

  return NULL;
}

typedef struct p4est3_qvt_non_const_wrapper
{
  const p4est3_quadrant_vtable_t *qvt;
} p4est3_qvt_non_const_wrapper_t;

#ifdef P4EST_ENABLE_DEBUG
static sc3_error_t *
p4est3_array_is_sorted (const void *q1, const void *q2, void *qvt_wrapper, int *j)
{
  const p4est3_quadrant_vtable_t *qvtable = ((p4est3_qvt_non_const_wrapper_t *) qvt_wrapper)->qvt;
  SC3A_IS (p4est3_quadrant_vtable_is_valid, qvtable);

  SC3A_IS (qvtable->quadrant_is_valid, q1);
  SC3A_IS (qvtable->quadrant_is_valid, q2);
  SC3A_CHECK (j != NULL);

  SC3E (qvtable->quadrant_compare (q1, q2, j));
  return NULL;
}
#endif

static int
p4est3_get_children_face_nb_id (p4est3_search_area_t * sa,
                                int child_id, int face)
{
  return sa->children_face_neighbors[child_id * sa->nfaces + face];
}

static int
p4est3_get_dual_face (p4est3_search_area_t * sa, int face)
{
  return sa->face_dual[face];
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
  SC3A_IS3 (sc3_array_is_sorted, array, p4est3_array_is_sorted, qvtw);

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
  SC3E (sc3_array_split (array, indices, p4est3_quadrant_num_children (qvt),
                         p4est3_array_split_ancestor_id, &data));
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

sc3_error_t        *
p4est3_iterate_volume (p4est3_t * p3,
                       p4est3_iterate_volume_t cvolume, void *user_data)
{
  p4est3_topidx       ntree;
  p4est3_tree_t      *tree;
  p4est3_iterate_volume_info_t info;
  p4est3_locidx       si;

  SC3A_IS (p4est3_is_setup, p3);
  if (p3->fltree < 0 || cvolume == NULL) {
    return NULL;
  }

  info.p3 = p3;
  info.user_data = user_data;
  for (ntree = p3->fltree; ntree <= p3->lltree; ++ntree) {
    info.ntree = ntree;
    SC3E (p4est3_tree_index (p3, ntree, &tree));
    for (si = 0; si < tree->num_quads; ++si) {
      info.quadrant = tree->tquads + si * p3->qvt->quadrant_size;
      info.nquad = si + tree->first_tquad;
      SC3E (cvolume (&info));
    }
  }
  return NULL;
}

static sc3_error_t *
p4est3_iterate_face_bound_init (p4est3_t * p3,
                                p4est3_search_area_t * search_area,
                                p4est3_topidx tree,
                                int face,
                                p4est3_iterate_face_side_t * fside,
                                int *is_lower)
{
  int                 orient;
  int                *Level_face = search_area->Level_face;
  int                *is_refine = search_area->is_refine;
  p4est3_topidx       tree_neighbor = tree;
  sc3_array_t       **idx_f_stack = search_area->idx_face_stack;
  sc3_array_t        *arr;
  p4est3_locidx     **b_f = search_area->begin_face,
    **e_f = search_area->end_face;
  int                 side;

  is_refine[0] = is_refine[1] = 0;
  fside[0].ntree = tree;
  fside[0].nface = face;
  SC3E (p4est3_connectivity_get_face
        (p3->conn, &tree_neighbor, &face, &orient));
  if (tree_neighbor > p3->lltree || tree_neighbor < p3->fltree) {
    *is_lower = 1;
    return NULL;
  }
  if (tree_neighbor < tree) {
    *is_lower = 1;
    return NULL;
  }
  search_area->finfo->orientation = orient;

  if (tree_neighbor == tree) {
    /* physical boundary */
    search_area->nsides = 1;
  }
  else {
    search_area->nsides = 2;
    fside[1].ntree = tree_neighbor;
    SC3E (p4est3_tree_index (p3, tree_neighbor, &search_area->tree_face[1]));
    fside[1].nface = face;
  }
  SC3E (sc3_array_resize (search_area->finfo->sides, search_area->nsides));
  for (side = 0; side < search_area->nsides; ++side) {
    Level_face[side] = 0;
    is_refine[side] = 1;
    SC3E (sc3_array_index (idx_f_stack[side], 0, &arr));
#ifdef P4EST_ENABLE_DEBUG
    int                 ecount;
    SC3E (sc3_array_get_elem_count (idx_f_stack[side], &ecount));
    SC3A_CHECK (ecount == 1);
#endif
    SC3E (sc3_array_index (*(sc3_array_t **) arr, 0, &b_f[side]));
    *(b_f[side]) = 0;
    SC3E (sc3_array_index (*(sc3_array_t **) arr, 1, &e_f[side]));
    *(e_f[side]) = search_area->tree_face[side]->num_quads;
  }
  *is_lower = 0;
  return NULL;
}

static sc3_error_t *
p4est3_internal_iterate_face (p4est3_t * p3,
                              p4est3_iterate_face_t cface,
                              p4est3_iterate_codim_t ccodim,
                              p4est3_search_area_t * search_area)
{
  const int           max_children = p3->num_children;
  const int           half_ch = max_children / 2;
  p4est3_tree_t     **trees = search_area->tree_face;
  p4est3_locidx     **b_f = search_area->begin_face;
  p4est3_locidx     **e_f = search_area->end_face;
  void               *stack_it[2];
  p4est3_locidx      *arr_it;
  sc3_array_t        *view_q = search_area->view_quads;
  sc3_array_t       **idx_face_stack = search_area->idx_face_stack;
  p4est3_iterate_face_side_t *fside;
  int                *is_refine = search_area->is_refine;
  int                *Level = search_area->Level_face;
  int                 i, side, level, idx, child_id;
  int                 ori = search_area->finfo->orientation;
  void               *first_quad;

  /* Check if both sides belong to the same process (at least, partly).
     If not, we ignore this face. */
  for (side = 0; side < search_area->nsides; ++side) {
    if (*b_f[side] == *e_f[side]) {
      return NULL;
    }
  }
  /* first check if the whole quadrant passes */
  SC3E (sc3_array_index (search_area->finfo->sides, 0, &fside));
  for (side = 0; side < search_area->nsides; ++side) {
    if (!is_refine[side]) {
      continue;
    }
    first_quad =
      (void *) (trees[side]->tquads +
                p3->qvt->quadrant_size * (*(b_f[side])));
    SC3E (p4est3_quadrant_level (p3->qvt, first_quad, &level));
    if (level == Level[side]) {
      is_refine[side] = 0;
      fside[side].nquad = *(b_f[side]) + trees[side]->first_tquad;
      fside[side].quadrant = first_quad;
    }
  }
  if (!is_refine[0] && !is_refine[1]) {
    if (cface != NULL) {
      SC3E (cface (search_area->finfo));
    }
    for (side = 0; side < search_area->nsides; ++side) {
      is_refine[side] = 1;
    }
    return NULL;
  }
  for (side = 0; side < search_area->nsides; ++side) {
    if (!is_refine[side]) {
      continue;
    }
    SC3E (sc3_array_push (idx_face_stack[side], &(stack_it[side])));
#ifdef P4EST_ENABLE_DEBUG
    SC3E (p4est3_array_set_zero (*(sc3_array_t **) (stack_it[side])));
#endif
    SC3E (sc3_array_renew_data (&view_q, trees[side]->tquads,
                                p3->qvt->quadrant_size, *(b_f[side]),
                                *(e_f[side]) - *(b_f[side])));
    SC3E (p4est3_quadrant_array_split
          (p3->qvt, view_q, Level[side], *(sc3_array_t **) (stack_it[side])));

    /* since array_split doesn't count shift from the beinning of quadrants
       in a tree, we shift result indices at the loop below */
    for (i = 0; i < max_children + 1; ++i) {
      SC3E (sc3_array_index (*(sc3_array_t **) (stack_it[side]), i, &arr_it));
      *arr_it += *(b_f[side]);
    }
  }
  for (i = 0; i < half_ch; ++i) {
    for (side = 0; side < search_area->nsides; ++side) {
      idx = i;
      SC3E (sc3_array_index (*(sc3_array_t **) (stack_it[side]), 0, &arr_it));
      if (!is_refine[side]) {
        continue;
      }
      if (side == 1) {
        SC3E (p4est3_connectivity_get_neighbor_face_corner
              (p3->conn, fside[0].nface, fside[1].nface, ori, &idx));
      }
      SC3E (p4est3_connectivity_get_face_child_id
            (p3->conn, fside[side].nface, idx, &child_id));
      b_f[side] = arr_it + child_id;
      e_f[side] = arr_it + child_id + 1;
    }
    Level[0]++;
    Level[1]++;
    SC3E (p4est3_internal_iterate_face (p3, cface, ccodim, search_area));
    Level[0]--;
    Level[1]--;
  }
  for (side = 0; side < search_area->nsides; ++side) {
    if (!is_refine[side]) {
      continue;
    }
    SC3E (sc3_array_pop (idx_face_stack[side]));
  }

  return NULL;
}

static sc3_error_t *
p4est3_iterate_face_inner_init (p4est3_t * p3,
                                p4est3_search_area_t * search_area,
                                p4est3_locidx * arr_it,
                                int child, int neighbor)
{
  int                *Level_face = search_area->Level_face;
  int                *is_refine = search_area->is_refine;

  sc3_array_t       **idx_f_stack = search_area->idx_face_stack;
  void               *arr;
  p4est3_locidx     **b_f = search_area->begin_face,
    **e_f = search_area->end_face;

  search_area->nsides = 2;
  Level_face[0] = Level_face[1] = search_area->Level;
  is_refine[0] = is_refine[1] = 1;
#ifdef P4EST_ENABLE_DEBUG
  int                 ecount, side;
  for (side = 0; side < search_area->nsides; ++side) {
    SC3E (sc3_array_get_elem_count (idx_f_stack[side], &ecount));
    SC3A_CHECK (ecount == 1);
  }
#endif

  SC3E (sc3_array_index (idx_f_stack[0], 0, &arr));
  SC3E (sc3_array_index (*(sc3_array_t **) arr, 0, &(b_f[0])));
  *(b_f[0]) = *(arr_it + child);
  SC3E (sc3_array_index (*(sc3_array_t **) arr, 1, &(e_f[0])));
  *(e_f[0]) = *(arr_it + child + 1);

  SC3E (sc3_array_index (idx_f_stack[1], 0, &arr));
  SC3E (sc3_array_index (*(sc3_array_t **) arr, 0, &(b_f[1])));
  *(b_f[1]) = *(arr_it + neighbor);
  SC3E (sc3_array_index (*(sc3_array_t **) arr, 1, &(e_f[1])));
  *(e_f[1]) = *(arr_it + neighbor + 1);

  return NULL;
}

static sc3_error_t *
p4est3_iterate_face_inner (p4est3_t * p3,
                           p4est3_iterate_face_t cface,
                           p4est3_iterate_codim_t ccodim,
                           p4est3_search_area_t * search_area,
                           p4est3_locidx * arr_it)
{

  int                 child, face, nb_id;
  p4est3_iterate_face_side_t *fside;
  search_area->tree_face[0] = search_area->tree_face[1] = search_area->tree;
  SC3E (sc3_array_index (search_area->finfo->sides, 0, &fside));
  for (child = 0; child < search_area->max_children; ++child) {
    for (face = 0; face < search_area->nfaces; ++face) {
      nb_id = p4est3_get_children_face_nb_id (search_area, child, face);
      if (nb_id < child) {
        continue;
      }
      fside[0].nface = face;
      fside[1].nface = p4est3_get_dual_face (search_area, face);
      SC3E (p4est3_iterate_face_inner_init
            (p3, search_area, arr_it, child, nb_id));
      SC3E (p4est3_internal_iterate_face (p3, cface, ccodim, search_area));
    }
  }
  return NULL;
}

static sc3_error_t *
p4est3_iterate_volume_rec_init (p4est3_t * p3,
                                p4est3_search_area_t * sa, p4est3_topidx tree)
{
  void               *arr;
  p4est3_iterate_face_side_t *fside;

  SC3E (p4est3_tree_index (p3, tree, &sa->tree));
  sa->vinfo->ntree = tree;
  sa->finfo->orientation = 0;
  sa->finfo->tree_boundary = 0;
  SC3E (sc3_array_resize (sa->finfo->sides, 2));

  memset (sa->level2nchildren, 0, sizeof (int) * p3->qvt->max_level);
  SC3A_CHECK (sa->Level == 0);

  SC3E (sc3_array_index (sa->idx_vol_stack, 0, &arr));
#ifdef P4EST_ENABLE_DEBUG
  int                 ecount;
  SC3E (sc3_array_get_elem_count (*(sc3_array_t **) arr, &ecount));
  SC3A_CHECK (ecount == 2);
  SC3E (sc3_array_get_elem_count (sa->idx_vol_stack, &ecount));
  SC3A_CHECK (ecount == 1);
#endif

  SC3E (sc3_array_index (*(sc3_array_t **) arr, 0, &sa->begin));
  *(sa->begin) = 0;
  sa->end = sa->begin + 1;
  *(sa->end) = sa->tree->num_quads;

  SC3E (sc3_array_index (sa->finfo->sides, 0, &fside));
  for (int side = 0; side < 2; ++side) {
    fside[side].ntree = tree;
    fside[side].is_ghost = 0;
  }

  return NULL;
}

static sc3_error_t *
p4est3_iterate_volume_rec (p4est3_t * p3,
                           p4est3_iterate_volume_t cvolume,
                           p4est3_iterate_face_t cface,
                           p4est3_iterate_codim_t ccodim,
                           p4est3_search_area_t * search_area)
{
  int                 i;
  void               *first_quad;       /*first quadrant in this search area */
  int                 level;
  void               *stack_it;
  p4est3_locidx      *arr_it;

  const int           max_children = p3->num_children;
  int                *l2nch = search_area->level2nchildren;
  int                *Level = &search_area->Level;
  const p4est3_locidx begin = *(search_area->begin);
  const p4est3_locidx end = *(search_area->end);
  p4est3_tree_t      *tree = search_area->tree;
  sc3_array_t        *view_q = search_area->view_quads;
  sc3_array_t        *idx_vol_stack = search_area->idx_vol_stack;
  p4est3_iterate_volume_info_t *vinfo = search_area->vinfo;

  /* Check if the considered search area intersect
     the area of the local process. If not, then skip it. */
  if (begin == end) {
    l2nch[*Level]++;
    return NULL;
  }
/*Are trees with no quads possible?*/
  first_quad = (void *) (tree->tquads + p3->qvt->quadrant_size * begin);
  SC3E (p4est3_quadrant_level (p3->qvt, first_quad, &level));
  if (level == *Level) {
    if (cvolume != NULL) {
      vinfo->quadrant = first_quad;
      vinfo->nquad = begin + tree->first_tquad;
      SC3E (cvolume (vinfo));
    }
    l2nch[*Level]++;
    return NULL;
  }

  if (l2nch[search_area->start_level] > 0) {
    return NULL;
  }

  SC3E (sc3_array_push (idx_vol_stack, &stack_it));
#ifdef P4EST_ENABLE_DEBUG
  SC3E (p4est3_array_set_zero (*(sc3_array_t **) stack_it));
#endif
  SC3E (sc3_array_renew_data (&view_q, tree->tquads,
                              p3->qvt->quadrant_size, begin, end - begin));
  SC3E (p4est3_quadrant_array_split (p3->qvt, view_q, *Level,
                                     *(sc3_array_t **) stack_it));
  l2nch[++(*Level)] = 0;

  /* since array_split doesn't count shift from the beinning of quadrants
     in a tree, we shift result indices at the loop below */
  for (i = 0; i < max_children + 1; ++i) {
    SC3E (sc3_array_index (*(sc3_array_t **) stack_it, i, &arr_it));
    *arr_it += begin;
  }
  SC3E (sc3_array_index (*(sc3_array_t **) stack_it, 0, &arr_it));
  for (i = 0; i < max_children; ++i) {
    search_area->begin = arr_it + i;
    search_area->end = arr_it + i + 1;
    /* if (*(search_area->begin) == *(search_area->end) - 1) {
       l2nch[*Level]++;
       continue;
       } */
    SC3E (p4est3_iterate_volume_rec
          (p3, cvolume, cface, ccodim, search_area));
  }
  SC3A_CHECK (l2nch[*Level] == max_children);
  SC3E (p4est3_iterate_face_inner (p3, cface, ccodim, search_area, arr_it));
  l2nch[--(*Level)]++;
  SC3E (sc3_array_pop (idx_vol_stack));
  return NULL;
}

sc3_error_t        *
p4est3_iterate_codim (p4est3_t * p3, int codims,
                      p4est3_iterate_volume_t cvolume,
                      p4est3_iterate_face_t cface,
                      p4est3_iterate_codim_t ccodim, void *user_data)
{
  p4est3_search_area_t ssa, *search_area = &ssa;
  p4est3_topidx       tree;
  p4est3_iterate_face_side_t *fside;
  int                 face, is_lower;

  /* This iteration is w/o ghost layer and for volumes only */
  if (codims < 0 || codims >= P4EST3_ITERATE_LAST) {
    return NULL;
  }
  SC3A_IS (p4est3_is_setup, p3);
  if (p3->fltree < 0 || (cvolume == NULL && cface == NULL && ccodim == NULL)) {
    return NULL;
  }

  if (codims == P4EST3_ITERATE_VOLUME) {
    SC3E (p4est3_iterate_volume (p3, cvolume, user_data));
    return NULL;
  }

  SC3E (p4est3_set_outer_data (p3, search_area, user_data));
  SC3E (sc3_array_index (search_area->finfo->sides, 0, &fside));
  for (tree = p3->fltree; tree <= p3->lltree; ++tree) {

    SC3E (p4est3_iterate_volume_rec_init (p3, search_area, tree));
    SC3E (p4est3_iterate_volume_rec
          (p3, cvolume, cface, ccodim, search_area));

    /* frame faces part */
    search_area->finfo->tree_boundary = 1;
    search_area->tree_face[0] = search_area->tree;
    for (face = 0; face < search_area->nfaces; ++face) {
      SC3E (p4est3_iterate_face_bound_init
            (p3, search_area, tree, face, fside, &is_lower));
      if (is_lower) {
        /* we iterate over such trees that tree_neighbor < tree */
        SC3E (sc3_array_pop (search_area->finfo->sides));
        continue;
      }
      SC3E (p4est3_internal_iterate_face (p3, cface, ccodim, search_area));
      if (is_lower) {
        SC3E (sc3_array_push (search_area->finfo->sides, NULL));
      }
    }
  }
  SC3E (p4est3_destroy_outer_data (p3, search_area));
  return NULL;
}

#ifdef __cplusplus
#if 0
{
#endif
}
#endif
