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
#include <p4est3_search.h>

#ifdef __cplusplus
extern              "C"
{
#if 0
}
#endif
#endif

typedef struct p4est3_search_area
{
  /* general section */
  int                 max_children;
  int                 nfaces;   /*Global number of faces */
  sc3_array_t        *view_quads;       /* array ptr to pass tree's quads
                                           into array_split */
  int                *children_face_neighbors;
  int                *face_dual;

  /* volume section */
  p4est3_tree_t      *tree;
  int                 child_id; /* Child id of the area under consideration */
  int                 Level;
  int                *level2nchildren;  /* Array specifing the number n
                                           of children processed on
                                           the particular level; n can't be
                                           greater than p3->num_children */
  sc3_array_t        *idx_vol_stack;    /* 2D stack storing arrays of indices,
                                           that are output of split_array */
  p4est3_gloidx       local_begin;      /* Id of the first local quadrant in
                                           the tree under iteration */
  p4est3_gloidx       local_end;        /* Id of the last local quadrant in
                                            the tree under iteration */
  p4est3_iterate_volume_info_t *vinfo;

  /* face section */
  int                 nsides;
  p4est3_topidx       treeid_face[2];
  int                 child_id_face[2];
  int                 Level_face[2];    /* array of 2, storing the level of
                                           current size */
  int                 is_refine[2];     /* array of 2, storing the information
                                           about neccesity of refenement. Must
                                           be initiated by 0 */
  sc3_array_t        *idx_face_stack[2];        /* 2D stacks storing arrays of indices,
                                                   that are output of split_array */
  p4est3_gloidx       remote_first[2];  /* Id of the first remote proc that
                                           potentially contains the neighbor
                                           quadrants */
  p4est3_gloidx       remote_last[2];   /* Id of the last remote proc that
                                           potentially contains the neighbor
                                           quadrants */
  p4est3_gloidx       local_begin_face[2]; /* Ids of the first local quadrant in
                                              the tree under iteration */
  p4est3_gloidx       local_end_face[2];   /* Ids of the last local quadrant in
                                              the tree under iteration */
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
  size_t              esize, ecount;
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
  int                 i, side, nodesize, noderank;
  void               *arr;
  /*set general section of sa */
  sa->max_children = p3->num_children;
  sa->nfaces = 2 * p3->qvt->dim;
  SC3E (p4est3_set_children_face_neighbors (p3, sa));
  SC3E (p4est3_set_face_dual (p3, sa));

  /*set volume section of sa */
  SC3E (p4est3_tree_index (p3, p3->fltree, &sa->tree));
  sa->child_id = -1;
  sa->Level = 0;
  SC3E (sc3_allocator_calloc (p3->alloc, p3->qvt->max_level, sizeof (int),
                              &sa->level2nchildren));
  memset (sa->level2nchildren, 0, sizeof (int) * p3->qvt->max_level);

  /* set 2d stack (array of arrays) */
  SC3E (p4est3_array_new (p3->alloc, sizeof (sc3_array_t *),
                          p3->qvt->max_level, p3->qvt->max_level, 1,
                          &sa->idx_vol_stack));
  SC3E (sc3_array_index (sa->idx_vol_stack, 0, &arr));
  SC3E (p4est3_array_new (p3->alloc, sizeof (p4est3_gloidx), 2, 2, 1,
                          (sc3_array_t **) arr));
  for (i = 1; i < p3->qvt->max_level; ++i) {
    SC3E (sc3_array_index (sa->idx_vol_stack, i, &arr));
    SC3E (p4est3_array_new (p3->alloc, sizeof (p4est3_gloidx), ntypes, ntypes,
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
  SC3E (sc3_mpienv_get_nodesize (p3->split_info, &nodesize));
  SC3E (sc3_mpienv_get_noderank (p3->split_info, &noderank));
  sa->finfo->p3 = p3;
  sa->finfo->user_data = user_data;
  for (side = 0; side < 2; ++side) {
    sa->treeid_face[side] = -1;
    sa->child_id_face[side] = -1;
    sa->Level_face[side] = 0;
    sa->is_refine[side] = 1;

    /* set 2d stack (array of arrays) */
    SC3E (p4est3_array_new (p3->alloc, sizeof (sc3_array_t *),
                            p3->qvt->max_level, p3->qvt->max_level, 1,
                            &sa->idx_face_stack[side]));
    SC3E (sc3_array_index (sa->idx_face_stack[side], 0, &arr));
    SC3E (p4est3_array_new
          (p3->alloc, sizeof (p4est3_gloidx), 2, 2, 1, (sc3_array_t **) arr));
    for (i = 1; i < p3->qvt->max_level; ++i) {
      SC3E (sc3_array_index (sa->idx_face_stack[side], i, &arr));
      SC3E (p4est3_array_new
            (p3->alloc, sizeof (p4est3_gloidx), ntypes, ntypes, 1,
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
                                p4est3_search_area_t * sa,
                                p4est3_topidx tree,
                                int face,
                                p4est3_iterate_face_side_t * fside,
                                int *is_iterate)
{
  int                 orient;
  int                *Level_face = sa->Level_face;
  int                *is_refine = sa->is_refine;
  p4est3_topidx       tree_neighbor = tree, tree_ids[2];
  sc3_array_t       **idx_f_stack = sa->idx_face_stack;
  sc3_array_t        *arr;
  p4est3_locidx      *begin, *end;
  int                 s;

  SC3E_RETVAL (is_iterate, 1);
  is_refine[0] = is_refine[1] = 0;
  fside[0].ntree = tree;
  fside[0].nface = face;
  SC3E (p4est3_connectivity_get_face
        (p3->conn, &tree_neighbor, &face, &orient));
  /*if (tree_neighbor > p3->lltree || tree_neighbor < p3->fltree) {
   *is_lower = 1;
   return NULL;
   }
   if (tree_neighbor < tree) {
   *is_lower = 1;
   return NULL;
   }*/
  if (p3->fltree <= tree_neighbor && tree_neighbor <= p3->lltree
      && tree_neighbor > tree) {
    *is_iterate = 0;
    return NULL;
  }
  sa->finfo->orientation = orient;

  if (tree_neighbor == tree) {
    /* physical boundary */
    sa->nsides = 1;
  }
  else {
    sa->nsides = 2;
    fside[1].ntree = tree_neighbor;
    sa->treeid_face[1] = tree_neighbor;
    fside[1].nface = face;
  }
  SC3E (sc3_array_resize (sa->finfo->sides, sa->nsides));
  tree_ids[0] = tree;
  tree_ids[1] = tree_neighbor;
  /* Find here two ranges of procs that:
     0th -- shares 0th side
     1st -- shares 1st side */
  for (s = 0; s < sa->nsides; ++s) {
    sa->remote_first[s] = -1;
    sa->remote_last[s] = -2;
  }
  if (sa->nsides == 2) {
    /* only if we probe not a physical boundary */
    for (s = 0; s < sa->nsides; ++s) {
      SC3E (p4est3_find_partition
            (p3->alloc, p3->mpisize, p3->goffset, p3->gtroffset[tree_ids[s]],
             p3->gtroffset[tree_ids[s] + 1] - 1,
             &sa->remote_first[s], &sa->remote_last[s]));
      if (p3->goffset[sa->remote_first[s]] > p3->gtroffset[tree_ids[s]]) {
        sa->remote_first[s]--;
      }
      if (p3->goffset[sa->remote_last[s]] >
          p3->gtroffset[tree_ids[s] + 1] - 1) {
        sa->remote_last[s]--;
      }
    }
  }

  for (s = 0; s < sa->nsides; ++s) {
    Level_face[s] = 0;
    is_refine[s] = 1;
    sa->child_id_face[s] = 0;
    SC3E (sc3_array_resize (idx_f_stack[s], 1));
    SC3E (sc3_array_index (idx_f_stack[s], 0, &arr));
    SC3E (sc3_array_index (*(sc3_array_t **) arr, 0, &begin));
    SC3E (sc3_array_index (*(sc3_array_t **) arr, 1, &end));

    sa->local_begin_face[s] =
      SC3_MAX (p3->gtroffset[tree_ids[s]], p3->goffset[p3->mpirank]);
    sa->local_end_face[s] =
      SC3_MIN (p3->gtroffset[tree_ids[s] + 1], p3->goffset[p3->mpirank + 1]);

    *begin = p3->gtroffset[tree_ids[s]];
    *end = p3->gtroffset[tree_ids[s] + 1];
  }

  return NULL;
}

static sc3_error_t *
p4est3_internal_iterate_face (p4est3_t * p3,
                              p4est3_iterate_face_t cface,
                              p4est3_iterate_codim_t ccodim,
                              p4est3_search_area_t * sa)
{
  const int           max_children = p3->num_children;
  const int           half_ch = max_children / 2;
  p4est3_topidx      *trees = sa->treeid_face;
  p4est3_gloidx      *b_f[2], *e_f[2];
  p4est3_gloidx       proc_owner;
  void               *stack_it[2];
  p4est3_gloidx      *arr_it;
  sc3_array_t        *view_q = sa->view_quads;
  sc3_array_t       **idx_face_stack = sa->idx_face_stack;
  p4est3_iterate_face_side_t *fside;
  int                *is_refine = sa->is_refine;
  int                *Level = sa->Level_face;
  int                 i, side, level, idx, child_id;
  int                 ori = sa->finfo->orientation;
  int                 is_lvl_increased[2] = { 0, 0 };
  void               *first_quad;

  for (side = 0; side < sa->nsides; ++side) {
    SC3E (sc3_array_index
          (idx_face_stack[side], Level[side], &(stack_it[side])));
    SC3E (sc3_array_index
          (*(sc3_array_t **) stack_it[side],
           sa->child_id_face[side], &(b_f[side])));
    SC3E (sc3_array_index
          (*(sc3_array_t **) stack_it[side],
           sa->child_id_face[side] + 1, &(e_f[side])));
  }

  /* Check if both sides belong to remote process(es).
     If so, we ignore this face. */
  if ((sa->nsides == 1 && 
        (*(b_f[0]) >= sa->local_end_face[0] || *(e_f[0]) <= sa->local_begin_face[0]))
        ||
        (sa->nsides == 2 && 
        (*(b_f[0]) >= sa->local_end_face[0] || *(e_f[0]) <= sa->local_begin_face[0]) &&
        (*(b_f[1]) >= sa->local_end_face[1] || *(e_f[1]) <= sa->local_begin_face[1]))) {
    return NULL;
  }

  /* first check if the whole quadrant passes */
  SC3E (sc3_array_index (sa->finfo->sides, 0, &fside));
  for (side = 0; side < sa->nsides; ++side) {
    if (!is_refine[side]) {
      continue;
    }

    /* TODO: make it through p4est3_find_partition and last_goffsets */
    /* Find process that owns this quadrant */
    /* Start binary search from a local process */
    proc_owner = p3->mpirank;
    SC3E (p4est3_search_lower_bound64
            (*(b_f[side]), p3->goffset, p3->mpisize + 1, &proc_owner));
    if (p3->goffset[proc_owner] > *(b_f[side])) {
      SC3A_CHECK (proc_owner > 0);
      proc_owner--;
    }
    SC3A_CHECK (proc_owner >= 0 && proc_owner < p3->mpisize);

    /* Get quadrant from the owning process's window */
    first_quad = (void *) (p3->nodequads[proc_owner] +
                          p3->qsize * (*(b_f[side]) - p3->goffset[proc_owner]));

    SC3E (p4est3_quadrant_level (p3->qvt, first_quad, &level));
    if (level == Level[side]) {
      is_refine[side] = 0;
      fside[side].nquad = *(b_f[side]) - p3->gtroffset[trees[side]];
      fside[side].quadrant = first_quad;
      fside[side].is_ghost = proc_owner == p3->mpirank ? 0 : 1;
    }
  }

  if (!is_refine[0] && !is_refine[1]) {
    if (cface != NULL) {
      SC3E (cface (sa->finfo));
    }
    /*for (side = 0; side < sa->nsides; ++side) {
       is_refine[side] = 1;
       } */
    return NULL;
  }
  for (side = 0; side < sa->nsides; ++side) {
    if (!is_refine[side]) {
      continue;
    }
    /* we split array taht is unite for and local and remote procs quads */

    SC3E (sc3_array_push (idx_face_stack[side], &(stack_it[side])));
#ifdef P4EST_ENABLE_DEBUG
    SC3E (p4est3_array_set_zero (*(sc3_array_t **) (stack_it[side])));
#endif
    /* here we start with the very beginning of not necessary local node
       quadrants, because of our specialized array_split_noncontig function */
    SC3E (sc3_array_renew_data (&view_q, p3->nodequads[0], p3->qsize,
                                0, *(e_f[side]) - *(b_f[side])));
    SC3E (p4est3_quadrant_array_split_noncontig
          (p3, view_q, Level[side], *(b_f[side]),
           *(sc3_array_t **) (stack_it[side])));

    /* since array_split doesn't count shift from the beinning of quadrants
       in a proc, we shift result indices at the loop below */
    for (i = 0; i < max_children + 1; ++i) {
      SC3E (sc3_array_index (*(sc3_array_t **) (stack_it[side]), i, &arr_it));
      *arr_it += *(b_f[side]);
    }
  }
  for (i = 0; i < half_ch; ++i) {
    is_lvl_increased[0] = is_lvl_increased[1] = 0;
    for (side = 0; side < sa->nsides; ++side) {
      /* Some service computations to derive direction of face-face connection */
      if (!is_refine[side]) {
        continue;
      }
      idx = i;
      if (side == 1) {
        SC3E (p4est3_connectivity_get_neighbor_face_corner
              (p3->conn, fside[0].nface, fside[1].nface, ori, &idx));
      }
      SC3E (p4est3_connectivity_get_face_child_id
            (p3->conn, fside[side].nface, idx, &child_id));
      sa->child_id_face[side] = child_id;
      Level[side]++;
      is_lvl_increased[side] = 1;
    }
    SC3E (p4est3_internal_iterate_face (p3, cface, ccodim, sa));
    for (side = 0; side < sa->nsides; ++side) {
      Level[side] = is_lvl_increased[side] ? Level[side] - 1 : Level[side];
      if (i != half_ch - 1) {
        is_refine[side] = is_lvl_increased[side];
      }
      else if (is_lvl_increased[side]) {
        SC3E (sc3_array_pop (idx_face_stack[side]));
      }
    }
  }

  return NULL;
}

static sc3_error_t *
p4est3_iterate_face_inner_init (p4est3_t * p3, p4est3_search_area_t * sa,
                                int child, int neighbor)
{
  int                *Level_face = sa->Level_face;
  int                *is_refine = sa->is_refine;
  int                 ch_neigh[2] = { child, neighbor }, s;

  sc3_array_t       **idx_f_stack = sa->idx_face_stack;
  void               *top;      /* generic top of a stack */
  p4est3_gloidx      *begin, *end;
  p4est3_gloidx      *arr_vol_it;

  sa->nsides = 2;
  for (s = 0; s < sa->nsides; ++s) {
    Level_face[s] = sa->Level;
    is_refine[s] = 1;
    sa->child_id_face[s] = ch_neigh[s];
    SC3E (sc3_array_resize (idx_f_stack[s], sa->Level + 1));
    SC3E (sc3_array_index (idx_f_stack[s], Level_face[s], &top));
    SC3E (sc3_array_index (*(sc3_array_t **) top, ch_neigh[s], &begin));
    SC3E (sc3_array_index (*(sc3_array_t **) top, ch_neigh[s] + 1, &end));
    SC3E (sc3_array_index (sa->idx_vol_stack, Level_face[s], &top));
    SC3E (sc3_array_index (*(sc3_array_t **) top, 0, &arr_vol_it));
    *(begin) = *(arr_vol_it + ch_neigh[s]);
    *(end) = *(arr_vol_it + ch_neigh[s] + 1);
    /* Since we iterate inner faces, we inherit the local boundaries for volumes */
    sa->local_begin_face[s] = sa->local_begin;
    sa->local_end_face[s] = sa->local_end;
  }

  return NULL;
}

static sc3_error_t *
p4est3_iterate_face_inner (p4est3_t * p3,
                           p4est3_iterate_face_t cface,
                           p4est3_iterate_codim_t ccodim,
                           p4est3_search_area_t * search_area)
{

  int                 child, face, nb_id;
  p4est3_iterate_face_side_t *fside;
  search_area->treeid_face[0] = search_area->treeid_face[1]
    = search_area->tree->treeid;
  SC3E (sc3_array_index (search_area->finfo->sides, 0, &fside));
  for (child = 0; child < search_area->max_children; ++child) {
    for (face = 0; face < search_area->nfaces; ++face) {
      nb_id = p4est3_get_children_face_nb_id (search_area, child, face);
      if (nb_id < child) {
        continue;
      }
      fside[0].nface = face;
      fside[1].nface = p4est3_get_dual_face (search_area, face);
      SC3E (p4est3_iterate_face_inner_init (p3, search_area, child, nb_id));
      SC3E (p4est3_internal_iterate_face (p3, cface, ccodim, search_area));
    }
  }
  return NULL;
}

static sc3_error_t *
p4est3_iterate_volume_rec_init (p4est3_t * p3,
                                p4est3_search_area_t * sa, p4est3_topidx tree)
{
  int                 side;
  void               *arr;
  p4est3_iterate_face_side_t *fside;
  p4est3_gloidx      *begin, *end;

  SC3E (p4est3_tree_index (p3, tree, &sa->tree));
  sa->vinfo->ntree = tree;
  sa->finfo->orientation = 0;
  sa->finfo->tree_boundary = 0;
  SC3E (sc3_array_resize (sa->finfo->sides, 2));

  memset (sa->level2nchildren, 0, sizeof (int) * p3->qvt->max_level);
  SC3A_CHECK (sa->Level == 0);

  SC3E (sc3_array_index (sa->idx_vol_stack, 0, &arr));
#ifdef P4EST_ENABLE_DEBUG
  size_t              ecount;
  SC3E (sc3_array_get_elem_count (*(sc3_array_t **) arr, &ecount));
  SC3A_CHECK (ecount == 2);
  SC3E (sc3_array_get_elem_count (sa->idx_vol_stack, &ecount));
  SC3A_CHECK (ecount == 1);
#endif

  sa->child_id = 0;
  SC3E (sc3_array_index (*(sc3_array_t **) arr, sa->child_id, &begin));
  SC3E (sc3_array_index (*(sc3_array_t **) arr, sa->child_id + 1, &end));

  *begin = p3->gtroffset[tree];
  *end = p3->gtroffset[tree + 1];

  sa->local_begin = SC3_MAX (p3->gtroffset[tree], p3->goffset[p3->mpirank]);
  sa->local_end = SC3_MIN (p3->gtroffset[tree + 1], p3->goffset[p3->mpirank + 1]);

  SC3E (sc3_array_index (sa->finfo->sides, 0, &fside));
  for (int side = 0; side < 2; ++side) {
    fside[side].ntree = tree;
    fside[side].is_ghost = -1;
  }

  /* Find here the range of procs that share the same tree as noderank */
  /* Since we process inner faces here:
     The procs < mpirank are responsible for 0th sides,
     the procs > mpirank are responsoble for 1st sides. Proof?*/
  for (side = 0; side < 2; ++side) {
    sa->remote_first[side] = -1;
    sa->remote_last[side] = -2;
  }
  if (tree == p3->fltree && p3->gtroffset[tree] != p3->goffset[p3->mpirank]) {
    SC3E (p4est3_find_partition
          (p3->alloc, p3->mpisize, p3->goffset, p3->gtroffset[tree],
           p3->goffset[p3->mpirank] - 1,
           &sa->remote_first[0], &sa->remote_last[0]));
    if (p3->goffset[sa->remote_first[0]] > p3->gtroffset[tree]) {
      sa->remote_first[0]--;
    }
    if (p3->goffset[sa->remote_last[0]] > p3->goffset[p3->mpirank] - 1) {
      sa->remote_last[0]--;
    }
    SC3A_CHECK (sa->remote_last[0] < p3->mpirank);
  }
  if (tree == p3->lltree
      && p3->gtroffset[tree + 1] != p3->goffset[p3->mpirank + 1]) {
    SC3E (p4est3_find_partition
          (p3->alloc, p3->mpisize, p3->goffset, p3->goffset[p3->mpirank + 1],
           p3->gtroffset[tree + 1] - 1,
           &sa->remote_first[1], &sa->remote_last[1]));
    if (p3->goffset[sa->remote_first[1]] > p3->goffset[p3->mpirank + 1]) {
      sa->remote_first[1]--;
    }
    if (p3->goffset[sa->remote_last[1]] > p3->gtroffset[tree + 1] - 1) {
      sa->remote_last[1]--;
    }
    SC3A_CHECK (sa->remote_first[1] > p3->mpirank);
  }

  return NULL;
}

static sc3_error_t *
p4est3_iterate_volume_rec (p4est3_t * p3,
                           p4est3_iterate_volume_t cvolume,
                           p4est3_iterate_face_t cface,
                           p4est3_iterate_codim_t ccodim,
                           p4est3_search_area_t * sa)
{
  int                 i;
  void               *first_quad;       /*first quadrant in this search area */
  int                 level;
  void               *stack_it;
  p4est3_gloidx      *arr_it;
  p4est3_gloidx       proc_owner;

  const int           max_children = p3->num_children;
  int                *l2nch = sa->level2nchildren;
  int                *Level = &sa->Level;
  p4est3_gloidx      *begin, *end;
  p4est3_tree_t      *tree = sa->tree;
  sc3_array_t        *view_q = sa->view_quads;
  sc3_array_t        *idx_vol_stack = sa->idx_vol_stack;
  p4est3_iterate_volume_info_t *vinfo = sa->vinfo;

  SC3E (sc3_array_index (idx_vol_stack, *Level, &stack_it));
  SC3E (sc3_array_index (*(sc3_array_t **) stack_it, sa->child_id, &begin));
  SC3E (sc3_array_index (*(sc3_array_t **) stack_it, sa->child_id + 1, &end));

  /* Check if the considered search area intersect
     the area of the local process. If not, then skip it. */
  if (*begin >= sa->local_end || *end <= sa->local_begin) {
    l2nch[*Level]++;
    return NULL;
  }

  /* TODO: make it through p4est3_find_partition and last_goffsets */
  /* Find process that owns this quadrant */
  /* Start binary search from a local process, we do this for volume, too,
     because we need to fill the metadata to proceed the recursion. */
  proc_owner = p3->mpirank;
  SC3E (p4est3_search_lower_bound64
         (*begin, p3->goffset, p3->mpisize + 1, &proc_owner));
  if (p3->goffset[proc_owner] > *begin) {
    SC3A_CHECK (proc_owner > 0);
    proc_owner--;
  }
  SC3A_CHECK (proc_owner >= 0 && proc_owner < p3->mpisize);

  /* Get quadrant from the owning process's window */
  first_quad = (void *) (p3->nodequads[proc_owner] +
                        p3->qsize * (*begin - p3->goffset[proc_owner]));

  SC3E (p4est3_quadrant_level (p3->qvt, first_quad, &level));
  if (level == *Level) {
    if (cvolume != NULL) {
      /* Check that this state is never reached by a remote process */
      SC3A_CHECK (proc_owner == p3->mpirank);

      vinfo->quadrant = first_quad;
      /* vinfo->nquad does not really make sense, since it is the quadrant index,
         calculated from the beginning of the tree, counting the quadrants from
         the other processes. Moreover, it cannot be local, since it can be the
         othe processes sharing the same tree. So far leave it as is, since 
         the tests are tailored for this value. */
         /** TODO: Maybe fix later. Possibly either by making it global,
          * or by making it the id/count of the local quadrant within the tree,
          * counting only local quadrant in this tree before it.
          */
      /* vinfo->nquad =  THE PREVIOUS VERSION
        *begin + p3->goffset[p3->mpirank] - p3->gtroffset[tree->treeid]; */
      vinfo->nquad = *begin - p3->gtroffset[tree->treeid];
      SC3E (cvolume (vinfo));
    }
    l2nch[*Level]++;
    return NULL;
  }

  SC3E (sc3_array_push (idx_vol_stack, &stack_it));
#ifdef P4EST_ENABLE_DEBUG
  SC3E (p4est3_array_set_zero (*(sc3_array_t **) stack_it));
#endif
  /* here we start with the very beginning of not necessary local node
     quadrants, because of our specialized array_split_noncontig function */
  SC3E (sc3_array_renew_data
        (&view_q, p3->nodequads[0], p3->qsize, 0, *end - *begin));

  SC3E (p4est3_quadrant_array_split_noncontig
        (p3, view_q, *Level, *begin, *(sc3_array_t **) stack_it));
  l2nch[++(*Level)] = 0;

  /* since array_split doesn't count shift from the beinning of quadrants
     in a node, we shift result indices at the loop below */
  for (i = 0; i < max_children + 1; ++i) {
    SC3E (sc3_array_index (*(sc3_array_t **) stack_it, i, &arr_it));
    *arr_it += *begin;
  }
  for (i = 0; i < max_children; ++i) {
    sa->child_id = i;
    SC3E (p4est3_iterate_volume_rec (p3, cvolume, cface, ccodim, sa));
  }
  SC3A_CHECK (l2nch[*Level] == max_children);
  SC3E (p4est3_iterate_face_inner (p3, cface, ccodim, sa));
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
  int                 face, is_iterate = 0;

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
    search_area->treeid_face[0] = search_area->tree->treeid;
    for (face = 0; face < search_area->nfaces; ++face) {
      SC3E (p4est3_iterate_face_bound_init
            (p3, search_area, tree, face, fside, &is_iterate));
      if (is_iterate) {
        /* we iterate over such trees that tree_neighbor < tree */
        SC3E (p4est3_internal_iterate_face (p3, cface, ccodim, search_area));
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
