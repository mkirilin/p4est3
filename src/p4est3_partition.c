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

#include <p4est3_partition.h>
#include <p4est3_internal.h>
#include <p4est3_search.h>

#ifdef __cplusplus
extern              "C"
{
#if 0
}
#endif
#endif

static sc3_error_t *
p4est3_part_array_new (sc3_allocator_t * alloc, size_t esize, int ealloc,
                       int ecount, sc3_array_t ** arr)
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

static sc3_error_t *
p4est3_partition_allocations (const p4est3_t * p3,
                              char ***precv_buf, char ***psend_buf,
                              p4est3_locidx ** pnum_recv_from,
                              p4est3_locidx ** pnum_send_to,
                              p4est3_locidx ** pnum_per_tree_local,
                              p4est3_locidx ** pnew_local_tree_elem_count,
                              p4est3_gloidx ** plast_goffsets,
                              p4est3_gloidx ** pbegin_send_to,
                              p4est3_gloidx ** pnew_last_goffsets)
{
  const p4est3_topidx num_send_trees
    = p3->gftree[p3->mpirank + 1] - p3->gftree[p3->mpirank] + 1;
  char              **char_pprt;
  p4est3_locidx      *locidx_prt;
  p4est3_gloidx      *gloidx_prt;
  p4est3_topidx       total_num_trees;

  SC3E (p4est3_connectivity_get_num_trees (p3->conn, &total_num_trees));

  SC3E_RETVAL (&char_pprt, NULL);
  SC3E (sc3_allocator_malloc
        (p3->alloc, p3->mpisize * sizeof (char *), &char_pprt));
  *precv_buf = char_pprt;

  SC3E_RETVAL (&char_pprt, NULL);
  SC3E (sc3_allocator_malloc
        (p3->alloc, p3->mpisize * sizeof (char *), &char_pprt));
  *psend_buf = char_pprt;

  SC3E_RETVAL (&locidx_prt, NULL);
  SC3E (sc3_allocator_calloc
        (p3->alloc, p3->mpisize, sizeof (p4est3_locidx), &locidx_prt));
  *pnum_recv_from = locidx_prt;

  SC3E_RETVAL (&locidx_prt, NULL);
  SC3E (sc3_allocator_calloc
        (p3->alloc, p3->mpisize, sizeof (p4est3_locidx), &locidx_prt));
  *pnum_send_to = locidx_prt;

  SC3E_RETVAL (&locidx_prt, NULL);
  SC3E (sc3_allocator_malloc
        (p3->alloc, num_send_trees * sizeof (p4est3_locidx)
          + 2 * sizeof (p4est3_gloidx), &locidx_prt));
  *pnum_per_tree_local = locidx_prt;

  SC3E_RETVAL (&locidx_prt, NULL);
  SC3E (sc3_allocator_malloc
        (p3->alloc, total_num_trees * sizeof (p4est3_locidx), &locidx_prt));
  *pnew_local_tree_elem_count = locidx_prt;

  SC3E_RETVAL (&gloidx_prt, NULL);
  SC3E (sc3_allocator_malloc
        (p3->alloc, p3->mpisize * sizeof (p4est3_gloidx), &gloidx_prt));
  *plast_goffsets = gloidx_prt;

  SC3E_RETVAL (&gloidx_prt, NULL);
  SC3E (sc3_allocator_malloc
        (p3->alloc, p3->mpisize * sizeof (p4est3_gloidx), &gloidx_prt));
#ifdef P4EST_ENABLE_DEBUG
  /* For checking the window of relevant processes. */
  memset (gloidx_prt, -1, p3->mpisize * sizeof (p4est3_gloidx));
#endif
  *pbegin_send_to = gloidx_prt;

  SC3E_RETVAL (&gloidx_prt, NULL);
  SC3E (sc3_allocator_malloc
        (p3->alloc, p3->mpisize * sizeof (p4est3_gloidx), &gloidx_prt));
  *pnew_last_goffsets = gloidx_prt;

  return NULL;
}

/* p4est3::goffsets should be updated before calling this function */
static sc3_error_t *
p4est3_procs_recv_from (const p4est3_t * p3,
                        p4est3_gloidx * last_goffsets,
                        p4est3_locidx * num_recv_from,
                        p4est3_gloidx * from_begin, p4est3_gloidx * from_end,
                        int *num_proc_recv_from)
{
  int                 from_proc;
  p4est3_gloidx       my_begin, my_end, lower_bound;

  my_begin = p3->goffset[p3->mpirank];
  my_end = p3->goffset[p3->mpirank + 1] - (p4est3_gloidx) 1;
  *num_proc_recv_from = 0;

  if (my_begin > my_end) {
    /* my_begin == my_end requires a search is legal for find_partition */
    *from_begin = p3->mpirank;
    *from_end = p3->mpirank;
  }
  else {
    SC3E (p4est3_find_partition
          (p3->alloc, p3->mpisize, last_goffsets,
           my_begin, my_end, from_begin, from_end));
    for (from_proc = *from_begin; from_proc <= *from_end; ++from_proc) {
      lower_bound = p3->old->goffset[from_proc];

      if ((lower_bound > my_end) || (last_goffsets[from_proc] < my_begin)) {
        continue;
      }
      /*TODO: make num_recv_from true/false only, sice we do not send quads */
      num_recv_from[from_proc] = SC3_MIN (my_end, last_goffsets[from_proc])
        - SC3_MAX (my_begin, lower_bound) + 1;
      SC3A_CHECK (num_recv_from[from_proc] >= 0);
      if (from_proc == p3->mpirank) {
        continue;
      }
      (*num_proc_recv_from)++;
    }
  }
  return NULL;
}

/* p4est3::goffsets should be updated before calling this function */
static sc3_error_t *
p4est3_procs_send_to (const p4est3_t * p3,
                      p4est3_gloidx * new_last_goffsets,
                      p4est3_locidx * num_send_to,
                      p4est3_gloidx * begin_send_to, p4est3_gloidx * to_begin,
                      p4est3_gloidx * to_end, int *num_proc_send_to)
{
  int                 to_proc;
  p4est3_gloidx       my_begin, my_end, lower_bound;

  my_begin = p3->old->goffset[p3->mpirank];
  my_end = p3->old->goffset[p3->mpirank + 1] - (p4est3_gloidx) 1;
  *num_proc_send_to = 0;

  if (my_begin > my_end) {
    /* my_begin == my_end requires a search is legal for find_partition */
    *to_begin = p3->mpirank;
    *to_end = p3->mpirank;
    memset (begin_send_to, -1, p3->mpisize * sizeof (p4est3_gloidx));
  }
  else {
    p4est3_find_partition (p3->alloc, p3->mpisize, new_last_goffsets,
                           my_begin, my_end, to_begin, to_end);
    for (to_proc = *to_begin; to_proc <= *to_end; ++to_proc) {
      /* I send to to_proc which may be empty */
      lower_bound = p3->goffset[to_proc];

      if ((lower_bound > my_end) || (new_last_goffsets[to_proc] < my_begin)) {
        continue;
      }
      /*TODO: make num_send_to true/false only, sice we do not send quads */
      num_send_to[to_proc] = SC3_MIN (my_end, new_last_goffsets[to_proc])
        - SC3_MAX (my_begin, lower_bound) + 1;
      begin_send_to[to_proc] = SC3_MAX (my_begin, lower_bound);
      SC3A_CHECK (num_send_to[to_proc] >= 0);
      if (to_proc == p3->mpisize) {
        continue;
      }
      (*num_proc_send_to)++;
    }
  }
  return NULL;
}

static sc3_error_t *
p4est3_trees_send_to (const p4est3_gloidx * begin_send_to,
                      const p4est3_locidx * num_send_to, p4est3_t * p3,
                      p4est3_topidx num_send_trees, int to_proc,
                      p4est3_locidx * num_per_tree_send_buf)
{
  const p4est3_gloidx my_base = p3->old->goffset[p3->mpirank];
  const p4est3_gloidx my_end
    = begin_send_to[to_proc] + (p4est3_gloidx) num_send_to[to_proc]
    - (p4est3_gloidx) 1 - my_base;
  p4est3_gloidx       my_begin = begin_send_to[to_proc] - my_base;
  p4est3_locidx       num_copy;
  p4est3_gloidx       from_begin, from_end, tree_from_begin, tree_from_end,
    num_copy_global;
  p4est3_topidx       which_tree;
  p4est3_tree_t      *tree;
  p4est3_gloidx      *first_last_tquad
    = (p4est3_gloidx *) (&num_per_tree_send_buf[num_send_trees]);
  first_last_tquad[0] = -1;
  first_last_tquad[1] = -1;

  /* Pack in the data to be sent */
  for (which_tree = p3->fltree; which_tree <= p3->lltree; ++which_tree) {
    SC3E (p4est3_tree_index (p3, which_tree, &tree));
    SC3A_CHECK (tree->num_quads >= 1);
    from_begin = tree->quad_offset;
    from_end = tree->quad_offset + tree->num_quads - 1;

    if (from_begin > my_end || from_end < my_begin) {
      continue;
    }
    /* Need to copy from tree which_tree */
    tree_from_begin = SC3_MAX (my_begin, from_begin) - from_begin;
    tree_from_end = SC3_MIN (my_end, from_end) - from_begin;
    num_copy_global = tree_from_end - tree_from_begin + 1;
    SC3A_CHECK (num_copy_global >= 0);
    SC3A_CHECK (num_copy_global <= (p4est3_gloidx) P4EST3_LOCIDX_MAX);
    num_copy = (p4est3_locidx) num_copy_global;
    num_per_tree_send_buf[which_tree - p3->fltree] = num_copy;
    if (num_copy > 0) {
      if (first_last_tquad[0] == -1) {
        first_last_tquad[0] = tree_from_begin;
      }
      first_last_tquad[1] = tree_from_end - 1;
    }

    if (to_proc == p3->mpirank) {
      continue;
    }
    /* move pointer to beginning of quads that need to be copied */
    my_begin += num_copy;
  }
  return NULL;
}

static sc3_error_t *
p4est3_trees_new_boundaries (const p4est3_t * p3,
                             const p4est3_locidx * num_recv_from,
                             char **recv_buf,
                             p4est3_locidx * num_per_tree_local,
                             p4est3_gloidx from_begin_global_quad,
                             p4est3_gloidx from_end_global_quad,
                             p4est3_locidx * new_local_tree_elem_count,
                             p4est3_gloidx * new_first_tquad,
                             p4est3_gloidx * new_last_tquad,
                             p4est3_topidx * new_first_local_tree,
                             p4est3_topidx * new_last_local_tree)
{
  int                 from_proc;
  p4est3_locidx      *num_per_tree_recv_buf;
  p4est3_topidx       first_from_tree, last_from_tree,
    num_recv_trees, it, from_tree;

#ifdef P4EST_ENABLE_DEBUG
  p4est3_topidx       total_num_trees;
  SC3E (p4est3_connectivity_get_num_trees (p3->conn, &total_num_trees));
#endif

  SC3E_RETVAL (new_first_tquad, -1);
  SC3E_RETVAL (new_last_tquad, -1);
  SC3E_RETVAL (new_first_local_tree, (p4est3_topidx) P4EST3_TOPIDX_MAX);
  SC3E_RETVAL (new_last_local_tree, 0);

  for (from_proc = from_begin_global_quad; from_proc <= from_end_global_quad;
       ++from_proc) {
    SC3A_CHECK (num_recv_from[from_proc] >= 0);
    if (num_recv_from[from_proc] == 0) {
      continue;
    }
    first_from_tree = p3->gftree[from_proc];
    last_from_tree = p3->gftree[from_proc + 1];
    num_recv_trees = last_from_tree - first_from_tree + (p4est3_topidx) 1;

    num_per_tree_recv_buf = (from_proc == p3->mpirank) ?
      num_per_tree_local : (p4est3_locidx *) recv_buf[from_proc];

    for (it = 0; it < num_recv_trees; ++it) {

      if (num_per_tree_recv_buf[it] <= 0) {
        continue;
      }
      from_tree = first_from_tree + it;

      SC3A_CHECK (from_tree >= 0 && from_tree < total_num_trees);
      /*P4EST_LDEBUGF ("partition recv %lld [%lld,%lld] quadrants"
         " from tree %lld from proc %d\n",
         (long long) num_per_tree_recv_buf[it],
         (long long) new_local_tree_elem_count[from_tree],
         (long long) new_local_tree_elem_count[from_tree]
         + num_per_tree_recv_buf[it], (long long) from_tree,
         from_proc); */
      *new_first_local_tree = SC3_MIN (*new_first_local_tree, from_tree);
      *new_last_local_tree = SC3_MAX (*new_last_local_tree, from_tree);
      new_local_tree_elem_count[from_tree] += num_per_tree_recv_buf[it];
    }
    if (*new_first_tquad == -1) {
      *new_first_tquad = num_per_tree_recv_buf[num_recv_trees];
    }
    *new_last_tquad = num_per_tree_recv_buf[num_recv_trees + 1];
  }
  if (*new_first_local_tree > *new_last_local_tree) {
    *new_first_local_tree = -1;
    *new_last_local_tree = -2;
  }
  /*P4EST_VERBOSEF ("partition new forest [%lld,%lld]\n",
     (long long) *new_first_local_tree,
     (long long) *new_last_local_tree); */
  return NULL;
}

static sc3_error_t *
p4est3_trees_local_reproduce (const p4est3_locidx * num_send_to,
                              const p4est3_locidx * new_local_tree_elem_count,
                              p4est3_t * p3,
                              p4est3_topidx new_first_local_tree,
                              p4est3_topidx new_last_local_tree,
                              p4est3_gloidx new_first_tquad,
                              p4est3_gloidx new_last_tquad)
{
  p4est3_topidx       tt;
  sc3_array_t        *trees;
  p4est3_tree_t      *tree, *prev_tree;

  p3->fltree = new_first_local_tree;
  p3->lltree = new_last_local_tree;
  p3->nltrees = p3->lltree - p3->fltree + 1;
  SC3E (p4est3_part_array_new
        (p3->alloc, sizeof (p4est3_tree_t), p3->nltrees, p3->nltrees,
         &trees));
  SC3E (sc3_array_unref (&p3->trees));
  p3->trees = trees;

  if (p3->nltrees == 0) {
    return NULL;
  }
  /* process the first local tree */
  SC3E (p4est3_tree_index (p3, p3->fltree, &tree));
  tree->treeid = p3->fltree;
  tree->num_quads = new_local_tree_elem_count[p3->fltree];
  tree->first_tquad = new_first_tquad;
  tree->end_tquad = new_first_tquad + tree->num_quads;
  tree->last_tquad = tree->end_tquad - 1;
  tree->quad_offset = 0;
  tree->tquads = p3->quads;
  /* processing middle trees */
  for (tt = p3->fltree + 1, prev_tree = tree; tt < p3->lltree;
       ++tt, prev_tree = tree) {
    SC3E (p4est3_tree_index (p3, tt, &tree));
    tree->treeid = p3->fltree + tt;
    tree->num_quads = new_local_tree_elem_count[p3->fltree + tt];
    tree->first_tquad = 0;
    tree->end_tquad = tree->num_quads;
    tree->last_tquad = tree->end_tquad - 1;
    tree->quad_offset = prev_tree->quad_offset + prev_tree->num_quads;
    tree->tquads = p3->quads + p3->qsize * tree->quad_offset;
  }
  /* processing the last local tree */
  SC3E (p4est3_tree_index (p3, p3->lltree, &tree));
  tree->treeid = p3->lltree;
  tree->num_quads = new_local_tree_elem_count[p3->lltree];
  tree->end_tquad = new_last_tquad;
  tree->last_tquad = new_last_tquad - 1;
  tree->first_tquad = tree->end_tquad - tree->num_quads;
  tree->quad_offset = prev_tree->quad_offset + prev_tree->quad_offset;
  tree->tquads = p3->quads + p3->qsize * tree->quad_offset;

  return NULL;
}

static sc3_error_t *
p4est3_partition_cleanup (const p4est3_t * p3,
                          int num_proc_recv_from, int num_proc_send_to,
                          p4est3_locidx from_begin_global_quad,
                          p4est3_locidx from_end_global_quad,
                          p4est3_locidx to_begin_global_quad,
                          p4est3_locidx to_end_global_quad,
                          MPI_Request * recv_request,
                          MPI_Request * send_request,
                          char **recv_buf, char **send_buf,
                          p4est3_locidx * num_recv_from,
                          p4est3_locidx * num_send_to,
                          p4est3_locidx * num_per_tree_local,
                          p4est3_locidx * new_local_tree_elem_count,
                          p4est3_gloidx * last_goffsets,
                          p4est3_gloidx * begin_send_to,
                          p4est3_gloidx * new_last_goffsets)
{
  int i;
#ifdef P4EST_ENABLE_MPI
#ifdef P4EST_ENABLE_DEBUG
    for (i = 0; i < num_proc_recv_from; ++i) {
      SC3A_CHECK (recv_request[i] == MPI_REQUEST_NULL);
    }
    for (i = 0; i < num_proc_send_to; ++i) {
      SC3A_CHECK (send_request[i] == MPI_REQUEST_NULL);
    }
#endif
    SC3E (sc3_allocator_free (p3->alloc, recv_request));
    SC3E (sc3_allocator_free (p3->alloc, send_request));
#endif

  for (i = from_begin_global_quad; i <= from_end_global_quad; ++i) {
    if (i != p3->mpirank && num_recv_from[i])
      SC3E (sc3_allocator_free (p3->alloc, recv_buf[i]));
  }
  SC3E (sc3_allocator_free (p3->alloc, recv_buf));

  for (i = to_begin_global_quad; i <= to_end_global_quad; ++i) {
    if (i != p3->mpirank && num_send_to[i])
      SC3E (sc3_allocator_free (p3->alloc, send_buf[i]));
  }
  SC3E (sc3_allocator_free (p3->alloc, send_buf));

  SC3E (sc3_allocator_free (p3->alloc, num_recv_from));
  SC3E (sc3_allocator_free (p3->alloc, num_send_to));
  SC3E (sc3_allocator_free (p3->alloc, num_per_tree_local));
  SC3E (sc3_allocator_free (p3->alloc, new_local_tree_elem_count));
  SC3E (sc3_allocator_free (p3->alloc, last_goffsets));
  SC3E (sc3_allocator_free (p3->alloc, begin_send_to));
  SC3E (sc3_allocator_free (p3->alloc, new_last_goffsets));

  return NULL;
}

static sc3_error_t *
p4est3_weighted_new_boundaries (p4est3_t * p3, int nodesize,
                                int node_offset, int noderank,
                                sc3_MPI_Comm_t nodecomm)
{
  int i;
  int64_t weight, weight_sum, cut;
  int64_t *local_weights, *global_weight_sums;
  ssize_t new_left_border;
  char *quad;
  p4est3_locidx kl = 0, lz;
  p4est3_topidx nt;
  p4est3_tree_t *tree;
  p4est3_quadrant_weight_info_t swi, *wi = &swi;
  wi->p3 = p3;
  wi->qvt = p3->qvt;
  wi->user_data = NULL; /**< Keep user_data NULL so far, may be we will use it in the future. */
  SC3E (sc3_allocator_malloc
        (p3->alloc, sizeof (int64_t) * p3->local_num_quads, &local_weights));
  SC3E (sc3_allocator_malloc
        (p3->alloc, sizeof (int64_t) * (nodesize + 1), &global_weight_sums));

  /* linearly sum weights across all trees */
  local_weights[0] = 0;
  for (nt = p3->fltree; nt <= p3->lltree; ++nt) {
    wi->ntree = nt;
    SC3E (p4est3_tree_index (p3, nt, &tree));
    for (lz = 0; lz < tree->num_quads; ++lz, ++kl) {
      quad = tree->tquads + p3->qsize * lz;
      wi->quadrant = quad;
      SC3E (p3->cweight (wi, &weight));
      SC3A_CHECK (weight >= 0);
      local_weights[kl + 1] = local_weights[kl] + weight;
    }
  }

  SC3A_CHECK (kl == p3->local_num_quads);
  weight_sum = local_weights[p3->local_num_quads];

  /* distribute local weight sums */
  global_weight_sums[0] = 0;
  SC3E (sc3_MPI_Allgather (&weight_sum, 1, MPI_LONG_LONG_INT,
                           &global_weight_sums[1], 1, MPI_LONG_LONG_INT,
                           nodecomm));
  /* adjust all arrays to reflect the global weight */
  for (i = 0; i < nodesize; ++i) {
    global_weight_sums[i + 1] += global_weight_sums[i];
  }
  if (noderank > 0) {
    weight_sum = global_weight_sums[noderank];
    for (kl = 0; kl <= p3->local_num_quads; ++kl) {
      local_weights[kl] += weight_sum;
    }
  }

  SC3A_CHECK (local_weights[0] == global_weight_sums[noderank]);
  SC3A_CHECK (local_weights[p3->local_num_quads] ==
                global_weight_sums[noderank + 1]);
  weight_sum = global_weight_sums[nodesize];

  /* if all quadrants have zero weight we do nothing */
  if (weight_sum == 0) {
    SC3E (sc3_allocator_free (p3->alloc, &local_weights));
    SC3E (sc3_allocator_free (p3->alloc, &global_weight_sums));
    return NULL;
  }

  /* determine processor ids to edit their boundary and just do it */
  for (i = 0; i < nodesize; ++i) {
    cut = p4est3_uint64cut (weight_sum, i, nodesize);
    if (global_weight_sums[noderank] > cut ||
        cut >= global_weight_sums[noderank + 1]) {
      continue;
    }
    SC3E (p4est3_search_lower_bound64
          (cut, local_weights, (ssize_t) p3->local_num_quads, &new_left_border));
    SC3A_CHECK (new_left_border > 0
                && (p4est3_locidx) new_left_border <= p3->local_num_quads);
    p3->goffset[node_offset + i] = new_left_border;
  }
  return NULL;
}

sc3_error_t        *
p4est3_partition (p4est3_t * p3)
{
  /* at this stage we have a completely setup refined forest */
  const p4est3_topidx num_send_trees
    = p3->gftree[p3->mpirank + 1] - p3->gftree[p3->mpirank] + 1;
  /* We are going to send p4est3_tree::num_quads for each sent tree,
     p4est3_tree::first_tquad for the first sent tree and
     p4est3_tree::last_tquad for the last sent tree. */
  const size_t        send_size = num_send_trees * sizeof (p4est3_locidx)
    + 2 * sizeof (p4est3_gloidx);
  int                 i, from_proc, sk, mpiret, to_proc;
  int                 nodesize, noderank, node_num, node_offset;
  const int          *node_offsets;
  int                 num_proc_recv_from, num_proc_send_to;
  char              **recv_buf, **send_buf;
  p4est3_topidx       num_recv_trees;
  p4est3_topidx       new_first_local_tree, new_last_local_tree;
  p4est3_locidx       from_begin_global_quad, from_end_global_quad;
  p4est3_locidx       to_begin_global_quad, to_end_global_quad;
  p4est3_locidx      *num_per_tree_send_buf, *new_local_tree_elem_count;
  p4est3_locidx      *num_send_to;
  p4est3_locidx      *num_recv_from;
                                /**< Numbers of quadrants coming from the i-th process */
  p4est3_locidx      *num_per_tree_local;
  p4est3_gloidx       from_begin, from_end, to_begin, to_end;
  p4est3_gloidx      *begin_send_to;
                                /**< Begin (quadrant) of proc where we want to send to
                                     or begin of current proc (MAX of them). From which quad we start
                                     sending info. */
  p4est3_gloidx       qcount_node;
  p4est3_gloidx       new_first_tquad, new_last_tquad;
  p4est3_gloidx      *last_goffsets, *new_last_goffsets;
                                                    /**< Offsets of last quadrant in a process */
  size_t              recv_size;
  sc3_MPI_Comm_t      nodecomm;
  MPI_Request        *recv_request, *send_request;

  /* We suppose to call this function after setting up routine */
  SC3A_CHECK (p3->old != NULL);
  SC3A_IS (p4est3_is_setup, p3->old);

  /* this function does nothing for processes without shared memory */
  /** TODO: Make sence to limit it for nodesize == 1, w/o sh.mem it leads to extra work. 
   * Keep it so far for testing purposes.
  */
  /*if (nodesize == 1) {
      return NULL;
    }*/
  if (p3->cweight == NULL) {
    /* Divide up the quadrants equally */
    SC3E (sc3_mpienv_get_nodesize (p3->split_info, &nodesize));
    SC3E (sc3_mpienv_get_node_num (p3->split_info, &node_num));
    SC3E (sc3_mpienv_get_node_offsets (p3->split_info, &node_offsets));
    SC3E (sc3_mpienv_get_noderank (p3->split_info, &noderank));
    SC3E (sc3_mpienv_get_nodecomm (p3->split_info, &nodecomm));

  node_offset = node_offsets[node_num];
  /* Adjust quadrant partition information of the forest */
  SC3E (sc3_MPI_Win_lock (SC3_MPI_LOCK_SHARED, 0, SC3_MPI_MODE_NOCHECK,
                          p3->goffsets->goffsetwin));

  if (p3->cweight == NULL) {
    qcount_node = p3->goffset[node_offsets[node_num + 1]]
      - p3->goffset[node_offset];
#ifdef P4EST_ENABLE_DEBUG
      if (noderank == 0) {
        SC3A_CHECK (p3->goffset[p3->mpirank] ==
                    p4est3_glocut (qcount_node, nodesize, noderank));
      }
      SC3E (sc3_MPI_Barrier (nodecomm));
#endif
      /** TODO: figure out mpisize vs nodesize stuff */
      p3->goffset[p3->mpirank] = p4est3_glocut (qcount_node, nodesize, noderank);
    }
    else {
      SC3E (p4est3_weighted_new_boundaries
            (p3, nodesize, node_offset, noderank, nodecomm));
    }
    SC3E (sc3_MPI_Win_sync (p3->goffsets->goffsetwin));
    SC3E (sc3_MPI_Win_unlock (0, p3->goffsets->goffsetwin));
    SC3E (sc3_MPI_Barrier (nodecomm));

    p3->local_num_quads = p3->goffset[p3->mpirank + 1] - p3->goffset[p3->mpirank];

    SC3E (sc3_MPI_Win_lock (SC3_MPI_LOCK_SHARED, 0, SC3_MPI_MODE_NOCHECK,
                            p3->gposition->gfposwin));
    SC3E (p4est3_quadrant_first_descendant
          (p3->qvt, p3->nodequads[0] + p3->qsize * p3->goffset[p3->mpirank],
            p3->qmaxlevel, p3->gfpos + p3->mpirank * p3->qsize));

    if (noderank == 0) {
      SC3E (p4est3_quadrant_translate
            (p3->old->qvt,
             p3->old->gfpos + node_offsets[node_num + 1] * p3->old->qsize,
             p3->qvt, p3->gfpos + node_offsets[node_num + 1] * p3->qsize));
    }
    else {
      p3->nodequads[noderank] =
        p3->nodequads[0] + p3->qsize * p3->goffset[p3->mpirank];
    }
    SC3E (sc3_MPI_Win_unlock (0, p3->gposition->gfposwin));

    p3->quads = p3->nodequads[noderank];

    /* Adjust trees partition information of the forest */
    /** TODO: Put some of the allocations after isend/irecv */
    SC3E (p4est3_partition_allocations
          (p3, &recv_buf, &send_buf, &num_recv_from, &num_send_to,
           &num_per_tree_local, &new_local_tree_elem_count, &last_goffsets,
           &begin_send_to, &new_last_goffsets));

    for (i = 0; i < p3->mpisize; ++i) {
      last_goffsets[i] = p3->old->goffset[i + 1] - 1;
    }
    SC3E (p4est3_procs_recv_from
          (p3, last_goffsets, num_recv_from,
           &from_begin, &from_end, &num_proc_recv_from));

    from_begin_global_quad = from_begin;
    from_end_global_quad = from_end;
    /* Post receives for the trees */
#ifdef P4EST_ENABLE_MPI
    SC3E (sc3_allocator_malloc (p3->alloc,
                                num_proc_recv_from * sizeof (MPI_Request),
                                &recv_request));
#endif

    /* Allocate space for receiving trees */
    for (from_proc = from_begin_global_quad, sk = 0;
         from_proc <= from_end_global_quad; ++from_proc) {
      if (from_proc != p3->mpirank && num_recv_from[from_proc]) {
        num_recv_trees =
          p3->gftree[from_proc + 1] - p3->gftree[from_proc] + 1;
        /* We are going to recv p4est3_tree::num_quads for each received tree,
           p4est3_tree::first_tquad for the first received tree and
           p4est3_tree::last_tquad for the last received tree. */
        recv_size = num_recv_trees * sizeof (p4est3_locidx)
          + 2 * sizeof (p4est3_gloidx);
        SC3E (sc3_allocator_malloc
              (p3->alloc, recv_size * sizeof (char), &recv_buf[from_proc]));
        /* Post receives for the quadrants and their data */
#ifdef P4EST_ENABLE_MPI
        mpiret = MPI_Irecv (recv_buf[from_proc], (int) recv_size, MPI_BYTE,
                            from_proc, 0, p3->mpicomm, recv_request + sk);
        SC3A_CHECK (mpiret == SC3_MPI_SUCCESS);
        ++sk;
#endif
      }
    }
#ifdef P4EST_ENABLE_MPI
    for (; sk < num_proc_recv_from; ++sk) {
      /* for empty processors in receiving range */
      recv_request[sk] = MPI_REQUEST_NULL;
    }
#endif

    /* For each processor calculate the number of quadrants sent */
    for (i = 0; i < p3->mpisize; ++i) {
      new_last_goffsets[i] = p3->goffset[i + 1] - 1;
    }
    SC3E (p4est3_procs_send_to
          (p3, new_last_goffsets, num_send_to,
           begin_send_to, &to_begin, &to_end, &num_proc_send_to));
    to_begin_global_quad = to_begin;
    to_end_global_quad = to_end;

    /* Communicate the trees */
#ifdef P4EST_ENABLE_MPI
    SC3E (sc3_allocator_malloc (p3->alloc,
                                num_proc_send_to * sizeof (MPI_Request),
                                &send_request));
#endif
    /* Set the num_per_tree_local */
    memset (num_per_tree_local, 0, send_size);
    SC3E (p4est3_trees_send_to
          (begin_send_to, num_send_to, p3,
           num_send_trees, p3->mpirank, num_per_tree_local));

    for (to_proc = to_begin_global_quad, sk = 0;
         to_proc <= to_end_global_quad; ++to_proc) {
      if (to_proc != p3->mpirank && num_send_to[to_proc]) {
        SC3E (sc3_allocator_malloc
              (p3->alloc, send_size, &send_buf[to_proc]));
        num_per_tree_send_buf = (p4est3_locidx *) send_buf[to_proc];
        memset (num_per_tree_send_buf, 0, send_size);
        SC3E (p4est3_trees_send_to
              (begin_send_to, num_send_to, p3,
               num_send_trees, to_proc, num_per_tree_send_buf));
        /* Post send operation for the quadrants and their data */
#ifdef P4EST_ENABLE_MPI
        mpiret = MPI_Isend (send_buf[to_proc], (int) send_size, MPI_BYTE,
                            to_proc, 0, p3->mpicomm, send_request + sk);
        SC3A_CHECK (mpiret == SC3_MPI_SUCCESS);
        ++sk;
#endif
      }
    }
#ifdef P4EST_ENABLE_MPI
    for (; sk < num_proc_send_to; ++sk) {
      send_request[sk] = MPI_REQUEST_NULL;
    }
    /* Fill in forest */
    mpiret =
      MPI_Waitall (num_proc_recv_from, recv_request, MPI_STATUSES_IGNORE);
    SC3A_CHECK (mpiret == SC3_MPI_SUCCESS);
#endif
    /* Calculate the local index of the end of each tree in the repartition */
    SC3E (p4est3_trees_new_boundaries (p3, num_recv_from, recv_buf,
                                       num_per_tree_local,
                                       from_begin_global_quad,
                                       from_end_global_quad,
                                       new_local_tree_elem_count,
                                       &new_first_tquad, &new_last_tquad,
                                       &new_first_local_tree,
                                       &new_last_local_tree));
    /* Restore local trees */
    SC3E (p4est3_trees_local_reproduce
          (num_send_to, new_local_tree_elem_count, p3, new_first_local_tree,
           new_last_local_tree, new_first_tquad, new_last_tquad));

    /* Free allocations */
#ifdef P4EST_ENABLE_MPI
    mpiret = MPI_Waitall (num_proc_send_to, send_request, MPI_STATUSES_IGNORE);
    SC3A_CHECK (mpiret == SC3_MPI_SUCCESS);
#endif
    SC3E (p4est3_partition_cleanup
          (p3, num_proc_recv_from, num_proc_send_to, from_begin_global_quad,
           from_end_global_quad, to_begin_global_quad, to_end_global_quad,
           recv_request, send_request, recv_buf, send_buf, num_recv_from,
           num_send_to, num_per_tree_local, new_local_tree_elem_count,
           last_goffsets, begin_send_to, new_last_goffsets));
  }
  return NULL;
}

#ifdef __cplusplus
#if 0
{
#endif
}
#endif
