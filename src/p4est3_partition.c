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
#include <p4est3_search.h>

#ifdef __cplusplus
extern              "C"
{
#if 0
}
#endif
#endif

/* disable code responsible for sending trees info,
   since we are in shared memory now */
#define SENDRECV_TREES 0

static sc3_error_t *
p4est3_part_array_new (sc3_allocator_t *alloc, size_t esize, int ealloc,
                       int ecount, sc3_array_t **arr)
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

/*static sc3_error_t *
p4est3_partition_allocations_prerecv (const p4est3_t *p3,
                                      p4est3_locidx **pnum_recv_from,
                                      p4est3_gloidx **plast_goffsets,
                                      p4est3_gloidx **plast_gtree_offsets,
                                      p4est3_gloidx **ploc_offsets)
{
  p4est3_locidx      *locidx_prt;
  p4est3_gloidx      *gloidx_prt;

  if (!p3->contiguous) {
    locidx_prt = NULL;
    SC3E (sc3_allocator_calloc
          (p3->alloc, p3->mpisize, sizeof (p4est3_locidx), &locidx_prt));
    *pnum_recv_from = locidx_prt;
  }

  gloidx_prt = NULL;
  SC3E (sc3_allocator_malloc
        (p3->alloc, p3->mpisize * sizeof (p4est3_gloidx), &gloidx_prt));
  *plast_goffsets = gloidx_prt;

  gloidx_prt = NULL;
  SC3E (sc3_allocator_malloc
        (p3->alloc, p3->num_trees * sizeof (p4est3_gloidx), &gloidx_prt));
  *plast_gtree_offsets = gloidx_prt;

  gloidx_prt = NULL;
  SC3E (sc3_allocator_malloc
        (p3->alloc, (p3->mpisize + 1) * sizeof (p4est3_gloidx), &gloidx_prt));
  *ploc_offsets = gloidx_prt;

  return NULL;
}*/

static sc3_error_t *
p4est3_procs_recv_from (const p4est3_t *p3,
                        p4est3_gloidx *last_goffsets,
                        p4est3_locidx *num_recv_from,
                        p4est3_gloidx *from_begin, p4est3_gloidx *from_end,
                        p4est3_gloidx *loc_offsets)
{
  p4est3_gloidx       from_proc;
  p4est3_gloidx       my_begin, my_end, lower_bound;

  SC3A_CHECK (num_recv_from != NULL);

  my_begin = loc_offsets[p3->mpirank];
  my_end = loc_offsets[p3->mpirank + 1] - (p4est3_gloidx) 1;

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
      num_recv_from[from_proc] = SC3_MIN (my_end, last_goffsets[from_proc])
        - SC3_MAX (my_begin, lower_bound) + 1;
      SC3A_CHECK (num_recv_from[from_proc] >= 0);
      if (from_proc == p3->mpirank) {
        continue;
      }
    }
  }
  return NULL;
}

static sc3_error_t *
p4est3_quads_copy_from (const p4est3_t *p3, p4est3_locidx *num_recv_from,
                        p4est3_gloidx from_begin, p4est3_gloidx from_end,
                        p4est3_gloidx *loc_offsets, int noderank)
{
  const p4est3_gloidx my_begin = loc_offsets[p3->mpirank];
  p4est3_gloidx       from_proc;
  p4est3_locidx       from_begin_copy, to_qid = 0, li;
  p4est3_gloidx       lower_bound;

  SC3A_CHECK (num_recv_from != NULL);

  for (from_proc = from_begin; from_proc <= from_end; ++from_proc) {
    if (num_recv_from[from_proc] == 0) {
      continue;
    }
    lower_bound = p3->old->goffset[from_proc];
    from_begin_copy =
      (p4est3_locidx) (SC3_MAX (my_begin, lower_bound) - lower_bound);
    //from_end_copy = from_begin_copy + num_recv_from[from_proc];
    if (p3->qvt == p3->old->qvt) {
      memcpy (p3->quads + (ptrdiff_t) to_qid * (ptrdiff_t) p3->qsize,
              p3->old->nodequads[from_proc] +
              (ptrdiff_t) from_begin_copy * (ptrdiff_t) p3->qsize,
              (ptrdiff_t) num_recv_from[from_proc] * (ptrdiff_t) p3->qsize);
    }
    else {
      for (li = 0; li < num_recv_from[from_proc]; ++li) {
        SC3E
          (p4est3_quadrant_translate
           (p3->old->qvt, p3->old->nodequads[from_proc]
            + (ptrdiff_t) (from_begin_copy + li) * (ptrdiff_t) p3->old->qsize,
            p3->qvt,
            p3->quads + (ptrdiff_t) (to_qid + li) * (ptrdiff_t) p3->qsize));
      }
    }
    to_qid += num_recv_from[from_proc];
    /*for (from_qid = from_begin_copy; from_qid < from_end_copy; ++from_qid) {
       SC3E (p4est3_quadrant_copy
       (p3->qvt,
       p3->old->nodequads[from_proc] + from_qid * p3->qsize,
       new_nodequad + (to_qid++) * p3->qsize));
       } */
  }
  return NULL;
}

/* p4est3::goffsets should be updated before calling this function */
/** Compute correction of partition for a process.
 *
 * The correction denotes how many quadrants the process with id `rank` takes
 * from (if correction positive) or gives to (if correction negative) the
 * previous process with id `rank-1` in order to assign a family of quadrants
 * to one process.
 * The process with the highest number of quadrants of a family gets all
 * quadrants belonging to this family from other processes. If this applies to
 * several processes, then the process with the lowest id gets the quadrants.
 * A process can give more quadrants than it owns, if it passes quadrants from
 * other processes.
 */
static sc3_error_t *
p4est3_partition_correction (const p4est3_t *p3,
                             int node_offset, int node_offset_next,
                             p4est3_gloidx *last_goffsets,
                             p4est3_gloidx *loc_offsets,
                             p4est3_locidx *correction)
{
  int                 level, is_parent;
  int                 i, rk, rank_with_max_quads;
  char               *cut_border_quad;  /* first quadrant in the new uncorrected cut */
  char               *family_parent = p3->temp_quad[0];
  p4est3_locidx       qid, quad_begin, quad_end;
  p4est3_gloidx       my_begin, my_end, lower_bound, from_begin, from_end;
  p4est3_gloidx       from_proc, min_quadrant_id, max_quadrant_id, h;
  p4est3_gloidx       max_num_quadrants;

  for (rk = 0; rk < p3->mpisize; ++rk) {
    rank_with_max_quads = rk;
    /* Determine first and last global ids of a family on the process border */
    /* set minimum possible borders containing the target family inside */
    my_begin = loc_offsets[rk] - (p4est3_gloidx) (p3->num_children - 1);
    my_end = loc_offsets[rk] + (p4est3_gloidx) (p3->num_children - 1);

    /* correct the borders w.r.t shared memory node boundaries */
    my_begin = SC3_MAX (my_begin, loc_offsets[node_offset]);
    my_end = SC3_MIN (my_end, loc_offsets[node_offset_next] - 1);

    /* find to which process the beginning of the (uncorrected) cut belongs */
    SC3E (p4est3_find_partition
          (p3->alloc, p3->mpisize, last_goffsets,
           loc_offsets[rk], loc_offsets[rk], &from_begin, &from_end));
    SC3A_CHECK (from_begin == from_end);
    qid = (p4est3_locidx) (loc_offsets[rk] - p3->old->goffset[from_begin]);
    /* if border quadrant is a tree, do nothing */
    cut_border_quad =
      p3->old->nodequads[from_begin - node_offset] +
      (ptrdiff_t) qid *(ptrdiff_t) p3->qsize;
    SC3E (p4est3_quadrant_level (p3->qvt, cut_border_quad, &level));
    if (level == 0) {
      continue;
    }

    /* find which processes new (but uncorrected) borders belong to */
    SC3E (p4est3_find_partition
          (p3->alloc, p3->mpisize, last_goffsets,
           my_begin, my_end, &from_begin, &from_end));

    min_quadrant_id = max_quadrant_id = -1;
    /* the parent of the target family, we find min and max global ids to */
    SC3E (p4est3_quadrant_parent (p3->qvt, cut_border_quad, family_parent));

    for (from_proc = from_begin; from_proc <= from_end; ++from_proc) {
      lower_bound = p3->old->goffset[from_proc];
      /* go from global ids to local ones w.r.t shared memory node boundaries */
      quad_begin = SC3_MAX (my_begin, lower_bound) - lower_bound;
      quad_end = SC3_MIN (my_end, last_goffsets[from_proc]) - lower_bound;
      for (qid = quad_begin; qid <= quad_end; ++qid) {
        SC3E (p4est3_quadrant_is_parent
              (p3->qvt, family_parent,
               p3->old->nodequads[from_proc - node_offset] +
               (ptrdiff_t) qid * (ptrdiff_t) p3->qsize, &is_parent));
        if (is_parent) {
          max_quadrant_id = (p4est3_gloidx) qid + lower_bound;
          min_quadrant_id =
            min_quadrant_id == -1 ? max_quadrant_id : min_quadrant_id;
        }
      }
    }

    /* Compute correction */
    /* no correction if num quadrants not sufficient for family */
    if (max_quadrant_id - min_quadrant_id + 1 != p3->num_children) {
      continue;
    }

    max_num_quadrants = SC3_MIN (max_quadrant_id, loc_offsets[rk + 1] - 1)
      - loc_offsets[rk] + 1;
    /* decreasing search for process with highest amount of quadrants */
    i = rank_with_max_quads - 1;
    while (min_quadrant_id < loc_offsets[i + 1]) {
      h = loc_offsets[i + 1] - SC3_MAX (min_quadrant_id, loc_offsets[i]);
      if (max_num_quadrants <= h) {
        max_num_quadrants = h;
        rank_with_max_quads = i;
      }
      i--;
    }

    /* increasing search for process with highest amount of quadrants */
    i = rank_with_max_quads + 1;
    while (loc_offsets[i] <= max_quadrant_id) {
      h =
        SC3_MIN (max_quadrant_id,
                 loc_offsets[i + 1] - 1) - loc_offsets[i] + 1;
      if (max_num_quadrants < h) {
        max_num_quadrants = h;
        rank_with_max_quads = i;
      }
      i++;
    }

    /* compute correction */
    if (rank_with_max_quads < rk) {
      correction[rk] =
        (p4est3_locidx) (loc_offsets[rk] - max_quadrant_id - 1);
    }
    else {
      correction[rk] = (p4est3_locidx) (loc_offsets[rk] - min_quadrant_id);
    }
  }
  return NULL;
}

/* A shortcut to allocate shared memory for local quadrants storage */
static sc3_error_t *
p4est3_quadrants_allocate (p4est3_t *p3, int nodesize)
{
  int                 i, dispunit;
  sc3_MPI_Aint_t      tempbytes;

  SC3E (p4est3_quadrants_new (p3->alloc, &p3->quadrants));
  SC3E (p4est3_glopartition_set_mpienv
        (NULL, NULL, NULL, NULL, p3->quadrants, p3->split_info));
  SC3E (p4est3_quadrants_set_local_num_quads
        (p3->quadrants, p3->local_num_quads));
  SC3E (p4est3_quadrants_set_global_alloc_quads
        (p3->quadrants, p3->global_num_quads));
  SC3E (p4est3_quadrants_set_qsize (p3->quadrants, p3->qsize));
  SC3E (p4est3_glopartition_setup (NULL, NULL, NULL, NULL, p3->quadrants));
  p3->quads = p3->quadrants->quads;
  for (i = 0; i < nodesize; ++i) {
    /* TODO: change functionality for non-cont, wrt pre-allocation functional
       of quadrant magic array */
    SC3E (sc3_MPI_Win_shared_query (p3->quadrants->meta->win, i,
                                    &tempbytes, &dispunit,
                                    &p3->nodequads[i]));
    SC3A_CHECK (dispunit == p3->qsize);
    SC3A_CHECK (p3->nodequads[i] != NULL || tempbytes == 0);
  }
  return NULL;
}

static sc3_error_t *
p4est3_weighted_new_boundaries (p4est3_t *p3, int nodesize,
                                int node_offset, int noderank,
                                sc3_MPI_Comm_t nodecomm)
{
  int                 i;
  int64_t             weight, weight_sum, cut;
  int64_t            *local_weights, *global_weight_sums;
  ssize_t             new_left_border = 0;
  char               *quad;
  p4est3_locidx       kl = 0, lz;
  p4est3_gloidx       goffset_current_rank = p3->old->goffset[p3->mpirank];
  p4est3_topidx       nt;
  p4est3_tree_t      *tree;
  p4est3_quadrant_weight_info_t swi, *wi = &swi;
  wi->p3 = p3;
  wi->qvt = p3->qvt;
  wi->user_data = NULL; /**< Keep user_data NULL so far, may be we will use it in the future. */
  SC3E (sc3_allocator_malloc
        (p3->alloc, sizeof (int64_t) * (p3->old->local_num_quads + 1),
         &local_weights));
  SC3E (sc3_allocator_malloc
        (p3->alloc, sizeof (int64_t) * (nodesize + 1), &global_weight_sums));

  /* linearly sum weights across all trees */
  local_weights[0] = 0;
  for (nt = p3->old->fltree; nt <= p3->old->lltree; ++nt) {
    wi->ntree = nt;
    SC3E (p4est3_tree_index (p3->old, nt, &tree));
    for (lz = 0; lz < tree->num_quads; ++lz, ++kl) {
      quad = tree->tquads + (ptrdiff_t) p3->qsize * (ptrdiff_t) lz;
      wi->quadrant = quad;
      SC3E (p3->cweight (wi, &weight));
      SC3A_CHECK (weight >= 0);
      local_weights[kl + 1] = local_weights[kl] + weight;
    }
  }

  SC3A_CHECK (kl == p3->old->local_num_quads);
  weight_sum = local_weights[p3->old->local_num_quads];

  /* distribute local weight sums */
  global_weight_sums[0] = 0;
  SC3E (sc3_MPI_Allgather (&weight_sum, 1, SC3_MPI_LONG_LONG,
                           &global_weight_sums[1], 1, SC3_MPI_LONG_LONG,
                           nodecomm));
  /* adjust all arrays to reflect the global weight */
  for (i = 0; i < nodesize; ++i) {
    global_weight_sums[i + 1] += global_weight_sums[i];
  }
  if (noderank > 0) {
    weight_sum = global_weight_sums[noderank];
    for (kl = 0; kl <= p3->old->local_num_quads; ++kl) {
      local_weights[kl] += weight_sum;
    }
  }

  SC3A_CHECK (local_weights[0] == global_weight_sums[noderank]);
  SC3A_CHECK (local_weights[p3->old->local_num_quads] ==
              global_weight_sums[noderank + 1]);
  weight_sum = global_weight_sums[nodesize];

  /* if all quadrants have zero weight we do nothing */
  if (weight_sum == 0) {
    SC3E (sc3_allocator_free (p3->alloc, local_weights));
    SC3E (sc3_allocator_free (p3->alloc, global_weight_sums));
    return NULL;
  }

  /* determine processor ids to edit their boundary and just do it */
  for (i = 0; i < nodesize; ++i) {
    cut = p4est3_uint64cut (weight_sum, nodesize, i);

    if (global_weight_sums[noderank] > cut ||
        cut >= global_weight_sums[noderank + 1]) {
      continue;
    }
    SC3E (p4est3_search_lower_bound64
          (cut, local_weights, (ssize_t) p3->old->local_num_quads + 1,
           &new_left_border));
    SC3A_CHECK (new_left_border >= 0
                && (p4est3_locidx) new_left_border <
                p3->old->local_num_quads);
    /* shift new left border by goffset since we searched it in local array */
    p3->goffset[node_offset + i] = new_left_border + goffset_current_rank;
  }
  SC3E (sc3_allocator_free (p3->alloc, local_weights));
  SC3E (sc3_allocator_free (p3->alloc, global_weight_sums));
  return NULL;
}

static sc3_error_t *
p4est3_find_first_last_local_trees (p4est3_t *p3,
                                    p4est3_gloidx *last_gtree_offsets,
                                    p4est3_gloidx *loc_offsets)
{
  sc3_array_t        *trees;
  p4est3_gloidx       first_quad_gloid, last_quad_gloid;
  p4est3_gloidx       from_begin, from_end;

  first_quad_gloid = loc_offsets[p3->mpirank];
  last_quad_gloid = loc_offsets[p3->mpirank + 1] - (p4est3_gloidx) 1;

  if (first_quad_gloid > last_quad_gloid) {
    /* we are in empty process */
    p3->fltree = -1;
    p3->lltree = -2;
    p3->nltrees = 0;
    SC3E (p4est3_part_array_new (p3->alloc, sizeof (p4est3_tree_t),
                                 p3->nltrees, p3->nltrees, &trees));
    p3->trees = trees;
    return NULL;
  }

  SC3E (p4est3_find_partition
        (p3->alloc, p3->num_trees, last_gtree_offsets,
         first_quad_gloid, last_quad_gloid, &from_begin, &from_end));

  p3->fltree = (p4est3_topidx) from_begin;
  p3->lltree = (p4est3_topidx) from_end;
  p3->nltrees = p3->lltree - p3->fltree + 1;

  /** TODO: Gather data from the other nodes */
  p3->gftree[p3->mpirank] = p3->fltree;
  if (p3->mpirank == 0) {
    p3->gftree[p3->mpisize] = p3->num_trees;
  }
  return NULL;
}

static sc3_error_t *
p4est3_local_trees_reproduce (p4est3_t *p3,
                              p4est3_gloidx *last_gtree_offsets,
                              p4est3_gloidx *loc_offsets)
{
  sc3_array_t        *trees;
  p4est3_tree_t      *tree, *prev_tree;
  p4est3_gloidx       first_quad_gloid, last_quad_gloid;
  p4est3_topidx       tt;

  first_quad_gloid = loc_offsets[p3->mpirank];
  last_quad_gloid = loc_offsets[p3->mpirank + 1] - (p4est3_gloidx) 1;

  if (first_quad_gloid > last_quad_gloid) {
    /* everything is already done by p4est3_find_first_last_local_trees */
    return NULL;
  }

  SC3E (p4est3_part_array_new
        (p3->alloc, sizeof (p4est3_tree_t), p3->nltrees, p3->nltrees,
         &trees));
  p3->trees = trees;

  SC3A_CHECK (p3->nltrees > 0);
  /* process the first local tree */
  SC3E (sc3_array_index (trees, 0, &tree));
  tree->treeid = p3->fltree;
  tree->first_tquad = first_quad_gloid - p3->gtroffset[p3->fltree];
  if (p3->nltrees == 1) {
    tree->last_tquad = last_quad_gloid - p3->gtroffset[p3->fltree];
  }
  else {
    tree->last_tquad
      = p3->gtroffset[p3->fltree + 1] - p3->gtroffset[p3->fltree] - 1;
  }
  tree->end_tquad = tree->last_tquad + 1;
  tree->num_quads = (p4est3_locidx) (tree->end_tquad - tree->first_tquad);
  tree->quad_offset = 0;
  tree->tquads = p3->quads;

  if (p3->nltrees == 1) {
    return NULL;
  }

  /* process middle local trees */
  for (tt = p3->fltree + 1, prev_tree = tree; tt < p3->lltree;
       ++tt, prev_tree = tree) {
    SC3E (sc3_array_index (trees, tt - p3->fltree, &tree));
    tree->treeid = tt;
    tree->first_tquad = 0;
    tree->num_quads =
      (p4est3_locidx) (p3->gtroffset[tt + 1] - p3->gtroffset[tt]);
    tree->end_tquad = tree->num_quads;
    tree->last_tquad = tree->end_tquad - 1;
    tree->quad_offset = prev_tree->quad_offset + prev_tree->num_quads;
    tree->tquads =
      p3->quads + (ptrdiff_t) p3->qsize * (ptrdiff_t) tree->quad_offset;
  }

  /* process the last local tree */
  SC3E (sc3_array_index (trees, p3->lltree - p3->fltree, &tree));
  tree->treeid = p3->lltree;
  tree->first_tquad = 0;
  tree->last_tquad = last_quad_gloid - p3->gtroffset[p3->lltree];
  tree->end_tquad = tree->last_tquad + 1;
  tree->num_quads = (p4est3_locidx) tree->end_tquad;
  tree->quad_offset = prev_tree->quad_offset + prev_tree->num_quads;
  tree->tquads =
    p3->quads + (ptrdiff_t) p3->qsize * (ptrdiff_t) tree->quad_offset;

  return NULL;
}

sc3_error_t        *
p4est3_partition (p4est3_t *p3)
{
  /* at this stage we have a completely setup refined forest */
  int                 i;
  int                 n, nodesize, noderank, node_num, node_offset,
    node_offset_next;
  const int          *node_offsets;
  p4est3_locidx      *correction;
  p4est3_locidx      *num_recv_from = NULL;
                                /**< Numbers of quadrants coming from the i-th process */
  p4est3_gloidx       from_begin = p3->mpirank, from_end = p3->mpirank;
                                /**< Range of processes to receive quadrants from */
  p4est3_gloidx       qcount_node;
  p4est3_gloidx      *last_goffsets = NULL;
                                  /**< Offsets of last quadrant in a process */
  p4est3_gloidx      *last_gtree_offsets = NULL;
                                  /**< Offsets of last quadrant in a tree */
  sc3_MPI_Comm_t      nodecomm;
  p4est3_gloidx      *loc_offsets = NULL;
  p4est3_locidx       li;
  p4est3_topidx       t;
  MPI_Request         req[2];
  MPI_Status          statuses[2];
  /* We suppose to call this function after setting up routine */
  SC3A_CHECK (p3->old != NULL);
  SC3A_IS (p4est3_is_setup, p3->old);

  if (!p3->shared) {
    /* so far it works for shared memory only */
    return NULL;
  }

  /* this function does nothing for processes without shared memory */
  /** TODO: Make sence to limit it for nodesize == 1, w/o sh.mem it leads to extra work.
   * Keep it so far for testing purposes.
  */
  /*if (nodesize == 1) {
     return NULL;
     } */
  /* Divide up the quadrants equally */
  SC3E (sc3_mpienv_get_nodesize (p3->split_info, &nodesize));
  SC3E (sc3_mpienv_get_node_num (p3->split_info, &node_num));
  SC3E (sc3_mpienv_get_node_offsets (p3->split_info, &node_offsets));
  SC3E (sc3_mpienv_get_noderank (p3->split_info, &noderank));
  SC3E (sc3_mpienv_get_nodecomm (p3->split_info, &nodecomm));
  SC3E (sc3_allocator_malloc
        (p3->alloc, (p3->mpisize + 1) * sizeof (p4est3_gloidx),
         &loc_offsets));
  p3->global_num_quads = p3->old->global_num_quads;

  SC3A_CHECK (loc_offsets != NULL);

  node_offset = node_offsets[node_num];
  node_offset_next = node_offsets[node_num + 1];
/* Adjust quadrant partition information of the forest */
  if (p3->cweight == NULL) {
    qcount_node =
      p3->old->goffset[node_offset_next] - p3->old->goffset[node_offset];
#ifdef P4EST_ENABLE_DEBUG
    if (noderank == 0) {
      SC3A_CHECK (p3->old->goffset[p3->mpirank] ==
                  p4est3_glocut (qcount_node, nodesize, noderank));
    }
#endif
    /** TODO: figure out mpisize vs nodesize stuff */
    for (i = 0; i < p3->mpisize; ++i) {
      loc_offsets[i] = p4est3_glocut (qcount_node, nodesize, i);
    }
  }
  else {
    /* TODO: Use local array for weighted correction instead of goffsets */
    SC3E (sc3_MPI_Win_lock (SC3_MPI_LOCK_SHARED, 0, SC3_MPI_MODE_NOCHECK,
                            p3->goffsets->meta->win));
    SC3E (p4est3_weighted_new_boundaries
          (p3, nodesize, node_offset, noderank, nodecomm));
    SC3E (sc3_MPI_Win_sync (p3->goffsets->meta->win));
    SC3E (sc3_MPI_Win_unlock (0, p3->goffsets->meta->win));
    SC3E (sc3_MPI_Barrier (nodecomm));
    for (i = 0; i < p3->mpisize; ++i) {
      loc_offsets[i] = p3->goffset[i];
    }
  }
  loc_offsets[p3->mpisize] = p3->global_num_quads;

  SC3E (sc3_allocator_malloc
        (p3->alloc, p3->mpisize * sizeof (p4est3_gloidx), &last_goffsets));
  if (p3->family) {
    for (i = 0; i < p3->mpisize; ++i) {
      last_goffsets[i] = p3->old->goffset[i + 1] - 1;
    }
    SC3E (sc3_allocator_calloc
          (p3->alloc, p3->mpisize, sizeof (p4est3_locidx), &correction));
    SC3E (p4est3_partition_correction
          (p3, node_offset, node_offset_next,
           last_goffsets, loc_offsets, correction));
    for (i = 0; i < p3->mpisize; ++i) {
      loc_offsets[i] -= correction[i];
    }
    SC3E (sc3_allocator_free (p3->alloc, correction));
  }

  SC3E (sc3_MPI_Win_lock (SC3_MPI_LOCK_SHARED, 0, SC3_MPI_MODE_NOCHECK,
                          p3->goffsets->meta->win));
  /* TODO: Move it as early as possible and use iBarrier */
  p3->goffset[p3->mpirank] = loc_offsets[p3->mpirank];
  if (p3->mpirank == 0) {
    p3->goffset[p3->mpisize] = p3->global_num_quads;
  }
  SC3E (sc3_MPI_Win_sync (p3->goffsets->meta->win));
  SC3E (sc3_MPI_Win_unlock (0, p3->goffsets->meta->win));
  MPI_Ibarrier (nodecomm, &req[0]);

  SC3E (sc3_allocator_malloc
        (p3->alloc, p3->num_trees * sizeof (p4est3_gloidx),
         &last_gtree_offsets));
  for (t = 0; t < p3->num_trees; ++t) {
    last_gtree_offsets[t] = p3->gtroffset[t + 1] - 1;
  }
  SC3E (sc3_MPI_Win_lock (SC3_MPI_LOCK_SHARED, 0, SC3_MPI_MODE_NOCHECK,
                          p3->gtreeoffsets->meta->win));
  SC3E (p4est3_find_first_last_local_trees
        (p3, last_gtree_offsets, loc_offsets));
  SC3E (sc3_MPI_Win_sync (p3->gtreeoffsets->meta->win));
  SC3E (sc3_MPI_Win_unlock (0, p3->gtreeoffsets->meta->win));
  MPI_Ibarrier (nodecomm, &req[1]);

  if (!p3->family) {
    for (i = 0; i < p3->mpisize; ++i) {
      last_goffsets[i] = p3->old->goffset[i + 1] - 1;
    }
  }

  p3->local_num_quads = (p4est3_locidx)
    (loc_offsets[p3->mpirank + 1] - loc_offsets[p3->mpirank]);

  SC3E (sc3_allocator_malloc
        (p3->alloc, nodesize * sizeof (char *), &p3->nodequads));

  if (p3->contiguous) {
    if (p3->qvt == p3->old->qvt) {
      /* Warning! This works only when qvt for p3 and p3->old are the same */
      /* In this case we don't allocate new memory,
         but simply ref to p3->old's one */
      SC3E (p4est3_quadrants_ref (p3->old->quadrants));
      p3->quadrants = p3->old->quadrants;
      p3->quads =
        p3->old->nodequads[0] + loc_offsets[p3->mpirank] * p3->old->qsize;
      for (n = 0; n < nodesize; ++n) {
        p3->nodequads[n]
          = p3->old->nodequads[0] + loc_offsets[n] * p3->old->qsize;
      }
      if (p3->old->_spin_quadrants != NULL) {
        SC3A_IS (sc3_refcount_is_last, &p3->old->_spin_quadrants->meta->rc);
        p3->_spin_quadrants = p3->old->_spin_quadrants;
        SC3E (p4est3_quadrants_ref (p3->old->_spin_quadrants));
      }
    }
    else {
      /* Allocate new shared memory to store quadrants */
      SC3E (p4est3_quadrants_allocate (p3, nodesize));
      SC3E (sc3_MPI_Win_lock (SC3_MPI_LOCK_SHARED, noderank,
                              SC3_MPI_MODE_NOCHECK,
                              p3->quadrants->meta->win));
      for (li = 0; li < p3->local_num_quads; ++li) {
        qcount_node = loc_offsets[p3->mpirank] + (p4est3_gloidx) li;
        SC3E (p4est3_quadrant_translate (p3->old->qvt, p3->old->nodequads[0]
                                         +
                                         qcount_node *
                                         (p4est3_gloidx) p3->old->qsize,
                                         p3->qvt,
                                         p3->quads +
                                         (ptrdiff_t) li *
                                         (ptrdiff_t) p3->qsize));
      }
    }
  }
  else {
    SC3E (sc3_allocator_calloc
          (p3->alloc, p3->mpisize, sizeof (p4est3_locidx), &num_recv_from));
    /* Allocate new shared memory to store quadrants */
    SC3E (p4est3_quadrants_allocate (p3, nodesize));
    SC3E (sc3_MPI_Win_lock (SC3_MPI_LOCK_SHARED, noderank,
                            SC3_MPI_MODE_NOCHECK, p3->quadrants->meta->win));
    SC3E (p4est3_procs_recv_from
          (p3, last_goffsets, num_recv_from,
           &from_begin, &from_end, loc_offsets));
    SC3E (p4est3_quads_copy_from
          (p3, num_recv_from, from_begin, from_end, loc_offsets, noderank));
    SC3E (sc3_allocator_free (p3->alloc, num_recv_from));
  }
  SC3A_CHECK (p3->nodequads[noderank] == p3->quads);

  SC3E (p4est3_local_trees_reproduce (p3, last_gtree_offsets, loc_offsets));

  if ((p3->contiguous && (p3->qvt != p3->old->qvt)) || !p3->contiguous) {
    SC3E (sc3_MPI_Win_sync (p3->quadrants->meta->win));
    SC3E (sc3_MPI_Win_unlock (noderank, p3->quadrants->meta->win));
  }

  SC3E (sc3_allocator_free (p3->alloc, loc_offsets));
  SC3E (sc3_allocator_free (p3->alloc, last_goffsets));
  SC3E (sc3_allocator_free (p3->alloc, last_gtree_offsets));
  MPI_Wait (&req[0], MPI_STATUS_IGNORE);
  MPI_Wait (&req[1], MPI_STATUS_IGNORE);
  return NULL;
}

#ifdef __cplusplus
#if 0
{
#endif
}
#endif
