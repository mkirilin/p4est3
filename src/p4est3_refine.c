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

#include <p4est3_refine.h>
#include <p4est3_iterate.h>
#include <p4est3_internal.h>

#ifdef __cplusplus
extern              "C"
{
#if 0
}
#endif
#endif

typedef struct coarse_callback_data
{
  sc3_array_t        *family;
  sc3_array_t        *pattern;
  int                 nsiblings;
}
coarse_callback_data_t;

/* A trivial callback that copies volumes' level into levels array */
static sc3_error_t *
p4est3_copy_volume_callback (p4est3_iterate_volume_info_t * vi)
{
  int                 level;
  int                *pattern_it;

  SC3E (p4est3_quadrant_level (vi->p3->qvt, vi->quadrant, &level));
  SC3E (sc3_array_push ((sc3_array_t *) vi->user_data, &pattern_it));
  *pattern_it = level;
  return NULL;
}

/* Specific volume iterator callback that processes volume info by refinement
   callback and returns (with user_data) population pattern
   (see refinement documentation) */
static sc3_error_t *
p4est3_refine_volume_callback (p4est3_iterate_volume_info_t * vi)
{
  SC3A_CHECK (vi->p3->crefine != NULL);
  int                 is_refine, level, i;
  int                *pattern_it;
  /* pack data for refinement callback input */
  p4est3_refine_callback_info_t ri = { vi->p3, vi->ntree, vi->quadrant };
  SC3E (p4est3_quadrant_level (vi->p3->qvt, vi->quadrant, &level));

  SC3E (vi->p3->crefine (&ri, &is_refine));
  if (!is_refine) {
    SC3E (sc3_array_push ((sc3_array_t *) vi->user_data, &pattern_it));
    *pattern_it = level;
  }
  else {
    for (i = 0; i < vi->p3->num_children; ++i) {
      SC3E (sc3_array_push ((sc3_array_t *) vi->user_data, &pattern_it));
      SC3A_CHECK (level + 1 <= vi->p3->qmaxlevel);
      *pattern_it = level + 1;
    }
  }
  return NULL;
}

/* Specific volume iterator callback that processes volume info by coarsining
   callback and returns (with user_data) population pattern
   (see coarsining documentation) */
static sc3_error_t *
p4est3_coarse_volume_callback (p4est3_iterate_volume_info_t * vi)
{
  SC3A_CHECK (vi->p3->ccoarse != NULL);
  int                 is_coarse, level, i, child_id;
  int                *pattern_it;
  void               *quad;
  coarse_callback_data_t *cdata = (coarse_callback_data_t *) vi->user_data;
  p4est3_coarse_callback_info_t ci;     /* = { vi->p3, vi->ntree, vi->quadrant }; */

  /* Decide if we call coarse callback.
     We do this only if we find a whole family. */
  SC3E (p4est3_quadrant_child_id (vi->p3->qvt, vi->quadrant, &child_id));
  if (cdata->nsiblings != child_id) {
    /* possibly the benning of a new family */
    if (child_id == 0) {
      /* the beginning indeed, start the family */
      SC3E (sc3_array_index (cdata->family, 0, &quad));
      quad = vi->quadrant;
      cdata->nsiblings = 1;
    }
    else {
      /* cannot be a part of a complete family */
      SC3E (sc3_array_push ((sc3_array_t *) vi->user_data, &pattern_it));
      SC3E (p4est3_quadrant_level (vi->p3->qvt, vi->quadrant, &level));
      *pattern_it = level;
    }
    return NULL;
  }

  /* cdata->nsiblings == child_id: the volume is a part of the family */
  SC3E (sc3_array_index (cdata->family, cdata->nsiblings, &quad));
  quad = vi->quadrant;
  cdata->nsiblings++;

  if (cdata->nsiblings == vi->p3->num_children) {
    /* we have complete family, pack data and call coarsening */
    ci.p3 = vi->p3;
    ci.ntree = vi->ntree;
    ci.family = cdata->family;
    SC3E (vi->p3->ccoarse (&ci, &is_coarse));
    cdata->nsiblings = 0;
  }
  else {
    /* nothing left to do here, go to the next volume */
    return NULL;
  }

  /* fill in the level information */
  SC3E (p4est3_quadrant_level (vi->p3->qvt, vi->quadrant, &level));
  if (is_coarse) {
    SC3E (sc3_array_push ((sc3_array_t *) vi->user_data, &pattern_it));
    SC3A_CHECK (level - 1 >= 0);
    *pattern_it = level - 1;
  }
  else {
    for (i = 0; i < vi->p3->num_children; ++i) {
      SC3E (sc3_array_push ((sc3_array_t *) vi->user_data, &pattern_it));
      *pattern_it = level;
    }
  }
  return NULL;
}

static sc3_error_t *
p4est3_refine_array_new (sc3_allocator_t * alloc, size_t esize, int ealloc,
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
p4est3_lowest_children_pattern (p4est3_t * p3, p4est3_tree_t * tree,
                                const void *q, sc3_array_t * levelq,
                                sc3_array_t * pattern, int *deep_idx)
{
  int                *deep;
  void               *child;
  int                 i, level;

  SC3E (p4est3_quadrant_level (p3->qvt, q, &level));
  SC3E (sc3_array_index (pattern, *deep_idx, &deep));
  SC3A_CHECK (level <= *deep && *deep <= p3->qmaxlevel);
  if (level < *deep) {
    SC3E (sc3_array_index (levelq, level + 1, &child));
    for (i = 0; i < p3->num_children; ++i) {
      SC3E (p4est3_quadrant_child (p3->qvt, q, i, child));
      SC3E (p4est3_lowest_children_pattern
            (p3, tree, child, levelq, pattern, deep_idx));
    }
  }
  else {
    SC3E (p4est3_quadrant_copy
          (p3->qvt, q, tree->tquads + *deep_idx * p3->qsize));
    (*deep_idx)++;
  }
  return NULL;
}

/** This functions does not support multithreading
   Initial value for variables level and deep_idx here must be 0.
   Initial values for range_begin and _end are 0
   and the number of quads at level of p3->qmaxlevel respectively.
   tree->first_tquad, last_tquad and end_tquad are supposed to
   consider uniform mesh of possible reached level of the forest.
   That means, that these values must be previously recalculated
   if necessary (i.e. via first_ last_descendat) */
static sc3_error_t *
p4est3_pattern_populate_tree_rec (p4est3_t * p3, int level,
                                  p4est3_tree_t * tree,
                                  sc3_array_t * pattern,
                                  sc3_array_t * levelq,
                                  p4est3_locidx range_begin,
                                  p4est3_locidx range_end, int *deep_idx)
{
  void               *q, *child;
  int                 i;
  p4est3_locidx       n_lowerq;

  if (range_begin >= tree->first_tquad && range_end <= tree->end_tquad) {
    SC3E (sc3_array_index (levelq, level, &q));
    SC3E (p4est3_lowest_children_pattern
          (p3, tree, q, levelq, pattern, deep_idx));
  }
  else if (range_begin > tree->last_tquad || range_end <= tree->first_tquad) {
    return NULL;
  }
  else {
    n_lowerq = (range_end - range_begin) / p3->num_children;
    range_end = range_begin + n_lowerq;
    for (i = 0; i < p3->num_children; ++i,
         range_begin += n_lowerq, range_end += n_lowerq) {
      SC3A_CHECK (level <= p3->qmaxlevel);
      SC3E (sc3_array_index (levelq, level, &q));
      SC3E (sc3_array_index (levelq, level + 1, &child));
      SC3E (p4est3_quadrant_child (p3->qvt, q, i, child));
      SC3E (p4est3_pattern_populate_tree_rec
            (p3, level + 1, tree, pattern, levelq,
             range_begin, range_end, deep_idx));
    }
  }
  return NULL;
}

static sc3_error_t *
p4est3_pattern_populate_tree (p4est3_t * p3, p4est3_tree_t * tree,
                              sc3_array_t * pattern, sc3_array_t * levelq,
                              int *deep_idx)
{
  void               *quad;
  p4est3_locidx       range_end;
  p4est3_tree_t      *oldtree;

  /* Prepare some temporaty values of the tree here */
  SC3E (p4est3_tree_index (p3->old, tree->treeid, &oldtree));
  quad = (void *) (oldtree->tquads + p3->old->qsize * oldtree->first_tquad);
  SC3E (p4est3_quadrant_first_descendant
        (p3->old->qvt, quad, p3->old->qmaxlevel, p3->old->temp_quad[0]));
  SC3E (p4est3_quadrant_linear_id
        (p3->old->qvt, p3->old->temp_quad[0],
         p3->old->qmaxlevel, &tree->first_tquad));

  quad = (void *) (oldtree->tquads + p3->old->qsize * oldtree->last_tquad);
  SC3E (p4est3_quadrant_last_descendant
        (p3->old->qvt, quad, p3->old->qmaxlevel, p3->old->temp_quad[0]));
  SC3E (p4est3_quadrant_linear_id
        (p3->old->qvt, p3->old->temp_quad[0],
         p3->old->qmaxlevel, &tree->last_tquad));
  tree->end_tquad = tree->last_tquad + 1;
  range_end = p4est3_quadrant_num_uniform (p3->qvt, p3->qmaxlevel);

  /* populate the tree */
  SC3E (p4est3_pattern_populate_tree_rec
        (p3, 0, tree, pattern, levelq, 0, range_end, deep_idx));

  return NULL;
}

sc3_error_t        *
p4est3_refine (p4est3_t * p3)
{
  int                 i, nodesize;
  int                 dispunit;
  int                 noderank;
  char               *quadmem, *nqmem;
  p4est3_locidx       lt_offset;
  p4est3_locidx      *local_num_quads;  /**< Array of the numbers of quadrants at every rank */
  p4est3_locidx      *first_tree_quads; /**< Array of the numbers of quadrants at the first local tree */
  sc3_MPI_Info_t      info_noncontig;
  sc3_MPI_Comm_t      nodecomm;
  sc3_MPI_Aint_t      tempbytes, goffsetbytes;
  p4est3_tree_t      *tree;
  sc3_array_t        *pattern; /**< Array of patterns that we will use
                              to run top-down forest creation */
  sc3_array_t        *levelq;
#ifdef P4EST_ENABLE_DEBUG
  int                 level, new_count;
#endif

  /* We suppose to call this function after setting up routine */
  SC3A_IS (p4est3_is_new, p3);
  SC3A_CHECK (p3->crefine != NULL);
  SC3A_CHECK (p3->old != NULL);
  SC3A_IS (p4est3_is_setup, p3->old);

  SC3E (sc3_mpienv_get_nodesize (p3->split_info, &nodesize));
  SC3E (sc3_mpienv_get_noderank (p3->split_info, &noderank));
  SC3E (sc3_mpienv_get_nodecomm (p3->split_info, &nodecomm));
  SC3E (sc3_mpienv_get_info_noncont (p3->split_info, &info_noncontig));

  /* Do some preliminary allocations */
  SC3E (sc3_allocator_calloc
        (p3->alloc, p3->mpisize, sizeof (p4est3_locidx), &local_num_quads));
  SC3E (sc3_allocator_calloc
        (p3->alloc, p3->mpisize, sizeof (p4est3_locidx), &first_tree_quads));

  SC3E (p4est3_refine_array_new
        (p3->alloc, sizeof (int),
         p3->local_num_quads * p3->num_children, 0, &pattern));

  SC3E (p4est3_refine_array_new (p3->alloc, p3->qsize, p3->qmaxlevel + 1,
                                 p3->qmaxlevel + 1, &levelq));

  SC3E (p4est3_refine_array_new (p3->alloc, sizeof (p4est3_tree_t),
                                 p3->nltrees, p3->nltrees, &p3->trees));

  /* Here we allocate shared p4est3_t::quadwin and p4est3_t::nodequads.
     We will fill in the latter later. */
  SC3E (sc3_allocator_malloc (p3->alloc, nodesize * sizeof (char *),
                              &p3->nodequads));
  SC3E (sc3_MPI_Win_allocate_shared
        (p3->local_num_quads * p3->qsize, p3->qsize,
         info_noncontig, nodecomm, &quadmem, &p3->quadwin));
  for (i = 0; i < nodesize; ++i) {
    SC3E (sc3_MPI_Win_shared_query (p3->quadwin, i,
                                    &tempbytes, &dispunit, &nqmem));
    SC3A_CHECK (dispunit == p3->qsize);
    SC3A_CHECK (nqmem != NULL || tempbytes == 0);
    p3->nodequads[i] = nqmem;
  }
  p3->quads = quadmem;
  SC3A_CHECK (p3->nodequads[noderank] == p3->quads);

  /* Allocate shared memory for global offsets */
  goffsetbytes = (p3->mpisize + 1) * sizeof (p4est3_gloidx);
  SC3E (sc3_MPI_Win_allocate_shared
        (noderank == 0 ? goffsetbytes : 0, p3->qsize,
         info_noncontig, nodecomm, &p3->goffset, &p3->goffsetwin));
  if (noderank > 0) {
    SC3E (sc3_MPI_Win_shared_query (p3->goffsetwin, 0,
                                    &tempbytes, &dispunit, &p3->goffset));
    SC3A_CHECK (tempbytes >= goffsetbytes);
    SC3A_CHECK (dispunit == sizeof (p4est3_gloidx));
    SC3A_CHECK (p3->goffset != NULL);
  }
  /* All the preparations are done. Begin with the algorithm. */
  /* A call to fill in level information, that is necessary for
     generation of a new forest's mesh */
  SC3E (p4est3_iterate_volume (p3, p4est3_refine_volume_callback, &pattern));
  SC3E (sc3_array_get_elem_count (pattern, &p3->local_num_quads));

#ifdef P4EST_ENABLE_DEBUG
  SC3E (sc3_array_get_elem_count (pattern, &new_count));
  for (i = 0; i < new_count; ++i) {
    SC3E (sc3_array_index (pattern, i, &level));
    /** TODO: Check on p3->qmaxlevel at populating stage and skip
     * if it is impossible to split*/
    SC3A_CHECK (level <= p3->qmaxlevel);
  }
#endif
  /* We got a pattern of population, and now we populate it
     and set p4est3_tree_t:: treeid, quad_offset and num_quads.
     We also initialize p4est3_tree_t::first_tquad by 0.
     We work on the process-local window onto the quadrants */
  SC3E (sc3_MPI_Win_lock (SC3_MPI_LOCK_EXCLUSIVE, noderank,
                          SC3_MPI_MODE_NOCHECK, p3->quadwin));
  lt_offset = 0;
  for (i = p3->fltree; i <= p3->lltree; ++i) {
    SC3E (p4est3_tree_index (p3, i, &tree));
    tree->treeid = i;
    tree->quad_offset = lt_offset;
    tree->tquads = p3->quads + p3->qsize * tree->quad_offset;
    SC3E (p4est3_pattern_populate_tree
          (p3, tree, pattern, levelq, &lt_offset));
    tree->num_quads = lt_offset - tree->quad_offset;
    tree->first_tquad = 0;
  }
  SC3E (sc3_MPI_Win_unlock (noderank, p3->quadwin));

  /* This check here is only to avoid creating a new mpi datatype. */
  SC3A_CHECK (sizeof (p4est3_locidx) == sizeof (int));
  SC3E (p4est3_tree_index (p3, p3->fltree, &tree));
  SC3E (sc3_MPI_Allgather
        (&p3->local_num_quads, 1, SC3_MPI_INT,
         local_num_quads, p3->mpisize, SC3_MPI_INT, p3->mpicomm));
  SC3E (sc3_MPI_Allgather
        (&tree->num_quads, 1, SC3_MPI_INT,
         first_tree_quads, p3->mpisize, SC3_MPI_INT, p3->mpicomm));

  if (noderank == 0) {
    p3->goffset[0] = 0;
    for (i = 1; i < p3->mpisize + 1; ++i) {
      p3->goffset[i] = p3->goffset[i - 1]
        + (p4est3_gloidx) local_num_quads[i - 1];
    }
  }
  if (p3->mpirank != 0) {
    for (i = p3->gftree[p3->mpirank - 1]; i == p3->fltree; --i) {
      tree->first_tquad += first_tree_quads[i];
    }
  }
  for (i = p3->fltree; i <= p3->lltree; ++i) {
    SC3E (p4est3_tree_index (p3, i, &tree));
    tree->end_tquad = tree->first_tquad + tree->num_quads;
    tree->last_tquad =
      (tree->end_tquad = tree->first_tquad + tree->num_quads) - 1;
  }

  SC3E (sc3_array_destroy (&pattern));
  SC3E (sc3_array_destroy (&levelq));
  SC3E (sc3_allocator_free (p3->alloc, &local_num_quads));
  SC3E (sc3_allocator_free (p3->alloc, &first_tree_quads));

  sc3_MPI_Barrier (nodecomm);
  p3->global_num_quads = p3->goffset[p3->mpisize];

  return NULL;
}

#ifdef __cplusplus
#if 0
{
#endif
}
#endif
