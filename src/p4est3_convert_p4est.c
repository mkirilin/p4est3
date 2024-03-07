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

#ifndef P4_TO_P8
#include <p4est3_convert_p4est.h>
#include <p4est3_p4est.h>
#else
#include <p4est3_convert_p8est.h>
#include <p4est3_p8est.h>
#endif
#include <p4est3_internal.h>
#include <sc3_omp.h>

#ifdef __cplusplus
extern              "C"
{
#if 0
}
#endif
#endif

static sc3_error_t *
p4est3_conv_array_new (sc3_allocator_t * alloc, size_t esize, int ealloc,
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

sc3_error_t        *
p4est3_convert_p4est (p4est_t * p, p4est3_t * p3)
{
  int n, is_comm_same;
  int                 dispunit;
  int nodesize, node_frank, noderank;
  p4est_topidx_t ti;
  p4est3_topidx ti3;
  p4est3_locidx locq_it;
  size_t i;
  p4est3_connectivity_t *conn;
  const p4est3_quadrant_vtable_t *qvt_standard;
  sc3_MPI_Aint_t tempbytes;
  sc3_MPI_Comm_t      nodecomm;
  p4est_tree_t *t;
  p4est3_tree_t *t3;
  SC3A_IS (p4est3_is_new, p3);
  SC3E_DEMAND (sc_MPI_Comm_compare (p3->mpicomm, p->mpicomm, &is_comm_same)
                == sc_MPI_SUCCESS, "Cannot compare p2 and p3 mpi comms");
  SC3E_DEMAND (is_comm_same, "p2 and p3 mpi communicators are different");

  /* inherit p4est_t's connectivity and set it up */
  p3->accessed_conn = 0;
  SC3E (p4est3_connectivity_new_p4est
        (p3->alloc, &conn, p->connectivity, 1));
  if (p3->conn != NULL) {
    SC3E (p4est3_connectivity_unref (p3->conn));
  }
  SC3E (p4est3_set_connectivity (p3, conn));

  /* pretend it never existed */
  SC3E_DEMAND (p3->pvt == NULL && p3->slf == NULL, "Forest vtable exists");

  /* if qvt is not seup in advance, we set up the one for classical quads */
  SC3E (p4est3_quadrant_vtable_p4est (&qvt_standard));
  if (p3->qvt == NULL) {
    /* this call also sets p4est3_t::sqvt */
    SC3E (p4est3_set_quadrant_vtable (p3, qvt_standard));
  }

/* the next two calls are not necessary, since we don't need these parameters,
  but leave them here for completeness */
  t = p4est_tree_array_index (p->trees, 0);
  SC3E (p4est3_set_level (p3, t->maxlevel));
  /* P4EST3_NEW_SUCCESSOR since it is used by p4est_t */
  SC3E (p4est3_set_setup_mode (p3, P4EST3_NEW_SUCCESSOR));

  /* variables populated during p4est3_setup: communicator related */
  /* we support only 1 shared memory node so far */
  SC3E_DEMAND (p3->shared == 1,
               "p3 is going to be non-shared, "
               "while we support 1 SM node so far");

  /* variables populated during p4est3_setup: partition related */
  /** TODO: why is it int type while p4est3_quadrant_size returs size_t? */
  p3->qsize = (int) p4est3_quadrant_size (p3->qvt);
  p3->qmaxlevel = p3->qvt->max_level;
  p3->num_children = p4est3_quadrant_num_children (p3->qvt);
  p3->max_threads = sc3_omp_max_threads ();

  /* setup mpi, this call also sets p4est3_t::commdup */
  SC3E (p4est3_set_comm (p3, p->mpicomm, 1));
  SC3E (p4est3_internal_setup_comm (p3));
  p3->mpisize = p->mpisize;
  p3->mpirank = p->mpirank;

  p3->num_trees = p->trees->elem_count;

  /* allocate magic structures */
  SC3E (p4est3_glotree_new (p3->alloc, &p3->gtrees));
  SC3E (p4est3_glopos_new (p3->alloc, &p3->gposition));
  SC3E (p4est3_glopos_set_qsize (p3->gposition, p3->qsize));
  SC3E (p4est3_glooffs_new (p3->alloc, &p3->goffsets));
  SC3E (p4est3_gtroffs_new (p3->alloc, &p3->gtreeoffsets));
  SC3E (p4est3_gtroffs_set_num_trees (p3->gtreeoffsets, p3->num_trees));
  SC3E (p4est3_quadrants_new (p3->alloc, &p3->quadrants));
  SC3E (p4est3_quadrants_set_local_num_quads (p3->quadrants,
                                              p->local_num_quadrants));
  SC3E (p4est3_quadrants_set_qsize (p3->quadrants, p3->qsize));
  SC3E (p4est3_glopartition_set_mpienv
        (p3->gtrees, p3->gposition, p3->goffsets,
         p3->gtreeoffsets, p3->quadrants, p3->split_info));
  SC3E (p4est3_glopartition_setup
        (p3->gtrees, p3->gposition, p3->goffsets,
         p3->gtreeoffsets, p3->quadrants));

  p3->gftree = p3->gtrees->gftree;
  p3->gfpos = p3->gposition->gfpos;
  p3->goffset = p3->goffsets->goffset;
  p3->gtroffset = p3->gtreeoffsets->gtreeoffset;
  p3->quads = p3->quadrants->quads;

  /* fill in gfpos */
  SC3E (sc3_MPI_Win_lock (SC3_MPI_LOCK_SHARED, 0, SC3_MPI_MODE_NOCHECK,
                          p3->gposition->meta->win));
  if (p3->qvt == qvt_standard) {
    /* just copy */
    SC3E (p4est3_quadrant_copy
          (p3->qvt, &p->global_first_position[p->mpirank],
           p3->gfpos + p3->mpirank * p3->qsize));
  }
  else {
    /* translate */
    SC3E (p4est3_quadrant_translate
          (qvt_standard, &p->global_first_position[p->mpirank],
          p3->qvt, p3->gfpos + p3->mpirank * p3->qsize));
  }
  if (p3->mpirank == 0) {
    /* same but for the (mpisize)'th element */
    if (p3->qvt == qvt_standard) {
    /* just copy */
    SC3E (p4est3_quadrant_copy
          (p3->qvt, &p->global_first_position[p->mpisize],
           p3->gfpos + p3->mpisize * p3->qsize));
    }
    else {
      /* translate */
      SC3E (p4est3_quadrant_translate
            (qvt_standard, &p->global_first_position[p->mpisize],
            p3->qvt, p3->gfpos + p3->mpisize * p3->qsize));
    }
  }
  SC3E (sc3_MPI_Win_sync (p3->gposition->meta->win));
  SC3E (sc3_MPI_Win_unlock (0, p3->gposition->meta->win));

  /* fill in goffset */
  SC3E (sc3_MPI_Win_lock (SC3_MPI_LOCK_SHARED, 0, SC3_MPI_MODE_NOCHECK,
                          p3->goffsets->meta->win));
  p3->goffset[p3->mpirank]
    = (p4est3_gloidx) p->global_first_quadrant[p->mpirank];
  if (p3->mpirank == 0) {
  /* same but for the (mpisize)'th element */
    p3->goffset[p3->mpisize]
      = (p4est3_gloidx) p->global_first_quadrant[p->mpisize];
  }
  SC3E (sc3_MPI_Win_sync (p3->goffsets->meta->win));
  SC3E (sc3_MPI_Win_unlock (0, p3->goffsets->meta->win));

  /* fill in gftree */
  SC3E (sc3_MPI_Win_lock (SC3_MPI_LOCK_SHARED, 0, SC3_MPI_MODE_NOCHECK,
                          p3->gtrees->meta->win));
  p3->gftree[p3->mpirank] = (p4est3_topidx) p->first_local_tree;
  if (p3->mpirank == 0) {
    p3->gftree[p3->mpisize] = p3->num_trees;
  }
  SC3E (sc3_MPI_Win_sync (p3->gtrees->meta->win));
  SC3E (sc3_MPI_Win_unlock (0, p3->gtrees->meta->win));

  /* allocate memory for node proc first quadrant pointers */
  p3->local_num_quads = p->local_num_quadrants;
  p3->global_num_quads = p->global_num_quadrants;

  SC3E (sc3_mpienv_get_nodesize (p3->split_info, &nodesize));
  SC3E (sc3_mpienv_get_node_frank (p3->split_info, &node_frank));
  SC3E (sc3_mpienv_get_noderank (p3->split_info, &noderank));

  SC3E (sc3_allocator_malloc
        (p3->alloc, nodesize * sizeof (char *), &p3->nodequads));
  for (n = 0; n < nodesize; ++n) {
    SC3E (sc3_MPI_Win_shared_query
          (p3->quadrants->meta->win, n, &tempbytes, &dispunit,
           &p3->nodequads[n]));
    SC3A_CHECK (tempbytes >= (sc3_MPI_Aint_t)
                ((p3->goffset[node_frank + n + 1] -
                  p3->goffset[node_frank + n]) * p3->qsize));
    SC3A_CHECK (dispunit == p3->qsize);
    SC3A_CHECK (p3->nodequads[n] != NULL || tempbytes == 0);
  }
  SC3A_CHECK (p3->nodequads[noderank] == p3->quads);

  /* copy quadrants in shared memory */
  if (p3->qvt == qvt_standard) {
    for (ti = p->first_local_tree, locq_it = 0; ti <= p->last_local_tree; ++ti)
    {
      t = p4est_tree_array_index (p->trees, ti - p->first_local_tree);
      for (i = 0; i < t->quadrants.elem_count; ++i) {
        SC3E (p4est3_quadrant_copy
              (p3->qvt, p4est_quadrant_array_index (&t->quadrants, i),
               p3->quads + p3->qsize * (locq_it++)));
      }
    }
  }
  else {
    for (ti = p->first_local_tree, locq_it = 0; ti <= p->last_local_tree; ++ti)
    {
      t = p4est_tree_array_index (p->trees, ti - p->first_local_tree);
      for (i = 0; i < t->quadrants.elem_count; ++i) {
        SC3E (p4est3_quadrant_translate
              (qvt_standard, p4est_quadrant_array_index (&t->quadrants, i),
               p3->qvt, p3->quads + p3->qsize * (locq_it++)));
      }
    }
  }
  SC3A_CHECK (locq_it == p3->local_num_quads);

  SC3E (sc3_allocator_malloc (p3->alloc, p3->max_threads * sizeof (char *),
                              &p3->temp_quad));
  for (n = 0; n < p3->max_threads; ++n) {
    SC3E (sc3_allocator_malloc (p3->alloc, p3->qsize, &p3->temp_quad[n]));
  }
  /* variables populated during p4est3_setup: tree and quadrant storage */
  p3->fltree = p->first_local_tree;
  p3->lltree = p->last_local_tree;
  p3->nltrees = p->last_local_tree - p->first_local_tree + 1;

  /* allocate and fill in trees and trees offsets */
  SC3E (sc3_mpienv_get_nodecomm (p3->split_info, &nodecomm));
  SC3E (p4est3_tree_offsets_communication (p3, noderank, nodecomm));

  /* fill in trees */
  SC3E (p4est3_conv_array_new (p3->alloc, sizeof (p4est3_tree_t),
                               p3->nltrees, p3->nltrees, &p3->trees));
  for (ti3 = p3->fltree; ti3 <= p3->lltree; ++ti3) {
    SC3E (p4est3_tree_index (p3, ti3, &t3));
    t = p4est_tree_array_index (p->trees, ti3 - p->first_local_tree);

    t3->treeid = ti3;
    t3->num_quads = t->quadrants.elem_count;
    t3->quad_offset = t->quadrants_offset;
    t3->tquads = p3->quads + p3->qsize * t3->quad_offset;
    t3->first_tquad = p3->goffset[p3->mpirank] -  p3->gtroffset[ti3];
    t3->end_tquad = t3->first_tquad + t3->num_quads;
    t3->last_tquad = t3->end_tquad - 1;
  }

  /* inherited user data stays as a sc_mempool_t, we leave it to user how to
   manage them */
  p3->user_data = p->user_data_pool;

  /* since refinement, coarsening and partitioning are not available upon
    convertion stage, we do not need these callback functions */
  p3->crefine = NULL;
  p3->ccoarse = NULL;
  p3->cweight = NULL;
  p3->old = NULL;

  /* the following parameters manipulates p4est3_t object on setup stage only,
     so we do not need them */
  p3->partition = 0;
  p3->family = 0;

  p3->setup = 1;
  SC3A_IS (p4est3_is_setup, p3);

  return NULL;
}

#ifdef __cplusplus
#if 0
{
#endif
}
#endif
