/*
  This file is part of p4est, version 3
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
#include <sc3_refcount_internal.h>

/* TODO add context information to sc3_error
        (i.e., which library is producing the error?) */

int
p4est3_is_valid (p4est3_t * p3, char *reason)
{
  SC3E_TEST (p3 != NULL, reason);
  SC3E_IS (sc3_refcount_is_valid, &p3->rc, reason);
  SC3E_IS (sc3_allocator_is_setup, p3->alloc, reason);

  SC3E_TEST (p3->mpicomm != SC3_MPI_COMM_NULL, reason);
  SC3E_TEST (p3->num_trees > 0, reason);
  SC3E_TEST (p3->level >= 0, reason);

  if (!p3->setup) {
    SC3E_TEST (p3->mpisize == 0 && p3->mpirank == 0, reason);
  }
  else {
    SC3E_TEST (p3->qvt != NULL, reason);

    SC3E_TEST (p3->nodesizewin != SC3_MPI_WIN_NULL, reason);
    SC3E_TEST (p3->headcomm != SC3_MPI_COMM_NULL || p3->noderank > 0, reason);
    SC3E_TEST (p3->nodecomm != SC3_MPI_COMM_NULL, reason);
  }

  SC3E_YES (reason);
}

int
p4est3_is_new (p4est3_t * p3, char *reason)
{
  SC3E_IS (p4est3_is_valid, p3, reason);
  SC3E_TEST (!p3->setup, reason);
  SC3E_YES (reason);
}

int
p4est3_is_setup (p4est3_t * p3, char *reason)
{
  SC3E_IS (p4est3_is_valid, p3, reason);
  SC3E_TEST (p3->setup, reason);
  SC3E_YES (reason);
}

sc3_error_t        *
p4est3_new (sc3_allocator_t * alloc, p4est3_t ** pp3)
{
  p4est3_t           *p3;

  SC3E_RETVAL (pp3, NULL);
  SC3A_IS (sc3_allocator_is_setup, alloc);

  SC3E (sc3_allocator_ref (alloc));
  SC3E_ALLOCATOR_CALLOC (alloc, p4est3_t, 1, p3);
  SC3E (sc3_refcount_init (&p3->rc));
  p3->alloc = alloc;
  p3->mpicomm = SC3_MPI_COMM_WORLD;
  p3->nodesizewin = SC3_MPI_WIN_NULL;
  p3->headcomm = SC3_MPI_COMM_NULL;
  p3->nodecomm = SC3_MPI_COMM_NULL;
  p3->num_trees = 1;
  SC3A_IS (p4est3_is_new, p3);

  *pp3 = p3;
  return NULL;
}

sc3_error_t        *
p4est3_set_comm (p4est3_t * p3, sc3_MPI_Comm_t comm, int dup)
{
  SC3A_IS (p4est3_is_new, p3);
  SC3A_CHECK (comm != SC3_MPI_COMM_NULL);

  /* remove previous communicator */
  if (p3->commdup) {
    SC3E (sc3_MPI_Comm_free (&p3->mpicomm));
  }

  /* register new communicator */
  if (dup) {
    SC3E (sc3_MPI_Comm_dup (comm, &p3->mpicomm));
  }
  else {
    p3->mpicomm = comm;
  }
  p3->commdup = dup;
  return NULL;
}

sc3_error_t        *
p4est3_set_vtable (p4est3_t * p3, p4est3_quadrant_vtable_t * qvt)
{
  SC3A_IS (p4est3_is_new, p3);
  SC3A_CHECK (qvt != NULL);

  p3->qvt = qvt;
  return NULL;
}

sc3_error_t        *
p4est3_set_num_trees (p4est3_t * p3, p4est3_topidx num_trees)
{
  SC3A_IS (p4est3_is_new, p3);
  SC3A_CHECK (num_trees > 0);

  p3->num_trees = num_trees;
  return NULL;
}

sc3_error_t        *
p4est3_set_level (p4est3_t * p3, int level)
{
  SC3A_IS (p4est3_is_new, p3);
  SC3A_CHECK (level >= 0);

  p3->level = level;
  return NULL;
}

sc3_error_t        *
p4est3_setup (p4est3_t * p3)
{
  int                 max_level, lev;
  int                 headsize, headrank;
  int                 p, next, *ofs;
  int                 dispunit;
  int                *nodesizemem;
  p4est3_gloidx       num_uniform, high_uniform;
  sc3_MPI_Aint_t      nodeabytes;

  /*
   * We will extend the functionality of p4est3_set_* and _setup in the future.
   * They will be used to create a new forest by refinement, partitioning, etc.
   * The forest created is generally immutable after being setup.
   * Currently, we are just creating a uniformly refined forest.
   */

  SC3A_IS (p4est3_is_new, p3);

  /* check conditions that arise due to omitting mandatory _set_ functions */
  SC3E_DEMAND (p3->qvt != NULL, "Quadrant virtual table must be set");

  /* query input communicator */
  SC3E (sc3_MPI_Comm_size (p3->mpicomm, &p3->mpisize));
  SC3E (sc3_MPI_Comm_rank (p3->mpicomm, &p3->mpirank));

  /* create one communicator on each shared-memory node */
  SC3E (sc3_MPI_Comm_split_type (p3->mpicomm, SC3_MPI_COMM_TYPE_SHARED,
                                 0, SC3_MPI_INFO_NULL, &p3->nodecomm));
  SC3E (sc3_MPI_Comm_size (p3->nodecomm, &p3->nodesize));
  SC3E (sc3_MPI_Comm_rank (p3->nodecomm, &p3->noderank));

  /* create communicator that contains the first rank on each node */
  SC3E (sc3_MPI_Comm_split (p3->mpicomm, p3->noderank == 0 ? 0 :
                            SC3_MPI_UNDEFINED, 0, &p3->headcomm));
  SC3A_CHECK ((p3->noderank != 0) == (p3->headcomm == SC3_MPI_COMM_NULL));
  if (p3->noderank == 0) {
    SC3E (sc3_MPI_Comm_size (p3->headcomm, &headsize));
    SC3E (sc3_MPI_Comm_rank (p3->headcomm, &headrank));
    nodeabytes = (2 + 2 * headsize + 1) * sizeof (int);
  }
  else {
    headsize = headrank = 0;
    nodeabytes = 0;
  }

  /* allocate shared memory for information on node and head communicators */
  SC3E (sc3_MPI_Win_allocate_shared
        (nodeabytes, sizeof (int), SC3_MPI_INFO_NULL,
         p3->nodecomm, &nodesizemem, &p3->nodesizewin));
  if (p3->noderank == 0) {
    SC3E (sc3_MPI_Win_lock (SC3_MPI_LOCK_EXCLUSIVE, 0, SC3_MPI_MODE_NOCHECK,
                            p3->nodesizewin));
    nodesizemem[0] = p3->num_nodes = headsize;
    nodesizemem[1] = p3->node_num = headrank;
    p3->node_sizes = &nodesizemem[2];

    /* allgather information about all nodes and compute offsets */
    SC3E (sc3_MPI_Allgather (&p3->nodesize, 1, SC3_MPI_INT,
                             p3->node_sizes, 1, SC3_MPI_INT, p3->headcomm));
    *(ofs = p3->node_offsets = &nodesizemem[2 + headsize]) = 0;
    for (p = 0; p < headsize; ++p) {
      next = *ofs + p3->node_sizes[p];
      *++ofs = next;
    }
    SC3A_CHECK (p3->node_offsets[p3->mpirank] == p3->mpirank);

    /* make sure shared memory contents are consistent */
    SC3E (sc3_MPI_Win_unlock (0, p3->nodesizewin));
    SC3E (sc3_MPI_Barrier (p3->nodecomm));
  }
  else {
    SC3E (sc3_MPI_Win_shared_query (p3->nodesizewin, 0,
                                    &nodeabytes, &dispunit, &nodesizemem));
    SC3A_CHECK (nodeabytes >= (sc3_MPI_Aint_t) sizeof (int));
    SC3A_CHECK (dispunit == (int) sizeof (int));
    SC3A_CHECK (nodesizemem != NULL);

    /* access shared memory written by other process */
    SC3E (sc3_MPI_Barrier (p3->nodecomm));
    SC3E (sc3_MPI_Win_lock (SC3_MPI_LOCK_SHARED, 0, SC3_MPI_MODE_NOCHECK,
                            p3->nodesizewin));
    p3->num_nodes = nodesizemem[0];
    SC3A_CHECK (nodeabytes ==
                (sc3_MPI_Aint_t) ((2 + 2 * p3->num_nodes + 1) *
                                  sizeof (int)));
    p3->node_num = nodesizemem[1];
    p3->node_sizes = &nodesizemem[2];
    p3->node_offsets = &nodesizemem[2 + p3->num_nodes];
    SC3E (sc3_MPI_Win_unlock (0, p3->nodesizewin));
  }

  /* determine uniform refinement level */
  max_level = p4est3_max_level (p3->qvt);
  SC3A_CHECK (max_level >= 0);
  max_level = SC3_MIN (p3->level, max_level);

  /* with number of children determine number of elements per tree */
  p3->num_children = p4est3_num_children (p3->qvt);
  SC3A_CHECK (p3->num_children > 0);
  high_uniform = P4EST3_GLOIDX_MAX / p3->num_children;
  for (num_uniform = 1, lev = 0;
       num_uniform <= high_uniform && lev < max_level; ++lev) {
    /* we iterate so we do not roll over the gloidx limit */
    num_uniform *= p3->num_children;
  }
  max_level = lev;
  SC3A_CHECK (p4est3_glopow (p3->num_children, max_level) == num_uniform);

  /* TODO create trees and quadrants */

  p3->setup = 1;
  SC3A_IS (p4est3_is_setup, p3);
  return NULL;
}

sc3_error_t        *
p4est3_ref (p4est3_t * p3)
{
  SC3A_IS (p4est3_is_setup, p3);
  SC3E (sc3_refcount_ref (&p3->rc));
  return NULL;
}

sc3_error_t        *
p4est3_unref (p4est3_t ** pp3)
{
  int                 waslast;
  sc3_allocator_t    *alloc;
  p4est3_t           *p3;

  SC3E_INOUTP (pp3, p3);
  SC3A_IS (p4est3_is_valid, p3);
  SC3E (sc3_refcount_unref (&p3->rc, &waslast));
  if (waslast) {
    *pp3 = NULL;

    alloc = p3->alloc;
    if (p3->setup) {
      SC3E (sc3_MPI_Win_free (&p3->nodesizewin));
      if (p3->noderank == 0) {
        SC3E (sc3_MPI_Comm_free (&p3->headcomm));
      }
      SC3E (sc3_MPI_Comm_free (&p3->nodecomm));

      /* deallocate internal storage */
    }
    if (p3->commdup) {
      SC3E (sc3_MPI_Comm_free (&p3->mpicomm));
    }
    SC3E_ALLOCATOR_FREE (alloc, p4est3_t, p3);
    SC3E (sc3_allocator_unref (&alloc));
  }
  return NULL;
}

sc3_error_t        *
p4est3_destroy (p4est3_t ** pp3)
{
  p4est3_t           *p3;

  SC3E_INULLP (pp3, p3);
  SC3E_DEMIS (sc3_refcount_is_last, &p3->rc);
  SC3E (p4est3_unref (&p3));

  SC3A_CHECK (p3 == NULL);
  return NULL;
}
