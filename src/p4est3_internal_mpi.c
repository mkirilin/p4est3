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

sc3_error_t        *
p4est3_internal_setup_comm (p4est3_t * p3)
{
  int                 headsize, headrank;
  int                 p, next, *ofs;
  int                 dispunit;
  int                *nodesizemem;
  sc3_MPI_Aint_t      nodeabytes;

  /* this is a special-purpose function to simplify p4est3_setup */
  SC3A_CHECK (p3 != NULL);

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

  /* create info structure to allow for per-rank allocation */
  SC3E (sc3_MPI_Info_create (&p3->info_noncontig));
  SC3E (sc3_MPI_Info_set
        (p3->info_noncontig, "alloc_shared_noncontig", "true"));

  /* allocate shared memory for information on node and head communicators */
  SC3E (sc3_MPI_Win_allocate_shared
        (nodeabytes, sizeof (int),
         p3->info_noncontig, p3->nodecomm, &nodesizemem, &p3->nodesizewin));
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
    p3->num_nodes = nodesizemem[0];
    SC3A_CHECK (nodeabytes ==
                (sc3_MPI_Aint_t) ((2 + 2 * p3->num_nodes + 1) *
                                  sizeof (int)));
    p3->node_num = nodesizemem[1];
    p3->node_sizes = &nodesizemem[2];
    p3->node_offsets = &nodesizemem[2 + p3->num_nodes];
  }

  return NULL;
}

sc3_error_t        *
p4est3_internal_setup_cut (p4est3_t * p3, p4est3_gloidx num_global, int qsize)
{
  int                 p;
  int                 dispunit;
  char               *gfposmem;
  p4est3_gloidx      *countmem;
  sc3_MPI_Aint_t      gfposbytes, countbytes, tempbytes;

  /* this is a special-purpose function to simplify p4est3_setup */
  SC3A_CHECK (p3 != NULL);
  SC3A_CHECK (0 <= p3->mpirank && p3->mpirank < p3->mpisize);
  SC3A_CHECK (0 <= p3->noderank && p3->noderank < p3->nodesize);
  SC3A_CHECK (num_global > 0);
  SC3A_CHECK (qsize > 0);

  /* create shared partition arrays */
  gfposbytes = (p3->mpisize + 1) * qsize;
  SC3E (sc3_MPI_Win_allocate_shared
        (p3->noderank == 0 ? gfposbytes : 0, qsize,
         p3->info_noncontig, p3->nodecomm, &gfposmem, &p3->gfposwin));
  countbytes = (p3->mpisize + 1) * sizeof (p4est3_gloidx);
  SC3E (sc3_MPI_Win_allocate_shared
        (p3->noderank == 0 ? countbytes : 0, sizeof (p4est3_gloidx),
         p3->info_noncontig, p3->nodecomm, &countmem, &p3->countwin));

  /* compute cuts for the whole program without communication */
  if (p3->noderank == 0) {
    SC3E (sc3_MPI_Win_lock (SC3_MPI_LOCK_EXCLUSIVE, 0, SC3_MPI_MODE_NOCHECK,
                            p3->countwin));
    for (p = 0; p <= p3->mpisize; ++p) {
      countmem[p] = p4est3_glocut (num_global, p3->mpisize, p);
    }
    SC3E (sc3_MPI_Win_unlock (0, p3->countwin));
    SC3E (sc3_MPI_Barrier (p3->nodecomm));
  }
  else {
    SC3E (sc3_MPI_Win_shared_query (p3->gfposwin, 0,
                                    &tempbytes, &dispunit, &gfposmem));
    SC3A_CHECK (gfposbytes == tempbytes);
    SC3A_CHECK (dispunit == qsize);
    SC3A_CHECK (gfposmem != NULL);
    SC3E (sc3_MPI_Win_shared_query (p3->countwin, 0,
                                    &tempbytes, &dispunit, &countmem));
    SC3A_CHECK (countbytes == tempbytes);
    SC3A_CHECK (dispunit == (int) sizeof (p4est3_gloidx));
    SC3A_CHECK (countmem != NULL);
    SC3E (sc3_MPI_Barrier (p3->nodecomm));
#ifdef P4EST_ENABLE_DEBUG
    for (p = 0; p <= p3->mpisize; ++p) {
      SC3A_CHECK (countmem[p] == p4est3_glocut (num_global, p3->mpisize, p));
    }
#endif
  }

  /* assign further object members */
  p3->qsize = qsize;
  p3->global_num_quads = num_global;
  p3->count = countmem;
  p3->gfpos = gfposmem;
  return NULL;
}
