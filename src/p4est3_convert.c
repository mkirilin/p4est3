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

#include <p4est3_convert.h>
#include <p4est3_internal.h>
#include <p4est3_p4est.h>

#ifdef __cplusplus
extern              "C"
{
#if 0
}
#endif
#endif

sc3_error_t        *
p4est3_convert_p4est (p4est_t * p, p4est3_t * p3)
{
  int is_comm_same;
  p4est3_connectivity_t *conn;
  const p4est3_quadrant_vtable_t *qvt;
  p4est_tree_t *t;
  SC3A_IS (p4est3_is_new, p3);
  SC3E_DEMAND (sc_MPI_Comm_compare (p3->mpicomm, p->mpicomm, &is_comm_same)
                == sc_MPI_SUCCESS, "Cannot compare p2 and p3 mpi comms");
  SC3E_DEMAND (is_comm_same, "p2 and p3 mpi communicators are different");

  /* inherit p4est_t's connectivity and set it up */
  p3->accessed_conn = 0;
  SC3E (p4est3_connectivity_new_p4est
        (p3->alloc, conn, p->connectivity, 1));
  if (p3->conn != NULL) {
    SC3E (p4est3_connectivity_unref (p3->conn));
  }
  SC3E (p4est3_set_connectivity (p3, conn));

  /* pretend it never existed */
  SC3E_DEMAND (p3->pvt == NULL && p3->slf == NULL, "Forest vtable exists");

  /* if qvt is not seup in advance, we set up the one for classical quads */
  if (p3->qvt == NULL) {
    /* this call also sets p4est3_t::sqvt */
    SC3E (p4est3_quadrant_vtable_p4est (&qvt));
    SC3E (p4est3_set_quadrant_vtable (p3, qvt));
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


  return NULL;
}

#ifdef __cplusplus
#if 0
{
#endif
}
#endif
