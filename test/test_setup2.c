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
#include <p4est_p4est3.h>
#else
#include <p8est_p4est3.h>
#endif

#define MAX_TEST_LEVEL 5
#define MAX_TEST_TREES 5

static sc3_error_t *
make_allocator (sc3_allocator_t * oa, sc3_allocator_t ** alloc)
{
  SC3A_IS (sc3_allocator_is_setup, oa);
  SC3E (sc3_allocator_new (oa, alloc));
  SC3E (sc3_allocator_setup (*alloc));
  return NULL;
}

static sc3_error_t *
make_new_p4est3 (p4est3_t ** p3, sc3_allocator_t * alloc,
                 p4est3_connectivity_t * conn, sc3_MPI_Comm_t mpicomm,
                 p4est3_quadrant_vtable_t * qvt, int level,
                 p4est3_setup_mode_t mode)
{
  SC3A_IS (sc3_allocator_is_setup, alloc);

  /* create p4est object with connectivity */
  SC3E (p4est3_new (alloc, p3));
  SC3E (p4est3_set_comm (*p3, mpicomm, 1));
  SC3E (p4est3_set_connectivity (*p3, conn));
  SC3E (p4est3_set_vtable (*p3, qvt));
  SC3E (p4est3_set_level (*p3, level));
  SC3E (p4est3_set_setup_mode (*p3, mode));
  SC3E (p4est3_setup (*p3));

  return NULL;
}

static sc3_error_t *
compare_p4est3_quadrants (const p4est3_t * lhs, const p4est3_t * rhs,
                          p4est3_quadrant_vtable_t * qvt)
{
  p4est3_gloidx       i;
  p4est3_gloidx       gln, grn;
  p4est3_locidx       lln, lrn;
  int                 is_eq;
  char               *lchar_q = NULL, *rchar_q = NULL;

  SC3E (p4est3_get_quadrants (lhs, &lchar_q));
  SC3E (p4est3_get_quadrants (rhs, &rchar_q));
  size_t              q_size = p4est3_quadrant_size (qvt);

  SC3E (p4est3_get_global_num_quads (lhs, &gln));
  SC3E (p4est3_get_global_num_quads (rhs, &grn));
  SC3E (p4est3_get_local_num_quads (lhs, &lln));
  SC3E (p4est3_get_local_num_quads (rhs, &lrn));
  SC3E_DEMAND (gln == grn, "different setup modes have different output");
  SC3E_DEMAND (lln == lrn, "different setup modes have different output");
  //lhs and rhs should be valid
  for (i = 0; i < lln; ++i, lchar_q += q_size, rchar_q += q_size) {
    SC3E (p4est3_quadrant_compare (qvt, lchar_q, rchar_q, &is_eq));
    SC3E_DEMAND (is_eq == 0, "setup mod's results differ");
  }

  return NULL;
}

static sc3_error_t *
free_allocator (sc3_allocator_t ** alloc)
{
  SC3A_IS (sc3_allocator_is_setup, *alloc);
  SC3E (sc3_allocator_destroy (alloc));
  return NULL;
}

static void
report_errors (sc3_allocator_t * mainalloc, sc3_error_t ** pe)
{
  char                eflat[SC3_BUFSIZE];
  char                reason[SC3_BUFSIZE];

  if (pe != NULL && *pe != NULL) {
    /* TODO print error messages in a nicer way */
    sc3_error_destroy_noerr (pe, eflat);
    fprintf (stderr, "Error: %s\n", eflat);
    SC_CHECK_ABORT (0, "Setup's tests failed\n");
  }

  if (!sc3_allocator_is_free (mainalloc, reason)) {
    fprintf (stderr, "Allocation error: %s\n", reason);
    SC_CHECK_ABORT (0, "Allocation's tests failed\n");
  }

#if 0
  e = sc3_error_destroy (pe);

  /* TODO synchronize e across MPI processes */

  if (e != NULL) {
    fprintf (stderr, "Errors remain\n");
    sc3_error_destroy (&e);
  }
#endif
}

int
main (int argc, char **argv)
{
  int                 level;
  p4est3_topidx       num_trees;
  sc3_allocator_t    *alloc, *mainalloc;
  sc3_error_t        *e;
  sc3_MPI_Comm_t      mpicomm;
  p4est3_quadrant_vtable_t vtable, *qvt = &vtable;
  p4est3_t           *p3m, *p3s, *p3r;
  p4est3_connectivity_t *conn;

  /* v3 standard procedure to isolate memory allocation contexts */
  mainalloc = sc3_allocator_nothread ();

  /* legacy wrapping for p4est quadrants */
  p4est_quadrant_vtable (qvt, 0);

  /* this is generally needed for MPI */
  SC3E_SET (e, sc3_MPI_Init (&argc, &argv));

  /* we don't need init calls for v3.  Just to check legacy wrapping */
  /* must not use SC3_MPI_COMM_WORLD due to incompatible non-mpi wrapping */
  mpicomm = SC3_MPI_COMM_WORLD;
  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);
  SC3E_NULL_SET (e, make_allocator (mainalloc, &alloc));
  for (level = 1; level < MAX_TEST_LEVEL; ++level) {
    for (num_trees = 1; num_trees < MAX_TEST_TREES; ++num_trees) {
      SC3E_NULL_SET (e, sc3_MPI_Barrier (mpicomm));

      SC3E_NULL_SET (e, p4est3_connectivity_new (alloc, &conn));
      SC3E_NULL_SET (e, p4est3_connectivity_set_num_trees (conn, num_trees));
      SC3E_NULL_SET (e, p4est3_connectivity_setup (conn));

      //P4EST3_NEW_MORTON
      SC3E_NULL_SET (e, make_new_p4est3 (&p3m, alloc, conn, mpicomm, qvt,
                                         level, P4EST3_NEW_MORTON));
      //P4EST3_NEW_SUCCESSOR
      SC3E_NULL_SET (e, make_new_p4est3 (&p3s, alloc, conn, mpicomm, qvt,
                                         level, P4EST3_NEW_SUCCESSOR));
      //P4EST3_NEW_RECURSIVE
      SC3E_NULL_SET (e, make_new_p4est3 (&p3r, alloc, conn, mpicomm, qvt,
                                         level, P4EST3_NEW_RECURSIVE));
      SC3E_NULL_SET (e, compare_p4est3_quadrants (p3m, p3s, qvt));
      SC3E_NULL_SET (e, compare_p4est3_quadrants (p3m, p3r, qvt));

      SC3E_NULL_SET (e, p4est3_destroy (&p3m));
      SC3E_NULL_SET (e, p4est3_destroy (&p3s));
      SC3E_NULL_SET (e, p4est3_destroy (&p3r));
      SC3E_NULL_SET (e, p4est3_connectivity_destroy (&conn));
    }
  }
  SC3E_NULL_SET (e, free_allocator (&alloc));

  /* again, just to check legacy wrapping */
  SC3E_NULL_REQ (e, !sc_finalize_noabort ());

  /* TODO: call finalize even with errors? */
  SC3E_NULL_SET (e, sc3_MPI_Finalize ());
  report_errors (mainalloc, &e);

  return 0;
}
