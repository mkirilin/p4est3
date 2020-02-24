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

static sc3_error_t *
make_allocator (sc3_allocator_t * oa, sc3_allocator_t ** alloc)
{
  SC3A_IS (sc3_allocator_is_setup, oa);
  SC3E (sc3_allocator_new (oa, alloc));
  SC3E (sc3_allocator_setup (*alloc));
  return NULL;
}

static sc3_error_t *
test_p4est_new (sc3_allocator_t * alloc,
                sc3_MPI_Comm_t mpicomm, p4est3_quadrant_vtable_t * qvt,
                p4est3_topidx num_trees, int level)
{
  int                 i;
  p4est_connectivity_t *c4;
  p4est3_connectivity_t *conn;
  p4est3_t           *p3;

  SC3A_IS (sc3_allocator_is_setup, alloc);

  for (i = 0; i < 4; ++i) {
    /* create connectivity structure */
    fprintf (stderr, "Trying %d\n", i);
    switch (i) {
    case 0:
      /* default connectivity with one tree */
      SC3E (p4est3_connectivity_new (alloc, &conn));
      SC3E (p4est3_connectivity_setup (conn));
      break;
    case 1:
      /* virtual connectivity with one tree */
      SC3E (p4est3_connectivity_new_num_trees (alloc, 1, &conn));
      break;
    case 2:
      /* at this stage this object is still dimension-independent */
      SC3E (p4est3_connectivity_new_unitcube (alloc, &conn));
      break;
    case 3:
      /* wrapping a p4est connectivity */
#ifndef P4_TO_P8
      c4 = p4est_connectivity_new_unitsquare ();
#else
      c4 = p8est_connectivity_new_unitcube ();
#endif
      SC3E (p4est3_connectivity_new_p4est (alloc, c4, 1, &conn));
      break;
    default:
      SC3E_UNREACH ("Invalid example counter");
    }

    /* create p4est object with connectivity */
    SC3E (p4est3_new (alloc, &p3));
    SC3E (p4est3_set_comm (p3, mpicomm, 1));
    SC3E (p4est3_set_connectivity (p3, conn));
    SC3E (p4est3_set_vtable (p3, qvt));
    SC3E (p4est3_set_level (p3, level));
    SC3E (p4est3_setup (p3));

    SC3E (p4est3_destroy (&p3));
    SC3E (p4est3_connectivity_destroy (&conn));
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
  }

  if (!sc3_allocator_is_free (mainalloc, reason)) {
    fprintf (stderr, "Allocation error: %s\n", reason);
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

  /* v3 standard procedure to isolate memory allocation contexts */
  mainalloc = sc3_allocator_nothread ();

  /* legacy wrapping for p4est quadrants */
  p4est_quadrant_vtable (qvt, 0);

  /* command line parameters */
  num_trees = 2;
  level = 3;

  /* this is generally needed for MPI */
  SC3E_SET (e, sc3_MPI_Init (&argc, &argv));

  /* we don't need init calls for v3.  Just to check legacy wrapping */
  mpicomm = SC3_MPI_COMM_WORLD;
  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  SC3E_NULL_SET (e, make_allocator (mainalloc, &alloc));

  SC3E_NULL_SET (e, test_p4est_new (alloc, mpicomm, qvt, num_trees, level));

  SC3E_NULL_SET (e, free_allocator (&alloc));

  /* again, just to check legacy wrapping */
  SC3E_NULL_REQ (e, !sc_finalize_noabort ());

  /* TODO: call finalize even with errors? */
  SC3E_NULL_SET (e, sc3_MPI_Finalize ());
  report_errors (mainalloc, &e);

  return 0;
}
