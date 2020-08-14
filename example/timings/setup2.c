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
#include <p4est3_quadrant_zyx.h>
#include <p4est_p4est3.h>
#else
#include <p8est3_quadrant_zyx.h>
#include <p8est_p4est3.h>
#endif

#include <time.h>

#define test(SETUP_MODE, qvt, t) do {                                       \
  SC3E_NULL_SET (e, make_new_p4est3 (&p3, alloc, conn, mpicomm, qvt,        \
                                     level, SETUP_MODE, &tb, &te));         \
  SC3E_NULL_SET (e, measure_setup (tb, te, mpicomm, mpirank, mpisize, &t)); \
  } while (0)

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
                 p4est3_setup_mode_t mode, clock_t * tb, clock_t * te)
{
  SC3A_IS (sc3_allocator_is_setup, alloc);

  /* create p4est object with connectivity */
  SC3E (p4est3_new (alloc, p3));
  SC3E (p4est3_set_comm (*p3, mpicomm, 1));
  SC3E (p4est3_set_connectivity (*p3, conn));
  SC3E (p4est3_set_vtable (*p3, qvt));
  SC3E (p4est3_set_level (*p3, level));
  SC3E (p4est3_set_setup_mode (*p3, mode));

  *tb = clock ();
  SC3E (p4est3_setup (*p3));
  *te = clock ();

  SC3E (p4est3_destroy (p3));
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

static sc3_error_t *
measure_setup (clock_t tb, clock_t te, sc3_MPI_Comm_t mpicomm,
               int mpirank, int mpisize, float *time)
{
  float              *times;
  int                 i;
  *time = (float) (te - tb) / CLOCKS_PER_SEC;
  times = (float *) malloc (sizeof (float) * mpisize);

  SC3E (sc3_MPI_Allgather (time, 1, SC3_MPI_FLOAT,
                           times, 1, SC3_MPI_FLOAT, mpicomm));
  for (i = 0, *time = -1.; i < mpisize; ++i) {
    if (*time < times[i]) {
      *time = times[i];
    }
  }

  free (times);
  return NULL;
}

void
print_stats (const char *name, float time, float time_avx, float min_rec,
             float min_rec_avx)
{
  printf ("\n%s: \n"
          "  Vectorized:            %f\n"
          "    Rec/Curr Ratio:      %f\n"
          "  Non-Vectorized:        %f\n"
          "    Rec/Curr Ratio:      %f\n"
          "  Vect/Non-Vect Ratio:   %f\n",
          name, time_avx, min_rec_avx / time_avx, time, min_rec / time,
          time_avx / time);
}

int
main (int argc, char **argv)
{
  p4est3_topidx       num_trees;
  sc3_allocator_t    *alloc, *mainalloc;
  sc3_error_t        *e;
  sc3_MPI_Comm_t      mpicomm;
  p4est3_quadrant_vtable_t vtable, *qvt = &vtable;
  p4est3_quadrant_vtable_t vtable_avx, *qvt_avx = &vtable_avx;
  p4est3_t           *p3;
  p4est3_connectivity_t *conn;
  clock_t             tb, te;
  int                 level, mpirank, mpisize;
  float               mtime, stime, rtime, rctime, rrtime,
    mtime_avx, stime_avx, rtime_avx, rctime_avx, rrtime_avx;
  float               min_rec, min_rec_avx;

  /* v3 standard procedure to isolate memory allocation contexts */
  mainalloc = sc3_allocator_nothread ();
  mpicomm = SC3_MPI_COMM_WORLD;

  /* legacy wrapping for p4est quadrants */
  p4est3_quadrant_zyx_vtable (qvt_avx);
  p4est_quadrant_vtable (qvt, 0);

  /* this is generally needed for MPI */
  SC3E_SET (e, sc3_MPI_Init (&argc, &argv));

  level = 1;
  num_trees = 1;
  if (argc == 2) {
    level = atoi (argv[1]);
    num_trees = 1;
  }
  else if (argc == 3) {
    level = atoi (argv[1]);
    num_trees = atoi (argv[2]);
  }
  else {
    sc_MPI_Abort (mpicomm, -1);
  }

  /* we don't need init calls for v3.  Just to check legacy wrapping */
  /* must not use SC3_MPI_COMM_WORLD due to incompatible non-mpi wrapping */
  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);
  SC3E_NULL_SET (e, make_allocator (mainalloc, &alloc));

  SC3E_NULL_SET (e, p4est3_connectivity_new (alloc, &conn));
  SC3E_NULL_SET (e, p4est3_connectivity_set_num_trees (conn, num_trees));
  SC3E_NULL_SET (e, p4est3_connectivity_setup (conn));

  SC3E_NULL_SET (e, sc3_MPI_Comm_rank (mpicomm, &mpirank));
  SC3E_NULL_SET (e, sc3_MPI_Comm_size (mpicomm, &mpisize));

  test (P4EST3_NEW_MORTON, qvt, mtime);
  test (P4EST3_NEW_SUCCESSOR, qvt, stime);
  test (P4EST3_NEW_RECURSIVE, qvt, rtime);
  test (P4EST3_NEW_RECURSIVE_CHILD, qvt, rctime);
  test (P4EST3_NEW_RECURSIVE_REGION, qvt, rrtime);

  //SIMD/AVX area
  test (P4EST3_NEW_MORTON, qvt_avx, mtime_avx);
  test (P4EST3_NEW_SUCCESSOR, qvt_avx, stime_avx);
  test (P4EST3_NEW_RECURSIVE, qvt_avx, rtime_avx);
  test (P4EST3_NEW_RECURSIVE_CHILD, qvt_avx, rctime_avx);
  test (P4EST3_NEW_RECURSIVE_REGION, qvt_avx, rrtime_avx);

  SC3E_NULL_SET (e, p4est3_connectivity_destroy (&conn));
  SC3E_NULL_SET (e, free_allocator (&alloc));

  if (mpirank == 0) {
    min_rec = SC_MIN (SC_MIN (rtime, rctime), rrtime);
    min_rec_avx = SC_MIN (SC_MIN (rtime_avx, rctime_avx), rrtime_avx);
    print_stats ("Morton", mtime, mtime_avx, min_rec, min_rec_avx);
    print_stats ("Successor", stime, stime_avx, min_rec, min_rec_avx);
    print_stats ("Recursive", rtime, rtime_avx, min_rec, min_rec_avx);
    print_stats ("Recursive_child", rctime, rctime_avx, min_rec, min_rec_avx);
    print_stats ("Recursive_region", rrtime, rrtime_avx, min_rec,
                 min_rec_avx);
  }
  /* again, just to check legacy wrapping */
  SC3E_NULL_REQ (e, !sc_finalize_noabort ());

  /* TODO: call finalize even with errors? */
  SC3E_NULL_SET (e, sc3_MPI_Finalize ());
  report_errors (mainalloc, &e);

  return 0;
}
