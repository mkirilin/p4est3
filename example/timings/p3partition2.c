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

#include <p4est3.h>
#ifndef P4_TO_P8
#include <p4est_extended.h>
#include <p4est_bits.h>
#include <p4est3_p4est.h>
#include <p4est3_quadrant_yx.h>
#include <p4est3_quadrant_mort2d.h>
#else
#include <p8est_extended.h>
#include <p8est_bits.h>
#include <p4est3_p8est.h>
#include <p4est3_quadrant_zyx.h>
#include <p4est3_quadrant_mort3d.h>
#endif

/* Generally speaking, it is not allowed to use _internal headers in
 applications. We use it only in exceptional cases as timings and tests. */
#include <p4est3_internal.h>
#include <sc_statistics.h>
#include <sc_flops.h>

#include <string.h>

static sc3_error_t *
check_quadrant_type (int argc, char **argv,
                     const p4est3_quadrant_vtable_t **qvt,
                     int mpirank, sc3_MPI_Comm_t mpicomm)
{
  if (argc <= 2) {
    return NULL;
  }
  if (strcmp (argv[1], "P4EST2") == 0) {
    return NULL;
  }
  if (strcmp (argv[2], "STANDARD") == 0) {
    SC3E (p4est3_quadrant_vtable_p4est (qvt));
  }
  else if (strcmp (argv[2], "AVX") == 0) {
    SC3E (p4est3_quadrant_yx_vtable (qvt));
  }
  else if (strcmp (argv[2], "MORT_ORD") == 0) {
    SC3E (p4est3_quadrant_mort2d_vtable (qvt));
  }
  else {
    if (mpirank == 0) {
      wrong_input (argv[2], 2);
      sc_MPI_Abort (mpicomm, -1);
    }
  }
  SC3E_DEMAND (*qvt != NULL, "AVX is not supported by hardware or"
               "p4est is not build neither in 2D nor 3D");
  return NULL;
}

int
main (int argc, char **argv)
{
  const p4est3_quadrant_vtable_t *qvt;
  p4est3_topidx       num_trees;
  sc3_allocator_t    *alloc, *mainalloc;
  sc3_error_t        *e;
  sc3_MPI_Comm_t      mpicomm;
  p4est3_t           *p3;
  p4est_t            *p;
  p4est3_connectivity_t *conn;
  p4est_connectivity_t *conn_old;
  int                 level, mpirank, mpisize;
  sc_flopinfo_t       fi, snapshot;
  sc_statinfo_t       stats;
  char               *heading;

  /* v3 standard procedure to isolate memory allocation contexts */
  mainalloc = sc3_allocator_nothread ();

  /* this is generally needed for MPI */
  SC3E_SET (e, sc3_MPI_Init (&argc, &argv));

  SC3E_NULL_SET (e, sc3_MPI_Comm_rank (mpicomm, &mpirank));
  SC3E_NULL_SET (e, sc3_MPI_Comm_size (mpicomm, &mpisize));

  /* default parameters */
  p4est3_setup_mode_t mode = P4EST3_NEW_MORTON;
  SC3E_NULL_SET (e, p4est3_quadrant_vtable_p4est (&qvt));
  level = 1;
  num_trees = 2;
  if (argc == 1 && mpirank == 0) {
    printf ("Execution without parameters. "
            "Default parameters are applied.\n"
            "Parameter's format: "
            "<QUADRANT TYPE> <#levels> <#trees>\n");
  }

  SC3E_NULL_SET (e, check_quadrant_type (argc, argv, &qvt, mpirank, mpicomm));
  if (argc > 2) {
    level = atoi (argv[3]);
    if (level == 0 && mpirank == 0) {
      wrong_input (argv[3], 3);
      sc_MPI_Abort (mpicomm, -1);
    }
  }
  if (argc > 3) {
    num_trees = atoi (argv[4]);
    if (num_trees == 0 && mpirank == 0) {
      wrong_input (argv[4], 4);
      sc_MPI_Abort (mpicomm, -1);
    }
  }
  sc3_MPI_Barrier (mpicomm);
  heading = set_heading (argc, argv, mpicomm);
}