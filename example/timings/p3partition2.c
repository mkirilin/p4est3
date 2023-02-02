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

static int          refine_level = 0;
static int          level_shift = 0;

static int
refine_fractal (p4est_t * p4est, p4est_topidx_t which_tree,
                p4est_quadrant_t * q)
{
  int                 qid;

  if ((int) q->level >= refine_level) {
    return 0;
  }
  if ((int) q->level < refine_level - level_shift) {
    return 1;
  }

  qid = p4est_quadrant_child_id (q);
  return (qid == 0 || qid == 3
#ifdef P4_TO_P8
          || qid == 5 || qid == 6
#endif
    );
}

static sc3_error_t *
refine_p3_normal_fn (p4est3_refine_callback_info_t * ri, int *is_refine)
{
  int                 level, child_id;

  SC3E (p4est3_quadrant_level (ri->qvt, ri->quadrant, &level));
  if (level >= refine_level) {
    *is_refine = 0;
    return NULL;
  }
  if (level < refine_level - level_shift) {
    *is_refine = 1;
    return NULL;
  }

  SC3E (p4est3_quadrant_child_id (ri->qvt, ri->quadrant, &child_id));
  *is_refine = (child_id == 0 || child_id == 3
#ifdef P4_TO_P8
          || child_id == 5 || child_id == 6
#endif
    );

  return NULL;
}

void
wrong_input (const char *name, int n)
{
  printf ("Wrong input parameter: \n");
  switch (n) {
  case 1:
    printf ("Quadrant type %s is not valid\n"
            "Valid quadrant type value: P4EST2 or"
            "STANDARD, AVX, MORT_ORD\n", name);
    break;
  case 2:
    printf ("The maximum number of levels %s is not valid\n"
            "Valid value: " "positiv int\n", name);
    break;
  case 3:
    printf ("The number of level shift %s is not valid\n"
            "Valid value: " "positiv int, "
            "less than maximum number of levels \n", name);
  default:
    break;
  }
  printf ("Parameter's format: "
          "<QUADRANT TYPE> <refine_level> <level_shift>\n");
}


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
  if (strcmp (argv[1], "STANDARD") == 0) {
    SC3E (p4est3_quadrant_vtable_p4est (qvt));
  }
  else if (strcmp (argv[1], "AVX") == 0) {
    SC3E (p4est3_quadrant_yx_vtable (qvt));
  }
  else if (strcmp (argv[1], "MORT_ORD") == 0) {
    SC3E (p4est3_quadrant_mort2d_vtable (qvt));
  }
  else {
    if (mpirank == 0) {
      wrong_input (argv[1], 1);
      sc_MPI_Abort (mpicomm, -1);
    }
  }
  SC3E_DEMAND (*qvt != NULL, "AVX is not supported by hardware or"
               "p4est is not build neither in 2D nor 3D");
  return NULL;
}

char *
set_heading (int argc, char **argv, sc3_MPI_Comm_t mpicomm)
{
  char               *heading;
  if (argc > 1) {
    heading = (char *) malloc (strlen (argv[1]) + 1);
    strcpy (heading, argv[1]);
  }
  else {
    heading = (char *) malloc (strlen ("P4EST2") + 1);
    strcpy (heading, "P4EST2");
  }
  if (heading == NULL) {
    sc_MPI_Abort (mpicomm, -1);
  }
  return heading;
}

int
main (int argc, char **argv)
{
  const p4est3_quadrant_vtable_t *qvt;
  sc3_allocator_t    *alloc, *mainalloc;
  sc3_error_t        *e;
  sc3_MPI_Comm_t      mpicomm;
  p4est3_t           *p3;
  p4est_t            *p;
  p4est3_connectivity_t *conn;
  p4est_connectivity_t *conn_old;
  int                 mpirank, mpisize;
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
  if (argc == 1 && mpirank == 0) {
    printf ("Execution without parameters. "
            "Default parameters are applied.\n"
            "Parameter's format: "
            "<QUADRANT TYPE> <refine_level> <level_shift>\n");
  }

  SC3E_NULL_SET (e, check_quadrant_type (argc, argv, &qvt, mpirank, mpicomm));
  if (argc > 2) {
    refine_level = atoi (argv[2]);
    if (refine_level == 0 && mpirank == 0) {
      wrong_input (argv[2], 2);
      sc_MPI_Abort (mpicomm, -1);
    }
  }
  if (argc > 3) {
    level_shift = atoi (argv[3]);
    if (level_shift == 0 && mpirank == 0
        && level_shift >= refine_level) {
      wrong_input (argv[3], 3);
      sc_MPI_Abort (mpicomm, -1);
    }
  }
  sc3_MPI_Barrier (mpicomm);
  heading = set_heading (argc, argv, mpicomm);
}