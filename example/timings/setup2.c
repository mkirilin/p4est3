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
#include <p4est3_quadrant_mort.h>
#include <p4est_p4est3.h>
#else
#include <p8est3_quadrant_zyx.h>
#include <p8est3_quadrant_mort.h>
#include <p8est_p4est3.h>
#endif

#include <sc_statistics.h>
#include <sc_flops.h>

#include <string.h>

static sc3_error_t *
make_allocator (sc3_allocator_t * oa, sc3_allocator_t ** alloc)
{
  SC3A_IS (sc3_allocator_is_setup, oa);
  SC3E (sc3_allocator_new (oa, alloc));
  SC3E (sc3_allocator_setup (*alloc));
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
    SC_CHECK_ABORT (0, "Setup's timings failed\n");
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

void
wrong_input (const char *name, int n)
{
  printf ("Wrong input parameter: \n");
  switch (n) {
  case 1:
    printf ("Setup mode %s is not valid\n"
            "Valid setup parameters: "
            "MORTON, SUCCESSOR, RECURSIVE, RECURSIVE_CHILD, RECURSIVE_REGION"
            "\n", name);
    break;
  case 2:
    printf ("Quadrant type %s is not valid\n"
            "Valid quadrant types: " "STANDART, AVX, MORT_ORD\n", name);
    break;
  case 3:
    printf ("The number of levels %s is not valid\n"
            "Valid value: " "positiv int\n", name);
    break;
  case 4:
    printf ("The number of trees %s is not valid\n"
            "Valid value: " "positiv int\n", name);
  default:
    break;
  }
  printf ("Parameter's format: "
          "<SETUP MODE> <QUADRANT TYPE> <#levels> <#trees>\n");
}

int
main (int argc, char **argv)
{
  p4est3_topidx       num_trees;
  sc3_allocator_t    *alloc, *mainalloc;
  sc3_error_t        *e;
  sc3_MPI_Comm_t      mpicomm;
  p4est3_quadrant_vtable_t vtable, *qvt = &vtable;
  p4est3_t           *p3;
  p4est3_connectivity_t *conn;
  int                 level, mpirank, mpisize;
  sc_flopinfo_t       fi, snapshot;
  sc_statinfo_t       stats;

  /* v3 standard procedure to isolate memory allocation contexts */
  mainalloc = sc3_allocator_nothread ();
  mpicomm = SC3_MPI_COMM_WORLD;

  /* this is generally needed for MPI */
  SC3E_SET (e, sc3_MPI_Init (&argc, &argv));

  /* default parameters */
  p4est3_setup_mode_t mode = P4EST3_NEW_MORTON;
  p4est_quadrant_vtable (qvt, 0);
  level = 1;
  num_trees = 1;

  SC3E_NULL_SET (e, sc3_MPI_Comm_rank (mpicomm, &mpirank));
  SC3E_NULL_SET (e, sc3_MPI_Comm_size (mpicomm, &mpisize));

  if (argc == 1 && mpirank == 0) {
    printf ("Execution without parameters. "
            "Default parameters are applied.\n"
            "Parameter's format: "
            "<SETUP MODE> <QUADRANT TYPE> <#levels> <#trees>\n");
  }
  if (argc > 1) {
    if (strcmp (argv[1], "MORTON") == 0) {
      mode = P4EST3_NEW_MORTON;
    }
    else if (strcmp (argv[1], "SUCCESSOR") == 0) {
      mode = P4EST3_NEW_SUCCESSOR;
    }
    else if (strcmp (argv[1], "RECURSIVE") == 0) {
      mode = P4EST3_NEW_RECURSIVE;
    }
    else if (strcmp (argv[1], "RECURSIVE_CHILD") == 0) {
      mode = P4EST3_NEW_RECURSIVE_CHILD;
    }
    else if (strcmp (argv[1], "RECURSIVE_REGION") == 0) {
      mode = P4EST3_NEW_RECURSIVE_REGION;
    }
    else {
      if (mpirank == 0) {
        wrong_input (argv[1], 1);
        sc_MPI_Abort (mpicomm, -1);
      }
    }
  }
  if (argc > 2) {
    if (strcmp (argv[2], "STANDART") == 0) {
      p4est_quadrant_vtable (qvt, 0);
    }
    else if (strcmp (argv[2], "AVX") == 0) {
      p4est3_quadrant_zyx_vtable (qvt);
    }
    else if (strcmp (argv[2], "MORT_ORD") == 0) {
      p4est3_quadrant_mort_vtable (qvt);
    }
    else {
      if (mpirank == 0) {
        wrong_input (argv[2], 2);
        sc_MPI_Abort (mpicomm, -1);
      }
    }
  }
  if (argc > 3) {
    level = atoi (argv[3]);
    if (level == 0 && mpirank == 0) {
      wrong_input (argv[3], 3);
      sc_MPI_Abort (mpicomm, -1);
    }
  }
  if (argc > 4) {
    num_trees = atoi (argv[4]);
    if (num_trees == 0 && mpirank == 0) {
      wrong_input (argv[4], 4);
      sc_MPI_Abort (mpicomm, -1);
    }
  }
  sc3_MPI_Barrier (mpicomm);

  char               *heading;
  if (argc > 2) {
    heading = (char *) malloc (strlen (argv[1]) + 1 + strlen (argv[2]) + 1);
    strcpy (heading, argv[1]);
    strcat (heading, " ");
    strcat (heading, argv[2]);
  }
  else if (argc == 2) {
    heading = (char *) malloc (strlen (argv[1]) + 1 + strlen ("STANDART") + 1);
    strcpy (heading, argv[1]);
    strcat (heading, " ");
    strcat (heading, "STANDART");
  }
  else {
    heading = (char *) malloc (strlen ("MORTON") + 1 + strlen ("STANDART") + 1);
    strcpy (heading, "MORTON");
    strcat (heading, " ");
    strcat (heading, "STANDART");
  }
  if (heading == NULL) {
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

  /* create p4est object with connectivity */
  SC3E_NULL_SET (e, p4est3_new (alloc, &p3));
  SC3E_NULL_SET (e, p4est3_set_comm (p3, mpicomm, 1));
  SC3E_NULL_SET (e, p4est3_set_connectivity (p3, conn));
  SC3E_NULL_SET (e, p4est3_set_vtable (p3, qvt));
  SC3E_NULL_SET (e, p4est3_set_level (p3, level));
  SC3E_NULL_SET (e, p4est3_set_setup_mode (p3, mode));

  sc_flops_snap (&fi, &snapshot);
  SC3E_NULL_SET (e, p4est3_setup (p3));
  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats, snapshot.iwtime, heading);

  sc_stats_compute (mpicomm, 1, &stats);
  sc_stats_print (p4est_package_id, SC_LP_ESSENTIAL, 1, &stats, 1, 1);

  free (heading);
  SC3E_NULL_SET (e, p4est3_destroy (&p3));

  SC3E_NULL_SET (e, p4est3_connectivity_destroy (&conn));
  SC3E_NULL_SET (e, free_allocator (&alloc));

  /* again, just to check legacy wrapping */
  SC3E_NULL_REQ (e, !sc_finalize_noabort ());

  /* TODO: call finalize even with errors? */
  SC3E_NULL_SET (e, sc3_MPI_Finalize ());
  report_errors (mainalloc, &e);

  return 0;
}
