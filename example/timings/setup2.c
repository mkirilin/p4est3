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
#include <p4est_extended.h>
#else
#include <p8est3_quadrant_zyx.h>
#include <p8est3_quadrant_mort.h>
#include <p8est_p4est3.h>
#include <p8est_extended.h>
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
            "P4EST2, "
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

void
check_setup_mode (int argc, char **argv, p4est3_setup_mode_t *mode,
                  int mpirank, sc3_MPI_Comm_t mpicomm)
{
  if (argc <= 1) {
    return;
  }
  if (strcmp (argv[1], "MORTON") == 0) {
    *mode = P4EST3_NEW_MORTON;
  }
  else if (strcmp (argv[1], "SUCCESSOR") == 0) {
    *mode = P4EST3_NEW_SUCCESSOR;
  }
  else if (strcmp (argv[1], "RECURSIVE") == 0) {
    *mode = P4EST3_NEW_RECURSIVE;
  }
  else if (strcmp (argv[1], "RECURSIVE_CHILD") == 0) {
    *mode = P4EST3_NEW_RECURSIVE_CHILD;
  }
  else if (strcmp (argv[1], "RECURSIVE_REGION") == 0) {
    *mode = P4EST3_NEW_RECURSIVE_REGION;
  }
  else if (strcmp (argv[1], "P4EST2") == 0) {
    *mode = P4EST3_NEW_MODE_LAST;
  }
  else {
    if (mpirank == 0) {
      wrong_input (argv[1], 1);
      sc_MPI_Abort (mpicomm, -1);
    }
  }
}

void
check_quadrant_type (int argc, char **argv,
                     p4est3_quadrant_vtable_t *qvt,
                     int mpirank, sc3_MPI_Comm_t mpicomm)
{
  if (argc <= 2) {
    return;
  }
  if (strcmp (argv[1], "P4EST2") == 0) {
    return;
  }
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

char *
set_heading (int argc, char **argv, sc3_MPI_Comm_t mpicomm)
{
  char               *heading;
  if (argc > 2) {
    heading = (char *) malloc (strlen (argv[1]) + 1 + strlen (argv[2]) + 1);
    strcpy (heading, argv[1]);
    strcat (heading, " ");
    strcat (heading, argv[2]);
  }
  else if (argc == 2) {
    heading =
      (char *) malloc (strlen (argv[1]) + 1 + strlen ("STANDART") + 1);
    strcpy (heading, argv[1]);
    strcat (heading, " ");
    strcat (heading, "STANDART");
  }
  else {
    heading =
      (char *) malloc (strlen ("MORTON") + 1 + strlen ("STANDART") + 1);
    strcpy (heading, "MORTON");
    strcat (heading, " ");
    strcat (heading, "STANDART");
  }
  if (heading == NULL) {
    sc_MPI_Abort (mpicomm, -1);
  }
  return heading;
}
#ifndef P4_TO_P8
p4est_connectivity_t *
p4est_connectivity_new_row (int trees)
{
  const p4est_topidx_t num_vertices = trees * 2 + 2;
  const p4est_topidx_t num_trees = trees;
  const p4est_topidx_t num_ctt = 0;
  double              *vertices;
  int                  i;
  vertices = (double *) malloc (sizeof (double) * num_vertices * 3);
  for (i = 0; i < num_trees; ++i) {
    vertices[i * 2 * 3] = i;
    vertices[i * 2 * 3 + 1] = 0;
    vertices[i * 2 * 3 + 2] = 0;

    vertices[i * 2 * 3 + 3] = i;
    vertices[i * 2 * 3 + 4] = 1;
    vertices[i * 2 * 3 + 5] = 0;
  }
  vertices[num_vertices * 3 - 6] = num_trees;
  vertices[num_vertices * 3 - 5] = 0;
  vertices[num_vertices * 3 - 4] = 0;

  vertices[num_vertices * 3 - 3] = num_trees;
  vertices[num_vertices * 3 - 2] = 1;
  vertices[num_vertices * 3 - 1] = 0;

  p4est_topidx_t *tree_to_vertex;
  tree_to_vertex =
    (p4est_topidx_t *) malloc (sizeof (p4est_topidx_t) * num_trees * 4);
  for (i = 0; i < num_trees; ++i) {
    tree_to_vertex[i * 4] = i * 2;
    tree_to_vertex[i * 4 + 1] = i * 2 + 2;
    tree_to_vertex[i * 4 + 2] = i * 2 + 1;
    tree_to_vertex[i * 4 + 3] = i * 2 + 3;
  }

  p4est_topidx_t *tree_to_tree;
  tree_to_tree =
    (p4est_topidx_t *) malloc (sizeof (p4est_topidx_t) * num_trees * 4);
  if (num_trees == 1) {
    tree_to_tree[0] = 0;
    tree_to_tree[1] = 0;
    tree_to_tree[2] = 0;
    tree_to_tree[3] = 0;
  }
  else {
    tree_to_tree[0] = 0;
    tree_to_tree[1] = 1;
    tree_to_tree[2] = 0;
    tree_to_tree[3] = 0;
    for (i = 1; i < num_trees - 1; ++i){
      tree_to_tree[i * 4] = i - 1;
      tree_to_tree[i * 4 + 1] = i + 1;
      tree_to_tree[i * 4 + 2] = i;
      tree_to_tree[i * 4 + 3] = i;
    }
    tree_to_tree[num_trees * 4 - 4] = num_trees - 2;
    tree_to_tree[num_trees * 4 - 3] = num_trees - 1;
    tree_to_tree[num_trees * 4 - 2] = num_trees - 1;
    tree_to_tree[num_trees * 4 - 1] = num_trees - 1;
  }

  int8_t *tree_to_face;
  tree_to_face =
    (int8_t *) malloc (sizeof (int8_t) * num_trees * 4);
  if (num_trees == 1) {
    tree_to_face[0] = 0;
    tree_to_face[1] = 1;
    tree_to_face[2] = 2;
    tree_to_face[3] = 3;
  }
  else {
    tree_to_face[0] = 0;
    tree_to_face[1] = 0;
    tree_to_face[2] = 2;
    tree_to_face[3] = 3;
    for (i = 1; i < num_trees - 1; ++i){
      tree_to_face[i * 4] = 1;
      tree_to_face[i * 4 + 1] = 0;
      tree_to_face[i * 4 + 2] = 2;
      tree_to_face[i * 4 + 3] = 3;
    }
    tree_to_face[num_trees * 4 - 4] = 1;
    tree_to_face[num_trees * 4 - 3] = 1;
    tree_to_face[num_trees * 4 - 2] = 2;
    tree_to_face[num_trees * 4 - 1] = 3;
  }
  return p4est_connectivity_new_copy (num_vertices, num_trees, 0,
                                      vertices, tree_to_vertex,
                                      tree_to_tree, tree_to_face,
                                      NULL, &num_ctt, NULL, NULL);
}
#endif /* !P4_TO_P8 */

int
main (int argc, char **argv)
{
  p4est3_topidx       num_trees;
  sc3_allocator_t    *alloc, *mainalloc;
  sc3_error_t        *e;
  sc3_MPI_Comm_t      mpicomm;
  p4est3_quadrant_vtable_t vtable, *qvt = &vtable;
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
  mpicomm = SC3_MPI_COMM_WORLD;

  /* this is generally needed for MPI */
  SC3E_SET (e, sc3_MPI_Init (&argc, &argv));

  /* default parameters */
  p4est3_setup_mode_t mode = P4EST3_NEW_MORTON;
  p4est_quadrant_vtable (qvt, 0);
  level = 1;
  num_trees = 2;

  SC3E_NULL_SET (e, sc3_MPI_Comm_rank (mpicomm, &mpirank));
  SC3E_NULL_SET (e, sc3_MPI_Comm_size (mpicomm, &mpisize));

  if (argc == 1 && mpirank == 0) {
    printf ("Execution without parameters. "
            "Default parameters are applied.\n"
            "Parameter's format: "
            "<SETUP MODE> <QUADRANT TYPE> <#levels> <#trees>\n");
  }
  check_setup_mode (argc, argv, &mode, mpirank, mpicomm);
  check_quadrant_type (argc, argv, qvt, mpirank, mpicomm);
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
  heading = set_heading (argc, argv, mpicomm);

  /* we don't need init calls for v3.  Just to check legacy wrapping */
  /* must not use SC3_MPI_COMM_WORLD due to incompatible non-mpi wrapping */
  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);
  SC3E_NULL_SET (e, make_allocator (mainalloc, &alloc));

  /* make the old style connectivity */
#ifdef P4_TO_P8
  conn_old = p8est_connectivity_new_unitcube ();
#else
  conn_old = p4est_connectivity_new_row (num_trees);
#endif /* P4_TO_P8 */
  /* make the p4est3 style connectivity */
  SC3E_NULL_SET (e, p4est3_connectivity_new (alloc, &conn));
  SC3E_NULL_SET (e, p4est3_connectivity_set_num_trees (conn, num_trees));
  SC3E_NULL_SET (e, p4est3_connectivity_setup (conn));

  /* create p4est object with connectivity */
  if (mode != P4EST3_NEW_MODE_LAST) {
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

    SC3E_NULL_SET (e, p4est3_destroy (&p3));
  }
  else {
    sc_flops_snap (&fi, &snapshot);
    p = p4est_new_ext (mpicomm, conn_old, 0, level, 1, 0, NULL, NULL);
    sc_flops_shot (&fi, &snapshot);
    sc_stats_set1 (&stats, snapshot.iwtime, heading);

    sc_stats_compute (mpicomm, 1, &stats);
    sc_stats_print (p4est_package_id, SC_LP_ESSENTIAL, 1, &stats, 1, 1);

    p4est_destroy (p);
  }

  free (heading);
  SC3E_NULL_SET (e, p4est3_connectivity_destroy (&conn));
  p4est_connectivity_destroy (conn_old);
  SC3E_NULL_SET (e, free_allocator (&alloc));

  /* again, just to check legacy wrapping */
  SC3E_NULL_REQ (e, !sc_finalize_noabort ());

  /* TODO: call finalize even with errors? */
  SC3E_NULL_SET (e, sc3_MPI_Finalize ());
  report_errors (mainalloc, &e);

  return 0;
}
