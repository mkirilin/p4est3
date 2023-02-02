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

static int          refine_level = 1;

static int
refine_pairs (p4est_t * p4est, p4est_topidx_t which_tree,
              p4est_quadrant_t * q)
{
  if ((which_tree & 1) == 0) {
    if (q->level < refine_level) {
      return 1;
    }
  }
  else {
    if (q->level < 1) {
      return 1;
    }
  }
  return 0;
}

static sc3_error_t *
refine_p3_pairs (p4est3_refine_callback_info_t * ri, int *is_refine)
{
  int                 level;
  *is_refine = 0;

  SC3E (p4est3_quadrant_level (ri->qvt, ri->quadrant, &level));
  if ((ri->ntree & 1) == 0) {
    if (level < refine_level) {
      *is_refine = 1;
      return NULL;
    }
  }
  else {
    if (level < 1) {
      *is_refine = 1;
      return NULL;
    }
  }
  return NULL;
}

static sc3_error_t *
array_new (sc3_allocator_t * alloc, size_t esize, int ealloc,
           int ecount, sc3_array_t ** arr)
{
  SC3E_RETVAL (arr, NULL);
  SC3A_IS (sc3_allocator_is_setup, alloc);
  SC3A_CHECK (ealloc >= 0);

  SC3E (sc3_array_new (alloc, arr));
  SC3E (sc3_array_set_elem_size (*arr, esize));
  SC3E (sc3_array_set_elem_alloc (*arr, ealloc));
  SC3E (sc3_array_set_elem_count (*arr, ecount));
  SC3E (sc3_array_set_initzero (*arr, 1));
  SC3E (sc3_array_setup (*arr));

  return NULL;
}

static sc3_error_t *
compare_results (sc3_allocator_t *alloc, p4est3_t * p3, p4est_t * p,
                 const p4est3_quadrant_vtable_t * qvt)
{

  char               *q3;
  int                *level, *p3level;
  p4est_quadrant_t   *q;
  p4est3_topidx       i, tt;
  size_t              nq;
  p4est3_topidx       fltree, lltree;
  p4est3_gloidx       num_glo_quads;
  p4est3_locidx       num_loc_quads;
  p4est3_locidx       processed_quads, processed_quads_p3;
  p4est_tree_t       *tree;
  sc3_array_t        *p3levels, *levels;

  SC3E (p4est3_get_local_num_trees (p3, &fltree, &lltree));
  SC3E (p4est3_get_global_num_quads (p3, &num_glo_quads));
  SC3E (p4est3_get_local_num_quads (p3, &num_loc_quads));
  SC3E_DEMAND (fltree == p->first_local_tree && lltree == p->last_local_tree,
               "Different trees at processor");
  SC3E_DEMAND (num_glo_quads == p->global_num_quadrants,
               "different #global quadrants");
  SC3E_DEMAND (num_loc_quads == p->local_num_quadrants,
               "different #local quadrants");
  SC3E (array_new (alloc, sizeof (int), num_loc_quads, 0, &p3levels));
  SC3E (array_new
        (alloc, sizeof (int), p->local_num_quadrants, 0, &levels));
  for (tt = p->first_local_tree; tt <= p->last_local_tree; ++tt) {
    tree = p4est_tree_array_index (p->trees, tt);
    for (nq = 0; nq < tree->quadrants.elem_count; ++nq) {
      q = (p4est_quadrant_t *) sc_array_index (&tree->quadrants, nq);
      SC3E (sc3_array_push (levels, &level));
      *level = q->level;
    }
  }

  SC3E (p4est3_get_quadrants (p3, &q3));
  if (num_loc_quads > 0) {
    SC3E (sc3_array_push (p3levels, &level));
    SC3E (p4est3_quadrant_level (qvt, q3, level));
  }
  for (i = 1; i < num_loc_quads; ++i) {
    q3 += qvt->quadrant_size;
    SC3E (sc3_array_push (p3levels, &level));
    SC3E (p4est3_quadrant_level (qvt, q3, level));
  }
  SC3E (sc3_array_get_elem_count (p3levels, &processed_quads_p3));
  SC3E (sc3_array_get_elem_count (levels, &processed_quads));
  SC3E_DEMAND (processed_quads_p3 == processed_quads, "wrong #p3levels");
  SC3E_DEMAND (processed_quads == num_loc_quads, "wrong #levels");
  for (i = 0; i < num_loc_quads; ++i) {
    SC3E (sc3_array_index (p3levels, i, &p3level));
    SC3E (sc3_array_index (levels, i, &level));
    SC3E_DEMAND (*level == *p3level, "levels mismatch");
  }

  SC3E (sc3_array_destroy (&p3levels));
  SC3E (sc3_array_destroy (&levels));
  return NULL;
}

static sc3_error_t *
refine (sc3_allocator_t *alloc, p4est3_t ** p3,
        const p4est3_quadrant_vtable_t * qvt,
        sc3_MPI_Comm_t mpicomm, p4est_connectivity_t *conn_old)
{
  int                 i;
  p4est_t * p;
  p4est3_t           *p3refined, *p3ptr = *p3;
  /* refine the old forest */
  p = p4est_new_ext (mpicomm, conn_old, 0, 1, 1, 0, NULL, NULL);
  p4est_refine (p, 1, refine_pairs, NULL);

  for (i = 0; i < refine_level; ++i) {
    SC3E (p4est3_new (alloc, &p3refined));
    SC3E (p4est3_set_quadrant_vtable (p3refined, qvt));
    SC3E (p4est3_set_refine (p3refined, refine_p3_pairs));
    SC3E (p4est3_set_source (p3refined, p3ptr));
    SC3E (p4est3_setup (p3refined));

    if (i != 0) {
      SC3E (p4est3_destroy (&p3ptr));
    }
    p3ptr = p3refined;
  }
  SC3E (compare_results (alloc, p3ptr, p, qvt));

  SC3E (p4est3_destroy (p3));
  *p3 = p3refined;
  p4est_destroy (p);
  return NULL;
}

static sc3_error_t *
partition (sc3_allocator_t *alloc, p4est3_t ** p3)
{
  p4est3_t *p3part, *p3ptr = *p3;

  SC3E (p4est3_new (alloc, &p3part));
  SC3E (p4est3_set_quadrant_vtable (p3part, p3ptr->qvt));
  SC3E (p4est3_set_source (p3part, p3ptr));
  SC3E (p4est3_set_partition (p3part, 1, NULL));
  SC3E (p4est3_set_shared (p3part, 1));
  SC3E (p4est3_set_contiguous (p3part, 0));
  SC3E (p4est3_setup (p3part));
  SC3E (p4est3_destroy (&p3ptr));
  *p3 = p3part;
  return NULL;
}

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
  default:
    break;
  }
  printf ("Parameter's format: "
          "<QUADRANT TYPE> <refine_level>\n");
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
  sc3_MPI_Comm_t      mpicomm = SC3_MPI_COMM_WORLD;
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
  SC3E_NULL_SET (e, p4est3_quadrant_vtable_p4est (&qvt));
  if (argc == 1 && mpirank == 0) {
    printf ("Execution without parameters. "
            "Default parameters are applied.\n"
            "Parameter's format: "
            "<QUADRANT TYPE> <refine_level>\n");
  }

  SC3E_NULL_SET (e, check_quadrant_type (argc, argv, &qvt, mpirank, mpicomm));
  if (argc > 2) {
    refine_level = atoi (argv[2]);
    if (refine_level == 0 && mpirank == 0) {
      wrong_input (argv[2], 2);
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
  conn_old = p8est_connectivity_new_brick (mpisize, 1, 1, 0, 0, 0);
#else
  conn_old = p4est_connectivity_new_brick (mpisize, 1, 0, 0);
#endif /* P4_TO_P8 */
  /* make the p4est3 style connectivity */
  SC3E_NULL_SET (e, p4est3_connectivity_new (alloc, &conn));
  SC3E_NULL_SET (e, p4est3_connectivity_set_dim (conn, P4EST_DIM));
  SC3E_NULL_SET (e, p4est3_connectivity_set_num_trees (conn, mpisize));
  SC3E_NULL_SET (e, p4est3_connectivity_setup (conn));

    /* create p4est object with connectivity */
  if (strcmp (argv[1], "P4EST2") != 0) {
    SC3E_NULL_SET (e, p4est3_new (alloc, &p3));
    SC3E_NULL_SET (e, p4est3_set_comm (p3, mpicomm, 1));
    SC3E_NULL_SET (e, p4est3_set_connectivity (p3, conn));
    SC3E_NULL_SET (e, p4est3_set_quadrant_vtable (p3, qvt));
    SC3E_NULL_SET (e, p4est3_set_level (p3, 1));
    SC3E_NULL_SET (e, p4est3_set_setup_mode (p3, P4EST3_NEW_RECURSIVE));
    SC3E_NULL_SET (e, p4est3_set_shared (p3, 0));

    SC3E_NULL_SET (e, p4est3_setup (p3));
    SC3E_NULL_SET (e, refine (alloc, &p3, qvt, mpicomm, conn_old));

    sc_flops_snap (&fi, &snapshot);
    SC3E_NULL_SET (e, partition (alloc, &p3));
    sc_flops_shot (&fi, &snapshot);
    sc_stats_set1 (&stats, snapshot.iwtime, heading);

    sc_stats_compute (mpicomm, 1, &stats);
    sc_stats_print (p4est_package_id, SC_LP_ESSENTIAL, 1, &stats, 1, 1);

    SC3E_NULL_SET (e, p4est3_destroy (&p3));
    SC3E_NULL_SET (e, p4est3_connectivity_destroy (&conn));
    SC3E_NULL_SET (e, sc3_allocator_destroy (&alloc));
  }
  else {
    p = p4est_new_ext (mpicomm, conn_old, 0, 1, 1, 0, NULL, NULL);
    p4est_refine (p, 1, refine_pairs, NULL);
    sc_flops_snap (&fi, &snapshot);
    p4est_partition (p, 0, NULL);
    sc_flops_shot (&fi, &snapshot);
    sc_stats_set1 (&stats, snapshot.iwtime, heading);

    sc_stats_compute (mpicomm, 1, &stats);
    sc_stats_print (p4est_package_id, SC_LP_ESSENTIAL, 1, &stats, 1, 1);

    p4est_destroy (p);
  }
  SC3E_NULL_SET (e, sc3_MPI_Finalize ());
  SC3X (e);
}