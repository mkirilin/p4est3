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
#ifndef P4_TO_P8
#include <p4est_extended.h>
#include <p4est_bits.h>
#include <p4est3_p4est.h>
#include <p4est3_quadrant_yx.h>
#include <p4est3_quadrant_mort2d.h>
#include <p4est_vtk.h>
#else
#include <p4est_to_p8est.h>
#include <p8est_extended.h>
#include <p8est_bits.h>
#include <p4est3_p8est.h>
#include <p4est3_quadrant_zyx.h>
#include <p4est3_quadrant_mort3d.h>
#include <p8est_vtk.h>
#endif

/* Generally speaking, it is not allowed to use _internal headers in
 applications. We use it only in exceptional cases as timings and tests. */
#include <p4est3_internal.h>
#include <sc_statistics.h>
#include <sc_options.h>
#include <sc_flops.h>

#include <stdlib.h>
#include <string.h>

static int          refine_level = 1;
#ifdef P4_TO_P8
static double       refinement_fraction = 1. / 7.;
#else
static double       refinement_fraction = 1. / 3.;
#endif

static int
refine_fractal (p4est_t * p, p4est_topidx_t which_tree,
                p4est_quadrant_t * q)
{
  /* Refine every 7th (3d) or 3rd (2d) global quadrant. */
  p4est_locidx_t     *quadrant_local_id = (p4est_locidx_t *) p->user_pointer;

  return
    (((p->global_first_quadrant[p->mpirank] + (*quadrant_local_id)++) %
#ifdef P4_TO_P8
      7)
#else
      3)
#endif
     == 0);
}

static sc3_error_t *
refine_p3_fractal (p4est3_refine_callback_info_t * ri, int *is_refine)
{
  /* Refine every 7th (3d) or 3rd (2d) global quadrant. */
  p4est3_locidx     *quadrant_local_id = (p4est3_locidx *) ri->user_data;
  SC3A_CHECK (is_refine != NULL);

  *is_refine =
    (((ri->p3->goffset[ri->p3->mpirank] + (*quadrant_local_id)++) %
#ifdef P4_TO_P8
      7)
#else
      3)
#endif
     == 0);
  return NULL;
}

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

static int
refine_fraction (p4est_t * p, p4est_topidx_t which_tree,
                 p4est_quadrant_t * q)
{
  /* The formula in the line below implies
   * quadrant_fraction = 7 * refinement_fraction + 1.
   * In particular we need for doubling the quadrants for a level increment
   * refinment_factor = 1 / 7.
  */
  const p4est_locidx_t quad_count_refinement_threshold =
    (p4est_locidx_t) (refinement_fraction * p->global_num_quadrants);
  p4est_locidx_t     *quadrant_local_id = (p4est_locidx_t *) p->user_pointer;

  return
    ((p->global_first_quadrant[p->mpirank] + (*quadrant_local_id)++) <=
      quad_count_refinement_threshold);
}

static sc3_error_t *
refine_p3_fraction (p4est3_refine_callback_info_t * ri, int *is_refine)
{
  /* The formula in the line below implies
   * quadrant_fraction = 7 * refinement_fraction + 1.
   * In particular we need for doubling the quadrants for a level increment
   * refinment_factor = 1 / 7.
   */
  /* p4est->global_num_quadrants is the old number of quadrants */
  SC3A_CHECK (is_refine != NULL);
  const p4est3_locidx quad_count_refinement_threshold =
    (p4est3_locidx) (refinement_fraction * ri->p3->global_num_quads);
  p4est3_locidx     *quadrant_local_id = (p4est3_locidx *) ri->user_data;
  *is_refine =
    ((ri->p3->goffset[ri->p3->mpirank] + (*quadrant_local_id)++) <=
      quad_count_refinement_threshold);

  return NULL;
}

#ifdef P4EST_ENABLE_DEBUG
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
#endif /* P4EST_ENABLE_DEBUG */

static sc3_error_t *
p4est3_new_shortcut (p4est3_t ** p3, sc3_allocator_t *alloc,
                     sc3_MPI_Comm_t mpicomm, p4est3_connectivity_t *conn,
                     const p4est3_quadrant_vtable_t * qvt, int start_level,
                     p4est3_t * src, p4est3_refine_callback_t p3crefine,
                     int is_partition, void *user_data)
{
  SC3E (p4est3_new (alloc, p3));
  SC3E (p4est3_set_comm (*p3, mpicomm, 1));
  SC3E (p4est3_set_connectivity (*p3, conn));
  SC3E (p4est3_set_quadrant_vtable (*p3, qvt));
  SC3E (p4est3_set_level (*p3, start_level));
  SC3E (p4est3_set_setup_mode (*p3, P4EST3_NEW_RECURSIVE));
  SC3E (p4est3_set_refine (*p3, p3crefine));
  SC3E (p4est3_set_source (*p3, src));
  SC3E (p4est3_set_shared (*p3, 1));
  SC3E (p4est3_set_contiguous (*p3, 1));
  SC3E (p4est3_set_partition (*p3, is_partition, NULL));
  /*SC3E (p4est3_set_user_data (*p3, user_data));*/
  (*p3)->user_data = user_data;

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
check_refinement_pattern (const char *opt_pattern,
                          p4est_refine_t *crefine,
                          p4est3_refine_callback_t *p3crefine)
{
  if (strcmp (opt_pattern, "PAIRS") == 0) {
    *crefine = refine_pairs;
    *p3crefine = refine_p3_pairs;
  }
  else if (strcmp (opt_pattern, "FRACTAL") == 0) {
    *crefine = refine_fractal;
    *p3crefine = refine_p3_fractal;
  }
  else if (strcmp (opt_pattern, "FRACTION") == 0) {
    *crefine = refine_fraction;
    *p3crefine = refine_p3_fraction;
  }
  else {
    SC3E_UNREACH ("unavailable refine mode");
  }
  return NULL;
}

static sc3_error_t *
check_quadrant_type (const char *opt_qtype,
                     const p4est3_quadrant_vtable_t **qvt)
{
  if (strcmp (opt_qtype, "P4EST2") == 0) {
    return NULL;
  }
  if (strcmp (opt_qtype, "STANDARD") == 0) {
    SC3E (p4est3_quadrant_vtable_p4est (qvt));
  }
  else if (strcmp (opt_qtype, "AVX") == 0) {
    SC3E (p4est3_quadrant_yx_vtable (qvt));
  }
  else if (strcmp (opt_qtype, "MORT_ORD") == 0) {
    SC3E (p4est3_quadrant_mort2d_vtable (qvt));
  }
  else {
    SC3E_UNREACH
      ("Invalid quadrant implementation name. "
       "Valid names are: P4EST2, STANDARD, AVX, MORT_ORD");
  }
  SC3E_DEMAND (*qvt != NULL, "AVX is not supported by hardware or"
               "p4est is not build neither in 2D nor 3D");
  return NULL;
}

int
main (int argc, char **argv)
{
  const p4est3_quadrant_vtable_t *qvt;
  int i, start_level, write_vtk, num_trees;
  sc3_allocator_t    *alloc, *mainalloc;
  sc3_error_t        *e;
  sc3_MPI_Comm_t      mpicomm = SC3_MPI_COMM_WORLD;
  p4est3_t           *p3, *p3refined;
  p4est_t            *p;
  p4est3_connectivity_t *conn;
  p4est_connectivity_t *conn_old;
  int                 mpirank, mpisize;
  sc_flopinfo_t       fi, snapshot;
  sc_statinfo_t       stats;
  p4est_refine_t crefine = refine_pairs;
  p4est3_refine_callback_t p3crefine = refine_p3_pairs;
  sc_options_t       *opt;
  const char         *opt_pattern, *opt_qtype;
  char heading[80], ref_lvl_string[10];
  p4est3_locidx     quadrant_local_id = 0;


  /* v3 standard procedure to isolate memory allocation contexts */
  mainalloc = sc3_allocator_nothread ();

  /* this is generally needed for MPI */
  SC3E_SET (e, sc3_MPI_Init (&argc, &argv));
  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  SC3E_NULL_SET (e, sc3_MPI_Comm_rank (mpicomm, &mpirank));
  SC3E_NULL_SET (e, sc3_MPI_Comm_size (mpicomm, &mpisize));

  /*** read command line parameters ***/
  opt = sc_options_new (argv[0]);
  sc_options_add_string
    (opt, 'P', "pattern", &opt_pattern, "PAIRS", "Refinement pattern");
  sc_options_add_string
    (opt, 'Q', "qtype", &opt_qtype, "STANDARD", "Quadrant implementation");
  sc_options_add_int
    (opt, 'l', "start level", &start_level, 4,
     "Level of the uniform starting grid");
  sc_options_add_int
    (opt, 'L', "refine level", &refine_level, 10, "Highest level");
  sc_options_add_int
    (opt, 'T', "trees", &num_trees, 1, "Number of trees in a forest");
  sc_options_add_switch (opt, 'V', "write-vtk", &write_vtk,
                         "write vtk output");
  sc_options_parse (p4est_package_id, SC_LP_DEFAULT, opt, argc, argv);

  SC3E_NULL_SET (e, check_refinement_pattern
                    (opt_pattern, &crefine, &p3crefine));
  SC3E_NULL_SET (e, check_quadrant_type (opt_qtype, &qvt));

  /*** set heading ***/
  sprintf(ref_lvl_string, "%d", refine_level);
  strcpy (heading, opt_pattern);
  strcat (heading, " ");
  strcat (heading, opt_qtype);
  strcat (heading, " ");
  strcat (heading, ref_lvl_string);

  /* we don't need init calls for v3.  Just to check legacy wrapping */
  /* must not use SC3_MPI_COMM_WORLD due to incompatible non-mpi wrapping */

  /* make the old style connectivity */
#ifdef P4_TO_P8
  conn_old = p8est_connectivity_new_brick (num_trees, 1, 1, 0, 0, 0);
#else
  conn_old = p4est_connectivity_new_brick (num_trees, 1, 0, 0);
#endif /* P4_TO_P8 */

  if (strcmp (opt_qtype, "P4EST2") != 0) {
    SC3E_NULL_SET (e, make_allocator (mainalloc, &alloc));
    /* make the p4est3 style connectivity */
    SC3E_NULL_SET (e, p4est3_connectivity_new (alloc, &conn));
    SC3E_NULL_SET (e, p4est3_connectivity_set_dim (conn, P4EST_DIM));
    SC3E_NULL_SET (e, p4est3_connectivity_set_num_trees (conn, num_trees));
    SC3E_NULL_SET (e, p4est3_connectivity_setup (conn));
    SC3E_NULL_SET (e, p4est3_new_shortcut
                      (&p3, alloc, mpicomm, conn, qvt, start_level,
                       NULL, NULL, 0, &quadrant_local_id));
    SC3E_NULL_SET (e, p4est3_setup (p3));
#ifdef P4EST_ENABLE_DEBUG
    p = p4est_new_ext
          (mpicomm, conn_old, 0, start_level, 1, 0, NULL, &quadrant_local_id);
#endif
 
    /*** refine in a loop ***/
    for (i = 0; i < refine_level; ++i, quadrant_local_id = 0) {
      SC3E_NULL_SET (e, p4est3_new_shortcut
                        (&p3refined, alloc, mpicomm, conn, qvt, 0,
                         p3, p3crefine, 0, &quadrant_local_id));
      SC3E_NULL_SET (e, p4est3_setup (p3refined));
      SC3E_NULL_SET (e, p4est3_destroy (&p3));
      SC3E_NULL_SET (e, p4est3_new_shortcut
                        (&p3, alloc, mpicomm, conn, qvt, 0,
                         p3refined, NULL, 1, &quadrant_local_id));
      if (i == refine_level - 1) {
        SC3E_NULL_SET (e, sc3_MPI_Barrier (mpicomm));
        sc_flops_snap (&fi, &snapshot);
        SC3E_NULL_SET (e, p4est3_setup (p3));
        sc_flops_shot (&fi, &snapshot);
        sc_stats_set1 (&stats, snapshot.iwtime, heading);
        sc_stats_compute (mpicomm, 1, &stats);
        sc_stats_print (p4est_package_id, SC_LP_ESSENTIAL, 1, &stats, 1, 1);
      }
      else {
        SC3E_NULL_SET (e, p4est3_setup (p3));
      }
#ifdef P4EST_ENABLE_DEBUG
      quadrant_local_id = 0;
      p4est_refine (p, 0, crefine, NULL);
      SC3E_NULL_SET (e, compare_results (alloc, p3, p, qvt));
#endif
      SC3E_NULL_SET (e, p4est3_destroy (&p3refined));
    }
#ifdef P4EST_ENABLE_DEBUG
    p4est_destroy (p);
#endif
    SC3E_NULL_SET (e, p4est3_destroy (&p3));
    SC3E_NULL_SET (e, p4est3_connectivity_destroy (&conn));
    SC3E_NULL_SET (e, sc3_allocator_destroy (&alloc));
  }
  else if ((strcmp (opt_qtype, "P4EST2") == 0) || write_vtk) {
    p = p4est_new_ext
          (mpicomm, conn_old, 0, start_level, 1, 0, NULL, &quadrant_local_id);
    for (i = 0; i < refine_level; ++i) {
      p4est_refine (p, 0, crefine, NULL);
      if ((i == refine_level - 1) && !write_vtk) {
        SC3E_NULL_SET (e, sc3_MPI_Barrier (mpicomm));
        sc_flops_snap (&fi, &snapshot);
        p4est_partition (p, 0, NULL);
        sc_flops_shot (&fi, &snapshot);
        sc_stats_set1 (&stats, snapshot.iwtime, heading);
        sc_stats_compute (mpicomm, 1, &stats);
        sc_stats_print (p4est_package_id, SC_LP_ESSENTIAL, 1, &stats, 1, 1);
      }
      else {
        p4est_partition (p, 0, NULL);
      }
    }
    if (write_vtk) {
      p4est_vtk_write_file (p, NULL, opt_pattern);
    }
    p4est_destroy (p);
  }
  p4est_connectivity_destroy (conn_old);
  sc_options_destroy (opt);
  SC3E_NULL_REQ (e, !sc_finalize_noabort ());
  SC3E_NULL_SET (e, sc3_MPI_Finalize ());
  SC3X (e);
  return 0;
}