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
#include <p4est3_internal.h>
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

#define MAX_TEST_LEVEL 5

typedef struct setup
{
  sc3_allocator_t    *alloc;
  sc3_allocator_t    *mainalloc;
  p4est_connectivity_t *conn2;
  p4est3_connectivity_t *conn3;
  sc3_MPI_Comm_t      mpicomm;
  int                 mpirank;
  int                 mpisize;
  int                 level;
}
setup_t;

static int
refine_fn (p4est_t * p4est, p4est_topidx_t which_tree,
           p4est_quadrant_t * quadrant)
{
  if ((int) quadrant->level >= p4est->mpirank) {
    return 0;
  }
  return 1;
}

static sc3_error_t *
refine_p3_fn (p4est3_refine_callback_info_t * ri, int *is_refine)
{
  int                 level;
  *is_refine = 1;

  SC3E (p4est3_quadrant_level (ri->qvt, ri->quadrant, &level));
  if (level >= ri->p3->mpirank) {
    *is_refine = 0;
    return NULL;
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
make_allocator (setup_t * t)
{
  SC3A_IS (sc3_allocator_is_setup, t->mainalloc);
  SC3E (sc3_allocator_new (t->mainalloc, &t->alloc));
  SC3E (sc3_allocator_setup (t->alloc));

  return NULL;
}

static sc3_error_t *
make_connectivity (setup_t * t)
{
  t->conn2 =
#ifdef P4_TO_P8
    p8est_connectivity_new_brick (t->mpisize, 1, 1, 0, 0, 0);
#else
    p4est_connectivity_new_brick (t->mpisize, 1, 0, 0);
#endif
  SC3E (p4est3_connectivity_new_p4est (t->alloc, &t->conn3, t->conn2, 1));
  return NULL;
}

static sc3_error_t *
make_new_p4est (p4est_t ** p, setup_t * t)
{
  *p = p4est_new_ext
    (t->mpicomm, t->conn2, 0, 0, 1, 0, NULL, NULL);

  return NULL;
}

static sc3_error_t *
make_new_p4est3 (p4est3_t ** p3, setup_t * t,
                 const p4est3_quadrant_vtable_t ** qvt)
{
  SC3A_IS (sc3_allocator_is_setup, t->alloc);

  /* create p4est object with connectivity */
  SC3E (p4est3_new (t->alloc, p3));
  SC3E (p4est3_set_comm (*p3, t->mpicomm, 1));
  SC3E (p4est3_set_connectivity (*p3, t->conn3));
  SC3E (p4est3_set_quadrant_vtable (*p3, *qvt));
  SC3E (p4est3_set_level (*p3, 0));
  SC3E (p4est3_set_shared (*p3, 1));
  SC3E (p4est3_set_contiguous (*p3, 1));
  SC3E (p4est3_setup (*p3));

  return NULL;
}

static sc3_error_t *
prepare_objects (p4est3_t ** p3, p4est_t ** p, setup_t * t,
                 const p4est3_quadrant_vtable_t ** qvt)
{
  t->mainalloc = sc3_allocator_nothread ();
  SC3E (make_allocator (t));
  SC3E (make_connectivity (t));
  SC3E (p4est3_quadrant_vtable_p4est (qvt));
  SC3E_DEMAND (*qvt != NULL, "p4est is not build neither in 2D nor 3D");
  SC3E (make_new_p4est (p, t));
  SC3E (make_new_p4est3 (p3, t, qvt));
  return NULL;
}

static sc3_error_t *
clean_up (p4est3_t * p3, p4est_t * p, setup_t * t)
{
  /*destroy forest, that was referenced for others */
  SC3E (p4est3_destroy (&p3));
  SC3E (p4est3_connectivity_destroy (&t->conn3));
  p4est_destroy (p);
  /* There is no need to destroy p4est2 connectivity,
      since it is destroyed at p4est3 conn destroying stage */
  /* p4est_connectivity_destroy (t->conn2); */
  SC3A_IS (sc3_allocator_is_setup, t->alloc);
  SC3E (sc3_allocator_destroy (&t->alloc));
  return NULL;
}

static sc3_error_t *
compare_results (setup_t * t, p4est3_t * p3, p4est_t * p,
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
  SC3E (array_new
        (t->alloc, qvt->quadrant_size, num_loc_quads, 0, &p3levels));
  SC3E (array_new
        (t->alloc, qvt->quadrant_size, p->local_num_quadrants, 0, &levels));
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

/**
 * 0 - Classic
 * 1 - AVX
 * 2 - Morton
*/
static sc3_error_t *
set_qvt (const p4est3_quadrant_vtable_t ** qvt, int i)
{
  SC3A_CHECK (0 <= i && i <= 2);
  switch (i) {
  case 0:
    SC3E (p4est3_quadrant_vtable_p4est (qvt));
    break;
  case 1:
    SC3E (p4est3_quadrant_yx_vtable (qvt));
    break;
  case 2:
    SC3E (p4est3_quadrant_mort2d_vtable (qvt));
    break;
  default:
    SC3E_UNREACH ("wrong qvt mode");
  }
  SC3E_DEMAND (*qvt != NULL, "AVX is not supported by hardware "
               "p4est is not build neither in 2D nor 3D");
  return NULL;
}

static sc3_error_t *
perform_test (p4est3_t * p3, p4est_t * p, setup_t * t,
              const p4est3_quadrant_vtable_t * qvt)
{
  const int refine_level = SC3_MIN (MAX_TEST_LEVEL, p3->mpisize);
  int i;
  p4est3_t           *p3refined, *p3ptr = p3;
  p4est_refine (p, 1, refine_fn, NULL);

  for (i = 0; i < refine_level; ++i) {
    SC3E (p4est3_new (t->alloc, &p3refined));
    SC3E (set_qvt (&qvt, i % 1));
    SC3E (p4est3_set_quadrant_vtable (p3refined, qvt));
    SC3E (p4est3_set_refine (p3refined, refine_p3_fn));
    SC3E (p4est3_set_source (p3refined, p3ptr));
    SC3E (p4est3_setup (p3refined));

    if (i != 0) {
      SC3E (p4est3_destroy (&p3ptr));
    }
    p3ptr = p3refined;
  }
  SC3E (compare_results (t, p3ptr, p, qvt));

  /* Test partition */
  SC3E (p4est3_new (t->alloc, &p3refined));
  SC3E (set_qvt (&qvt, 0));
  SC3E (p4est3_set_quadrant_vtable (p3refined, qvt));
  SC3E (p4est3_set_source (p3refined, p3ptr));
  SC3E (p4est3_set_partition (p3refined, 1));
  SC3E (p4est3_set_shared (p3refined, 1));
  SC3E (p4est3_setup (p3refined));

  SC3E (p4est3_destroy (&p3ptr));
  p3ptr = p3refined;
  SC3E (compare_results (t, p3ptr, p, qvt));

  return NULL;
}

int
main (int argc, char **argv)
{
  p4est3_t *p3;
  p4est_t *p;
  setup_t             st, *t = &st;
  const p4est3_quadrant_vtable_t *qvt;

  SC3X (sc3_MPI_Init (&argc, &argv));
  t->mpicomm = SC3_MPI_COMM_WORLD;
  SC3X (sc3_MPI_Comm_rank (t->mpicomm, &t->mpirank));
  SC3X (sc3_MPI_Comm_size (t->mpicomm, &t->mpisize));
  sc_init (t->mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  SC3X (prepare_objects (&p3, &p, t, &qvt));
  SC3X (perform_test (p3, p, t, qvt));
  SC3X (clean_up (p3, p, t));

  sc_finalize_noabort ();
  SC3X (sc3_MPI_Finalize ());
  return 0;
}