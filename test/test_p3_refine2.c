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

#define MAX_TEST_LEVEL 5
#define MAX_TEST_TREES 5
#define FOREST_START_LEVEL 1

/* Let's assume the overall max_level for every quadrant implementation */
#define TEST_MAX_LEVEL 15
#define TEST_ROOT_LEN ((int32_t) 1 << TEST_MAX_LEVEL)
#define TEST_QUADRANT_LEN(l) ((int32_t) 1 << (TEST_MAX_LEVEL - (l)))
#define TEST_LAST_OFFSET(l) (TEST_ROOT_LEN - TEST_QUADRANT_LEN (l))

static int          refine_level = 0;

typedef struct setup
{
  sc3_allocator_t    *alloc;
  sc3_allocator_t    *mainalloc;
  p4est_connectivity_t *conn2;
  p4est3_connectivity_t *conn3;
  sc3_MPI_Comm_t      mpicomm;
  int                 mpirank;
  int                 level;
  p4est3_topidx       num_trees;
}
setup_t;

static int
refine_normal_fn (p4est_t * p4est, p4est_topidx_t which_tree,
                  p4est_quadrant_t * quadrant)
{
  if ((int) quadrant->level >= (refine_level - (int) (which_tree % 3))) {
    return 0;
  }
  if (quadrant->level == 1 && p4est_quadrant_child_id (quadrant) == 3) {
    return 1;
  }
  if (quadrant->x == TEST_LAST_OFFSET (2) &&
      quadrant->y == TEST_LAST_OFFSET (2)) {
    return 1;
  }
  if (quadrant->x >= TEST_QUADRANT_LEN (2)) {
    return 0;
  }

  return 1;
}

static int
coarsen_normal_fn (p4est_t * p4est, p4est_topidx_t which_tree,
                   p4est_quadrant_t * quadrants[])
{
  const int           condition = refine_level - 2 < 1 ? 1 : refine_level - 2;
  if ((int) quadrants[0]->level > condition) {
    return 1;
  }
  return 0;
}

static sc3_error_t *
refine_p3_normal_fn (p4est3_refine_callback_info_t * ri, int *is_refine)
{
  int                 level, child_id, coords[P4EST_DIM];
  *is_refine = 1;

  SC3E (p4est3_quadrant_level (ri->qvt, ri->quadrant, &level));
  if (level >= (refine_level - (ri->ntree % 3))) {
    *is_refine = 0;
    return NULL;
  }

  SC3E (p4est3_quadrant_child_id (ri->qvt, ri->quadrant, &child_id));
  if (level == 1 && child_id == 3) {
    *is_refine = 1;
    return NULL;
  }

  SC3E (p4est3_quadrant_coordinates
        (ri->qvt, ri->quadrant, ri->qvt->dim, coords));
  if (coords[0] == TEST_LAST_OFFSET (2) && coords[1] == TEST_LAST_OFFSET (2)) {
    *is_refine = 1;
    return NULL;
  }

  if (coords[0] >= TEST_QUADRANT_LEN (2)) {
    *is_refine = 0;
    return NULL;
  }
  return NULL;
}

static sc3_error_t *
coarsen_p3_normal_fn (p4est3_coarsen_callback_info_t * ci, int *is_coarse)
{
  const int           condition = refine_level - 2 < 1 ? 1 : refine_level - 2;
  int                 level;
  void              **q;

  SC3E_RETVAL (is_coarse, 0);
  SC3E (sc3_array_index (ci->family, 0, &q));
  SC3E (p4est3_quadrant_level (ci->qvt, *(void **) q, &level));
  if (level > condition) {
    *is_coarse = 1;
  }
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
make_connectivity (setup_t * t, int dim)
{
  t->conn2 =
#ifdef P4_TO_P8
    p8est_connectivity_new_brick (t->num_trees, 1, 1, 0, 0, 0);
#else
    p4est_connectivity_new_brick (t->num_trees, 1, 0, 0);
#endif
  SC3E (p4est3_connectivity_new_p4est (t->alloc, &t->conn3, t->conn2, 1));
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
set_parameters (setup_t * t, const p4est3_quadrant_vtable_t ** qvt)
{
  t->mainalloc = sc3_allocator_nothread ();
  SC3E (make_allocator (t));
  SC3E (p4est3_quadrant_vtable_p4est (qvt));
  SC3E_DEMAND (*qvt != NULL, "p4est is not build neither in 2D nor 3D");
  return NULL;
}

static sc3_error_t *
make_new_p4est (p4est_t ** p, setup_t * t)
{
  *p = p4est_new_ext
    (t->mpicomm, t->conn2, 0, FOREST_START_LEVEL, 1, 0, NULL, NULL);

  return NULL;
}

static sc3_error_t *
make_new_p4est3 (p4est3_t ** p3, setup_t * t, const p4est3_quadrant_vtable_t ** qvt)
{
  SC3A_IS (sc3_allocator_is_setup, t->alloc);

  /* create p4est object with connectivity */
  SC3E (p4est3_new (t->alloc, p3));
  SC3E (p4est3_set_comm (*p3, t->mpicomm, 1));
  SC3E (p4est3_set_connectivity (*p3, t->conn3));
  SC3E (p4est3_set_quadrant_vtable (*p3, *qvt));
  SC3E (p4est3_set_level (*p3, FOREST_START_LEVEL));
  SC3E (p4est3_setup (*p3));

  return NULL;
}

static sc3_error_t *
free_allocator (sc3_allocator_t ** alloc)
{
  SC3A_IS (sc3_allocator_is_setup, *alloc);
  SC3E (sc3_allocator_destroy (alloc));
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
  SC3E (array_new (t->alloc, qvt->quadrant_size, num_loc_quads, 0, &levels));
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
perform_test (setup_t * t, p4est3_t * p3, p4est_t * p)
{
  int                 i;
  p4est3_t           *p3refined, *p3ptr = p3;
  const p4est3_quadrant_vtable_t *qvt;
  /* refine the old forest */
  refine_level = t->level;
  p4est_refine (p, 1, refine_normal_fn, NULL);

  for (i = 0; i < refine_level; ++i) {
    SC3E (p4est3_new (t->alloc, &p3refined));
    SC3E (set_qvt (&qvt, i % 3));
    SC3E (p4est3_set_quadrant_vtable (p3refined, qvt));
    SC3E (p4est3_set_refine (p3refined, refine_p3_normal_fn));
    SC3E (p4est3_set_source (p3refined, p3ptr));
    SC3E (p4est3_setup (p3refined));

    if (i != 0) {
      SC3E (p4est3_destroy (&p3ptr));
    }
    p3ptr = p3refined;
  }
  SC3E (compare_results (t, p3ptr, p, qvt));

  p4est_coarsen (p, 1, coarsen_normal_fn, NULL);
  for (i = 0; i < refine_level; ++i) {
    SC3E (p4est3_new (t->alloc, &p3refined));
    SC3E (set_qvt (&qvt, i % 3));
    SC3E (p4est3_set_quadrant_vtable (p3refined, qvt));
    SC3E (p4est3_set_coarsen (p3refined, coarsen_p3_normal_fn));
    SC3E (p4est3_set_source (p3refined, p3ptr));
    SC3E (p4est3_setup (p3refined));

    SC3E (p4est3_destroy (&p3ptr));
    p3ptr = p3refined;
  }
  SC3E (compare_results (t, p3ptr, p, qvt));

  SC3E (p4est3_new (t->alloc, &p3refined));
  SC3E (set_qvt (&qvt, refine_level % 3));
  SC3E (p4est3_set_quadrant_vtable (p3refined, qvt));
  SC3E (p4est3_set_source (p3refined, p3ptr));
  SC3E (p4est3_setup (p3refined));
  SC3E (p4est3_destroy (&p3ptr));
  p3ptr = p3refined;

  SC3E (compare_results (t, p3ptr, p, qvt));

  SC3E (p4est3_destroy (&p3ptr));
  return NULL;
}

static sc3_error_t *
perform_tests (setup_t * t, const p4est3_quadrant_vtable_t ** qvt)
{
  p4est_t            *p;
  p4est3_t           *p3;

  for (t->level = 1; t->level <= MAX_TEST_LEVEL; ++t->level) {
    for (t->num_trees = 1; t->num_trees <= MAX_TEST_TREES; ++t->num_trees) {
#ifdef P4EST_ENABLE_DEBUG
      if (t->mpirank == 0) {
        printf ("l = %d, t = %d\n", t->level, t->num_trees);
      }
#endif /* P4EST_ENABLE_DEBUG */
      SC3E (make_connectivity (t, (*qvt)->dim));
      SC3E (make_new_p4est (&p, t));
      SC3E (make_new_p4est3 (&p3, t, qvt));

      SC3E (perform_test (t, p3, p));

      /*destroy forest, that was referenced for others */
      SC3E (p4est3_destroy (&p3));
      SC3E (p4est3_connectivity_destroy (&t->conn3));
      p4est_destroy (p);
      /* There is no need to destroy p4est2 connectivity,
         since it is destroyed at p4est3 conn destroying stage */
      /* p4est_connectivity_destroy (t->conn2); */
    }
  }
  return NULL;
}

int
main (int argc, char **argv)
{
  setup_t             st, *t = &st;
  const p4est3_quadrant_vtable_t *qvt;

  SC3X (sc3_MPI_Init (&argc, &argv));
  t->mpicomm = SC3_MPI_COMM_WORLD;
  SC3X (sc3_MPI_Comm_rank (t->mpicomm, &t->mpirank));
  sc_init (t->mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  SC3X (set_parameters (t, &qvt));
  SC3X (perform_tests (t, &qvt));
  SC3X (free_allocator (&t->alloc));
  sc_finalize_noabort ();
  SC3X (sc3_MPI_Finalize ());
  return 0;
}
