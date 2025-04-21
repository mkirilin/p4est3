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
#include <p4est_ghost.h>
#include <p4est_vtk.h>
#include <p4est3_p4est.h>
#include <p4est3_convert_p4est.h>
#include <p4est3_ghost_p4est.h>

#else

#include <p8est_extended.h>
#include <p8est_bits.h>
#include <p8est_ghost.h>
#include <p8est_vtk.h>
#include <p4est3_p8est.h>
#include <p4est3_convert_p8est.h>
#include <p4est3_ghost_p8est.h>

#endif

#ifndef P4_TO_P8
static int          refine_level = 1;
#else
static int          refine_level = 4;
#endif

#define TEST_SIMPLE 1

#if TEST_SIMPLE == 1

static int
refine_fractal (p4est_t *p, p4est_topidx_t which_tree, p4est_quadrant_t *q)
{
  /* Refine every 7th (3d) or 3rd (2d) global quadrant. */
  p4est_locidx_t     *quadrant_local_id = (p4est_locidx_t *) p->user_pointer;

  return (((p->global_first_quadrant[p->mpirank] + (*quadrant_local_id)++) %
#ifdef P4_TO_P8
           7)
#else
           3)
#endif
          == 0);
}

#else

static int
refine_fn (p4est_t *p4est, p4est_topidx_t which_tree,
           p4est_quadrant_t *quadrant)
{
  int                 cid;
  p4est = (void *) p4est;

  if (which_tree == 2 || which_tree == 3) {
    return 0;
  }

  cid = p4est_quadrant_child_id (quadrant);

  if (cid == P4EST_CHILDREN - 1 ||
      (quadrant->x >= P4EST_LAST_OFFSET (P4EST_MAXLEVEL - 2) &&
       quadrant->y >= P4EST_LAST_OFFSET (P4EST_MAXLEVEL - 2)
#ifdef P4_TO_P8
       && quadrant->z >= P4EST_LAST_OFFSET (P4EST_MAXLEVEL - 2)
#endif
      )) {
    if (quadrant->level >= refine_level) {
      return 0;
    }
    return 1;
  }
  if ((int) quadrant->level >= (refine_level - (int) (which_tree % 3))) {
    return 0;
  }
  if (quadrant->level == 1 && cid == 2) {
    return 1;
  }
  if (quadrant->x == P4EST_QUADRANT_LEN (2) &&
      quadrant->y == P4EST_LAST_OFFSET (2)) {
    return 1;
  }
  if (quadrant->y >= P4EST_QUADRANT_LEN (2)) {
    return 0;
  }

  return 1;
}

#endif

static p4est_ghost_t *
init_ghost_layer (p4est_t *p)
{
  const p4est_topidx_t num_trees = p->connectivity->num_trees;
  p4est_ghost_t      *gl;

  gl = P4EST_ALLOC (p4est_ghost_t, 1);
  sc_array_init (&gl->ghosts, sizeof (p4est_quadrant_t));
  gl->tree_offsets = P4EST_ALLOC (p4est_locidx_t, num_trees + 1);
  gl->proc_offsets = P4EST_ALLOC (p4est_locidx_t, p->mpisize + 1);

  sc_array_init (&gl->mirrors, sizeof (p4est_quadrant_t));
  gl->mirror_tree_offsets = P4EST_ALLOC (p4est_locidx_t, num_trees + 1);
  gl->mirror_proc_mirrors = NULL;
  gl->mirror_proc_offsets = P4EST_ALLOC (p4est_locidx_t, p->mpisize + 1);

  return gl;
}

static sc3_error_t *
make_allocator (sc3_allocator_t *oa, sc3_allocator_t **alloc)
{
  SC3A_IS (sc3_allocator_is_setup, oa);
  SC3E (sc3_allocator_new (oa, alloc));
  SC3E (sc3_allocator_setup (*alloc));
  return NULL;
}

static sc3_error_t *
compare_ghost_results (p4est_ghost_t *ghost_p4est,
                       p4est_ghost_t *ghost_p4est3)
{
  int                 i;
  p4est_quadrant_t   *q_p4est, *q_p3;
  size_t              ghost_count_p4est, ghost_count_p3;
  size_t              mirror_count_p4est, mirror_count_p3;
  p4est_topidx_t      num_trees;
  p4est_locidx_t      n_offsetx;

  /* Check basic properties */
  SC3E_DEMAND (ghost_p4est->mpisize == ghost_p4est3->mpisize,
               "Ghost mpisize mismatch");
  SC3E_DEMAND (ghost_p4est->num_trees == ghost_p4est3->num_trees,
               "Ghost num_trees mismatch");
  SC3E_DEMAND (ghost_p4est->btype == ghost_p4est3->btype,
               "Ghost btype mismatch");

  /* Check ghost counts */
  ghost_count_p4est = ghost_p4est->ghosts.elem_count;
  ghost_count_p3 = ghost_p4est3->ghosts.elem_count;
  SC3E_DEMAND (ghost_count_p4est == ghost_count_p3,
               "Number of ghost elements mismatch");

  /* Check mirror counts */
  mirror_count_p4est = ghost_p4est->mirrors.elem_count;
  mirror_count_p3 = ghost_p4est3->mirrors.elem_count;
  SC3E_DEMAND (mirror_count_p4est == mirror_count_p3,
               "Number of mirror elements mismatch");

  /* Check tree offsets */
  num_trees = ghost_p4est->num_trees;
  for (i = 0; i <= num_trees; ++i) {
    SC3E_DEMAND (ghost_p4est->tree_offsets[i] ==
                 ghost_p4est3->tree_offsets[i],
                 "Ghost tree offsets mismatch");

    SC3E_DEMAND (ghost_p4est->mirror_tree_offsets[i] ==
                 ghost_p4est3->mirror_tree_offsets[i],
                 "Mirror tree offsets mismatch");
  }

  /* Check proc offsets */
  SC3E_DEMAND (ghost_p4est->mirror_proc_offsets != NULL &&
               ghost_p4est3->mirror_proc_offsets != NULL,
               "Ghost mirror_proc_offsets array is NULL");

  SC3E_DEMAND (ghost_p4est->proc_offsets != NULL &&
               ghost_p4est3->proc_offsets != NULL,
               "Ghost proc_offsets array is NULL");

  for (i = 0; i <= ghost_p4est->mpisize; ++i) {
    SC3E_DEMAND (ghost_p4est->proc_offsets[i] ==
                 ghost_p4est3->proc_offsets[i],
                 "Ghost proc offsets mismatch");

    SC3E_DEMAND (ghost_p4est->mirror_proc_offsets[i] ==
                 ghost_p4est3->mirror_proc_offsets[i],
                 "Mirror proc offsets mismatch");
  }

  /* Compare individual ghost quadrants */
  for (i = 0; i < (int) ghost_count_p4est; ++i) {
    q_p4est = p4est_quadrant_array_index (&ghost_p4est->ghosts, i);
    q_p3 = p4est_quadrant_array_index (&ghost_p4est3->ghosts, i);

    /* Compare position and level */
    SC3E_DEMAND (q_p4est->level == q_p3->level,
                 "Ghost quadrant level mismatch");
    SC3E_DEMAND (q_p4est->x == q_p3->x && q_p4est->y == q_p3->y,
                 "Ghost quadrant position mismatch");
#ifdef P4_TO_P8
    SC3E_DEMAND (q_p4est->z == q_p3->z, "Ghost quadrant z position mismatch");
#endif
    /* Compare piggy3 data */
    SC3E_DEMAND (q_p4est->p.piggy3.which_tree == q_p3->p.piggy3.which_tree,
                 "Ghost quadrant tree id mismatch");
    SC3E_DEMAND (q_p4est->p.piggy3.local_num == q_p3->p.piggy3.local_num,
                 "Ghost quadrant local_num mismatch");
  }

  /* Compare individual mirror quadrants */
  for (i = 0; i < (int) mirror_count_p4est; ++i) {
    q_p4est = p4est_quadrant_array_index (&ghost_p4est->mirrors, i);
    q_p3 = p4est_quadrant_array_index (&ghost_p4est3->mirrors, i);

    /* Compare position and level */
    SC3E_DEMAND (q_p4est->level == q_p3->level,
                 "Mirror quadrant level mismatch");
    SC3E_DEMAND (q_p4est->x == q_p3->x && q_p4est->y == q_p3->y,
                 "Mirror quadrant position mismatch");
#ifdef P4_TO_P8
    SC3E_DEMAND (q_p4est->z == q_p3->z,
                 "Mirror quadrant z position mismatch");
#endif
    /* Compare piggy3 data */
    SC3E_DEMAND (q_p4est->p.piggy3.which_tree == q_p3->p.piggy3.which_tree,
                 "Mirror quadrant tree id mismatch");
    SC3E_DEMAND (q_p4est->p.piggy3.local_num == q_p3->p.piggy3.local_num,
                 "Mirror quadrant local_num mismatch");
  }

  /* Check mirror_proc_mirrors */
  SC3E_DEMAND (ghost_p4est->mirror_proc_mirrors != NULL &&
               ghost_p4est3->mirror_proc_mirrors != NULL,
               "Ghost mirror_proc_mirrors array is NULL");

  n_offsetx = ghost_p4est->mirror_proc_offsets[ghost_p4est->mpisize];
  for (i = 0; i < n_offsetx; ++i) {
    SC3E_DEMAND (ghost_p4est->mirror_proc_mirrors[i] ==
                 ghost_p4est3->mirror_proc_mirrors[i],
                 "Ghost mirror_proc_mirrors array element mismatch");
  }

  return NULL;
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  sc_MPI_Comm         mpicomm;
  p4est_t            *p4est;
  p4est3_t           *p4est3;
  p4est_connectivity_t *conn;
  p4est3_connectivity_t *conn3;
  p4est_ghost_t      *ghost_p4est, *ghost_p4est3;
  sc3_allocator_t    *alloc, *mainalloc;
  sc3_error_t        *e = NULL;
#if TEST_SIMPLE == 1
  p4est3_locidx       quadrant_local_id = 0;
#endif

  /* initialize MPI */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpicomm = sc_MPI_COMM_WORLD;

  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  /*--------------------------------------------------------------*/
  /************************** P4EST 2 *****************************/
  /*--------------------------------------------------------------*/

#ifndef P4_TO_P8

#if TEST_SIMPLE == 1
  conn = p4est_connectivity_new_twotrees (1, 0, 0);
#else
  conn = p4est_connectivity_new_moebius ();
#endif

#else
  conn = p8est_connectivity_new_rotcubes ();
#endif

#if TEST_SIMPLE == 1
  p4est = p4est_new_ext (mpicomm, conn, 0, 0, 1, 0, NULL, &quadrant_local_id);
#else
  p4est = p4est_new (mpicomm, conn, 0, NULL, NULL);
#endif

  /* refine to make the number of elements interesting */
#if TEST_SIMPLE == 1
  /*** refine in a loop ***/
  for (int i = 0; i < refine_level; ++i) {
    quadrant_local_id = 0;
    p4est_refine (p4est, 0, refine_fractal, NULL);
  }
#else
  p4est_refine (p4est, 1, refine_fn, NULL);
#endif

  /* balance the forest */
  p4est_balance (p4est, P4EST_CONNECT_FULL, NULL);

  /* do a uniform partition */
  p4est_partition (p4est, 0, NULL);

  /* create the ghost layer for p4est */
  ghost_p4est = p4est_ghost_new (p4est, P4EST_CONNECT_FACE);
  printf ("p4est_ghost_new done\n");

  /*--------------------------------------------------------------*/
  /************************** P4EST 3 *****************************/
  /*--------------------------------------------------------------*/

  mainalloc = sc3_allocator_nothread ();
  SC3E_NULL_SET (e, make_allocator (mainalloc, &alloc));
  SC3E_NULL_SET (e, p4est3_new (alloc, &p4est3));
  SC3E_NULL_SET (e, p4est3_set_shared (p4est3, 1));
  SC3E_NULL_SET (e, p4est3_set_contiguous (p4est3, 1));

  SC3E_NULL_SET (e, p4est3_convert_p4est (p4est, p4est3, &conn3));
  printf ("p4est3_convert_p4est done\n");

  printf ("writing p4est3 file...\n");
  p4est_vtk_write_file (p4est, NULL, "test_ghost_p3");

  printf ("making p4est ghost layer...\n");
  ghost_p4est3 = init_ghost_layer (p4est);
  printf ("making p4est3 ghost layer...\n");
  SC3E_NULL_SET (e, p4est3_ghost_fill_p4est (p4est3, ghost_p4est3));

  printf ("comparing ghost layers...\n");
  SC3E_NULL_SET (e, compare_ghost_results (ghost_p4est, ghost_p4est3));

  /* clean up */
  p4est_ghost_destroy (ghost_p4est);
  p4est_ghost_destroy (ghost_p4est3);
  p4est_destroy (p4est);
  p4est_connectivity_destroy (conn);

  SC3E_NULL_SET (e, p4est3_destroy (&p4est3));
  SC3E_NULL_SET (e, p4est3_connectivity_destroy (&conn3));
  SC3E_NULL_SET (e, sc3_allocator_destroy (&alloc));

  /* exit */
  SC3E_NULL_REQ (e, !sc_finalize_noabort ());

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  SC3X (e);
  return 0;
}
