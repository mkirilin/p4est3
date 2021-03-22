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
#include <p4est3_quadrant_yx.h>
#include <p4est3_quadrant_mort2d.h>
#include <p4est3_p4est.h>
#else
#include <p4est3_quadrant_zyx.h>
#include <p4est3_quadrant_mort3d.h>
#include <p4est3_p8est.h>
#endif

#define MAX_TEST_LEVEL 5
#define MAX_TEST_TREES 5

typedef struct setup
{
  sc3_allocator_t    *alloc;
  sc3_allocator_t    *mainalloc;
  p4est3_connectivity_t *conn;
  sc3_MPI_Comm_t      mpicomm;
  int                 mpirank;
  int                 level;
  p4est3_topidx       num_trees;
}
setup_t;

static sc3_error_t *
set_vtables (p4est3_quadrant_vtable_t * q, p4est3_quadrant_vtable_t * qmort,
             p4est3_quadrant_vtable_t * qavx, sc3_error_t ** e)
{
  SC3X (p4est3_quadrant_vtable_p4est (q, 0));
  SC3X (p4est3_quadrant_mort2d_vtable (qmort));
  /* the AVX virtual table can only be set with hardware support */
  SC3F (p4est3_quadrant_yx_vtable (qavx), *e);

  return NULL;
}

static sc3_error_t *
make_allocator (setup_t *t)
{
  SC3A_IS (sc3_allocator_is_setup, t->mainalloc);
  SC3E (sc3_allocator_new (t->mainalloc, &t->alloc));
  SC3E (sc3_allocator_setup (t->alloc));

  return NULL;
}

static sc3_error_t *
make_connectivity (setup_t * t)
{
  SC3E (p4est3_connectivity_new (t->alloc, &t->conn));
  SC3E (p4est3_connectivity_set_num_trees (t->conn, t->num_trees));
  SC3E (p4est3_connectivity_setup (t->conn));

  return NULL;
}

static sc3_error_t *
make_new_p4est3 (p4est3_t ** p3, sc3_allocator_t * alloc,
                 p4est3_connectivity_t * conn, sc3_MPI_Comm_t mpicomm,
                 p4est3_quadrant_vtable_t * qvt, int level,
                 p4est3_setup_mode_t mode)
{
  SC3A_IS (sc3_allocator_is_setup, alloc);

  /* create p4est object with connectivity */
  SC3E (p4est3_new (alloc, p3));
  SC3E (p4est3_set_comm (*p3, mpicomm, 1));
  SC3E (p4est3_set_connectivity (*p3, conn));
  SC3E (p4est3_set_quadrant_vtable (*p3, qvt));
  SC3E (p4est3_set_level (*p3, level));
  SC3E (p4est3_set_setup_mode (*p3, mode));
  SC3E (p4est3_setup (*p3));

  return NULL;
}

static sc3_error_t *
setup_forests (setup_t * t, p4est3_quadrant_vtable_t * q,
               p4est3_t ** m, p4est3_t ** s, p4est3_t ** rc)
{
  SC3E (make_new_p4est3 (m, t->alloc, t->conn, t->mpicomm, q, t->level,
        P4EST3_NEW_MORTON));
  SC3E (make_new_p4est3 (s, t->alloc, t->conn, t->mpicomm, q, t->level,
      P4EST3_NEW_SUCCESSOR));
  SC3E (make_new_p4est3 (rc, t->alloc, t->conn, t->mpicomm, q, t->level,
      P4EST3_NEW_RECURSIVE_CHILD));

  return NULL;
}

static sc3_error_t *
compare_p4est3_quadrants (const p4est3_t * lhs, const p4est3_t * rhs,
                          p4est3_quadrant_vtable_t * lqvt,
                          p4est3_quadrant_vtable_t * rqvt)
{
  p4est3_gloidx       i;
  p4est3_gloidx       gln, grn;
  p4est3_locidx       lln, lrn;
  int                 is_eq, ll, rl;
  char               *lchar_q = NULL;
  char               *rchar_q = NULL;
  int                 lc[P4EST_DIM];
  int                 rc[P4EST_DIM];

  SC3E (p4est3_get_quadrants (lhs, &lchar_q));
  SC3E (p4est3_get_quadrants (rhs, &rchar_q));

  SC3E (p4est3_get_global_num_quads (lhs, &gln));
  SC3E (p4est3_get_global_num_quads (rhs, &grn));
  SC3E (p4est3_get_local_num_quads (lhs, &lln));
  SC3E (p4est3_get_local_num_quads (rhs, &lrn));
  SC3E_DEMAND (gln == grn, "different setup modes have different output");
  SC3E_DEMAND (lln == lrn, "different setup modes have different output");

  size_t              lq_size = p4est3_quadrant_size (lqvt);
  size_t              rq_size = p4est3_quadrant_size (rqvt);
  //lhs and rhs should be valid
  for (i = 0; i < lln; ++i, lchar_q += lq_size, rchar_q += rq_size) {
    SC3E (p4est3_quadrant_level (lqvt, lchar_q, &ll));
    SC3E (p4est3_quadrant_level (rqvt, rchar_q, &rl));
    SC3E (p4est3_quadrant_coordinates (lqvt, lchar_q, P4EST_DIM, lc));
    SC3E (p4est3_quadrant_coordinates (rqvt, rchar_q, P4EST_DIM, rc));
    is_eq = (int) (ll == rl && lc[0] == rc[0] && lc[1] == rc[1] &&
#ifdef P4_TO_P8
                   lc[2] == rc[2] &&
#endif
                   1);
    SC3E_DEMAND (is_eq == 1, "setup mod's results differ");
  }

  return NULL;
}

static sc3_error_t *
perform_test_old (setup_t * t,
                  p4est3_t ** m, p4est3_t ** s, p4est3_t ** rc,
                  p4est3_quadrant_vtable_t * qref,
                  p4est3_quadrant_vtable_t * q)
{
  SC3E (setup_forests (t, q, m, s, rc));
  SC3E (compare_p4est3_quadrants (*m, *s, q, q));
  SC3E (compare_p4est3_quadrants (*m, *rc, q, q));

  /*do not destroy m fprest, due to it is referenced for others*/
  SC3E (p4est3_destroy (s));
  SC3E (p4est3_destroy (rc));

  return NULL;
}

static sc3_error_t *
perform_test (setup_t * t,
              p4est3_t * ref, p4est3_t ** m, p4est3_t ** s, p4est3_t ** rc,
              p4est3_quadrant_vtable_t * qref, p4est3_quadrant_vtable_t * q)
{
  SC3E (setup_forests (t, q, m, s, rc));
  SC3E (compare_p4est3_quadrants (ref, *m, qref, q));
  SC3E (compare_p4est3_quadrants (*m, *s, q, q));
  SC3E (compare_p4est3_quadrants (*m, *rc, q, q));

  SC3E (p4est3_destroy (m));
  SC3E (p4est3_destroy (s));
  SC3E (p4est3_destroy (rc));

  return NULL;
}

static sc3_error_t *
free_allocator (sc3_allocator_t ** alloc)
{
  SC3A_IS (sc3_allocator_is_setup, *alloc);
  SC3E (sc3_allocator_destroy (alloc));
  return NULL;
}

int
main (int argc, char **argv)
{
  setup_t             st, *t = &st;
  sc3_error_t        *e_avx;
  p4est3_quadrant_vtable_t vtable, *qvt = &vtable;
  p4est3_quadrant_vtable_t vtavx, *qvtavx = &vtavx;
  p4est3_quadrant_vtable_t vtmort, *qvtmort = &vtmort;
  p4est3_t           *p3m, *p3s, *p3rc;
  p4est3_t           *p3m_avx, *p3s_avx, *p3rc_avx;
  p4est3_t           *p3m_mort, *p3s_mort, *p3rc_mort;

  /* v3 standard procedure to isolate memory allocation contexts */
  t->mainalloc = sc3_allocator_nothread ();

  /* this is generally needed for MPI */
  SC3X (sc3_MPI_Init (&argc, &argv));
  t->mpicomm = SC3_MPI_COMM_WORLD;
  SC3X (sc3_MPI_Comm_rank (t->mpicomm, &t->mpirank));

  SC3X (set_vtables (qvt, qvtmort, qvtavx, &e_avx));

  SC3X (make_allocator (t));
  for (t->level = 1; t->level < MAX_TEST_LEVEL; ++(t->level)) {
    for (t->num_trees = 1; t->num_trees < MAX_TEST_TREES; ++(t->num_trees)) {

#ifdef P4EST_ENABLE_DEBUG
      if (t->mpirank == 0) {
        printf ("l = %d, t = %d\n", t->level, t->num_trees);
      }
#endif /* P4EST_ENABLE_DEBUG */
      SC3X (sc3_MPI_Barrier (t->mpicomm));
      SC3X (make_connectivity (t));

      SC3X (perform_test_old (t, &p3m, &p3s, &p3rc, qvt, qvt));
#ifdef P4EST_ENABLE_DEBUG
      if (t->mpirank == 0) {
        printf ("Old-quadrants are equal for every setup option\n");
      }
#endif /* P4EST_ENABLE_DEBUG */

      SC3X (perform_test (t, p3m, &p3m_mort, &p3s_mort, &p3rc_mort, qvt,
                          qvtmort));
#ifdef P4EST_ENABLE_DEBUG
      if (t->mpirank == 0) {
        printf ("Mort-quadrants are equal for every setup option\n");
      }
#endif /* P4EST_ENABLE_DEBUG */

      if (!sc3_error_is2_kind (e_avx, SC3_ERROR_RUNTIME, NULL)) {
        SC3X (perform_test (t, p3m, &p3m_avx, &p3s_avx, &p3rc_avx, qvt, qvtavx));
#ifdef P4EST_ENABLE_DEBUG
        if (t->mpirank == 0) {
          printf ("AVX-quadrants are equal for every setup option\n");
        }
#endif /* P4EST_ENABLE_DEBUG */
      }

      /*destroy forest, that was referenced for others*/
      SC3X (p4est3_destroy (&p3m));
      SC3X (p4est3_connectivity_destroy (&t->conn));
    }
  }
  if (e_avx != NULL) {
    SC3X (sc3_error_unref (&e_avx));
  }
  SC3X (free_allocator (&t->alloc));

  SC3X (sc3_MPI_Finalize ());
  return 0;
}
