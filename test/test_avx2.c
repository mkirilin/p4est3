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
#include <p4est3_p4est.h>
#else
#include <p4est3_quadrant_zyx.h>
#include <p4est3_p8est.h>
#endif

#define N_QUADS_2_TEST 6000000

sc3_error_t        *
create_quadrants (p4est3_quadrant_vtable_t * qvt_avx,
                  p4est3_quadrant_vtable_t * qvt, sc3_array_t * v,
                  sc3_array_t * q, int n_quad)
{
  int                 quad, put_ind, child;
  void               *in, *out;

  for (quad = 0, put_ind = 1; P4EST_CHILDREN * (quad + 1) < n_quad;
       ++quad, put_ind += P4EST_CHILDREN) {
    for (child = 0; child < P4EST_CHILDREN; ++child) {
      SC3E (sc3_array_index (v, quad, &in));
      SC3E (sc3_array_index (v, put_ind + child, &out));
      SC3E (p4est3_quadrant_child (qvt_avx, in, child, out));

      SC3E (sc3_array_index (q, quad, &in));
      SC3E (sc3_array_index (q, put_ind + child, &out));
      SC3E (p4est3_quadrant_child (qvt, in, child, out));
    }
  }

  for (child = 0; child < n_quad - put_ind; ++child) {
    SC3E (sc3_array_index (v, quad, &in));
    SC3E (sc3_array_index (v, put_ind + child, &out));
    SC3E (p4est3_quadrant_child (qvt_avx, in, child, out));

    SC3E (sc3_array_index (q, quad, &in));
    SC3E (sc3_array_index (q, put_ind + child, &out));
    SC3E (p4est3_quadrant_child (qvt, in, child, out));
  }
  return NULL;
}

sc3_error_t        *
is_equal_child (p4est3_quadrant_vtable_t * qvt_avx,
                p4est3_quadrant_vtable_t * qvt, sc3_array_t * v,
                sc3_array_t * q, int32_t n_quad)
{
  int32_t             is_equal, i;
  int                 xy[P4EST_DIM], l;
  void               *p;
  void               *v_p;
  p4est_quadrant_t   *q_row;

  SC3E (create_quadrants (qvt_avx, qvt, v, q, n_quad));

  SC3E (sc3_array_index (q, 0, &p));
  q_row = (p4est_quadrant_t *) p;

  for (i = 0; i < n_quad; ++i) {
    SC3E (sc3_array_index (v, i, &v_p));
    SC3E (p4est3_quadrant_level (qvt_avx, v_p, &l));
    SC3E (p4est3_quadrant_coordinates (qvt_avx, v_p, P4EST_DIM, xy));
    is_equal =
      (int32_t) q_row[i].level == l && q_row[i].x == xy[0] && q_row[i].y == xy[1] &&
#ifdef P4_TO_P8
      q_row[i].z == xy[2] &&
#endif
      1;
    SC3E_DEMAND (is_equal, "Comparing of quadrant_child results");
  }
  return NULL;
}

sc3_error_t        *
is_equal_parent (p4est3_quadrant_vtable_t * qvt_avx,
                 p4est3_quadrant_vtable_t * qvt, sc3_array_t * v,
                 sc3_array_t * q, int n_quad)
{
  int32_t             is_equal, i;
  int                 xy[P4EST_DIM], l;
  void               *v_p;
  p4est_quadrant_t   *q_p;
  void               *in, *p;

  SC3E (sc3_array_index (v, 0, &v_p));

  SC3E (sc3_array_index (q, 0, &p));
  q_p = (p4est_quadrant_t *) p;

  for (i = 1; i < n_quad; ++i) {
    SC3E (sc3_array_index (v, i, &in));
    SC3E (p4est3_quadrant_parent (qvt_avx, in, v_p));

    SC3E (sc3_array_index (q, i, &in));
    SC3E (p4est3_quadrant_parent (qvt, in, q_p));

    SC3E (p4est3_quadrant_level (qvt_avx, v_p, &l));
    SC3E (p4est3_quadrant_coordinates (qvt_avx, v_p, P4EST_DIM, xy));
    is_equal = (int32_t) q_p->level == l && q_p->x == xy[0] && q_p->y == xy[1] &&
#ifdef P4_TO_P8
      q_p->z == xy[2] &&
#endif
      1;
    SC3E_DEMAND (is_equal, "Comparing of quadrant_parent results");
  }
  return NULL;
}

sc3_error_t        *
is_equal_compare (p4est3_quadrant_vtable_t * qvt_avx,
                  p4est3_quadrant_vtable_t * qvt, sc3_array_t * v,
                  sc3_array_t * q, int n_quad)
{
  int                 res, res_avx, i;
  void               *lhs, *rhs;
  for (i = 0; i < n_quad; ++i) {
    SC3E (sc3_array_index (v, i, &lhs));
    SC3E (sc3_array_index (v, n_quad - i - 1, &rhs));
    SC3E (p4est3_quadrant_compare (qvt_avx, lhs, rhs, &res_avx));

    SC3E (sc3_array_index (q, i, &lhs));
    SC3E (sc3_array_index (q, n_quad - i - 1, &rhs));
    SC3E (p4est3_quadrant_compare (qvt, lhs, rhs, &res));
    SC3E_DEMAND (res_avx == res, "Comparing of quadrant_compare results");
  }
  return NULL;
}

sc3_error_t        *
is_equal_successor (p4est3_quadrant_vtable_t * qvt_avx,
                    p4est3_quadrant_vtable_t * qvt, sc3_array_t * v,
                    sc3_array_t * q, int n_quad)
{
  int32_t             is_equal, i, j;
  int32_t             xy[P4EST_DIM], l;
  void               *v_p;
  p4est_quadrant_t   *q_p;
  void               *in, *p;

  SC3E (sc3_array_index (v, 0, &v_p));

  SC3E (sc3_array_index (q, 0, &p));
  q_p = (p4est_quadrant_t *) p;

  for (i = 1; i < n_quad; i += P4EST_CHILDREN) {
    for (j = 0; j < P4EST_CHILDREN - 1 && i + j < n_quad; ++j) {
      SC3E (sc3_array_index (v, i + j, &in));
      SC3E (p4est3_quadrant_successor (qvt_avx, in, v_p));

      SC3E (sc3_array_index (q, i + j, &in));
      SC3E (p4est3_quadrant_successor (qvt, in, q_p));

      SC3E (p4est3_quadrant_level (qvt_avx, v_p, &l));
      SC3E (p4est3_quadrant_coordinates (qvt_avx, v_p, P4EST_DIM, xy));
      is_equal = (int32_t) q_p->level == l && q_p->x == xy[0] && q_p->y == xy[1] &&
#ifdef P4_TO_P8
        q_p->z == xy[2] &&
#endif
        1;
      if (!is_equal) {
        return NULL;
      }
      SC3E_DEMAND (is_equal, "Comparing of quadrant_successor results");
    }
  }
  return NULL;
}

static void
report_errors (sc3_error_t ** pe)
{
  char                eflat[SC3_BUFSIZE];

  if (pe != NULL && *pe != NULL) {
    sc3_error_destroy_noerr (pe, eflat);
    fprintf (stderr, "Error: %s\n", eflat);
  }
}

int
main (int argc, char **argv)
{
  const int32_t       n_quads = N_QUADS_2_TEST;
  sc3_error_t        *e;
  sc3_error_kind_t    kind = SC3_ERROR_KIND_LAST;
  p4est3_quadrant_vtable_t sqvt_avx, *qvt_avx = &sqvt_avx;
  p4est3_quadrant_vtable_t sqvt, *qvt = &sqvt;
  sc3_array_t        *qarr_avx, *qarr;
  void               *p;

  e = p4est3_quadrant_yx_vtable (qvt_avx);
  if (e != NULL) {
    sc3_error_t        *e_;
    SC3E_SET (e_, sc3_error_get_kind (e, &kind));
    report_errors (&e_);
    if (kind == SC3_ERROR_RUNTIME) {
      report_errors (&e);
      return 0;
    }
  }

  SC3E_SET (e, sc3_MPI_Init (&argc, &argv));

  SC3E_NULL_SET (e, p4est3_quadrant_vtable_p4est (qvt, 0));

  SC3E_NULL_SET (e, p4est3_quadrant_array_new (sc3_allocator_nocount (),
                                               qvt_avx, n_quads, &qarr_avx));
  SC3E_NULL_SET (e, p4est3_quadrant_array_new (sc3_allocator_nocount (),
                                               qvt, n_quads, &qarr));

  SC3E_NULL_SET (e, sc3_array_index (qarr_avx, 0, &p));
  SC3E_NULL_SET (e, p4est3_quadrant_root (qvt_avx, p));

  SC3E_NULL_SET (e, sc3_array_index (qarr, 0, &p));
  SC3E_NULL_SET (e, p4est3_quadrant_root (qvt, p));

  SC3E_NULL_SET (e, is_equal_child (qvt_avx, qvt, qarr_avx, qarr, n_quads));
  SC_CHECK_ABORT (e == NULL, "quadrant_child: non-vec != vec");

  SC3E_NULL_SET (e, is_equal_parent (qvt_avx, qvt, qarr_avx, qarr, n_quads));
  SC_CHECK_ABORT (e == NULL, "quadrant_parent: non-vec != vec");

  SC3E_NULL_SET (e, is_equal_compare (qvt_avx, qvt, qarr_avx, qarr, n_quads));
  SC_CHECK_ABORT (e == NULL, "quadrant_compare: non-vec != vec");

  SC3E_NULL_SET (e,
                 is_equal_successor (qvt_avx, qvt, qarr_avx, qarr, n_quads));
  SC_CHECK_ABORT (e == NULL, "quadrant_successor: non-vec != vec");

  SC3E_NULL_SET (e, sc3_array_destroy (&qarr_avx));
  SC3E_NULL_SET (e, sc3_array_destroy (&qarr));

  SC3E_NULL_REQ (e, !sc_finalize_noabort ());
  SC3E_NULL_SET (e, sc3_MPI_Finalize ());
  report_errors (&e);

  return 0;
}
