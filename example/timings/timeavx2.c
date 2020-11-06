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
#include <p4est_p4est3.h>
#else
#include <p8est3_quadrant_zyx.h>
#include <p8est_p4est3.h>
#endif

static sc3_error_t *
test_child (sc3_array_t * pull, p4est3_quadrant_vtable * qvt,
            p4est3_locidx n_quads, double *exec_time)
{
  p4est3_locidx       quad, put_ind;
  int                 child;
  void               *p, *q;
  double              t_b, t_e;

  t_b = sc3_MPI_Wtime ();
  for (quad = 0, put_ind = 1; P4EST_CHILDREN * (quad + 1) < n_quads;
       ++quad, put_ind += P4EST_CHILDREN) {
    for (child = 0; child < P4EST_CHILDREN; ++child) {
      SC3E (sc3_array_index (pull, quad, &p));
      SC3E (sc3_array_index (pull, put_ind + child, &q));
      SC3E (p4est3_quadrant_child (qvt, p, child, q));
    }
  }

  for (child = 0; child < n_quads - put_ind; ++child) {
    SC3E (sc3_array_index (pull, quad, &p));
    SC3E (sc3_array_index (pull, put_ind + child, &q));
    SC3E (p4est3_quadrant_child (qvt, p, child, q));
  }
  t_e = sc3_MPI_Wtime ();

  *exec_time = t_e - t_b;
  return NULL;
}

#define test_parent(q, pull, qvt, n_quads, exec_time) do {        \
  int                 quad;                                       \
  void               *p;                                          \
  double              t_b, t_e;                                   \
                                                                  \
  t_b = sc3_MPI_Wtime ();                                         \
  for (quad = 1; quad < n_quads; ++quad) {                        \
    SC3E (sc3_array_index (pull, quad, &p));                      \
    SC3E (p4est3_quadrant_parent (qvt, p, q));                    \
  }                                                               \
  t_e = sc3_MPI_Wtime ();                                         \
  exec_time = t_e - t_b;                                          \
  } while (0)

#define test_compare(pull, qvt, n_quads, exec_time) do {          \
  void               *p, *q;                                      \
  int                 j, quad;                                    \
  double              t_b, t_e;                                   \
                                                                  \
  t_b = sc3_MPI_Wtime ();                                         \
  for (quad = 0; quad < n_quads; ++quad) {                        \
    SC3E (sc3_array_index (pull, quad, &p));                      \
    SC3E (sc3_array_index (pull, n_quads - quad - 1, &q));        \
    SC3E (p4est3_quadrant_compare (qvt, p, q, &j));               \
  }                                                               \
  t_e = sc3_MPI_Wtime ();                                         \
  exec_time = t_e - t_b;                                          \
  } while (0)

#define test_successor(q, pull, qvt, n_quads, exec_time) do {               \
  int                 i, quad;                                              \
  void               *p;                                                    \
  double              t_b, t_e;                                             \
                                                                            \
  t_b = sc3_MPI_Wtime ();                                                   \
  for (quad = 1; quad < n_quads-P4EST_CHILDREN; quad += P4EST_CHILDREN) {   \
    for (i = 0; i < P4EST_CHILDREN - 1; ++i) {                              \
      SC3E (sc3_array_index (pull, quad + i, &p));                          \
      SC3E (p4est3_quadrant_successor (qvt, p, q));                         \
    }                                                                       \
  }                                                                         \
  t_e = sc3_MPI_Wtime ();                                                   \
  exec_time = t_e - t_b;                                                    \
  } while (0)

static void
print_time_info (double exec_avx, double exec_nonavx, const char *name)
{
  printf ("  %s: \n"
          "    Vectorized:        %g\n"
          "    Non-Vectorized:    %g\n"
          "    Vect/Non-Vect Ratio:  %g\n",
          name, exec_avx, exec_nonavx,
          exec_nonavx == 0. ? 0. : exec_avx / exec_nonavx);
}

static sc3_error_t *
measure_child (sc3_array_t * v_pull2check, sc3_array_t * q_pull2check,
               p4est3_quadrant_vtable_t * qvt_avx,
               p4est3_quadrant_vtable_t * qvt, p4est3_locidx n_quads)
{
  double              exec_avx, exec_nonavx;

  SC3E (test_child (v_pull2check, qvt_avx, n_quads, &exec_avx));
  SC3E (test_child (q_pull2check, qvt, n_quads, &exec_nonavx));

  print_time_info (exec_avx, exec_nonavx, "Child");
  return NULL;
}

static sc3_error_t *
measure_parent (sc3_array_t * v_pull2check, sc3_array_t * q_pull2check,
                p4est3_quadrant_vtable_t * qvt_avx,
                p4est3_quadrant_vtable_t * qvt, p4est3_locidx n_quads)
{
  double              exec_avx, exec_nonavx;
  void               *v, *q;

  SC3E (sc3_array_index (v_pull2check, 0, &v));
  SC3E (sc3_array_index (q_pull2check, 0, &q));
  test_parent (v, v_pull2check, qvt_avx, n_quads, exec_avx);
  test_parent (q, q_pull2check, qvt, n_quads, exec_nonavx);

  print_time_info (exec_avx, exec_nonavx, "Parent");
  return NULL;
}

static sc3_error_t *
measure_compare (sc3_array_t * v_pull2check, sc3_array_t * q_pull2check,
                 p4est3_quadrant_vtable_t * qvt_avx,
                 p4est3_quadrant_vtable_t * qvt, p4est3_locidx n_quads)
{
  double              exec_avx, exec_nonavx;

  test_compare (v_pull2check, qvt_avx, n_quads, exec_avx);
  test_compare (q_pull2check, qvt, n_quads, exec_nonavx);

  print_time_info (exec_avx, exec_nonavx, "Compare");
  return NULL;
}

static sc3_error_t *
measure_successor (sc3_array_t * v_pull2check, sc3_array_t * q_pull2check,
                   p4est3_quadrant_vtable_t * qvt_avx,
                   p4est3_quadrant_vtable_t * qvt, p4est3_locidx n_quads)
{
  float               exec_avx, exec_nonavx;
  void               *v, *q;

  SC3E (sc3_array_index (v_pull2check, 0, &v));
  SC3E (sc3_array_index (q_pull2check, 0, &q));
  test_successor (v, v_pull2check, qvt_avx, n_quads, exec_avx);
  test_successor (q, q_pull2check, qvt, n_quads, exec_nonavx);

  print_time_info (exec_avx, exec_nonavx, "Successor");
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
  sc3_error_t        *e;
  sc3_error_kind_t    kind = SC3_ERROR_KIND_LAST;
  p4est3_locidx       n_quads;
  p4est3_quadrant_vtable_t sqvt_avx, *qvt_avx = &sqvt_avx;
  p4est3_quadrant_vtable_t sqvt, *qvt = &sqvt;
  sc3_array_t        *qarr_avx, *qarr;
  void               *p;

  e = p4est3_quadrant_zyx_vtable (qvt_avx);
  if (e != NULL) {
    sc3_error_t        *e_;
    SC3E_SET (e_, sc3_error_get_kind (e, &kind));
    report_errors (&e_);
    if (kind == SC3_ERROR_RUNTIME) {
      report_errors (&e);
      return 0;
    }
  }
  p4est_quadrant_vtable (qvt, 0);

  SC3E_SET (e, sc3_MPI_Init (&argc, &argv));

  if (argc == 1) {
    n_quads = 0;
    for (int i = 0; i < 10; ++i) {
      n_quads += (1 << (i * P4EST_DIM));
    }
  }
  else {
    n_quads = atoll (argv[1]);
  }

  /* TODO create a dedicated allocator for this program */
  SC3E_NULL_SET (e, p4est3_quadrant_array_new (sc3_allocator_nocount (),
                                               qvt_avx, n_quads, &qarr_avx));
  SC3E_NULL_SET (e, p4est3_quadrant_array_new (sc3_allocator_nocount (),
                                               qvt, n_quads, &qarr));

  SC3E_NULL_SET (e, sc3_array_index (qarr_avx, 0, &p));
  SC3E_NULL_SET (e, p4est3_quadrant_root (qvt_avx, p));

  SC3E_NULL_SET (e, sc3_array_index (qarr, 0, &p));
  SC3E_NULL_SET (e, p4est3_quadrant_root (qvt, p));

  if (e == NULL) {
    printf ("Executing time: \n");
    SC3E_SET (e, measure_child (qarr_avx, qarr, qvt_avx, qvt, n_quads));
    report_errors (&e);

    SC3E_SET (e, measure_parent (qarr_avx, qarr, qvt_avx, qvt, n_quads));
    report_errors (&e);

    SC3E_SET (e, measure_compare (qarr_avx, qarr, qvt_avx, qvt, n_quads));
    report_errors (&e);

    SC3E_SET (e, measure_successor (qarr_avx, qarr, qvt_avx, qvt, n_quads));
    report_errors (&e);
  }

  SC3E_NULL_SET (e, sc3_array_destroy (&qarr_avx));
  SC3E_NULL_SET (e, sc3_array_destroy (&qarr));

  SC3E_NULL_REQ (e, !sc_finalize_noabort ());
  SC3E_NULL_SET (e, sc3_MPI_Finalize ());
  report_errors (&e);
  return 0;
}
