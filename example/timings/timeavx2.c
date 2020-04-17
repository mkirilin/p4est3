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

#include <immintrin.h>
#include <smmintrin.h>
#include <emmintrin.h>

#include <time.h>

#define test_child(pull, qvt, n_quads, exec_time)                 \
({                                                                \
  int quad, put_ind;                                              \
  clock_t t_b = clock ();                                         \
                                                                  \
  for (quad = 0, put_ind = 1                                      \
       ; P4EST_CHILDREN * (quad + 1) < n_quads                    \
       ; ++quad,  put_ind += P4EST_CHILDREN                       \
      ){                                                          \
    for (int child = 0; child < P4EST_CHILDREN; ++child){         \
      SC3E (p4est3_quadrant_child (qvt, &pull[quad]               \
                          , child, &pull[put_ind + child]));      \
    }                                                             \
  }                                                               \
                                                                  \
  for (int child = 0; child < n_quads - put_ind; ++child) {       \
    SC3E (p4est3_quadrant_child (qvt, &pull[quad]                 \
                        , child, &pull[put_ind + child]));        \
  }                                                               \
  clock_t t_e = clock ();                                         \
                                                                  \
  exec_time = (float) (t_e - t_b) / CLOCKS_PER_SEC;               \
});

#define test_parent(p, pull, qvt, n_quads, exec_time)             \
({                                                                \
  clock_t t_b = clock();                                          \
  for (int quad = 1; quad < n_quads; ++quad) {                    \
    SC3E (p4est3_quadrant_parent (qvt, &pull[quad], p));          \
  }                                                               \
  clock_t t_e = clock();                                          \
  exec_time = (float) (t_e - t_b) / CLOCKS_PER_SEC;               \
});

#define test_compare(pull, qvt, n_quads, exec_time)               \
({                                                                \
  clock_t t_b = clock();                                          \
  int j;                                                          \
  for (int quad = 0; quad < n_quads; ++quad) {                    \
    SC3E (p4est3_quadrant_compare (qvt, &pull[quad]               \
                        , &pull[n_quads - quad - 1], &j));        \
  }                                                               \
  clock_t t_e = clock();                                          \
  exec_time = (float) (t_e - t_b) / CLOCKS_PER_SEC;               \
});

#define test_successor(p, pull, qvt, n_quads, exec_time)                      \
({                                                                            \
  clock_t t_b = clock();                                                      \
  for (int quad = 1; quad < n_quads-P4EST_CHILDREN; quad += P4EST_CHILDREN) { \
    for (int i = 0; i < P4EST_CHILDREN - 1; ++i) {                            \
      SC3E (p4est3_quadrant_successor (qvt, &pull[quad + i], p));             \
    }                                                                         \
  }                                                                           \
  clock_t t_e = clock();                                                      \
  exec_time = (float) (t_e - t_b) / CLOCKS_PER_SEC;                           \
});

static inline void
print_time_info (float exec_avx, float exec_nonavx, const char * name)
{
    printf ("  %s: \n"
            "    Vectorized:        %f\n"
            "    Non-Vectorized:    %f\n"
            "    Vect/Non-Vect Ratio:  %f\n"
            , name
            , exec_avx, exec_nonavx
            , exec_nonavx == 0. ? 0. : exec_avx / exec_nonavx);
}

static sc3_error_t *
measure_child (__m128i * v_pull2check, p4est_quadrant_t * q_pull2check,
               p4est3_quadrant_vtable_t * qvt_avx,
               p4est3_quadrant_vtable_t * qvt, int32_t n_quads)
{
  float               exec_avx, exec_nonavx;

  test_child (v_pull2check, qvt_avx, n_quads, exec_avx);
  test_child (q_pull2check, qvt, n_quads, exec_nonavx);

  print_time_info (exec_avx, exec_nonavx, "Child");
  return NULL;
}

static sc3_error_t *
measure_parent (__m128i * v_pull2check, p4est_quadrant_t * q_pull2check,
                p4est3_quadrant_vtable_t * qvt_avx,
                p4est3_quadrant_vtable_t * qvt, int32_t n_quads)
{
  float               exec_avx, exec_nonavx;
  __m128i             v_r;
  p4est_quadrant_t    q_r;

  test_parent (&v_r, v_pull2check, qvt_avx, n_quads, exec_avx);
  test_parent (&q_r, q_pull2check, qvt, n_quads, exec_nonavx);

  print_time_info (exec_avx, exec_nonavx, "Parent");
  return NULL;
}

static sc3_error_t *
measure_compare (__m128i * v_pull2check, p4est_quadrant_t * q_pull2check,
                 p4est3_quadrant_vtable_t * qvt_avx,
                 p4est3_quadrant_vtable_t * qvt, int32_t n_quads)
{
  float               exec_avx, exec_nonavx;

  test_compare (v_pull2check, qvt_avx, n_quads, exec_avx);
  test_compare (q_pull2check, qvt,  n_quads, exec_nonavx);

  print_time_info (exec_avx, exec_nonavx, "Compare");
  return NULL;
}

static sc3_error_t *
measure_successor(__m128i * v_pull2check, p4est_quadrant_t * q_pull2check,
                  p4est3_quadrant_vtable_t * qvt_avx,
                  p4est3_quadrant_vtable_t * qvt, int32_t n_quads)
{
  float               exec_avx, exec_nonavx;
  __m128i             v_r;
  p4est_quadrant_t    q_r;

  test_successor (&v_r, v_pull2check, qvt_avx, n_quads, exec_avx);
  test_successor (&q_r, q_pull2check, qvt, n_quads, exec_nonavx);

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
  p4est3_quadrant_vtable_t sqvt_avx, *qvt_avx = &sqvt_avx;
  p4est3_quadrant_vtable_t sqvt, *qvt = &sqvt;
  sc3_array_t        *qarr_avx, *qarr;

  p4est3_quadrant_zyx_vtable (qvt_avx);
  p4est_quadrant_vtable (qvt, 0);

  SC3E_SET (e, sc3_MPI_Init (&argc, &argv));

  int32_t n_quads;
  if(argc == 1) {
    n_quads = 0;
    for (int i = 0; i < 10; ++i){
      n_quads += (1 << (i * P4EST_DIM));
    }
  } else {
    n_quads = atoll(argv[1]);
  }

  /* TODO create a dedicated allocator for this program */
  /* TODO use quadrant array instead of manual allocation.
          This will allow us to not include immintrin.h in this program. */
  SC3E_SET (e, p4est3_quadrant_array_new (sc3_allocator_nocount (), qvt_avx,
                                          n_quads, &qarr_avx));
  SC3E_SET (e, p4est3_quadrant_array_new (sc3_allocator_nocount (), qvt,
                                          n_quads, &qarr));

  __m128i * v_pull2check =
         (__m128i *)malloc (p4est3_quadrant_size (qvt_avx) * n_quads);
  p4est_quadrant_t * q_pull2check =
          (p4est_quadrant_t *)malloc (p4est3_quadrant_size (qvt) * n_quads);

  v_pull2check[0] = _mm_setzero_si128 ();
  ((void) memset ((&q_pull2check[0]), 0, sizeof (p4est_quadrant_t)));

  if (e == NULL) {
    printf("Executing time: \n");
    SC3E_SET (e,
      measure_child (v_pull2check, q_pull2check, qvt_avx, qvt, n_quads));
    report_errors (&e);

    SC3E_SET (e,
      measure_parent (v_pull2check, q_pull2check, qvt_avx, qvt, n_quads));
    report_errors (&e);

    SC3E_SET (e,
      measure_compare (v_pull2check, q_pull2check, qvt_avx, qvt, n_quads));
    report_errors (&e);

    SC3E_SET (e,
      measure_successor (v_pull2check, q_pull2check, qvt_avx, qvt, n_quads));
    report_errors (&e);
  }

  free (v_pull2check);
  free (q_pull2check);
  SC3E_SET (e, sc3_array_destroy (&qarr_avx));
  SC3E_SET (e, sc3_array_destroy (&qarr));

  SC3E_NULL_REQ (e, !sc_finalize_noabort ());
  SC3E_NULL_SET (e, sc3_MPI_Finalize ());
  report_errors (&e);
  return 0;
}
