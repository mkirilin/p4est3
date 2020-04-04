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

#include <time.h>

#ifndef SC3E_TERR
/* copied here until we're clear on whether we're keeping it */
#define SC3E_TERR(f,r) do {                                             \
  sc3_error_t * _e = (f);                                               \
  if (_e != NULL) {                                                     \
    sc3_error_destroy_noerr (&_e, r); return 0; }} while (0)
#endif

#define test_child(e, pull, qvt, n_quads, exec_time)              \
({                                                                \
  int quad, put_ind;                                              \
  clock_t t_b = clock ();                                         \
                                                                  \
  for (quad = 0, put_ind = 1                                      \
       ; P4EST_CHILDREN * (quad + 1) < n_quads                    \
       ; ++quad,  put_ind += P4EST_CHILDREN                       \
      ){                                                          \
    for (int child = 0; child < P4EST_CHILDREN; ++child){         \
      SC3E_NULL_SET (e, p4est3_quadrant_child (qvt, &pull[quad]   \
                          , child, &pull[put_ind + child]));      \
    }                                                             \
  }                                                               \
                                                                  \
  for (int child = 0; child < n_quads - put_ind; ++child) {       \
    SC3E_NULL_SET (e, p4est3_quadrant_child (qvt, &pull[quad]     \
                        , child, &pull[put_ind + child]));        \
  }                                                               \
  clock_t t_e = clock ();                                         \
                                                                  \
  exec_time = (float) (t_e - t_b) / CLOCKS_PER_SEC;               \
});

#define test_parent(e, p, pull, qvt, n_quads, exec_time)          \
({                                                                \
  clock_t t_b = clock();                                          \
  for (int quad = 1; quad < n_quads; ++quad) {                    \
    SC3E_NULL_SET (e, p4est3_quadrant_parent (qvt                 \
                        , &pull[quad], p));                       \
  }                                                               \
  clock_t t_e = clock();                                          \
  exec_time = (float) (t_e - t_b) / CLOCKS_PER_SEC;               \
});

#define test_compare(e, pull, qvt, n_quads, exec_time)            \
({                                                                \
  clock_t t_b = clock();                                          \
  int j;                                                          \
  for (int quad = 0; quad < n_quads; ++quad) {                    \
    SC3E_NULL_SET (e, p4est3_quadrant_compare (qvt, &pull[quad]   \
                        , &pull[n_quads - quad - 1], &j));        \
  }                                                               \
  clock_t t_e = clock();                                          \
  exec_time = (float) (t_e - t_b) / CLOCKS_PER_SEC;               \
});

#define test_successor(e, p, pull, qvt, n_quads, exec_time)                   \
({                                                                            \
  clock_t t_b = clock();                                                      \
  for (int quad = 1; quad < n_quads-P4EST_CHILDREN; quad += P4EST_CHILDREN) { \
    for (int i = 0; i < P4EST_CHILDREN - 1; ++i) {                            \
      SC3E_NULL_SET (e, p4est3_quadrant_successor (qvt, &pull[quad + i], p)); \
    }                                                                         \
  }                                                                           \
  clock_t t_e = clock();                                                      \
  exec_time = (float) (t_e - t_b) / CLOCKS_PER_SEC;                           \
});

int
main (int argc, char **argv)
{
  sc3_error_t        *e;
  char                reason[SC3_BUFSIZE];

  p4est3_quadrant_vtable_t *qvt_avx
      = (p4est3_quadrant_vtable_t *)malloc (sizeof (p4est3_quadrant_vtable_t));
  p4est3_quadrant_vtable_t *qvt
      = (p4est3_quadrant_vtable_t *)malloc (sizeof (p4est3_quadrant_vtable_t));
  
  p4est3_quadrant_zyx_vtable (qvt_avx);
  p4est_quadrant_vtable (qvt, 0);

  SC3E_SET (e, sc3_MPI_Init (&argc, &argv));

  __m128i v_r;
  p4est_quadrant_t q_r; 

  int32_t n_quads;
  if(argc == 1) {
    n_quads = 0;
    for (int i = 0; i < 10; ++i){
      n_quads += (1 << (i * P4EST_DIM));
    }
  } else {
    n_quads = atoll(argv[1]);
  }

  __m128i * v_pull2check =
         (__m128i *)malloc (p4est3_quadrant_size (qvt_avx) * n_quads);
  p4est_quadrant_t * q_pull2check =
          (p4est_quadrant_t *)malloc (p4est3_quadrant_size (qvt) * n_quads);
          
  v_pull2check[0] = _mm_setzero_si128 ();
  ((void) memset ((&q_pull2check[0]), 0, sizeof (p4est_quadrant_t)));

  float exec_avx = 0., exec_nonavx = 0.;
  test_child (e, v_pull2check, qvt_avx, n_quads, exec_avx);
  test_child (e, q_pull2check, qvt, n_quads, exec_nonavx);
  printf ("Executing time: \n"
          "  Child: \n"
          "    Vectorized:        %f\n"
          "    Non-Vectorized:    %f\n"
          "    Vect/Non-Vect Ratio:  %f\n"
          , exec_avx, exec_nonavx
          , exec_nonavx == 0. ? 0. : exec_avx / exec_nonavx);
  SC3E_TERR (e, reason);
  
  test_parent (e, &v_r, v_pull2check, qvt_avx, n_quads, exec_avx);
  test_parent (e, &q_r, q_pull2check, qvt, n_quads, exec_nonavx);
  printf ("  Parent: \n"
          "    Vectorized:        %f\n"
          "    Non-Vectorized:    %f\n"
          "    Vect/Non-Vect Ratio:  %f\n"
          , exec_avx, exec_nonavx
          , exec_nonavx == 0. ? 0. : exec_avx / exec_nonavx);
  SC3E_TERR (e, reason);

  test_compare (e, v_pull2check, qvt_avx, n_quads, exec_avx);
  test_compare (e, q_pull2check, qvt,  n_quads, exec_nonavx);
  printf ("  Compare: \n"
          "    Vectorized:        %f\n"
          "    Non-Vectorized:    %f\n"
          "    Vect/Non-Vect Ratio:  %f\n"
          , exec_avx, exec_nonavx
          , exec_nonavx == 0. ? 0. : exec_avx / exec_nonavx);
  SC3E_TERR (e, reason);

  test_successor (e, &v_r, v_pull2check, qvt_avx, n_quads, exec_avx);
  test_successor (e, &q_r, q_pull2check, qvt, n_quads, exec_nonavx);
  printf ("  Successor: \n"
          "    Vectorized:        %f\n"
          "    Non-Vectorized:    %f\n"
          "    Vect/Non-Vect Ratio:  %f\n"
          , exec_avx, exec_nonavx
          , exec_nonavx == 0. ? 0. : exec_avx / exec_nonavx);
  SC3E_TERR (e, reason);

  free (qvt_avx);
  free (qvt);
  free (v_pull2check);
  free (q_pull2check);

  SC3E_NULL_REQ (e, !sc_finalize_noabort ());
  SC3E_NULL_SET (e, sc3_MPI_Finalize ());
  SC3E_TERR (e, reason);
  return 0;
}
