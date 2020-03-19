/*
  This file is part of the SC Library, version 3.
  The SC Library provides support for parallel scientific applications.

  Copyright (C) 2019 individual authors

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:

  1. Redistributions of source code must retain the above copyright notice,
  this list of conditions and the following disclaimer.

  2. Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
  POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef P4_TO_P8
#include <p4est3_quadrant_zyx.h>
#include <p4est_p4est3.h>
#else
#include <p8est3_quadrant_zyx.h>
#include <p8est_p4est3.h>
#endif

#define N_QUADS_2_TEST 6000000

sc3_error_t *
create_quadrants (p4est3_quadrant_vtable_t * qvt_avx
                , p4est3_quadrant_vtable_t * qvt
                , __m128i * v
                , p4est_quadrant_t * q, int n_quad)
{
  int quad, put_ind;
  for (quad = 0, put_ind = 1
       ; P4EST_CHILDREN * (quad + 1) < n_quad
       ; ++quad,  put_ind += P4EST_CHILDREN
      ){
    for (int child = 0; child < P4EST_CHILDREN; ++child){
      SC3E (p4est3_quadrant_child (qvt_avx, &v[quad], child
                           , &v[put_ind + child]));
      SC3E (p4est3_quadrant_child (qvt, &q[quad], child
                           , &q[put_ind + child]));
    }
  }

  for (int child = 0; child < n_quad - put_ind; ++child) {
    SC3E (p4est3_quadrant_child (qvt_avx, &v[quad], child
                         , &v[put_ind + child]));
    SC3E (p4est3_quadrant_child (qvt, &q[quad], child
                         , &q[put_ind + child]));
  }
  return NULL;
}

sc3_error_t *
is_equal_child (p4est3_quadrant_vtable_t * qvt_avx
              , p4est3_quadrant_vtable_t * qvt
              , __m128i * v, p4est_quadrant_t * q, int32_t n_quad)
{
  int32_t is_equal;
  SC3E (create_quadrants (qvt_avx, qvt, v, q, n_quad));
  for (int i = 0; i < n_quad; ++i) {
    is_equal =
      (int32_t) q[i].level == _mm_extract_epi32 (v[i], 0) &&
      q[i].x == _mm_extract_epi32 (v[i], 3) &&
      q[i].y == _mm_extract_epi32 (v[i], 2) &&
#ifdef P4_TO_P8
      q[i].z == _mm_extract_epi32 (v[i], 1) &&
#endif
      1;
    SC3E_DEMAND (is_equal, "Comparing of quadrant_child results");
  }
  return NULL;
}

sc3_error_t *
is_equal_parent (p4est3_quadrant_vtable_t * qvt_avx
               , p4est3_quadrant_vtable_t * qvt, const __m128i * v
               , const p4est_quadrant_t * q, int n_quad, __m128i * v_p
               , p4est_quadrant_t * q_p)
{
  int32_t is_equal;
  for (int i = 1; i < n_quad; ++i) {
    SC3E (p4est3_quadrant_parent (qvt_avx, &v[i], v_p));
    SC3E (p4est3_quadrant_parent (qvt, &q[i], q_p));
    is_equal =
      (int32_t) q_p->level == _mm_extract_epi32 (*v_p, 0) &&
      q_p->x == _mm_extract_epi32 (*v_p, 3) &&
      q_p->y == _mm_extract_epi32 (*v_p, 2) &&
#ifdef P4_TO_P8
      q_p->z == _mm_extract_epi32 (*v_p, 1) &&
#endif
      1;
    SC3E_DEMAND (is_equal, "Comparing of quadrant_parent results");
  }
  return NULL;
}

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

  __m128i v;
  p4est_quadrant_t q; 

  const int32_t n_quads = N_QUADS_2_TEST;
  __m128i * v_quads =
          (__m128i *)malloc (p4est3_quadrant_size (qvt_avx) * n_quads);
  p4est_quadrant_t * q_quads =
          (p4est_quadrant_t *)malloc (p4est3_quadrant_size (qvt) * n_quads);

  v_quads[0] = _mm_setzero_si128 ();
  ((void) memset ((q_quads), 0, sizeof (p4est_quadrant_t)));

  SC3E_NULL_SET (
    e, is_equal_child (qvt_avx, qvt, v_quads, q_quads, n_quads));
  SC_CHECK_ABORT (e == NULL, "quadrant_child: non-vec != vec");
  
  SC3E_NULL_SET (
    e, is_equal_parent (qvt_avx, qvt, v_quads, q_quads, n_quads, &v, &q));
  SC_CHECK_ABORT (e == NULL, "quadrant_parent: non-vec != vec");

  free (qvt_avx);
  free (qvt);
  free (v_quads);
  free (q_quads);

  SC3E_NULL_REQ (e, !sc_finalize_noabort ());
  SC3E_NULL_SET (e, sc3_MPI_Finalize ());
  SC3E_TERR (e, reason);
  return 0;
}