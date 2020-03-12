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

#ifdef P4_TO_P8
#include <p8est3_quadrant_zyx.h>
#else
#include <p4est3_quadrant_zyx.h>
#endif /* !P4_TO_P8 */

static int
p4est3_quadrant_zyx_is_inside_root (const __m128i * q, char *reason)
{

  SC3E_TEST (
    _mm_test_all_ones (
      _mm_cmpgt_epi32 (*q, _mm_set1_epi32(-1)) ) == 1
  , reason
  );
  SC3E_TEST (
    _mm_test_all_ones (
      _mm_cmplt_epi32 (*q, _mm_set1_epi32 (P4EST_ROOT_LEN)) ) == 1
  , reason
  );
  SC3E_YES (reason);
}

static int
p4est3_quadrant_zyx_is_valid (const __m128i * q, char *reason)
{
  int32_t  level = _mm_extract_epi32 (*q, 0);
  SC3E_TEST (level >= 0, reason);
  SC3E_TEST (level <= P4EST_MAXLEVEL, reason);
  SC3E_TEST (
    _mm_testz_si128 (*q, _mm_set_epi32 (P4EST_QUADRANT_LEN (level) - 1
                                      , P4EST_QUADRANT_LEN (level) - 1
                                      , P4EST_QUADRANT_LEN (level) - 1
                                      , 0)) == 1
  , reason
  );
  SC3E_IS (p4est3_quadrant_zyx_is_inside_root, q, reason);
  SC3E_YES (reason);
}

static int
p4est3_quadrant_zyx_is_parent (const __m128i * q, const __m128i * r
                            , char *reason)
{
  SC3E_IS (p4est3_quadrant_zyx_is_valid, q, reason);
  SC3E_IS (p4est3_quadrant_zyx_is_valid, r, reason);

  int32_t r_level = _mm_extract_epi32 (*r, 0);

  SC3E_TEST (_mm_extract_epi32 (*q, 0) + 1 == r_level, reason);

  __m128i lhs = _mm_add_epi32 (*q, _mm_set_epi32 (0, 0, 0, 1));
  __m128i rhs = _mm_and_si128 (*r
                , _mm_set_epi32 ( //Optimized by compiler?
                                 ~P4EST_QUADRANT_LEN (r_level)
                               , ~P4EST_QUADRANT_LEN (r_level)
#ifdef P4_TO_P8
                               , ~P4EST_QUADRANT_LEN (r_level)
#else
                               , 0
#endif
                               , 0xFFFFFFFF)
      );
  SC3E_TEST (_mm_testc_si128 (rhs, lhs) == 1, reason);
  SC3E_YES (reason);
}

static int
p4est3_quadrant_zyx_is_node (const __m128i * q, int inside, char *reason)
{
  SC3E_TEST (_mm_extract_epi32 (*q, 0) == P4EST_MAXLEVEL, reason);
  SC3E_TEST (_mm_test_all_zeros (_mm_set_epi32 ((P4EST_ROOT_LEN) << 1
                                              , (P4EST_ROOT_LEN) << 1
                                              , (P4EST_ROOT_LEN) << 1
                                              , 0), *q) == 1
            , reason);
  SC3E_TEST (_mm_test_all_zeros (_mm_set_epi32 (P4EST_MAXLEVEL
                                              , P4EST_MAXLEVEL
                                              , P4EST_MAXLEVEL
                                              , 0), *q) == 1
            , reason);
  SC3E_TEST (
    (inside ? 1 : _mm_test_all_zeros (_mm_set_epi32 ((P4EST_MAXLEVEL) - 1 
                                                  , (P4EST_MAXLEVEL) - 1 
                                                  , (P4EST_MAXLEVEL) - 1
                                                  , 0), *q)) == 1
    , reason);
  __m128i isNodeCell =
    _mm_or_si128 (
      _mm_andnot_si128 (*q
      , _mm_set1_epi32 ((1 << (P4EST_MAXLEVEL - P4EST_QMAXLEVEL)) - 1))

      , _mm_and_si128 (
            _mm_set_epi32 (inside, inside, inside, 0xFFFFFFFF)
          , _mm_cmpeq_epi32 (*q, _mm_set_epi32 (P4EST_ROOT_LEN - 1
                                              , P4EST_ROOT_LEN - 1
                                              , P4EST_ROOT_LEN - 1
                                              , _mm_extract_epi32 (*q, 0)))
        )
    );
  SC3E_TEST (_mm_test_all_zeros (
                _mm_cmpeq_epi32 (isNodeCell, _mm_setzero_si128 ())
              , _mm_set1_epi32 (0xFFFFFFFF)) == 1
          , reason);
  SC3E_YES (reason);
}