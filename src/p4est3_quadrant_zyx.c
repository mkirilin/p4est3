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

static sc3_error_t *
p4est3_quadrant_zyx_child (const __m128i * q, int child_id, __m128i * r)
{
  const p4est_qcoord_t level = _mm_extract_epi32 (*q, 0);

  SC3A_IS (p4est3_quadrant_zyx_is_valid, q);
  SC3A_CHECK (level < P4EST_QMAXLEVEL);
  SC3A_CHECK (child_id >= 0 && child_id < P4EST_CHILDREN);

  *r = 
  _mm_or_si128 (
    _mm_sllv_epi32 (
      _mm_srlv_epi32 (
        _mm_and_si128 (_mm_set_epi32 (child_id, child_id, child_id, 0x00)
                     , _mm_set_epi32 (0x01, 0x02, 0x04, 0x00))
      , _mm_set_epi32 (0, 1, 2, 0)
      )
    , _mm_set1_epi32 (P4EST_MAXLEVEL - (level + 1))
    )
  , *q 
  );
  *r = _mm_insert_epi32 (*r, level + 1, 0);
#ifndef P4_TO_P8
  *r = _mm_insert_epi32 (*r, 0, 1);
#endif
#ifdef SC_ENABLE_DEBUG
  { //Waiting for the corresponding SC3A_.. macro
    char _r[SC3_BUFSIZE];
    if (!(p4est3_quadrant_zyx_is_parent (q, r, _r))) {
      char _errmsg[SC3_BUFSIZE];
      snprintf (_errmsg, SC3_BUFSIZE, "%s(%s, %s): %s",
                 "p4est3_quadrant_zyx_is_parent", "q", "r", _r);
      return sc3_error_new_fatal (__FILE__, __LINE__, _errmsg);
    }
  }
#endif //SC_ENABLE_DEBUG
  return NULL;
}

static sc3_error_t *
p4est3_quadrant_zyx_parent (const __m128i * q, __m128i * r)
{
  int32_t level = _mm_extract_epi32 (*q, 0);
  SC3A_IS (p4est3_quadrant_zyx_is_valid, q);
  SC3A_CHECK (level > 0);

  *r = _mm_sub_epi32 (
          _mm_and_si128 (*q, _mm_set_epi32 (~P4EST_QUADRANT_LEN (level)
                                          , ~P4EST_QUADRANT_LEN (level)
                                          , ~P4EST_QUADRANT_LEN (level)
                                          , 0xFFFFFFFF))
        , _mm_set_epi32 (0, 0, 0, 1) 
      );

  SC3A_IS (p4est3_quadrant_zyx_is_valid, r);
#ifdef SC_ENABLE_DEBUG
  { //Waiting for the corresponding SC3A_.. macro
    //There is no such check in original function
    char _r[SC3_BUFSIZE];
    if (!(p4est3_quadrant_zyx_is_parent (r, q, _r))) {
      char _errmsg[SC3_BUFSIZE];
      snprintf (_errmsg, SC3_BUFSIZE, "%s(%s, %s): %s",
                 "p4est3_quadrant_zyx_is_parent", "r", "q", _r);
      return sc3_error_new_fatal (__FILE__, __LINE__, _errmsg);
    }
  }
#endif //SC_ENABLE_DEBUG
  return NULL;
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

static sc3_error_t *
p4est3_quadrant_zyx_compare (const __m128i * q1, const __m128i * q2, int * j)
{
  int64_t diff;
  __m128i p, q;

  { //Waiting for the corresponding SC3A_.. macro
    SC3A_CHECK (p4est3_quadrant_zyx_is_node (q1, 1, NULL) ||
                  p4est3_quadrant_zyx_is_valid (q1, NULL));
    SC3A_CHECK (p4est3_quadrant_zyx_is_node (q2, 1, NULL) ||
                  p4est3_quadrant_zyx_is_valid (q2, NULL));  
  }
  __m128i exclor_coords = _mm_xor_si128 (*q1, *q2);
  __m128i exclor = _mm_set1_epi32 (_mm_extract_epi32 (exclor_coords, 3)
                                 | _mm_extract_epi32 (exclor_coords, 2));
#ifdef P4_TO_P8
  exclor = _mm_or_si128 (_mm_set_epi32 (0
                                      , 0
                                      , _mm_extract_epi32 (exclor_coords, 1)
                                      , 0)
                      , exclor);
#endif

  if (!_mm_extract_epi32 (exclor, 1)) {
    *j = _mm_cvtsi128_si32 (_mm_sub_epi32 (*q1, *q2));
    return NULL;
  }

  __m128i cond = _mm_cmpgt_epi32 (exclor_coords, _mm_xor_si128 (exclor, exclor_coords));
#ifdef P4_TO_P8
  if (_mm_extract_epi32 (cond, 1)) {
    q = _mm_set_epi64x ((int64_t) _mm_extract_epi32 (*q2, 1)
                      , (int64_t) _mm_extract_epi32 (*q1, 1));
  } else 
#if 0
    ;
#endif
#endif
  if (_mm_extract_epi32 (cond, 2)) {
    q = _mm_set_epi64x ((int64_t) _mm_extract_epi32 (*q2, 2)
                      , (int64_t) _mm_extract_epi32 (*q1, 2));
  } else {
    q = _mm_set_epi64x ((int64_t) _mm_extract_epi32 (*q2, 3)
                      , (int64_t) _mm_extract_epi32 (*q1, 3));
  }
  
  p = _mm_or_si128 (_mm_cmpeq_epi64 (q, _mm_setzero_si128 ())
                  , _mm_cmpgt_epi64 (q, _mm_setzero_si128 ()));

  q = _mm_add_epi64 (q
      , _mm_set_epi64x (
            _mm_extract_epi64 (p, 0) ? 0 : ((int64_t) 1 << (P4EST_MAXLEVEL + 2))
          , _mm_extract_epi64 (p, 1) ? 0 : ((int64_t) 1 << (P4EST_MAXLEVEL + 2))
        )
      );
  diff = _mm_extract_epi64 (q, 0) - _mm_extract_epi64 (q, 1);
  *j = ((diff == 0) ? 0 : ((diff < 0) ? -1 : 1));
  return NULL;
}

static sc3_error_t *
p4est3_quadrant_zyx_ancestor_id (const __m128i * q, int level, int * j)
{
  SC3A_IS (p4est3_quadrant_zyx_is_valid, q);
  SC3A_CHECK (0 <= level && level <= P4EST_MAXLEVEL);
  SC3A_CHECK (_mm_extract_epi32 (*q, 0) >= level);

  if (level == 0) {
    *j = 0;
    return NULL;
  }

  __m128i t = _mm_and_si128 (*q, _mm_set1_epi32 (P4EST_QUADRANT_LEN (level)));
  *j |= _mm_extract_epi32 (t, 3) ? 0x01 : 0;
  *j |= _mm_extract_epi32 (t, 2) ? 0x02 : 0;
#ifdef P4_TO_P8
  *j |= _mm_extract_epi32 (t, 1) ? 0x04 : 0;
#endif

  return NULL;
}

static sc3_error_t *
p4est3_quadrant_zyx_ancestor (const __m128i * q, int level, __m128i * r)
{
  SC3A_IS (p4est3_quadrant_zyx_is_valid, q);
  SC3A_CHECK (_mm_extract_epi32 (*q, 0) > level && level >= 0);

  *r = _mm_insert_epi32 (
          _mm_and_si128 (*q
                       , _mm_set1_epi32 (~(P4EST_QUADRANT_LEN (level) - 1)))
        , level
        , 0
      );
  SC3A_IS (p4est3_quadrant_zyx_is_valid, r);
  return NULL;
}

static sc3_error_t *
p4est3_quadrant_zyx_sibling (const __m128i * q, __m128i * r, int sibling_id)
{
  const p4est_qcoord_t q_level = _mm_extract_epi32 (*q, 0);
  const p4est_qcoord_t shift = P4EST_QUADRANT_LEN (q_level);

  const __m128i add = _mm_slli_epi32 (
                        _mm_srlv_epi32 (
                          _mm_and_si128 (_mm_set1_epi32 (sibling_id)
                                       , _mm_set_epi32 (0x01, 0x02, 0x04, 0))
                        , _mm_set_epi32 (0, 1, 2, 0))
                      , P4EST_MAXLEVEL - q_level
                      );

  SC3A_IS (p4est3_quadrant_zyx_is_valid, q);
  SC3A_CHECK (q_level > 0);
  SC3A_CHECK (sibling_id >= 0 && sibling_id < P4EST_CHILDREN);

  *r = _mm_or_si128 (add, _mm_andnot_si128 (_mm_set1_epi32 (shift), *q));
  *r = _mm_insert_epi32 (*r, q_level, 0);
  #ifndef P4_TO_P8
  *r = _mm_insert_epi32 (*r, 0, 1);
  #endif
  SC3A_IS (p4est3_quadrant_zyx_is_valid, r);
  return NULL;
}

static sc3_error_t *
p4est3_quadrant_zyx_first_descendant (const __m128i * q, int level
                                    , __m128i * fd)
{
  SC3A_IS (p4est3_quadrant_zyx_is_valid, q);
  SC3A_CHECK (_mm_extract_epi32 (*q, 0) <= level && level <= P4EST_MAXLEVEL);

  *fd = _mm_insert_epi32 (*q, level, 0);
  return NULL;
}

static sc3_error_t *
p4est3_quadrant_zyx_last_descendant (const __m128i * q, int level
                                    , __m128i * ld)
{
  p4est_qcoord_t shift;

  SC3A_IS (p4est3_quadrant_zyx_is_valid, q);
  SC3A_CHECK (_mm_extract_epi32 (*q, 0) <= level && level <= P4EST_QMAXLEVEL);

  shift = 
  P4EST_QUADRANT_LEN (_mm_extract_epi32 (*q, 0)) - P4EST_QUADRANT_LEN (level);

  *ld = _mm_insert_epi32 (
            _mm_add_epi32 (*q, _mm_set1_epi32 (shift))
          , level
          , 0
        );
  return NULL;
}
