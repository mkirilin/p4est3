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
#include <p4est.h>
#include <p4est3_quadrant_yx.h>
#else
#include <p8est.h>
#include <p4est3_quadrant_zyx.h>
#endif /* !P4_TO_P8 */

/* These have the same values for the 2D and 3D implementations */
#define P4EST3_YX_MAXLEVEL 30
#define P4EST3_YX_QMAXLEVEL 30

#ifdef P4EST_ENABLE_AVX2

#include <immintrin.h>
#include <smmintrin.h>
#include <emmintrin.h>

static              p4est3_gloidx
p4est3_quadrant_zyx_num_uniform (int level)
{
  if (level < 0) {
    return -1;
  }
  return p4est3_glopow (P4EST_CHILDREN, level);
}

static int
p4est3_quadrant_zyx_is_inside_root (const __m128i * q, char *reason)
{
/* *INDENT-OFF* */
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
/* *INDENT-ON* */
  SC3E_YES (reason);
}

static int
p4est3_quadrant_zyx_is_valid (const __m128i * q, char *reason)
{
  int32_t             level = _mm_extract_epi32 (*q, 0);
  SC3E_TEST (level >= 0, reason);
  SC3E_TEST (level <= P4EST3_YX_MAXLEVEL, reason);
/* *INDENT-OFF* */
  SC3E_TEST (
    _mm_testz_si128 (*q, _mm_set_epi32 (P4EST_QUADRANT_LEN (level) - 1
                                      , P4EST_QUADRANT_LEN (level) - 1
                                      , P4EST_QUADRANT_LEN (level) - 1
                                      , 0)) == 1
  , reason
  );
/* *INDENT-ON* */
  SC3E_IS (p4est3_quadrant_zyx_is_inside_root, q, reason);
  SC3E_YES (reason);
}

#ifdef P4EST_ENABLE_DEBUG

static int
p4est3_quadrant_zyx_is_parent (const __m128i * q, const __m128i * r,
                               char *reason)
{
  int32_t             r_level = _mm_extract_epi32 (*r, 0);

  SC3E_IS (p4est3_quadrant_zyx_is_valid, q, reason);
  SC3E_IS (p4est3_quadrant_zyx_is_valid, r, reason);

  SC3E_TEST (_mm_extract_epi32 (*q, 0) + 1 == r_level, reason);

/* *INDENT-OFF* */
  __m128i             lhs = _mm_add_epi32 (*q, _mm_set_epi32 (0, 0, 0, 1));
  __m128i             rhs = _mm_and_si128 (*r
                            , _mm_set_epi32 ( /* Optimized by compiler? */
                                             ~P4EST_QUADRANT_LEN (r_level)
                                           , ~P4EST_QUADRANT_LEN (r_level)
#ifdef P4_TO_P8
                                           , ~P4EST_QUADRANT_LEN (r_level)
#else
                                           , 0
#endif
                                           , 0xFFFFFFFF)
                            );
/* *INDENT-ON* */
  SC3E_TEST (_mm_testc_si128 (rhs, lhs) == 1, reason);
  SC3E_YES (reason);
}

#endif

static sc3_error_t *
p4est3_quadrant_zyx_is_ancestor (const __m128i * q, const __m128i * r, int *j)
{
  int                 lq, lr;
  __m128i             exclor;

  SC3A_IS (p4est3_quadrant_zyx_is_valid, q);
  SC3A_IS (p4est3_quadrant_zyx_is_valid, r);

  lq = _mm_extract_epi32 (*q, 0);
  lr = _mm_extract_epi32 (*r, 0);
  if (lq >= lr) {
    *j = 0;
    return NULL;
  }
/* *INDENT-OFF* */
  exclor = _mm_srlv_epi32 (
             _mm_xor_si128 (*q, *r)
           , _mm_set_epi32 (P4EST3_YX_MAXLEVEL - lq
                          , P4EST3_YX_MAXLEVEL - lq
                          , P4EST3_YX_MAXLEVEL - lq
                          , 32)
           );
/* *INDENT-ON* */
  *j = _mm_test_all_zeros (exclor, _mm_set1_epi32 (0xFFFFFFFF));
  return NULL;
}

static sc3_error_t *
p4est3_quadrant_zyx_child (const __m128i * q, int child_id, __m128i * r)
{
  const p4est_qcoord_t level = _mm_extract_epi32 (*q, 0);

  SC3A_IS (p4est3_quadrant_zyx_is_valid, q);
  SC3A_CHECK (level < P4EST3_YX_QMAXLEVEL);
  SC3A_CHECK (0 <= child_id && child_id < P4EST_CHILDREN);

/* *INDENT-OFF* */
  *r =
   _mm_or_si128 (
     _mm_sllv_epi32 (
       _mm_srlv_epi32 (
         _mm_and_si128 (_mm_set_epi32 (child_id, child_id, child_id, 0x00)
                      , _mm_set_epi32 (0x01, 0x02, 0x04, 0x00))
       , _mm_set_epi32 (0, 1, 2, 0)
       )
     , _mm_set1_epi32 (P4EST3_YX_MAXLEVEL - (level + 1))
     )
   , *q
   );
/* *INDENT-ON* */
  *r = _mm_insert_epi32 (*r, level + 1, 0);
#ifndef P4_TO_P8
  *r = _mm_insert_epi32 (*r, 0, 1);
#endif
  SC3A_IS2 (p4est3_quadrant_zyx_is_parent, q, r);
  return NULL;
}

static sc3_error_t *
p4est3_quadrant_zyx_parent (const __m128i * q, __m128i * r)
{
  int32_t             level = _mm_extract_epi32 (*q, 0);
  SC3A_IS (p4est3_quadrant_zyx_is_valid, q);
  SC3A_CHECK (level > 0);

/* *INDENT-OFF* */
  *r =
   _mm_sub_epi32 (
     _mm_and_si128 (*q, _mm_set_epi32 (~P4EST_QUADRANT_LEN (level)
                                     , ~P4EST_QUADRANT_LEN (level)
                                     , ~P4EST_QUADRANT_LEN (level)
                                     , 0xFFFFFFFF))
   , _mm_set_epi32 (0, 0, 0, 1)
   );
/* *INDENT-ON* */

  SC3A_IS (p4est3_quadrant_zyx_is_valid, r);
  SC3A_IS2 (p4est3_quadrant_zyx_is_parent, r, q);
  return NULL;
}

static sc3_error_t *
p4est3_quadrant_zyx_copy (const __m128i * q, __m128i * copy)
{
  SC3A_IS (p4est3_quadrant_zyx_is_valid, q);
  *copy = *q;
  return NULL;
}

static sc3_error_t *
p4est3_quadrant_zyx_compare (const __m128i * q1, const __m128i * q2, int *j)
{
  int64_t             diff;
  __m128i             p, q;

/* *INDENT-OFF* */
  { /* Waiting for the corresponding SC3A_.. macro */
    SC3A_CHECK (p4est3_quadrant_zyx_is_valid (q1, NULL));
    SC3A_CHECK (p4est3_quadrant_zyx_is_valid (q2, NULL));
  }
  __m128i            exclor_coords = _mm_xor_si128 (*q1, *q2);
  __m128i            exclor =
                      _mm_set1_epi32 (_mm_extract_epi32 (exclor_coords, 3)
                                    | _mm_extract_epi32 (exclor_coords, 2));
#ifdef P4_TO_P8
  exclor = _mm_or_si128 (
             _mm_set_epi32 (0, 0, _mm_extract_epi32 (exclor_coords, 1), 0)
          , exclor);
#endif

  if (!_mm_extract_epi32 (exclor, 1)) {
    *j = _mm_cvtsi128_si32 (_mm_sub_epi32 (*q1, *q2));
    return NULL;
  }

  __m128i             cond = _mm_cmpgt_epi32 (
                               exclor_coords,
                              _mm_xor_si128 (exclor, exclor_coords)
                             );
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
            _mm_extract_epi64 (p, 0) ? 0 : ((int64_t) 1 << (P4EST3_YX_MAXLEVEL + 2))
          , _mm_extract_epi64 (p, 1) ? 0 : ((int64_t) 1 << (P4EST3_YX_MAXLEVEL + 2))
        )
      );
/* *INDENT-ON* */
  diff = _mm_extract_epi64 (q, 0) - _mm_extract_epi64 (q, 1);
  *j = ((diff == 0) ? 0 : ((diff < 0) ? -1 : 1));
  return NULL;
}

static sc3_error_t *
p4est3_quadrant_zyx_ancestor_id (const __m128i * q, int level, int *j)
{
  SC3A_IS (p4est3_quadrant_zyx_is_valid, q);
  SC3A_CHECK (0 <= level && level <= P4EST3_YX_MAXLEVEL);
  SC3A_CHECK (_mm_extract_epi32 (*q, 0) >= level);
  *j = 0;

  if (level == 0) {
    return NULL;
  }

  __m128i             t =
    _mm_and_si128 (*q, _mm_set1_epi32 (P4EST_QUADRANT_LEN (level)));
  *j |= _mm_extract_epi32 (t, 3) ? 0x01 : 0;
  *j |= _mm_extract_epi32 (t, 2) ? 0x02 : 0;
#ifdef P4_TO_P8
  *j |= _mm_extract_epi32 (t, 1) ? 0x04 : 0;
#endif

  return NULL;
}

static sc3_error_t *
p4est3_quadrant_zyx_child_id (const __m128i * q, int *j)
{
  int                 level = _mm_extract_epi32 (*q, 0);
  SC3E (p4est3_quadrant_zyx_ancestor_id (q, level, j));
  return NULL;
}

static sc3_error_t *
p4est3_quadrant_zyx_coordinates (const __m128i * q, int n, int *j)
{
  __m128i             r;
  int                 d = P4EST3_REF_MAXLEVEL - P4EST3_YX_MAXLEVEL;
  SC3A_IS (p4est3_quadrant_zyx_is_valid, q);
  SC3A_CHECK (n == P4EST_DIM);
  SC3A_CHECK (d >= 0);

  r = _mm_slli_epi32 (*q, d);
  j[0] = _mm_extract_epi32 (r, 3);
  j[1] = _mm_extract_epi32 (r, 2);
#ifdef P4_TO_P8
  j[2] = _mm_extract_epi32 (r, 1);
#endif
  return NULL;
}

static sc3_error_t *
p4est3_quadrant_zyx_level (const __m128i * q, int *l)
{
  SC3A_IS (p4est3_quadrant_zyx_is_valid, q);
  *l = _mm_extract_epi32 (*q, 0);
  return NULL;
}

static sc3_error_t *
p4est3_quadrant_zyx_ancestor (const __m128i * q, int level, __m128i * r)
{
  SC3A_IS (p4est3_quadrant_zyx_is_valid, q);
  SC3A_CHECK (_mm_extract_epi32 (*q, 0) > level && level >= 0);
/* *INDENT-OFF* */
  *r =
    _mm_insert_epi32 (
     _mm_and_si128 (*q, _mm_set1_epi32 (~(P4EST_QUADRANT_LEN (level) - 1))),
     level, 0);
/* *INDENT-ON* */
  SC3A_IS (p4est3_quadrant_zyx_is_valid, r);
  return NULL;
}

static sc3_error_t *
p4est3_quadrant_zyx_sibling (const __m128i * q, __m128i * r, int sibling_id)
{
  const p4est_qcoord_t q_level = _mm_extract_epi32 (*q, 0);
  const p4est_qcoord_t shift = P4EST_QUADRANT_LEN (q_level);
/* *INDENT-OFF* */
  const __m128i       add =
                      _mm_slli_epi32 (
                       _mm_srlv_epi32 (
                        _mm_and_si128 (_mm_set1_epi32 (sibling_id)
                                     , _mm_set_epi32 (0x01, 0x02, 0x04, 0))
                      , _mm_set_epi32 (0, 1, 2, 0))
                     , P4EST3_YX_MAXLEVEL - q_level);
/* *INDENT-ON* */
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
p4est3_quadrant_zyx_first_descendant (const __m128i * q, int level,
                                      __m128i * fd)
{
  SC3A_IS (p4est3_quadrant_zyx_is_valid, q);
  SC3A_CHECK (_mm_extract_epi32 (*q, 0) <= level &&
              level <= P4EST3_YX_MAXLEVEL);

  *fd = _mm_insert_epi32 (*q, level, 0);
  return NULL;
}

static sc3_error_t *
p4est3_quadrant_zyx_last_descendant (const __m128i * q, int level,
                                     __m128i * ld)
{
  p4est_qcoord_t      shift;

  SC3A_IS (p4est3_quadrant_zyx_is_valid, q);
  SC3A_CHECK (_mm_extract_epi32 (*q, 0) <= level &&
              level <= P4EST3_YX_QMAXLEVEL);

  shift = P4EST_QUADRANT_LEN (_mm_extract_epi32 (*q, 0))
    - P4EST_QUADRANT_LEN (level);

  *ld = _mm_insert_epi32 (_mm_add_epi32 (*q, _mm_set1_epi32 (shift)),
                          level, 0);
  return NULL;
}

static sc3_error_t *
p4est3_zyx_nearest_common_ancestor (const __m128i * q1, const __m128i * q2,
                                    __m128i * r)
{
  int                 maxlevel, min;
  __m128i             exclor;
  int32_t             maxclor;

  SC3A_IS (p4est3_quadrant_zyx_is_valid, q1);
  SC3A_IS (p4est3_quadrant_zyx_is_valid, q2);

  exclor = _mm_xor_si128 (*q1, *q2);
#ifdef P4_TO_P8
/* *INDENT-OFF* */
  maxclor = _mm_extract_epi32 (exclor, 3)
          | _mm_extract_epi32 (exclor, 2)
          | _mm_extract_epi32 (exclor, 1);
/* *INDENT-ON* */
#else
/* *INDENT-OFF* */
  maxclor = _mm_extract_epi32 (exclor, 3)
          | _mm_extract_epi32 (exclor, 2);
/* *INDENT-ON* */
#endif
  maxlevel = SC_LOG2_32 (maxclor) + 1;

  SC3A_CHECK (maxlevel <= P4EST3_YX_MAXLEVEL);
  *r = _mm_and_si128 (*q1, _mm_set1_epi32 (~((1 << maxlevel) - 1)));
/* *INDENT-OFF* */
  min = (int) SC_MIN (P4EST3_YX_MAXLEVEL - maxlevel,
                       SC_MIN (_mm_extract_epi32 (*q1, 0)
                             , _mm_extract_epi32 (*q2, 0)));
/* *INDENT-ON* */
  *r = _mm_insert_epi32 (*r, min, 0);
  SC3A_IS (p4est3_quadrant_zyx_is_valid, r);
  return NULL;
}

static sc3_error_t *
p4est3_quadrant_zyx_successor (const __m128i * q, __m128i * r)
{
  int                 level, q_level;
  int                 successor_id = 0;
  q_level = level = _mm_extract_epi32 (*q, 0);

  SC3A_IS (p4est3_quadrant_zyx_is_valid, q);
  SC3A_CHECK (q_level > 0);

  SC3E (p4est3_quadrant_zyx_ancestor_id (q, level, &successor_id));
  successor_id++;

  /* iterate until it is possible to increment the child/ancestor_id */
  while (successor_id == P4EST_CHILDREN) {
    SC3E (p4est3_quadrant_zyx_ancestor_id (q, --level, &successor_id));
    successor_id++;
    SC3A_CHECK (level > 0);
  }
  SC3A_CHECK (0 < successor_id && successor_id < P4EST_CHILDREN);

  /* compute result */
  if (level < q_level) {
    /* coarsen to level - 1 and add shifts according to the successor_id */
/* *INDENT-OFF* */
    *r =
    _mm_add_epi32 (
      _mm_sllv_epi32 (
        _mm_srlv_epi32 (
          _mm_and_si128 (_mm_set1_epi32 (successor_id)
                       , _mm_set_epi32  (0x01, 0x02, 0x04, 0x00))
        , _mm_set_epi32 (0, 1, 2, 0)
        )
      , _mm_set1_epi32 (P4EST3_YX_MAXLEVEL - level)
      )
    , _mm_and_si128 (*q, _mm_set_epi32 (~(P4EST_QUADRANT_LEN (level - 1) - 1)
                                      , ~(P4EST_QUADRANT_LEN (level - 1) - 1)
                                      , ~(P4EST_QUADRANT_LEN (level - 1) - 1)
                                      , 0xFFFFFFFF))
    );
/* *INDENT-ON* */
  }
  else {
    SC3E (p4est3_quadrant_zyx_sibling (q, r, successor_id));
  }
  SC3A_IS (p4est3_quadrant_zyx_is_valid, r);
  return NULL;
}

static sc3_error_t *
p4est3_quadrant_zyx_predecessor (const __m128i * q, __m128i * r)
{
  int                 level, q_level;
  int                 predecessor_id;
  q_level = level = _mm_extract_epi32 (*q, 0);

  SC3A_IS (p4est3_quadrant_zyx_is_valid, q);
  SC3A_CHECK (q_level > 0);

  SC3E (p4est3_quadrant_zyx_ancestor_id (q, level, &predecessor_id));
  predecessor_id--;

  /* iterate until it is possible to decrement the child/ancestor_id */
  while (predecessor_id == -1) {
    SC3E (p4est3_quadrant_zyx_ancestor_id (q, --level, &predecessor_id));
    predecessor_id--;
    SC3A_CHECK (level > 0);
  }
  SC3A_CHECK (0 <= predecessor_id && predecessor_id < P4EST_CHILDREN - 1);

  /* compute result */
  if (level < q_level) {
    /* coarsen to level - 1 and add shifts according to the predecessor_id */
/* *INDENT-OFF* */
    *r =
    _mm_add_epi32 (
      _mm_sllv_epi32 (
        _mm_srlv_epi32 (
          _mm_and_si128 (_mm_set1_epi32 (predecessor_id)
                       , _mm_set_epi32  (0x01, 0x02, 0x04, 0x00))
        , _mm_set_epi32 (0, 1, 2, 0)
        )
      , _mm_set1_epi32 (P4EST3_YX_MAXLEVEL - level)
      )
    , _mm_and_si128 (*q, _mm_set_epi32 (~(P4EST_QUADRANT_LEN (level - 1) - 1)
                                      , ~(P4EST_QUADRANT_LEN (level - 1) - 1)
                                      , ~(P4EST_QUADRANT_LEN (level - 1) - 1)
                                      , 0xFFFFFFFF))
    );
/* *INDENT-ON* */
  }
  else {
    SC3E (p4est3_quadrant_zyx_sibling (q, r, predecessor_id));
  }
  SC3A_IS (p4est3_quadrant_zyx_is_valid, r);
  return NULL;
}

static sc3_error_t *
p4est3_quadrant_zyx_linear_id (const __m128i * quadrant, int level,
                               p4est3_gloidx * id)
{
  int                 i;
  __m128i             shifted, res;

  SC3A_IS (p4est3_quadrant_zyx_is_valid, quadrant);
  SC3A_CHECK (0 <= level && level <= P4EST3_YX_MAXLEVEL);

  /* this preserves the high bits from negative numbers */
  shifted =
    _mm_srlv_epi32 (*quadrant, _mm_set1_epi32 (P4EST3_YX_MAXLEVEL - level));

  *id = (p4est3_gloidx) 0;
  for (i = 0; i < level; ++i) {
/* *INDENT-OFF* */
    res = _mm_sllv_epi32 (
            _mm_and_si128 (shifted, _mm_set1_epi32 ((uint32_t) 1 << i))
          , _mm_set_epi32 ((P4EST_DIM - 1) * i
                         , (P4EST_DIM - 1) * i + 1
                         , (P4EST_DIM - 1) * i + 2
                         , 0)
          );
/* *INDENT-ON* */
    *id |= _mm_extract_epi32 (res, 3);
    *id |= _mm_extract_epi32 (res, 2);
    *id |= _mm_extract_epi32 (res, 1);
  }
  SC3A_CHECK (*id >= 0L && *id < ((int64_t) 1 << P4EST_DIM * level));
  return NULL;
}

static sc3_error_t *
p4est3_quadrant_zyx_morton (int level, p4est3_gloidx id, __m128i * quadrant)
{
  int                 i;

  SC3A_CHECK (0 <= level && level <= P4EST3_YX_QMAXLEVEL);
  if (level < P4EST3_YX_QMAXLEVEL) {
    SC3A_CHECK (id < (((p4est3_gloidx) 1) << P4EST_DIM * level));
  }

  *quadrant = _mm_setzero_si128 ();

  /* this may set the sign bit to create negative numbers */
  for (i = 0; i < level; ++i) {
/* *INDENT-OFF* */
    __m128i xy_coord_id =
    _mm_srlv_epi64 (
      _mm_and_si128 (
        _mm_set1_epi64x (id)
      , _mm_sllv_epi64 (
          _mm_set1_epi64x (1ULL)
        , _mm_set_epi64x (P4EST_DIM * i
                        , P4EST_DIM * i + 1)
        )
      )
    , _mm_set_epi64x ((P4EST_DIM - 1) * i
                    , (P4EST_DIM - 1) * i + 1)
    );
    *quadrant
    = _mm_or_si128 (*quadrant
                  , _mm_set_epi32 (
                      (p4est_qcoord_t) _mm_extract_epi64 (xy_coord_id, 1)
                    , (p4est_qcoord_t) _mm_extract_epi64 (xy_coord_id, 0)
#ifdef P4_TO_P8
                    , (p4est_qcoord_t) ((id & (1ULL << (P4EST_DIM * i + 2)))
                        >> ((P4EST_DIM - 1) * i + 2))
#else
                    , 0
#endif /* P4_TO_P8 */
                    , 0
                  ));
/* *INDENT-ON* */
  }

  *quadrant = _mm_slli_epi32 (*quadrant, P4EST3_YX_MAXLEVEL - level);
  *quadrant = _mm_insert_epi32 (*quadrant, level, 0);
  SC3A_IS (p4est3_quadrant_zyx_is_valid, quadrant);
  return NULL;
}

static sc3_error_t *
p4est3_quadrant_zyx_root (__m128i * r)
{
  SC3E (p4est3_quadrant_zyx_morton (0, 0, r));
  return NULL;
}

#endif /* P4EST_ENABLE_AVX2 */

sc3_error_t        *
p4est3_quadrant_yx_vtable (p4est3_quadrant_vtable_t * qvt)
{
  SC3A_CHECK (qvt != NULL);
  memset (qvt, 0, sizeof (p4est3_quadrant_vtable_t));

#ifdef P4EST_ENABLE_AVX2
  qvt->dim = P4EST_DIM;
  qvt->max_level = P4EST3_YX_MAXLEVEL;
  qvt->max_children = P4EST_CHILDREN;
  qvt->quadrant_size = sizeof (__m128i);

  qvt->quadrant_num_uniform = p4est3_quadrant_zyx_num_uniform;

  qvt->quadrant_is_valid =
    (p4est3_quadrant_is_t) p4est3_quadrant_zyx_is_valid;

  qvt->quadrant_level = (p4est3_quadrant_level_t) p4est3_quadrant_zyx_level;

  qvt->quadrant_child_id =
    (p4est3_quadrant_child_id_t) p4est3_quadrant_zyx_child_id;

  qvt->quadrant_ancestor_id =
    (p4est3_quadrant_ancestor_id_t) p4est3_quadrant_zyx_ancestor_id;

  qvt->quadrant_coordinates =
    (p4est3_quadrant_in_i_out_t) p4est3_quadrant_zyx_coordinates;

  qvt->quadrant_compare =
    (p4est3_quadrant_compare_t) p4est3_quadrant_zyx_compare;

  qvt->quadrant_root = (p4est3_quadrant_root_t) p4est3_quadrant_zyx_root;

  qvt->quadrant_copy = (p4est3_quadrant_copy_t) p4est3_quadrant_zyx_copy;

  qvt->quadrant_parent =
    (p4est3_quadrant_parent_t) p4est3_quadrant_zyx_parent;

  qvt->quadrant_predecessor =
    (p4est3_quadrant_predecessor_t) p4est3_quadrant_zyx_predecessor;

  qvt->quadrant_successor =
    (p4est3_quadrant_successor_t) p4est3_quadrant_zyx_successor;

  qvt->quadrant_child = (p4est3_quadrant_child_t) p4est3_quadrant_zyx_child;

  qvt->quadrant_ancestor =
    (p4est3_quadrant_ancestor_t) p4est3_quadrant_zyx_ancestor;

  qvt->quadrant_first_descendant =
    (p4est3_quadrant_first_descendant_t) p4est3_quadrant_zyx_first_descendant;

  qvt->quadrant_last_descendant =
    (p4est3_quadrant_last_descendant_t) p4est3_quadrant_zyx_last_descendant;

  qvt->quadrant_morton =
    (p4est3_quadrant_morton_t) p4est3_quadrant_zyx_morton;

  qvt->nearest_common_ancestor =
    (p4est3_nearest_common_ancestor_t) p4est3_zyx_nearest_common_ancestor;

  qvt->quadrant_linear_id =
    (p4est3_quadrant_linear_id_t) p4est3_quadrant_zyx_linear_id;

  qvt->quadrant_is_ancestor =
    (p4est3_quadrant_is_ancestor_t) p4est3_quadrant_zyx_is_ancestor;

  SC3A_IS (p4est3_quadrant_vtable_is_valid, qvt);
  return NULL;
#else
  return sc3_error_new_kind (SC3_ERROR_RUNTIME, __FILE__, __LINE__, "Creation"
                             " of AVX2-based virtual table is denied since AVX2"
                             " is disabled or not found working");
#endif /* !P4EST_ENABLE_AVX2 */
}
