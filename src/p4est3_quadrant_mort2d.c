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

#define CRD_BITS 56
#ifndef P4_TO_P8
#include <p4est3_quadrant_mort2d.h>
 /*(CRD_BITS - 8) / 2*/
#define P4EST3_MORT_MAXLEVEL 28
#define P4EST3_MORT_QMAXLEVEL 28
/* This mask is used to extract a coordinate at 0 position.
  It's binary representation is 0..0|0101..01 (00 x8|01 x28)*/
#define P4EST3_MORT_COORD_MASK ((89478485UL << 28) | 89478485UL)
#else
#include <p4est3_quadrant_mort3d.h>
 /*(CRD_BITS - 8) / 3*/
#define P4EST3_MORT_MAXLEVEL 18
#define P4EST3_MORT_QMAXLEVEL 18
/* This mask is used to extract a coordinate at 0 position.
  It's binary representation is 0..0|001001001..001 (0 x10| 001 x18)*/
#define P4EST3_MORT_COORD_MASK ((((37449UL << 18) | 37449UL) << 18) | 37449UL)
#endif

#define P4EST3_QUADRANT_MORT_LEN(n, l) \
        ((uint64_t) (n) << P4EST_DIM * (P4EST3_MORT_MAXLEVEL - (l)))
#define P4EST3_ROOT_MORT_LEN ((p4est_qcoord_t) 1 << P4EST3_MORT_MAXLEVEL)


/*(64 - CRD_BITS) highest bits for level, CRD_BITS bits for coordinates*/
typedef uint64_t p4est3_quadrant_mort_t;
/* Extract level from compressed morton type */
#define P4EST3_MORT_EXT_LEVEL(q) (((uint64_t) q) >> 56)

static              p4est3_gloidx
p4est3_quadrant_mort_num_uniform (int level)
{
  if (level < 0) {
    return -1;
  }
  return p4est3_glopow (P4EST_CHILDREN, level);
}

static int
p4est3_quadrant_mort_is_valid (const p4est3_quadrant_mort_t * q, char *reason)
{
  const int level = P4EST3_MORT_EXT_LEVEL(*q);
  SC3E_TEST ((level >= 0 && level <= P4EST3_MORT_QMAXLEVEL) &&
             ((*q & (P4EST3_QUADRANT_MORT_LEN (0x01, level) - 1)) ==
              0), reason);
  SC3E_YES (reason);
}

static sc3_error_t *
p4est3_quadrant_mort_coords (const p4est3_quadrant_mort_t * q,
                             int n, p4est_qcoord_t * coords)
{
  const int level = P4EST3_MORT_EXT_LEVEL(*q);
  int                 i, forward, backward;
  int                 d = P4EST3_REF_MAXLEVEL - P4EST3_MORT_MAXLEVEL;

  SC3A_IS (p4est3_quadrant_mort_is_valid, q);
  SC3A_CHECK (n == P4EST_DIM);
  SC3A_CHECK (d >= 0);

  coords[0] = 0;
  coords[1] = 0;
#ifdef P4_TO_P8
  coords[2] = 0;
#endif /* P4_TO_P8 */
  for (i = 1; i < level + 2; ++i) {
    forward = P4EST_DIM * (P4EST3_MORT_MAXLEVEL - i);
    backward = forward - (P4EST3_MORT_MAXLEVEL - i);
    coords[0] |=
      (p4est_qcoord_t) ((*q & (1ULL << forward)) >> backward);
    coords[1] |=
      (p4est_qcoord_t) ((*q & (1ULL << (forward + 1))) >> (backward +
                                                                  1));
#ifdef P4_TO_P8
    coords[2] |=
      (p4est_qcoord_t) ((*q & (1ULL << (forward + 2))) >> (backward +
                                                                  2));
#endif /* P4_TO_P8 */
  }

  /**Check if coordinates belong to a root.
   * Use (coords[0] & P4EST3_ROOT_MORT_LEN) == 0 instead of comparison due to
   * the lack of bits in signed int.
  */
  SC3A_CHECK (0 <= coords[0] && (coords[0] & P4EST3_ROOT_MORT_LEN) == 0);
  SC3A_CHECK (0 <= coords[1] && (coords[1] & P4EST3_ROOT_MORT_LEN) == 0);
#ifdef P4_TO_P8
  SC3A_CHECK (0 <= coords[2] && (coords[2] & P4EST3_ROOT_MORT_LEN) == 0);
#endif /* P4_TO_P8 */
  coords[0] <<= d;
  coords[1] <<= d;
#ifdef P4_TO_P8
  coords[2] <<= d;
#endif
  return NULL;
}

static sc3_error_t *
p4est3_quadrant_mort_quadrant (const p4est_qcoord_t * c, int l,
                               p4est3_quadrant_mort_t * q)
{
  int                 i;
  p4est_qcoord_t      x, y;
#ifdef P4_TO_P8
  p4est_qcoord_t      z;
#endif
  const int           d = P4EST3_REF_MAXLEVEL - P4EST3_MORT_MAXLEVEL;
  SC3A_CHECK (c != NULL);
  SC3A_CHECK (q != NULL);
  SC3A_CHECK (0 <= l && l <= P4EST3_MORT_MAXLEVEL);
  /* Since the coordinates are normalized by the P4EST3_REF_MAXLEVEL,
     we shift them according to qvt maxlevel */
  x = (c[0] >> d) >> (P4EST3_MORT_MAXLEVEL - l);
  y = (c[1] >> d) >> (P4EST3_MORT_MAXLEVEL - l);
#ifdef P4_TO_P8
  z = (c[2] >> d) >> (P4EST3_MORT_MAXLEVEL - l);
#endif
  *q = 0;
  for (i = 0; i < l; ++i) {
    *q |= ((x & ((p4est_qcoord_t) 1 << i)) << ((P4EST_DIM - 1) * i));
    *q |= ((y & ((p4est_qcoord_t) 1 << i)) << ((P4EST_DIM - 1) * i + 1));
#ifdef P4_TO_P8
    *q |= ((z & ((p4est_qcoord_t) 1 << i)) << ((P4EST_DIM - 1) * i + 2));
#endif
  }
  *q = P4EST3_QUADRANT_MORT_LEN (*q, l);
  *q |= ((uint64_t) l) << CRD_BITS;
  SC3A_IS (p4est3_quadrant_mort_is_valid, q);
  return NULL;
}

#if 0

static              p4est_qcoord_t
p4est3_quadrant_mort_coord_noerr (const p4est3_quadrant_mort_t * q,
                                  int32_t dim, char *reason)
{
  int                 i, forward, backward;
  p4est_qcoord_t      coord;
  if (dim > P4EST_DIM) {
    return -1;
  }
  if (!p4est3_quadrant_mort_is_valid (q, reason)) {
    return -1;
  }

  coord = 0;
  for (i = 1; i < q->level + 2; ++i) {
    forward = P4EST_DIM * (P4EST3_MORT_MAXLEVEL - i) + dim;
    backward = forward - (P4EST3_MORT_MAXLEVEL - i);
    coord |= (p4est_qcoord_t) ((q->coords & (1ULL << forward)) >> backward);
  }

  /* Please explain the issue */
  /**???*/
  if (0 > coord || coord >= P4EST3_ROOT_MORT_LEN) {
    return -1;
  }
  return coord;
}

#endif

static sc3_error_t *
p4est3_quadrant_mort_level (const p4est3_quadrant_mort_t * q, int *l)
{
  SC3A_IS (p4est3_quadrant_mort_is_valid, q);
  *l = P4EST3_MORT_EXT_LEVEL(*q);
  return NULL;
}

#ifdef P4EST_ENABLE_DEBUG

static int
p4est3_quadrant_mort_is_parent (const p4est3_quadrant_mort_t * q,
                                const p4est3_quadrant_mort_t * r,
                                char *reason)
{
  const int r_level = P4EST3_MORT_EXT_LEVEL(*r);
  int32_t             mask;
  SC3E_IS (p4est3_quadrant_mort_is_valid, q, reason);
  SC3E_IS (p4est3_quadrant_mort_is_valid, r, reason);
#ifdef P4_TO_P8
  mask = 0x07;
#else
  mask = 0x03;
#endif
  SC3E_TEST ((*q + ((uint64_t) 1 << CRD_BITS) == *r & ~P4EST3_QUADRANT_MORT_LEN (mask, r_level)), reason);
  SC3E_YES (reason);
}

#endif

static sc3_error_t *
p4est3_quadrant_mort_is_ancestor (const p4est3_quadrant_mort_t * q,
                                  const p4est3_quadrant_mort_t * r,
                                  int32_t * j)
{
  const int q_level = P4EST3_MORT_EXT_LEVEL(*q);
  p4est_qcoord_t      exclor;

  SC3A_IS (p4est3_quadrant_mort_is_valid, q);
  SC3A_IS (p4est3_quadrant_mort_is_valid, r);

  if (q_level >= P4EST3_MORT_EXT_LEVEL(*r)) {
    *j = 0;
    return NULL;
  }

  exclor =
    ((*q ^ *r) << 8) >> (P4EST_DIM * (P4EST3_MORT_MAXLEVEL - q_level) + 8);
  *j = exclor == 0 ? 1 : 0;
  return NULL;
}

static sc3_error_t *
p4est3_quadrant_mort_get_tree_boundary (const p4est3_quadrant_mort_t * q,
                                        int face, int *j)
{
  const uint64_t      l_mask =
    ~(P4EST3_QUADRANT_MORT_LEN (0x01, P4EST3_MORT_EXT_LEVEL(*q)) - 1);
  const uint64_t      my_coord_mask =
    (P4EST3_MORT_COORD_MASK & l_mask) << (face / 2);
  const uint64_t      bound = face % 2 == 0 ? 0 : my_coord_mask;

  *j = ((*q & my_coord_mask) == bound);
  return NULL;
}

static sc3_error_t *
p4est3_quadrant_mort_tree_boundaries (const p4est3_quadrant_mort_t * q,
                                      sc3_array_t * nf)
{
  const int level = P4EST3_MORT_EXT_LEVEL(*q);
  SC3A_CHECK (nf != NULL);
  SC3A_IS (p4est3_quadrant_mort_is_valid, q);
  const uint64_t      l_mask =
    ~(P4EST3_QUADRANT_MORT_LEN (0x01, level) - 1);
  const uint64_t      x_mask = P4EST3_MORT_COORD_MASK & l_mask;
  const uint64_t      y_mask = x_mask << 1;
  const uint64_t      x_extracted = *q & x_mask;
  const uint64_t      y_extracted = *q & y_mask;
#ifdef P4_TO_P8
  const uint64_t      z_mask = y_mask << 1;
  const uint64_t      z_extracted = *q & z_mask;
#endif

  int                *x, *y;
#ifdef P4_TO_P8
  int                *z;
#endif

  SC3E (sc3_array_index (nf, 0, &x));
  SC3E (sc3_array_index (nf, 1, &y));
#ifdef P4_TO_P8
  SC3E (sc3_array_index (nf, 2, &z));
#endif

  if (level == 0) {
    *x = *y = -2;
#ifdef P4_TO_P8
    *z = -2;
#endif
    return NULL;
  }

  *x = x_extracted == 0 ? 0 : (x_extracted == x_mask) ? 1 : -1;
  *y = y_extracted == 0 ? 2 : (y_extracted == y_mask) ? 3 : -1;
#ifdef P4_TO_P8
  *z = z_extracted == 0 ? 4 : (z_extracted == z_mask) ? 5 : -1;
#endif

  return NULL;
}

static sc3_error_t *
p4est3_quadrant_mort_child (const p4est3_quadrant_mort_t * q,
                            int32_t child_id, p4est3_quadrant_mort_t * r)
{
  const int level = P4EST3_MORT_EXT_LEVEL(*q);
  const uint64_t      shift = P4EST3_QUADRANT_MORT_LEN (0x01, level + 1);

  SC3A_IS (p4est3_quadrant_mort_is_valid, q);
  SC3A_CHECK (level < P4EST3_MORT_QMAXLEVEL);
  SC3A_CHECK (child_id >= 0 && child_id < P4EST_CHILDREN);

  *r = child_id & 0x01 ? (*q | shift) : *q;
  *r = child_id & 0x02 ? (*r | (shift << 1)) : *r;
#ifdef P4_TO_P8
  *r = child_id & 0x04 ? (*r | (shift << 2)) : *r;
#endif
  *r += ((uint64_t) 1 << CRD_BITS);
  SC3A_IS2 (p4est3_quadrant_mort_is_parent, q, r);
  return NULL;
}

static sc3_error_t *
p4est3_quadrant_mort_parent (const p4est3_quadrant_mort_t * q,
                             p4est3_quadrant_mort_t * r)
{
  const int level = P4EST3_MORT_EXT_LEVEL(*q);
  SC3A_IS (p4est3_quadrant_mort_is_valid, q);
  SC3A_CHECK (level > 0);

  int32_t             mask;
#ifdef P4_TO_P8
  mask = 0x07;
#else
  mask = 0x03;
#endif
  *r = *q & ~P4EST3_QUADRANT_MORT_LEN (mask, level);
  *r -= ((uint64_t) 1 << CRD_BITS);
  SC3A_IS (p4est3_quadrant_mort_is_valid, r);
  return NULL;
}

static sc3_error_t *
p4est3_quadrant_mort_face_neighbor (const p4est3_quadrant_mort_t * q,
                                    int i, p4est3_quadrant_mort_t * r)
{
  SC3A_IS (p4est3_quadrant_mort_is_valid, q);
  SC3A_CHECK (0 <= i && i < P4EST_FACES);

  const int8_t        sign = i & 0x01 ? 1 : -1;
  const uint64_t      l_mask =
    ~(P4EST3_QUADRANT_MORT_LEN (0x01, P4EST3_MORT_EXT_LEVEL(*q)) - 1);
  const uint64_t      dir_mask = (P4EST3_MORT_COORD_MASK & l_mask) << (i / 2);

  if (sign == 1) {
    *r = (*q | ~dir_mask) + 1;
  }
  else {
    *r = (*q & dir_mask) - 1;
  }
  *r = (*r & dir_mask) | (*q & ~dir_mask);
  SC3A_IS (p4est3_quadrant_mort_is_valid, r);
  return NULL;
}

static sc3_error_t *
p4est3_quadrant_mort_tree_face_neighbor (const p4est3_quadrant_mort_t * q,
                                         sc3_array_t * transform,
                                         int face, p4est3_quadrant_mort_t * r)
{
  const uint64_t      l_mask =
    ~(P4EST3_QUADRANT_MORT_LEN (0x01, q->level) - 1);
  const uint64_t      dir_level_mask = P4EST3_MORT_COORD_MASK & l_mask;
  int                *my_axis;
  int                *target_axis;
  int                *edge_reverse;
  uint64_t            extracted, mh, Rmh;

  SC3A_IS (p4est3_quadrant_mort_is_valid, q);
  SC3A_CHECK (q != r);

  SC3E (sc3_array_index (transform, 0, &my_axis));
  SC3E (sc3_array_index (transform, 3, &target_axis));
  SC3E (sc3_array_index (transform, 6, &edge_reverse));

#ifdef P4EST_ENABLE_DEBUG
  int                 i;
  for (i = 0; i < 3; ++i) {
    SC3A_CHECK (0 <= my_axis[i] && my_axis[i] < P4EST_DIM);
    SC3A_CHECK (0 <= target_axis[i] && target_axis[i] < P4EST_DIM);
  }
  SC3A_CHECK (my_axis[0] != my_axis[2]);
  SC3A_CHECK (target_axis[0] != target_axis[2]);
  SC3A_CHECK (0 <= edge_reverse[0] && edge_reverse[0] < 2);
  SC3A_CHECK (0 <= edge_reverse[2] && edge_reverse[2] < 4);

#ifdef P4_TO_P8
  SC3A_CHECK (my_axis[0] != my_axis[1] && my_axis[1] != my_axis[2]);
  SC3A_CHECK (target_axis[0] != target_axis[1] &&
              target_axis[1] != target_axis[2]);
  SC3A_CHECK (0 <= edge_reverse[1] && edge_reverse[1] < 2);
#else
  SC3A_CHECK (my_axis[1] == 0 && target_axis[1] == 0);
  SC3A_CHECK (edge_reverse[1] == 0);
#endif /* P4_TO_P8 */
#endif

  r->level = q->level;
  r->coords = (uint64_t) 0;

  if (q->level == P4EST3_MORT_MAXLEVEL) {
    /* If (P4EST3_YX_MAXLEVEL == 31) Rmh will overflow. */
    mh = 0;
  }
  else {
    mh = P4EST3_QUADRANT_MORT_LEN (0x01, q->level);
  }
  Rmh = P4EST3_QUADRANT_MORT_LEN (0x01, 0) - mh;

  if (!edge_reverse[0]) {
    extracted = q->coords & (dir_level_mask << my_axis[0]);
  }
  else {
    extracted = (Rmh ^ q->coords) & (dir_level_mask << my_axis[0]);
  }
  r->coords |= ((extracted >> my_axis[0]) << target_axis[0]);

#ifdef P4_TO_P8
  if (!edge_reverse[1]) {
    extracted = q->coords & (dir_level_mask << my_axis[1]);
  }
  else {
    extracted = (Rmh ^ q->coords) & (dir_level_mask << my_axis[1]);
  }
  r->coords |= ((extracted >> my_axis[1]) << target_axis[1]);
#endif

  if (edge_reverse[2] == 1 || edge_reverse[2] == 3) {
    r->coords |= (Rmh & (dir_level_mask << target_axis[2]));
  }

  SC3A_IS (p4est3_quadrant_mort_is_valid, r);
  return NULL;
}

static sc3_error_t *
p4est3_quadrant_mort_copy (const p4est3_quadrant_mort_t * q,
                           p4est3_quadrant_mort_t * copy)
{
  SC3A_IS (p4est3_quadrant_mort_is_valid, q);
  copy->level = q->level;
  copy->coords = q->coords;
  return NULL;
}

static sc3_error_t *
p4est3_quadrant_mort_compare (const p4est3_quadrant_mort_t * q1,
                              const p4est3_quadrant_mort_t * q2, int32_t * j)
{
  SC3A_CHECK (p4est3_quadrant_mort_is_valid (q1, NULL));
  SC3A_CHECK (p4est3_quadrant_mort_is_valid (q2, NULL));

  if (q1->coords == q2->coords) {
    *j = 0;
    return NULL;
  }
  *j = q1->coords > q2->coords ? 1 : -1;
  return NULL;
}

static sc3_error_t *
p4est3_quadrant_mort_ancestor_id (const p4est3_quadrant_mort_t * q,
                                  int32_t level, int32_t * j)
{
  int32_t             mask;
#ifdef P4_TO_P8
  mask = 0x07;
#else
  mask = 0x03;
#endif

  SC3A_IS (p4est3_quadrant_mort_is_valid, q);
  SC3A_CHECK (0 <= level && level <= P4EST3_MORT_MAXLEVEL);
  SC3A_CHECK ((int32_t) q->level >= level);

  *j = 0;
  if (level == 0) {
    return NULL;
  }

  *j = (q->coords & P4EST3_QUADRANT_MORT_LEN (mask, level))
    >> P4EST_DIM * (P4EST3_MORT_MAXLEVEL - level);

  return NULL;
}

static sc3_error_t *
p4est3_quadrant_mort_child_id (const p4est3_quadrant_mort_t * q, int *j)
{
  SC3E (p4est3_quadrant_mort_ancestor_id (q, q->level, j));
  return NULL;
}

static sc3_error_t *
p4est3_quadrant_mort_ancestor (const p4est3_quadrant_mort_t * q,
                               int32_t level, p4est3_quadrant_mort_t * r)
{
  SC3A_IS (p4est3_quadrant_mort_is_valid, q);
  SC3A_CHECK (q->level > level && level >= 0);

  r->coords = q->coords & ~(P4EST3_QUADRANT_MORT_LEN (0x01, level) - 1);
  r->level = (int8_t) level;
  SC3A_IS (p4est3_quadrant_mort_is_valid, r);
  return NULL;
}

static sc3_error_t *
p4est3_quadrant_mort_first_descendant (const p4est3_quadrant_mort_t * q,
                                       int32_t level,
                                       p4est3_quadrant_mort_t * fd)
{
  SC3A_IS (p4est3_quadrant_mort_is_valid, q);
  SC3A_CHECK ((int32_t) q->level <= level && level <= P4EST3_MORT_QMAXLEVEL);

  fd->coords = q->coords;
  fd->level = (int8_t) level;
  SC3A_IS (p4est3_quadrant_mort_is_valid, fd);
  return NULL;
}

static sc3_error_t *
p4est3_quadrant_mort_last_descendant (const p4est3_quadrant_mort_t * q,
                                      int level, p4est3_quadrant_mort_t * ld)
{
  uint64_t            shift;

  SC3A_IS (p4est3_quadrant_mort_is_valid, q);
  SC3A_CHECK ((int) q->level <= level && level <= P4EST3_MORT_QMAXLEVEL);

  shift = P4EST3_QUADRANT_MORT_LEN (0x01, q->level)
    - P4EST3_QUADRANT_MORT_LEN (0x01, level);

  ld->coords = q->coords | shift;
  ld->level = (int8_t) level;
  return NULL;
}

static sc3_error_t *
p4est3_mort_nearest_common_ancestor (const p4est3_quadrant_mort_t * q1,
                                     const p4est3_quadrant_mort_t * q2,
                                     p4est3_quadrant_mort_t * r)
{
  int32_t             maxlevel, shift;
  uint64_t            maxclor;

  SC3A_IS (p4est3_quadrant_mort_is_valid, q1);
  SC3A_IS (p4est3_quadrant_mort_is_valid, q2);

  maxclor = q1->coords ^ q2->coords;
  shift = SC_LOG2_64 (maxclor);
  maxlevel = shift / P4EST_DIM + 1;

  SC3A_CHECK (maxlevel <= P4EST3_MORT_MAXLEVEL);

  r->coords = q1->coords & ~(((uint64_t) 1 << (maxlevel * P4EST_DIM)) - 1);
  r->level = (int8_t) SC_MIN (P4EST3_MORT_MAXLEVEL - maxlevel,
                              (int) SC_MIN (q1->level, q2->level));

  SC3A_IS (p4est3_quadrant_mort_is_valid, r);
  return NULL;
}

static sc3_error_t *
p4est3_quadrant_mort_linear_id (const p4est3_quadrant_mort_t * quadrant,
                                int level, p4est3_gloidx * id)
{
  SC3A_IS (p4est3_quadrant_mort_is_valid, quadrant);
  SC3A_CHECK (0 <= level && level <= P4EST3_MORT_MAXLEVEL);

  *id = quadrant->coords >> (P4EST_DIM * (P4EST3_MORT_MAXLEVEL - level));
  if (level < P4EST3_MORT_QMAXLEVEL) {
    SC3A_CHECK (0 <= *id && *id < ((p4est3_gloidx) 1 << P4EST_DIM * level));
  }
  return NULL;
}

static sc3_error_t *
p4est3_quadrant_mort_morton (int level, p4est3_gloidx id,
                             p4est3_quadrant_mort_t * quadrant)
{
  SC3A_CHECK (0 <= level && level <= P4EST3_MORT_QMAXLEVEL);
  if (level < P4EST3_MORT_QMAXLEVEL) {
    SC3A_CHECK (id < ((p4est3_gloidx) 1 << P4EST_DIM * level));
  }

  quadrant->level = (int8_t) level;
  quadrant->coords = P4EST3_QUADRANT_MORT_LEN (id, level);

  SC3A_IS (p4est3_quadrant_mort_is_valid, quadrant);
  return NULL;
}

static sc3_error_t *
p4est3_quadrant_mort_successor (const p4est3_quadrant_mort_t * q,
                                p4est3_quadrant_mort_t * r)
{
  SC3A_IS (p4est3_quadrant_mort_is_valid, q);
  SC3A_CHECK (q->coords < (1ULL << (P4EST3_MORT_MAXLEVEL * P4EST_DIM)) - 1);

  r->coords = q->coords + P4EST3_QUADRANT_MORT_LEN (0x01, q->level);
  r->level = (int8_t) q->level;

  SC3A_IS (p4est3_quadrant_mort_is_valid, r);
  return NULL;
}

static sc3_error_t *
p4est3_quadrant_mort_predecessor (const p4est3_quadrant_mort_t * q,
                                  p4est3_quadrant_mort_t * r)
{
  SC3A_IS (p4est3_quadrant_mort_is_valid, q);
  SC3A_CHECK (q->coords > (uint64_t) 0);

  r->coords = q->coords - P4EST3_QUADRANT_MORT_LEN (0x01, q->level);
  r->level = (int8_t) q->level;

  SC3A_IS (p4est3_quadrant_mort_is_valid, r);
  return NULL;
}

static sc3_error_t *
p4est3_quadrant_mort_root (p4est3_quadrant_mort_t * r)
{
  r->level = (int8_t) 0;
  r->coords = (uint64_t) 0;
  return NULL;
}

sc3_error_t        *
p4est3_quadrant_mort2d_vtable (p4est3_quadrant_vtable_t * qvt)
{
  SC3A_CHECK (qvt != NULL);
  memset (qvt, 0, sizeof (p4est3_quadrant_vtable_t));

#if (P4EST_DIM == 2 && defined(P4EST_ENABLE_BUILD_2D)) \
 || (P4EST_DIM == 3 && defined(P4EST_ENABLE_BUILD_3D))
  qvt->dim = P4EST_DIM;
  qvt->max_level = P4EST3_MORT_MAXLEVEL;
  qvt->max_children = P4EST_CHILDREN;
  qvt->quadrant_size = sizeof (p4est3_quadrant_mort_t);

  qvt->quadrant_num_uniform = p4est3_quadrant_mort_num_uniform;

  qvt->quadrant_is_valid =
    (p4est3_quadrant_is_t) p4est3_quadrant_mort_is_valid;

  qvt->quadrant_level = (p4est3_quadrant_level_t) p4est3_quadrant_mort_level;

  qvt->quadrant_child_id =
    (p4est3_quadrant_child_id_t) p4est3_quadrant_mort_child_id;

  qvt->quadrant_ancestor_id =
    (p4est3_quadrant_ancestor_id_t) p4est3_quadrant_mort_ancestor_id;

  qvt->quadrant_get_tree_boundary = (p4est3_quadrant_get_tree_boundary_t)
    p4est3_quadrant_mort_get_tree_boundary;

  qvt->quadrant_tree_boundaries =
    (p4est3_quadrant_tree_boundaries_t) p4est3_quadrant_mort_tree_boundaries;

  qvt->quadrant_coordinates =
    (p4est3_quadrant_in_i_out_t) p4est3_quadrant_mort_coords;

  qvt->quadrant_quadrant =
    (p4est3_quadrant_quadrant_t) p4est3_quadrant_mort_quadrant;

  qvt->quadrant_compare =
    (p4est3_quadrant_compare_t) p4est3_quadrant_mort_compare;

  qvt->quadrant_root = (p4est3_quadrant_root_t) p4est3_quadrant_mort_root;

  qvt->quadrant_copy = (p4est3_quadrant_copy_t) p4est3_quadrant_mort_copy;

  qvt->quadrant_parent =
    (p4est3_quadrant_parent_t) p4est3_quadrant_mort_parent;

  qvt->quadrant_face_neighbor =
    (p4est3_quadrant_face_neighbor_t) p4est3_quadrant_mort_face_neighbor;

  qvt->quadrant_tree_face_neighbor = (p4est3_quadrant_tree_face_neighbor_t)
    p4est3_quadrant_mort_tree_face_neighbor;

  qvt->quadrant_face_neighbor =
    (p4est3_quadrant_face_neighbor_t) p4est3_quadrant_mort_face_neighbor;

  qvt->quadrant_predecessor =
    (p4est3_quadrant_predecessor_t) p4est3_quadrant_mort_predecessor;

  qvt->quadrant_successor =
    (p4est3_quadrant_successor_t) p4est3_quadrant_mort_successor;

  qvt->quadrant_child = (p4est3_quadrant_child_t) p4est3_quadrant_mort_child;

  qvt->quadrant_ancestor =
    (p4est3_quadrant_ancestor_t) p4est3_quadrant_mort_ancestor;

  qvt->quadrant_first_descendant = (p4est3_quadrant_first_descendant_t)
    p4est3_quadrant_mort_first_descendant;

  qvt->quadrant_last_descendant =
    (p4est3_quadrant_last_descendant_t) p4est3_quadrant_mort_last_descendant;

  qvt->quadrant_morton =
    (p4est3_quadrant_morton_t) p4est3_quadrant_mort_morton;

  qvt->nearest_common_ancestor =
    (p4est3_nearest_common_ancestor_t) p4est3_mort_nearest_common_ancestor;

  qvt->quadrant_linear_id =
    (p4est3_quadrant_linear_id_t) p4est3_quadrant_mort_linear_id;

  qvt->quadrant_is_ancestor =
    (p4est3_quadrant_is_ancestor_t) p4est3_quadrant_mort_is_ancestor;

  /* verify correctness */
  SC3A_IS (p4est3_quadrant_vtable_is_valid, qvt);
  return NULL;
#else
  return sc3_error_new_kind (SC3_ERROR_RUNTIME, __FILE__, __LINE__, "Creation"
                             " of Morton-based virtual table is denied since"
                             " this dimension is disabled or not supported");
#endif /* !(P4EST_DIM == ? && defined(P4EST_ENABLE_BUILD_?D)) */
}
