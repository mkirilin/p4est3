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

#include <p4est3_connectivity.h>

/* *INDENT-OFF* */
static const int    face_corners_2d[4][2] =
{{ 0, 2 },
 { 1, 3 },
 { 0, 1 },
 { 2, 3 }};
static const int    face_corners_3d[6][4] =
{{ 0, 2, 4, 6 },
 { 1, 3, 5, 7 },
 { 0, 1, 4, 5 },
 { 2, 3, 6, 7 },
 { 0, 1, 2, 3 },
 { 4, 5, 6, 7 }};
/* *INDENT-ON* */

static sc3_error_t *
make_connectivity (sc3_allocator_t * alloc,
                   p4est3_connectivity_t ** conn, int dim)
{
  SC3E (p4est3_connectivity_new (alloc, conn));
  SC3E (p4est3_connectivity_set_dim (*conn, dim));
  SC3E (p4est3_connectivity_set_num_trees (*conn, 1));
  SC3E (p4est3_connectivity_setup (*conn));

  return NULL;
}

/** Checks for each face corner if the corner indices match on both
 * sides.
 * Let face corner fci, corresponding to corner index ci, be
 * adjacent to face corner fcj, corresponding to corner index cj. We
 * test if fci is seen from fcj and vice versa.
 */
static sc3_error_t *
test_face_neighbor_face_corner (sc3_allocator_t * alloc,
                                p4est3_connectivity_t * conn, int dim)
{
  const int           numFaces = 2 * dim;
  const int           nFaceCorners = 1 << (dim - 1);
  int                 l_face;   /* left face index */
  int                 r_face;   /* right face index */
  int                 ori;      /* the orientation that has been set */
  int                 c0, c1, lowerFaceIdx, higherFaceIdx;

  SC3A_IS (p4est3_connectivity_is_setup, conn);

  for (l_face = 0; l_face < numFaces; ++l_face) {       /* set l_face */
    for (r_face = 0; r_face < numFaces; ++r_face) {     /* set r_face */
      for (ori = 0; ori < nFaceCorners; ++ori) {        /* set orientation */
        /* swap face indices if necessary */
        if (l_face <= r_face) {
          lowerFaceIdx = l_face;
          higherFaceIdx = r_face;
        }
        else {
          lowerFaceIdx = r_face;
          higherFaceIdx = l_face;
        }
        /* verify bijectivity of transformation */
        for (c0 = 0; c0 < nFaceCorners; ++c0) {
          c1 = c0;
          SC3E (p4est3_connectivity_face_neighbor_face_corner
                (conn, &c1, lowerFaceIdx, higherFaceIdx, ori));
          SC3E (p4est3_connectivity_face_neighbor_face_corner
                (conn, &c1, higherFaceIdx, lowerFaceIdx, ori));
          SC3E_DEMAND (c0 == c1, "Face <-> neighbor face corner "
                       "transformation is not bijective");
        }
      }
    }
  }
  return NULL;
}

static sc3_error_t *
test_face_child (sc3_allocator_t * alloc,
                 p4est3_connectivity_t * conn, int dim)
{
  const int           numFaces = 2 * dim;
  const int           nFaceCorners = 1 << (dim - 1);
  int                 face;
  int                 c0, c1;

  SC3A_IS (p4est3_connectivity_is_setup, conn);
  for (face = 0; face < numFaces; ++face) {
    for (c0 = 0; c0 < nFaceCorners; ++c0) {
      c1 = c0;
      SC3E (p4est3_connectivity_get_face_child_id (conn, face, &c1));
      if (dim == 2) {
        SC3E_DEMAND (c1 == face_corners_2d[face][c0],
                     "Connectivity provides a wrong child id for a corner");
      }
      else if (dim == 3) {
        SC3E_DEMAND (c1 == face_corners_3d[face][c0],
                     "Connectivity provides a wrong child id for a corner");
      }
      else {
        SC3E_UNREACH ("Wrong dimension");
      }
    }
  }
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
  sc3_allocator_t    *alloc;
  p4est3_connectivity_t *conn;

  /* make allocator */
  SC3X (sc3_allocator_new (sc3_allocator_nothread (), &alloc));
  SC3X (sc3_allocator_setup (alloc));

  for (int dim = 2; dim <= 3; ++dim) {
    SC3X (make_connectivity (alloc, &conn, dim));
    SC3X (test_face_neighbor_face_corner (alloc, conn, dim));
    SC3X (test_face_child (alloc, conn, dim));
    SC3X (p4est3_connectivity_destroy (&conn));
  }
  SC3X (free_allocator (&alloc));

  return 0;
}
