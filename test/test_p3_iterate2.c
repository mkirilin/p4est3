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

#include <p4est3_iterate.h>
#include <p4est3_internal.h>
#ifndef P4_TO_P8
#include <p4est3_p4est.h>
#include <p4est3_quadrant_yx.h>
#include <p4est3_quadrant_mort2d.h>

#else
#include <p4est3_p8est.h>
#include <p4est3_quadrant_zyx.h>
#include <p4est3_quadrant_mort3d.h>
#endif

#define MAX_TEST_LEVEL 2
#define TEST_QUADRANT_LEN(l) ((int32_t) 1 << (MAX_TEST_LEVEL - (l)))
#define DIM2 2

//#define DISABLE_TEST_FACES

/* *INDENT-OFF* */
static const int    nface_predef_2d[4][2] =
{{1, 0}, {3, 2}, {3, 2}, {1, 0}};

static const int    nquad_predef_2d[4][2] =
{{0, 1}, {0, 2}, {1, 3}, {2, 3}};

static const int    nface_predef_3d[12][2] =
{{1, 0}, {3, 2},
 {5, 4}, {3, 2},
 {5, 4}, {1, 0},
 {5, 4}, {5, 4},
 {1, 0}, {3, 2},
 {3, 2}, {1, 0}};

static const int    nquad_predef_3d[12][2] =
{{0, 1}, {0, 2},
 {0, 4}, {1, 3},
 {1, 5}, {2, 3},
 {2, 6}, {3, 7},
 {4, 5}, {4, 6},
 {5, 7}, {6, 7}};

static const int    bound_dir_2d[4] =
{0, 1, 1, 0};

static const int    bound_dir_3d[12] =
{0, 1, 2, 1, 2, 0, 2, 2, 0, 1, 1, 0};
/* *INDENT-ON* */

sc3_error_t        *
p4est3_connectivity_new_p4est_twotrees (sc3_allocator_t * alloc,
                                        p4est3_connectivity_t ** pc,
                                        int l_face, int r_face,
                                        int orientation)
{
  p4est_connectivity_t *c4;

  /* verify arguments */
  SC3E_RETVAL (pc, NULL);
  SC3A_IS (sc3_allocator_is_valid, alloc);
#ifdef P4_TO_P8
  SC3A_CHECK (0 <= l_face && l_face < 6);
  SC3A_CHECK (0 <= r_face && r_face < 6);
  SC3A_CHECK (0 <= orientation && orientation < 4);
#else
  SC3A_CHECK (0 <= l_face && l_face < 4);
  SC3A_CHECK (0 <= r_face && r_face < 4);
  SC3A_CHECK (0 <= orientation && orientation < 2);
#endif

  /* create two trees connectivity */
  c4 = p4est_connectivity_new_twotrees (l_face, r_face, orientation);

  /* wrap two trees into p4est3 connectivity */
  SC3E (p4est3_connectivity_new_p4est (alloc, pc, c4, 1));
  return NULL;
}

typedef struct setup
{
  sc3_allocator_t    *alloc;
  sc3_allocator_t    *mainalloc;
  p4est3_connectivity_t *conn;
  sc3_MPI_Comm_t      mpicomm;
  int                 mpirank;
  int                 level;
  p4est3_topidx       num_trees;
  sc3_array_t        *transform;
  sc3_array_t        *nf;
}
setup_t;

typedef struct callback_data
{
  sc3_array_t        *volumes;
  sc3_array_t        *faces;
}
callback_data_t;

static sc3_error_t *
quadrant_2d_mid (const p4est_quadrant_t * quadrant, int level,
                 uint64_t * mid)
{
  int                 i;
  uint64_t            id;
  uint64_t            x, y;

  SC3A_CHECK (0 <= level && level <= MAX_TEST_LEVEL);

  /* this preserves the high bits from negative numbers */
  x = quadrant->x >> (MAX_TEST_LEVEL - level);
  y = quadrant->y >> (MAX_TEST_LEVEL - level);

  id = 0;
  for (i = 0; i < level + 2; ++i) {
    id |= ((x & ((uint64_t) 1 << i)) << ((DIM2 - 1) * i));
    id |= ((y & ((uint64_t) 1 << i)) << ((DIM2 - 1) * i + 1));
  }

  *mid = id;
}

static sc3_error_t *
make_allocator (setup_t * t)
{
  SC3A_IS (sc3_allocator_is_setup, t->mainalloc);
  SC3E (sc3_allocator_new (t->mainalloc, &t->alloc));
  SC3E (sc3_allocator_setup (t->alloc));

  return NULL;
}

static sc3_error_t *
array_new (sc3_allocator_t * alloc, size_t esize, int ealloc,
           int ecount, sc3_array_t ** arr)
{
  SC3E_RETVAL (arr, NULL);
  SC3A_IS (sc3_allocator_is_setup, alloc);
  SC3A_CHECK (ealloc >= 0);

  SC3E (sc3_array_new (alloc, arr));
  SC3E (sc3_array_set_elem_size (*arr, esize));
  SC3E (sc3_array_set_elem_alloc (*arr, ealloc));
  SC3E (sc3_array_set_elem_count (*arr, ecount));
  SC3E (sc3_array_set_initzero (*arr, 1));
  SC3E (sc3_array_setup (*arr));

  return NULL;
}

static sc3_error_t *
make_connectivity (setup_t * t, int dim, int l_face, int r_face, int ori)
{
/* SC3E (p4est3_connectivity_new (t->alloc, &t->conn));
  SC3E (p4est3_connectivity_set_dim (t->conn, dim));
  SC3E (p4est3_connectivity_set_num_trees (t->conn, t->num_trees));
  SC3E (p4est3_connectivity_setup (t->conn)); */
  SC3E (p4est3_connectivity_new_p4est_twotrees
        (t->alloc, &t->conn, l_face, r_face, ori));

  return NULL;
}

static sc3_error_t *
volume_callback (p4est3_iterate_volume_info_t * vi)
{
  int i, nfaces = 1 << vi->p3->qvt->dim;
  p4est3_tree_t *tree;
  p4est3_locidx qid;
  int qlevel, face_area;
  int8_t **qinfo_array = (int8_t **) vi->user_data;

  SC3E (p4est3_quadrant_level (vi->p3->qvt, vi->quadrant, &qlevel));
  face_area = sc3_intpow (TEST_QUADRANT_LEN (qlevel), vi->p3->qvt->dim - 1);
  SC3E (p4est3_tree_index (vi->p3, vi->ntree, &tree));
  qid = tree->quad_offset + vi->nquad;
  SC3A_CHECK (qid < vi->p3->local_num_quads);
  for (i = 0; i < nfaces; ++i) {
    SC3E_DEMAND
      (qinfo_array[qid * nfaces + i] == NULL,
       "volume is visited the second time");
    SC3E (sc3_allocator_calloc
          (vi->p3->alloc, sizeof (int8_t), face_area,
           &qinfo_array[qid * nfaces + i]));
  }

  return NULL;
}

static int
check_q_in_proc (const p4est3_t * p3, const p4est3_topidx t,
                 const p4est3_locidx nquad)
{
  if (t < p3->fltree || p3->lltree < t) {
    return 0;
  }
  if (p3->fltree < t && t < p3->lltree) {
    return 1;
  }
  if (t == p3->fltree) {
    return (nquad + p3->gtroffset[p3->fltree] >= p3->goffset[p3->mpirank]);
  }
  if (t == p3->lltree) {
    return (nquad + p3->gtroffset[p3->lltree] < p3->goffset[p3->mpirank]);
  }
  return 0;
}

static sc3_error_t *
convert_quad_to_mid (const p4est3_t * const p3, const void * q)
{
  /* convert a quadrant to d-1 morton index */
  return NULL;
}

static sc3_error_t *
fill_in_side_info_compl (const p4est3_t * p3, int8_t ** const qinfo_array,
                         const p4est3_iterate_face_side_t * const side,
                         const int nsides)
{
  int8_t *finfo_array;
  p4est3_locidx qid_loc;
  int i, qlevel, face_area, nfaces = 1 << p3->qvt->dim;


  qid_loc = (p4est3_gloidx) side->nquad
            + p3->gtroffset[side->ntree] - p3->goffset[p3->mpirank];
  SC3E (p4est3_quadrant_level (p3->qvt, side->quadrant, &qlevel));
  face_area = sc3_intpow (TEST_QUADRANT_LEN (qlevel), p3->qvt->dim - 1);
  finfo_array = qinfo_array[qid_loc * nfaces + side->nface];
  SC3A_CHECK (finfo_array != NULL);
  for (i = 0; i < face_area; ++i) {
    SC3E_DEMAND (finfo_array[i] == 0,
                  "trying to write into non-empty face-info cell");
    finfo_array[i] = nsides;
  }
  return NULL;
}

static sc3_error_t *
face_callback (p4est3_iterate_face_info_t * fi)
{
  char *tempq[2];
  int i, qlevel, face_area, nfaces = 1 << fi->p3->qvt->dim;
  p4est3_iterate_face_side_t *sides[2], *side_small, *side_big;
  sc3_array_t *ftransform;
  int8_t **qinfo_array = (int8_t **) fi->user_data;
  int8_t *finfo_array;
  int ntree, nsides, levels[2], ss_id /* smaller side index */;
  SC3E (sc3_array_get_elem_count (fi->sides, &nsides));
  SC3E_DEMAND ((nsides == 2) || ((nsides == 1) && fi->tree_boundary),
               "one face's side not on a tree's boundary");

  /* test adjacency */
  for (i = 0; i < nsides; ++i) {
    SC3E (sc3_array_index (fi->sides, 0, &sides[i]));
    SC3E (p4est3_quadrant_level
          (fi->p3->qvt, sides[i]->quadrant, &levels[i]));
    SC3E (sc3_allocator_calloc_one (fi->p3->alloc, fi->p3->qsize, &tempq[i]));
  }
  if (nsides == 2) {
    ss_id = levels[0] > levels[1] ? 0 : 1;
    side_small = sides[ss_id];
    side_big = sides[1 - ss_id];
    SC3E (p4est3_quadrant_ancestor
          (fi->p3->qvt, side_small->quadrant, levels[1 - ss_id], &tempq[0]));
    if (side_small->ntree == side_big->ntree) {
      SC3E (p4est3_quadrant_face_neighbor
            (fi->p3->qvt, tempq[0], side_small->nface, &tempq[1]));
    }
    else {
      ntree = side_small->ntree;
      SC3E (array_new (fi->p3->alloc, sizeof (int), 9, 9, &ftransform));
      SC3E (p4est3_connectivity_get_face_transform
            (fi->p3->conn, side_small->nface, &ntree, ftransform));
      SC3E_DEMAND (ntree == side_big->ntree,
                   "face transform trees mismatch");
      SC3E (p4est3_quadrant_tree_face_neighbor
            (fi->p3->qvt, tempq[0], ftransform,
             side_small->nface, &tempq[1]));
    }
    SC3E_DEMIS3 (
      p4est3_quadrant_is3_equal, fi->p3->qvt, tempq[1], side_big->quadrant);
  } else { /* nsides == 1, nothing else is possible */
    side_big = side_small = sides[0];
  }
  /* fill testing arrays */
  /* for the smaller and/or one-sided quad we fill
     the whole face-related (part of) array */
  if (check_q_in_proc (fi->p3, side_small->ntree, side_small->nquad)) {
    SC3E (fill_in_side_info_compl (fi->p3, qinfo_array, side_small, nsides));
  }
  /* if quadrants are of equal size, we full the bigger one completely, too */
  if (nsides == 2 && levels[0] == levels[1]) {
    if (check_q_in_proc (fi->p3, side_big->ntree, side_big->nquad)) {
      SC3E (fill_in_side_info_compl (fi->p3, qinfo_array, side_big, nsides));
    }
  }
  return NULL;
}

static sc3_error_t *
allocate_test_tracking_array (p4est3_t *p3, void **ptr_qinfo_array)
{
  const size_t out_array_size = (1 << p3->qvt->dim) * p3->local_num_quads;
  char * out_array;

  SC3E (sc3_allocator_calloc
        (p3->alloc, sizeof (int8_t *), out_array_size, &out_array));
  return NULL;
}

static sc3_error_t *
make_new_p4est3 (p4est3_t ** p3, setup_t * t, const p4est3_quadrant_vtable_t * qvt)
{
  SC3A_IS (sc3_allocator_is_setup, t->alloc);

  /* create p4est object with connectivity */
  SC3E (p4est3_new (t->alloc, p3));
  SC3E (p4est3_set_comm (*p3, t->mpicomm, 1));
  SC3E (p4est3_set_connectivity (*p3, t->conn));
  SC3E (p4est3_set_quadrant_vtable (*p3, qvt));
  SC3E (p4est3_set_level (*p3, t->level));
  SC3E (p4est3_set_shared (*p3, 1));
  SC3E (p4est3_setup (*p3));

  return NULL;
}

static sc3_error_t *
free_allocator (sc3_allocator_t ** alloc)
{
  SC3A_IS (sc3_allocator_is_setup, *alloc);
  SC3E (sc3_allocator_destroy (alloc));
  return NULL;
}

static sc3_error_t *
set_parameters (setup_t * t, const p4est3_quadrant_vtable_t ** qvt)
{
  t->mainalloc = sc3_allocator_nothread ();
  SC3E (make_allocator (t));
  SC3E (p4est3_quadrant_vtable_p4est (qvt));
  SC3E (array_new (t->alloc, sizeof (int), 9, 9, &t->transform));
  SC3E (array_new (t->alloc, sizeof (int), (*qvt)->dim, (*qvt)->dim, &t->nf));

  t->level = 1;
  t->num_trees = 2;

  return NULL;
}

static sc3_error_t *
clean_up (setup_t * t)
{
  SC3E (sc3_array_destroy (&t->transform));
  SC3E (sc3_array_destroy (&t->nf));
  SC3E (free_allocator (&t->alloc));
  return NULL;
}

int
main (int argc, char **argv)
{
  setup_t             st, *t = &st;
  const p4est3_quadrant_vtable_t *qvt;

  SC3X (sc3_MPI_Init (&argc, &argv));
  t->mpicomm = SC3_MPI_COMM_WORLD;
  SC3X (sc3_MPI_Comm_rank (t->mpicomm, &t->mpirank));
  SC3X (set_parameters (t, &qvt));

  SC3X (clean_up (t));

  SC3X (sc3_MPI_Finalize ());
  return 0;
}
