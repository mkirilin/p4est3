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
make_allocator (setup_t * t)
{
  SC3A_IS (sc3_allocator_is_setup, t->mainalloc);
  SC3E (sc3_allocator_new (t->mainalloc, &t->alloc));
  SC3E (sc3_allocator_setup (t->alloc));

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
  int qlevel, qlen;
  int8_t *qinfo_array = (int8_t *) vi->user_data;

  SC3E (p4est3_quadrant_level (vi->p3->qvt, vi->quadrant, &qlevel));
  qlen = TEST_QUADRANT_LEN (qlevel);
  SC3E (p4est3_tree_index (vi->p3, vi->ntree, &tree));
  qid = tree->quad_offset + vi->nquad;
  SC3A_CHECK (qid < vi->p3->local_num_quads);
  for (i = 0; i < nfaces; ++i) {
    SC3E_DEMAND
      (qinfo_array[qid + i] == NULL, "volume is visited the second time");
    SC3E (sc3_allocator_calloc
          (vi->p3->alloc, sizeof (int8_t), qlen, &qinfo_array[qid + i]));
  }

  return NULL;
}

static sc3_error_t *
face_callback (p4est3_iterate_face_info_t * fi)
{
  return NULL;
}

static sc3_error_t *
allocate_test_tracking_array (p4est3_t *p3, void **ptr_qinfo_array)
{
  const size_t out_array_size = (1 << p3->qvt->dim) * p3->local_num_quads;
  char * out_array;

  SC3E (sc3_allocator_calloc
        (p3->alloc, sizeof (int8_t *), out_array_size, &out_array));
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
  const p4est3_quadrant_vtable_t *qvt, *qvt_avx, *qvt_mrt;

  SC3X (sc3_MPI_Init (&argc, &argv));
  t->mpicomm = SC3_MPI_COMM_WORLD;
  SC3X (sc3_MPI_Comm_rank (t->mpicomm, &t->mpirank));
  SC3X (set_parameters (t, &qvt, &qvt_avx, &qvt_mrt));

  SC3X (clean_up (t));

  SC3X (sc3_MPI_Finalize ());
  return 0;
}
