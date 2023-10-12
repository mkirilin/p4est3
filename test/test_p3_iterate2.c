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

#ifdef P4EST_ENABLE_DEBUG
#include <p4est_extended.h>
#include <p4est_vtk.h>
#endif /* P4EST_ENABLE_DEBUG */

#include <p4est3_p4est.h>
#include <p4est3_quadrant_yx.h>
#include <p4est3_quadrant_mort2d.h>

#else
#ifdef P4EST_ENABLE_DEBUG
#include <p4est_to_p8est.h>
#include <p8est_extended.h>
#include <p8est_vtk.h>
#endif /* P4EST_ENABLE_DEBUG */

#include <p4est3_p8est.h>
#include <p4est3_quadrant_zyx.h>
#include <p4est3_quadrant_mort3d.h>
#endif

#define MAX_TEST_LEVEL 10
#define TEST_QUADRANT_LEN(l) ((int32_t) 1 << (MAX_TEST_LEVEL - (l)))

static int          refine_level = 5;

#if 0
#ifdef P4_TO_P8
static double       refinement_fraction = 1. / 7.;
#else
static double       refinement_fraction = 1. / 3.;
#endif
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

static int
refine_fraction (p4est_t * p, p4est_topidx_t which_tree, p4est_quadrant_t * q)
{
  /* The formula in the line below implies
   * quadrant_fraction = 7 * refinement_fraction + 1.
   * In particular we need for doubling the quadrants for a level increment
   * refinment_factor = 1 / 7.
   */
  const p4est_locidx_t quad_count_refinement_threshold =
    (p4est_locidx_t) (refinement_fraction * p->global_num_quadrants);
  p4est_locidx_t     *quadrant_local_id = (p4est_locidx_t *) p->user_pointer;

  return
    ((p->global_first_quadrant[p->mpirank] + (*quadrant_local_id)++) <=
     quad_count_refinement_threshold);
}

static sc3_error_t *
refine_p3_fraction (p4est3_refine_callback_info_t * ri, int *is_refine)
{
  /* The formula in the line below implies
   * quadrant_fraction = 7 * refinement_fraction + 1.
   * In particular we need for doubling the quadrants for a level increment
   * refinment_factor = 1 / 7.
   */
  /* p4est->global_num_quadrants is the old number of quadrants */
  SC3A_CHECK (is_refine != NULL);
  const p4est3_locidx quad_count_refinement_threshold =
    (p4est3_locidx) (refinement_fraction * ri->p3->global_num_quads);
  p4est3_locidx      *quadrant_local_id = (p4est3_locidx *) ri->user_data;
  *is_refine =
    ((ri->p3->goffset[ri->p3->mpirank] + (*quadrant_local_id)++) <=
     quad_count_refinement_threshold);

  return NULL;
}
#endif

#ifdef P4EST_ENABLE_DEBUG
static int
refine_fractal (p4est_t * p, p4est_topidx_t which_tree, p4est_quadrant_t * q)
{
  /* Refine every 7th (3d) or 3rd (2d) global quadrant. */
  p4est_locidx_t     *quadrant_local_id = (p4est_locidx_t *) p->user_pointer;

  return (((p->global_first_quadrant[p->mpirank] + (*quadrant_local_id)++) %
#ifdef P4_TO_P8
           7)
#else
           3)
#endif
          == 0);
}
#endif /* P4EST_ENABLE_DEBUG */

static sc3_error_t *
refine_p3_fractal (p4est3_refine_callback_info_t * ri, int *is_refine)
{
  /* Refine every 7th (3d) or 3rd (2d) global quadrant. */
  p4est3_locidx      *quadrant_local_id = (p4est3_locidx *) ri->user_data;
  SC3A_CHECK (is_refine != NULL);

  *is_refine = (((ri->p3->goffset[ri->p3->mpirank] + (*quadrant_local_id)++) %
#ifdef P4_TO_P8
                 7)
#else
                 3)
#endif
                == 0);
  return NULL;
}

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
  int                 level;
}
setup_t;

static sc3_error_t *
quadrant_to_mid (const p4est_qcoord_t coords[P4EST_DIM - 1], int level,
                 uint64_t * const mid)
{
  int                 i;
  uint64_t            id;
  uint64_t            x;
#ifdef P4_TO_P8
  uint64_t            y;
#endif

  SC3A_CHECK (0 <= level && level <= MAX_TEST_LEVEL);

  /* this preserves the high bits from negative numbers */
  x = coords[0] >> (MAX_TEST_LEVEL - level);
#ifdef P4_TO_P8
  y = coords[1] >> (MAX_TEST_LEVEL - level);
#endif

  id = 0;
  for (i = 0; i <= level; ++i) {
    id |= ((x & ((uint64_t) 1 << i)) << (((P4EST_DIM - 1) - 1) * i));
#ifdef P4_TO_P8
    id |= ((y & ((uint64_t) 1 << i)) << (((P4EST_DIM - 1) - 1) * i + 1));
#endif
  }

  *mid = id * (1 << (MAX_TEST_LEVEL - level));
  return NULL;
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
volume_callback (p4est3_iterate_volume_info_t * vi)
{
  /* for every face of every volume we allocate
     an array by the length of face's area */

  int                 i;
  p4est3_tree_t      *tree;
  p4est3_locidx       qid;
  int                 qlevel, face_area;
  int8_t            **qinfo_array = (int8_t **) vi->user_data;

  SC3E (p4est3_quadrant_level (vi->p3->qvt, vi->quadrant, &qlevel));
  face_area = sc3_intpow (TEST_QUADRANT_LEN (qlevel), vi->p3->qvt->dim - 1);
  SC3E (p4est3_tree_index (vi->p3, vi->ntree, &tree));
  SC3A_CHECK (vi->ntree == tree->treeid);
  qid = (p4est3_locidx)
    ((p4est3_gloidx) vi->nquad + vi->p3->gtroffset[tree->treeid])
    - vi->p3->goffset[vi->p3->mpirank];
  SC3A_CHECK (qid < vi->p3->local_num_quads);
  for (i = 0; i < P4EST_FACES; ++i) {
    SC3E_DEMAND
      (qinfo_array[qid * P4EST_FACES + i] == NULL,
       "volume is visited the second time");
    SC3E (sc3_allocator_calloc
          (vi->p3->alloc, sizeof (int8_t), face_area,
           &qinfo_array[qid * P4EST_FACES + i]));
  }

  return NULL;
}

static int
check_q_in_proc (const p4est3_t * p3, const p4est3_topidx t,
                 const p4est3_locidx nquad)
{
  p4est3_gloidx       nquad_glo = (p4est3_gloidx) nquad + p3->gtroffset[t];
  return (p3->goffset[p3->mpirank] <= nquad_glo
          && nquad_glo < p3->goffset[p3->mpirank + 1]);
}

static sc3_error_t *
convert_quad_to_mid (const p4est3_t * const p3,
                     const p4est3_iterate_face_side_t * const side,
                     const p4est3_iterate_face_side_t * const patch,
                     uint64_t * const mid)
{
  /* convert a quadrant to d-1 morton index */
  const int           axis = patch->nface / P4EST_DIM;
  p4est_qcoord_t      patch_crdDIM[P4EST_DIM], side_crdDIM[P4EST_DIM],
    coords[P4EST_DIM - 1];
  int                 i, c, patch_lvl;
#ifdef P4EST_ENABLE_DEBUG
  /*int32_t patch_len, side_len; */
  int                 side_lvl;
#endif

  SC3E (p4est3_quadrant_coordinates (p3->qvt, patch->quadrant, patch_crdDIM));
  SC3E (p4est3_quadrant_coordinates (p3->qvt, side->quadrant, side_crdDIM));

  SC3E (p4est3_quadrant_level (p3->qvt, patch->quadrant, &patch_lvl));
#ifdef P4EST_ENABLE_DEBUG
  SC3E (p4est3_quadrant_level (p3->qvt, side->quadrant, &side_lvl));
  SC3A_CHECK (side_lvl <= patch_lvl);

  /* is not correct for inter tree border */
  /*if (patch_crdDIM[axis] < side_crdDIM[axis]) {
     patch_len = (int32_t) 1 << (p3->qvt->max_level - patch_lvl);
     SC3A_CHECK (patch_crdDIM[axis] + patch_len == side_crdDIM[axis]);
     }
     else if (patch_crdDIM[axis] > side_crdDIM[axis]) {
     side_len = (int32_t) 1 << (p3->qvt->max_level - side_lvl);
     SC3A_CHECK (side_crdDIM[axis] + side_len == patch_crdDIM[axis]);
     } */
#endif

  memset (coords, 0, sizeof (p4est_qcoord_t) * (P4EST_DIM - 1));
  for (i = 0; i < p3->qvt->dim - 1; ++i) {
    coords[i] = 0;
  }
  /* DIM x d coords -> DIM-1 x d coords */
  for (i = 0, c = 0; i < p3->qvt->dim; ++i) {
    if (i == axis) {
      continue;
    }
    /* the ckech is invalid for inter-tree iteration */
    /*SC3A_CHECK ((side_lvl == patch_lvl && patch_crdDIM[i] == side_crdDIM[i])
       || (side_lvl != patch_lvl && patch_crdDIM[i] >= side_crdDIM[i])); */
    coords[c++] =
      (patch_crdDIM[i] - side_crdDIM[i]) >> (P4EST3_REF_MAXLEVEL -
                                             MAX_TEST_LEVEL);
  }
  SC3E (quadrant_to_mid (coords, patch_lvl, mid));
  return NULL;
}

static sc3_error_t *
fill_in_side_info (const p4est3_t * p3, int8_t ** const qinfo_array,
                   const p4est3_iterate_face_side_t * const side,
                   const p4est3_iterate_face_side_t * const patch,
                   const int nsides)
{
  int8_t             *finfo_array;
  p4est3_locidx       side_qid_loc;
  int                 patch_qlevel, patch_area;
  uint64_t            i, begin_mid = 0;

#ifdef P4EST_ENABLE_DEBUG
  if (side != patch) {
    SC3A_CHECK (nsides == 2);
  }
#endif

  side_qid_loc = (p4est3_gloidx) side->nquad
    + p3->gtroffset[side->ntree] - p3->goffset[p3->mpirank];

  SC3E (p4est3_quadrant_level (p3->qvt, patch->quadrant, &patch_qlevel));
  patch_area =
    sc3_intpow (TEST_QUADRANT_LEN (patch_qlevel), p3->qvt->dim - 1);
  finfo_array = qinfo_array[side_qid_loc * P4EST_FACES + side->nface];
  SC3A_CHECK (finfo_array != NULL);

  SC3E (convert_quad_to_mid (p3, side, patch, &begin_mid));

  for (i = begin_mid; i < begin_mid + patch_area; ++i) {
    SC3E_DEMAND (finfo_array[i] == 0,
                 "trying to write into non-empty face-info cell");
    finfo_array[i] = nsides;
  }
  return NULL;
}

static sc3_error_t *
face_callback (p4est3_iterate_face_info_t * fi)
{
  /* Test adjacency (necessity) and fill in
     `test tracking array` to check later (sufficiency) */

  char               *tempq[2];
  p4est3_iterate_face_side_t *sides[2], *side_small, *side_big;
  sc3_array_t        *ftransform;
  int8_t            **qinfo_array = (int8_t **) fi->user_data;
  p4est3_topidx       ntree;
  int                 i, nsides, levels[2], ss_id /* smaller side index */ ;
  SC3E (sc3_array_get_elem_count (fi->sides, &nsides));
  SC3E_DEMAND ((nsides == 2) || ((nsides == 1) && fi->tree_boundary),
               "one face's side not on a tree's boundary");

  /* test adjacency */
  for (i = 0; i < nsides; ++i) {
    SC3E (sc3_array_index (fi->sides, i, &sides[i]));
    SC3E (p4est3_quadrant_level
          (fi->p3->qvt, sides[i]->quadrant, &levels[i]));
    SC3E (sc3_allocator_calloc_one (fi->p3->alloc, fi->p3->qsize, &tempq[i]));
  }
  if (nsides == 2) {
    ss_id = levels[0] > levels[1] ? 0 : 1;
    side_small = sides[ss_id];
    side_big = sides[1 - ss_id];
    if (levels[0] == levels[1]) {
      SC3E (p4est3_quadrant_copy
            (fi->p3->qvt, side_small->quadrant, tempq[0]));
    }
    else {
      SC3E (p4est3_quadrant_ancestor
            (fi->p3->qvt, side_small->quadrant, levels[1 - ss_id], tempq[0]));
    }
    if (side_small->ntree == side_big->ntree) {
      SC3E (p4est3_quadrant_face_neighbor
            (fi->p3->qvt, tempq[0], side_small->nface, tempq[1]));
    }
    else {
      ntree = side_small->ntree;
      SC3E (array_new (fi->p3->alloc, sizeof (int), 9, 9, &ftransform));
      SC3E (p4est3_connectivity_get_face_transform
            (fi->p3->conn, side_small->nface, &ntree, ftransform));
      SC3E_DEMAND (ntree == side_big->ntree, "face transform trees mismatch");
      SC3E (p4est3_quadrant_tree_face_neighbor
            (fi->p3->qvt, tempq[0], ftransform, side_small->nface, tempq[1]));
      SC3E (sc3_array_destroy (&ftransform));
    }
    SC3E_DEMIS3 (p4est3_quadrant_is3_equal, fi->p3->qvt, tempq[1],
                 side_big->quadrant);
  }
  else {                        /* nsides == 1, nothing else is possible */
    side_big = side_small = sides[0];
  }
  /* fill testing arrays */
  /* for the smaller and/or one-sided quad we fill
     the whole face-related (part of) array */
  if (check_q_in_proc (fi->p3, side_small->ntree, side_small->nquad)) {
    SC3E (fill_in_side_info
          (fi->p3, qinfo_array, side_small, side_small, nsides));
  }
  if (nsides == 2 && levels[0] == levels[1]) {
    /* if quadrants are of equal size, we full the bigger one completely, too */
    if (check_q_in_proc (fi->p3, side_big->ntree, side_big->nquad)) {
      SC3E (fill_in_side_info
            (fi->p3, qinfo_array, side_big, side_big, nsides));
    }
  }
  else if (nsides == 2) {
    /* if quads have different size, we fill the bigger one only partially */
    if (check_q_in_proc (fi->p3, side_big->ntree, side_big->nquad)) {
      SC3E (fill_in_side_info
            (fi->p3, qinfo_array, side_big, side_small, nsides));
    }
  }
  for (i = 0; i < nsides; ++i) {
    SC3E (sc3_allocator_free (fi->p3->alloc, tempq[i]));
  }
  return NULL;
}

static sc3_error_t *
test_tracking_array (p4est3_t * p3, int8_t ** const qinfo_array)
{
  int                 qlevel, face_area, is_boundary;
  int                 f, p /* patch number */ , nface, orient;
  char               *q;
  int8_t             *finfo_array;
  p4est3_topidx       t, which_tree;
  p4est3_locidx       qtid /* local number of quadrant within on a tree */ ;
  p4est3_tree_t      *tree;

  for (t = p3->fltree; t <= p3->lltree; ++t) {
    SC3E (p4est3_tree_index (p3, t, &tree));
    for (qtid = 0; qtid < tree->num_quads; ++qtid) {
      q = tree->tquads + qtid * p3->qvt->quadrant_size;
      SC3E (p4est3_quadrant_level (p3->qvt, q, &qlevel));
      face_area = sc3_intpow (TEST_QUADRANT_LEN (qlevel), p3->qvt->dim - 1);

      for (f = 0; f < P4EST_FACES; ++f) {
        nface = f;
        finfo_array =
          qinfo_array[(qtid + tree->quad_offset) * P4EST_FACES + f];
        SC3E_DEMAND (finfo_array != NULL, "face info array is not allocated");
        SC3E_DEMAND (finfo_array[0] == 1 || finfo_array[0] == 2,
                     "face info array has an illigal value");
        which_tree = t;
        SC3E (p4est3_connectivity_get_face
              (p3->conn, &which_tree, &nface, &orient));
        SC3E (p4est3_quadrant_get_tree_boundary
              (p3->qvt, q, f, &is_boundary));

        for (p = 0; p < face_area; ++p) {
          SC3E_DEMAND (finfo_array[p] == finfo_array[0],
                       "face info array values differ");
          if (finfo_array[p] == 1) {
            /* iff it's a physical boundary */
            SC3E_DEMAND (which_tree == t,
                         "Neighbor tree exists => not a physical boundary");
            SC3E_DEMAND (is_boundary,
                         "Not a tree's boundary => not a physical boundary");
          }
          else if (finfo_array[p] == 2) {
            /* iff it's not a physical boundary */
            SC3E_DEMAND (!is_boundary || (is_boundary && which_tree != t),
                         "A physical boundary");
          }
        }
      }
    }
  }
  return NULL;
}

static sc3_error_t *
allocate_test_tracking_array (p4est3_t * p3, int8_t *** ptr_qinfo_array)
{
  const size_t        out_array_size = P4EST_FACES * p3->local_num_quads;
  int8_t            **out_array;

  SC3E (sc3_allocator_calloc
        (p3->alloc, sizeof (int8_t *), out_array_size, &out_array));
  *ptr_qinfo_array = out_array;
  return NULL;
}

static sc3_error_t *
p4est3_new_shortcut (p4est3_t ** p3, const setup_t * t, int start_level,
                     const p4est3_quadrant_vtable_t * qvt,
                     p4est3_t * src, p4est3_refine_callback_t p3crefine,
                     int is_partition, void *user_data)
{
  SC3E (p4est3_new (t->alloc, p3));
  SC3E (p4est3_set_comm (*p3, SC3_MPI_COMM_WORLD, 1));
  SC3E (p4est3_set_connectivity (*p3, t->conn));
  SC3E (p4est3_set_quadrant_vtable (*p3, qvt));
  SC3E (p4est3_set_level (*p3, start_level));
  SC3E (p4est3_set_setup_mode (*p3, P4EST3_NEW_RECURSIVE));
  SC3E (p4est3_set_refine (*p3, p3crefine));
  SC3E (p4est3_set_source (*p3, src));
  SC3E (p4est3_set_shared (*p3, 1));
  SC3E (p4est3_set_contiguous (*p3, 1));
  SC3E (p4est3_set_family (*p3, 0));
  SC3E (p4est3_set_partition (*p3, is_partition, NULL));
  /*SC3E (p4est3_set_user_data (*p3, user_data)); */
  if ((*p3)->old != NULL) {
    (*p3)->old->user_data = user_data;
  }

  return NULL;
}

static sc3_error_t *
make_forest_for_test (p4est3_t ** p3, setup_t * t,
                      p4est3_quadrant_vtable_t * qvt)
{
  int                 i;
  p4est3_locidx       quadrant_local_id = 0;
  p4est3_t           *p3refined;
  p4est3_refine_callback_t p3crefine = refine_p3_fractal;
#ifdef P4EST_ENABLE_DEBUG
  p4est_t            *p;
  p4est_connectivity_t *conn_old;
  p4est_refine_t      crefine = refine_fractal;
#endif

  SC3E (p4est3_new_shortcut (p3, t, t->level, qvt, NULL, NULL, 0, NULL));
  SC3E (p4est3_setup (*p3));
#ifdef P4EST_ENABLE_DEBUG
  conn_old = p4est_connectivity_new_twotrees (1, 0, 0);
  p = p4est_new_ext
    (sc_MPI_COMM_WORLD, conn_old, 0, t->level,
     1, 0, NULL, &quadrant_local_id);
#endif

  /*** refine in a loop ***/
  for (i = 0; i < refine_level; ++i, quadrant_local_id = 0) {
    SC3E (p4est3_new_shortcut
          (&p3refined, t, 0, qvt, *p3, p3crefine, 0, &quadrant_local_id));
    SC3E (p4est3_setup (p3refined));
#ifdef P4EST_ENABLE_DEBUG
    quadrant_local_id = 0;
    p4est_refine (p, 0, crefine, NULL);
#endif

    SC3E (p4est3_destroy (p3));
    SC3E (p4est3_new_shortcut (p3, t, 0, qvt, p3refined, NULL, 1, NULL));
    SC3E (p4est3_setup (*p3));
#ifdef P4EST_ENABLE_DEBUG
    p4est_partition (p, 0, NULL);
    if (i == refine_level - 1) {
      p4est_vtk_write_file (p, NULL, "p3_iter_test");
    }
#endif

    SC3E (p4est3_destroy (&p3refined));
  }

#ifdef P4EST_ENABLE_DEBUG
  p4est_destroy (p);
  p4est_connectivity_destroy (conn_old);
#endif

  return NULL;
}

static sc3_error_t *
set_parameters (setup_t * t, const p4est3_quadrant_vtable_t ** qvt)
{
  t->mainalloc = sc3_allocator_nothread ();
  SC3E (make_allocator (t));
  SC3E (p4est3_quadrant_vtable_p4est (qvt));
  SC3E (p4est3_connectivity_new_p4est_twotrees (t->alloc, &t->conn, 1, 0, 0));

  t->level = 2;

  return NULL;
}

static sc3_error_t *
clean_up (setup_t * t)
{
  SC3E (p4est3_connectivity_destroy (&t->conn));
  SC3E (sc3_allocator_destroy (&t->alloc));

  return NULL;
}

static sc3_error_t *
perform_test (setup_t * t, p4est3_quadrant_vtable_t * qvt)
{
  size_t              i, out_array_size;
  int8_t            **qinfo_array = NULL;
  p4est3_t           *p3;

  /* preparations */
  SC3E (make_forest_for_test (&p3, t, qvt));
  SC3E (allocate_test_tracking_array (p3, &qinfo_array));

  /* face iterator test run */
  SC3E (p4est3_iterate_face
        (p3, volume_callback, face_callback, qinfo_array));

  /* test arrays filled during iteration */
  SC3E (test_tracking_array (p3, qinfo_array));

  /* cleaning up */
  out_array_size = P4EST_FACES * p3->local_num_quads;
  for (i = 0; i < out_array_size; ++i) {
    SC3E_DEMAND (qinfo_array[i] != NULL,
                 "Cleaning TRA error: array's cell is empty");
    SC3E (sc3_allocator_free (p3->alloc, qinfo_array[i]));
  }
  SC3E (sc3_allocator_free (p3->alloc, qinfo_array));
  SC3E (p4est3_destroy (&p3));
  return NULL;
}

int
main (int argc, char **argv)
{
  setup_t             st, *t = &st;
  const p4est3_quadrant_vtable_t *qvt;
  sc3_error_t        *e = NULL;

  SC3E_NULL_SET (e, sc3_MPI_Init (&argc, &argv));
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  p4est_init (NULL, SC_LP_ESSENTIAL);

  SC3E_NULL_SET (e, set_parameters (t, &qvt));
  SC3E_NULL_SET (e, perform_test (t, qvt));
  SC3E_NULL_SET (e, clean_up (t));

  SC3E_NULL_REQ (e, !sc_finalize_noabort ());
  SC3E_NULL_SET (e, sc3_MPI_Finalize ());
  SC3X (e);
  return 0;
}
