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
#else
#include <p4est3_p8est.h>
#endif

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
  p4est3_iterate_volume_info_t *idx;
  sc3_array_t        *out = ((callback_data_t *) vi->user_data)->volumes;
  SC3A_IS (sc3_array_is_setup, out);
  SC3E (sc3_array_push (out, &idx));
  *idx = *vi;

  return NULL;
}

static sc3_error_t *
face_callback (p4est3_iterate_face_info_t * fi)
{
  p4est3_iterate_face_info_t *idx;
  p4est3_iterate_face_side_t *in_side, *out_side;
  sc3_array_t        *out = ((callback_data_t *) fi->user_data)->faces;
  int                 nsides, side;

  SC3A_IS (sc3_array_is_setup, out);
  SC3A_IS (sc3_array_is_setup, fi->sides);
  SC3E (sc3_array_push (out, &idx));
  SC3A_IS (sc3_array_is_setup, idx->sides);
  idx->p3 = fi->p3;
  idx->orientation = fi->orientation;
  idx->tree_boundary = fi->tree_boundary;

  SC3E (sc3_array_get_elem_count (fi->sides, &nsides));
  SC3E (sc3_array_index (fi->sides, 0, &in_side));
  for (side = 0; side < nsides; ++side, ++in_side) {
    SC3E (sc3_array_push (idx->sides, &out_side));
    *out_side = *in_side;
  }
  return NULL;
}

static sc3_error_t *
make_new_p4est3 (p4est3_t ** p3, sc3_allocator_t * alloc,
                 p4est3_connectivity_t * conn,
                 p4est3_quadrant_vtable_t * qvt, int level)
{
  SC3A_IS (sc3_allocator_is_setup, alloc);

  /* create p4est object with connectivity */
  SC3E (p4est3_new (alloc, p3));
  /* SC3E (p4est3_set_comm (*p3, mpicomm, 1)); */
  SC3E (p4est3_set_connectivity (*p3, conn));
  SC3E (p4est3_set_quadrant_vtable (*p3, qvt));
  SC3E (p4est3_set_level (*p3, level));
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
make_ref_array_volume (setup_t * t, p4est3_t * p3,
                       p4est3_quadrant_vtable_t * qvt, sc3_array_t ** vpredef)
{
  const int           nquad = (1 << (qvt->dim * t->level)) * t->num_trees;
  p4est3_iterate_volume_info_t *vit;
  p4est3_tree_t      *tree;
  int                 i, ntree;
  char               *qit;

  SC3E (array_new (t->alloc, sizeof (p4est3_iterate_volume_info_t), nquad,
                   nquad, vpredef));

  SC3E (sc3_array_index (*vpredef, 0, &vit));
  for (ntree = 0; ntree < t->num_trees; ++ntree) {
    SC3E (p4est3_tree_index (p3, ntree, &tree));
    qit = tree->tquads;
    for (i = 0; i < tree->num_quads; ++i, ++vit, qit += qvt->quadrant_size) {
      vit->ntree = ntree;
      vit->nquad = i;
      vit->quadrant = (void *) qit;
    }
  }
  return NULL;
}

static sc3_error_t *
iterate_unimesh_inner_simple_children (setup_t * t, p4est3_t * p3,
                                       p4est3_quadrant_vtable_t * qvt,
                                       p4est3_gloidx first_qid,
                                       p4est3_tree_t * tree,
                                       sc3_array_t * fpredef)
{
  p4est3_iterate_face_info_t *finfo;
  p4est3_iterate_face_side_t *sinfo;
  int                 i, nfaces;

  if (qvt->dim == 2) {
    nfaces = 2 * qvt->dim;
    int                 nface_predef[4][2] =
      { {1, 0}, {3, 2}, {3, 2}, {1, 0} };
    int                 nquad_predef[4][2] =
      { {0L, 1L}, {0L, 2L}, {1L, 3L}, {2L, 3L} };

    for (i = 0; i < nfaces; ++i) {
      SC3E (sc3_array_push (fpredef, &finfo));
      SC3E (array_new
            (t->alloc, sizeof (p4est3_iterate_face_side_t), 2, 0,
             &finfo->sides));
      finfo->orientation = 0;
      finfo->tree_boundary = 0;
      SC3E (sc3_array_push (finfo->sides, &sinfo));
      sinfo->ntree = tree->treeid;
      sinfo->is_ghost = 0;

      sinfo->nface = nface_predef[i][0];
      sinfo->nquad = first_qid + nquad_predef[i][0];
      sinfo->quadrant =
        (void *) (tree->tquads + qvt->quadrant_size * sinfo->nquad);

      SC3E (sc3_array_push (finfo->sides, &sinfo));
      sinfo->ntree = tree->treeid;
      sinfo->is_ghost = 0;

      sinfo->nface = nface_predef[i][1];
      sinfo->nquad = first_qid + nquad_predef[i][1];
      sinfo->quadrant =
        (void *) (tree->tquads + qvt->quadrant_size * sinfo->nquad);
    }
  }
  else {
    SC3E_UNREACH ("wrong dimension");
  }
  return NULL;
}

static sc3_error_t *
fill_face_info (sc3_array_t * fpredef, int nquad, int face, int face_neighbor,
                int nsides, int orientation, int is_tree_boundary,
                p4est3_tree_t * tree, p4est3_tree_t * tree_neighbor,
                p4est3_quadrant_vtable_t * qvt, void *q, void *r, setup_t * t)
{
  p4est3_iterate_face_info_t *finfo;
  p4est3_iterate_face_side_t *sinfo;
  p4est3_gloidx       linear_id;
  p4est3_topidx       ntree;

  SC3A_CHECK (nsides == 1 || nsides == 2);
  SC3E (sc3_array_push (fpredef, &finfo));
  SC3E (array_new (t->alloc, sizeof (p4est3_iterate_face_side_t), 2, 0,
                   &finfo->sides));
  finfo->orientation = orientation;
  finfo->tree_boundary = is_tree_boundary;
  SC3E (sc3_array_push (finfo->sides, &sinfo));
  sinfo->ntree = tree->treeid;
  sinfo->is_ghost = 0;

  sinfo->nface = face;
  sinfo->nquad = nquad;
  sinfo->quadrant =
    (void *) (tree->tquads + qvt->quadrant_size * sinfo->nquad);

  if (nsides == 1) {
    return NULL;
  }

  SC3E (sc3_array_push (finfo->sides, &sinfo));
  sinfo->ntree = tree_neighbor->treeid;
  sinfo->is_ghost = 0;

  sinfo->nface = face_neighbor;
  if (finfo->tree_boundary) {
    ntree = tree->treeid;
    SC3E (p4est3_quadrant_face_neighbor_extra
          (qvt, t->conn, r, face, &ntree, &face_neighbor, q));
    SC3A_CHECK (sinfo->nface == face_neighbor % (2 * qvt->dim));
  }
  else {
    SC3E (p4est3_quadrant_face_neighbor (qvt, r, face, q));
  }
  SC3E (p4est3_quadrant_linear_id (qvt, q, t->level, &linear_id));
  sinfo->nquad = (p4est3_locidx) linear_id;
  sinfo->quadrant =
    (void *) (tree_neighbor->tquads + qvt->quadrant_size * sinfo->nquad);

  return NULL;
}

static sc3_error_t *
iterate_unimesh_inner_face_compl (setup_t * t, p4est3_quadrant_vtable_t * qvt,
                                  p4est3_gloidx nquads_compl, int level,
                                  void *q, void *r, p4est3_tree_t * tree,
                                  sc3_array_t * fpredef)
{
  p4est3_gloidx       nquads_per_next_level = 1L << (qvt->dim * (level + 1));
  p4est3_gloidx       nquads_per_level = 1L << (qvt->dim * level);
  p4est3_gloidx       start_id = nquads_compl - nquads_per_next_level;
  p4est3_gloidx       bound_id;
  int                 i, face, nfaces = 2 * qvt->dim, nsides = 2;
  int                 is_tree_boundary = 0, orientation = 0;
  int                 nface_predef[4][2] = { {1, 0}, {3, 2}, {3, 2}, {1, 0} };

  if (qvt->dim == 2) {
    int                 bound_dir[4] = { 0, 1, 1, 0 };
    int                 start_ids[4];
    start_ids[0] = start_id;
    start_ids[1] = start_id;
    start_ids[2] = start_id + nquads_per_level;
    start_ids[3] = start_id + 2 * nquads_per_level;
    bound_id =
      start_id + nquads_per_next_level / (p4est3_gloidx) qvt->max_children;
    int                 bound_coords[2];
    int                 coords[2];
    SC3E (p4est3_quadrant_morton (qvt, t->level, bound_id - 1, q));
    SC3E (p4est3_quadrant_coordinates (qvt, q, qvt->dim, bound_coords));
    for (face = 0; face < nfaces; ++face) {
      SC3E (p4est3_quadrant_morton (qvt, t->level, start_ids[face], r));
      for (i = 0; i < nquads_per_level - 1; ++i) {
        SC3E (p4est3_quadrant_coordinates (qvt, r, qvt->dim, coords));
        if (coords[bound_dir[face]] == bound_coords[bound_dir[face]]) {
          SC3E (fill_face_info (fpredef, start_ids[face] + i,
                                nface_predef[face][0], nface_predef[face][1],
                                nsides, orientation, is_tree_boundary,
                                tree, tree, qvt, q, r, t));
        }
        SC3E (p4est3_quadrant_successor (qvt, r, r));
      }
      SC3E (p4est3_quadrant_coordinates (qvt, r, qvt->dim, coords));
      if (coords[bound_dir[face]] == bound_coords[bound_dir[face]]) {
        SC3E (fill_face_info (fpredef, start_ids[face] + nquads_per_level - 1,
                              nface_predef[face][0], nface_predef[face][1],
                              nsides, orientation, is_tree_boundary,
                              tree, tree, qvt, q, r, t));
      }
    }
  }
  else {
    SC3E_UNREACH ("wrong dimension");
  }
  return NULL;
}

static sc3_error_t *
iterate_unimesh_tree_boundary_face (setup_t * t, p4est3_t * p3,
                                    sc3_array_t * fpredef,
                                    p4est3_topidx ntree, void *q, void *r,
                                    p4est3_tree_t * tree,
                                    p4est3_quadrant_vtable_t * qvt)
{
  /* iterate over tree boundary */
  const int           is_tree_bound = 1;
  p4est3_topidx       ntree_neighbor;
  p4est3_tree_t      *tree_neighbor;
  p4est3_gloidx       nquads_per_level = 1L << (qvt->dim * t->level);
  int                 i, face, face_neighbor, ori, nsides;
  int                 nfaces = 2 * qvt->dim;
  if (qvt->dim == 2) {
    int                 nf[2];
    for (face = 0; face < nfaces; ++face) {
      nsides = 2;
      ntree_neighbor = ntree;
      face_neighbor = face;
      SC3E (p4est3_connectivity_get_face
            (t->conn, &ntree_neighbor, &face_neighbor, &ori));
      if (ntree_neighbor < ntree) {
        continue;
      }
      if (ntree_neighbor == ntree) {
        nsides = 1;
      }
      SC3E (p4est3_tree_index (p3, ntree_neighbor, &tree_neighbor));
      SC3E (p4est3_quadrant_morton (qvt, t->level, 0L, r));
      for (i = 0; i < nquads_per_level - 1; ++i) {
        SC3E (p4est_quadrant_tree_boundary (qvt, r, nf));
        if (nf[0] == face || nf[1] == face || nf[0] == -2 || nf[1] == -2) {
          SC3E (fill_face_info (fpredef, i, face, face_neighbor, nsides, ori,
                                is_tree_bound, tree, tree_neighbor, qvt, q, r,
                                t));
        }
        SC3E (p4est3_quadrant_successor (qvt, r, r));
      }
      SC3E (p4est_quadrant_tree_boundary (qvt, r, nf));
      if (nf[0] == face || nf[1] == face || nf[0] == -2 || nf[1] == -2) {
        SC3E (fill_face_info (fpredef, nquads_per_level - 1, face,
                              face_neighbor, nsides, ori, is_tree_bound,
                              tree, tree_neighbor, qvt, q, r, t));
      }
    }
  }
  else {
    SC3E_UNREACH ("wrong dimension");
  }
  return NULL;
}

static sc3_error_t *
iterate_unimesh_face (setup_t * t, p4est3_t * p3,
                      p4est3_quadrant_vtable_t * qvt, sc3_array_t * fpredef)
{
  p4est3_gloidx       nquads_compl = 0;
  p4est3_topidx       ntree;
  p4est3_tree_t      *tree;
  int                 level = 0;
  int                *level2nchildren;  /*how many processed quads of level level */
  void               *q, *r;
  SC3E (sc3_allocator_calloc
        (t->alloc, t->level + 1, sizeof (int), &level2nchildren));
  SC3E (sc3_allocator_calloc_one (t->alloc, qvt->quadrant_size, &q));
  SC3E (sc3_allocator_calloc_one (t->alloc, qvt->quadrant_size, &r));
  memset (level2nchildren, 0, t->level * sizeof (int));

  for (ntree = 0; ntree < t->num_trees; ++ntree) {
    SC3E (p4est3_tree_index (p3, ntree, &tree));
    while (level < t->level) {
      SC3A_CHECK (level2nchildren[level] <= qvt->max_children);
      if (level2nchildren[level] == qvt->max_children) {
        if (level == 0) {
          SC3E (iterate_unimesh_inner_simple_children
                (t, p3, qvt, nquads_compl, tree, fpredef));
          nquads_compl += qvt->max_children;
        }
        else {
          SC3E (iterate_unimesh_inner_face_compl
                (t, qvt, nquads_compl, level, q, r, tree, fpredef));
        }
        level2nchildren[level] = 0;
        level2nchildren[++level]++;
      }
      else {
        level = 0;
        level2nchildren[level] += 4;
      }
    }
    SC3E (iterate_unimesh_tree_boundary_face
          (t, p3, fpredef, ntree, q, r, tree, qvt));
  }
  SC3E (sc3_allocator_free (t->alloc, level2nchildren));
  SC3E (sc3_allocator_free (t->alloc, q));
  SC3E (sc3_allocator_free (t->alloc, r));
  return NULL;
}

static sc3_error_t *
make_result_arrays (setup_t * t, p4est3_t * p3,
                    p4est3_quadrant_vtable_t * qvt, sc3_array_t ** voutput,
                    sc3_array_t ** vpredef, sc3_array_t ** foutput,
                    sc3_array_t ** fpredef)
{
  p4est3_iterate_face_info_t *fit;
  const int           nquad = (1 << (qvt->dim * t->level)) * t->num_trees;
  int                 nface, i;

  SC3E (array_new (t->alloc, sizeof (p4est3_iterate_volume_info_t), nquad,
                   0, voutput));
  SC3E (make_ref_array_volume (t, p3, qvt, vpredef));

  SC3E (array_new
        (t->alloc, sizeof (p4est3_iterate_face_info_t), 0, 0, fpredef));
  SC3E (iterate_unimesh_face (t, p3, qvt, *fpredef));
  SC3E (sc3_array_get_elem_count (*fpredef, &nface));
  SC3E (array_new (t->alloc, sizeof (p4est3_iterate_face_info_t),
                   nface, nface, foutput));

  /* allocate memory for sides arrays at faces */
  /* side arrays for foutput have size 0, since we will extend them at
     corresponding callback */
  SC3E (sc3_array_index (*foutput, 0, &fit));
  for (i = 0; i < nface; ++i, ++fit) {
    SC3E (array_new (t->alloc, sizeof (p4est3_iterate_face_side_t), 2, 0,
                     &fit->sides));
  }
  SC3E (sc3_array_resize (*foutput, 0));

  return NULL;
}

static sc3_error_t *
perform_test (setup_t * t, p4est3_t * p3, p4est3_quadrant_vtable_t * qvt,
              callback_data_t * user_data)
{
  SC3E (p4est3_iterate_codim
        (p3, P4EST3_ITERATE_VOLUME || P4EST3_ITERATE_FACE,
         volume_callback, face_callback, NULL, user_data));
  return NULL;
}

static sc3_error_t *
compare_results (setup_t * t, p4est3_t * p3, p4est3_quadrant_vtable_t * qvt,
                 sc3_array_t * voutput, sc3_array_t * vpredef,
                 sc3_array_t * foutput, sc3_array_t * fpredef)
{
  const int           nquad = (1 << (qvt->dim * t->level));
  p4est3_iterate_volume_info_t *vit_out, *vit_pre;
  p4est3_iterate_face_info_t *fit_out, *fit_pre;
  p4est3_iterate_face_side_t *sit_out, *sit_pre;

  int                 nface_out, nface_pre;
  int                 i, side, nsides_out, nsides_pre;
  int                 out_coord[2], pre_coord[2];
  int                 out_l, pre_l;

  SC3E (sc3_array_index (voutput, 0, &vit_out));
  SC3E (sc3_array_index (vpredef, 0, &vit_pre));
  for (i = 0; i < nquad; ++i, ++vit_out, ++vit_pre) {
    SC3E_DEMAND (vit_out->ntree == vit_pre->ntree, "Volume's ntree differs");
    SC3E_DEMAND (vit_out->nquad == vit_pre->nquad, "Volume's nquad differs");

    /* degub purposes only */
    SC3E (p4est3_quadrant_coordinates
          (qvt, vit_out->quadrant, qvt->dim, out_coord));
    SC3E (p4est3_quadrant_coordinates
          (qvt, vit_pre->quadrant, qvt->dim, pre_coord));
    SC3E (p4est3_quadrant_level (qvt, vit_out->quadrant, &out_l));
    SC3E (p4est3_quadrant_level (qvt, vit_pre->quadrant, &pre_l));
    /* degub purposes only */

    SC3E_DEMAND (vit_out->quadrant == vit_pre->quadrant,
                 "Volumes poiner to various quadrants");
  }

  SC3E (sc3_array_get_elem_count (foutput, &nface_out));
  SC3E (sc3_array_get_elem_count (fpredef, &nface_pre));
  SC3E_DEMAND (nface_out == nface_pre, "Face: number of faces differs");
  SC3E (sc3_array_index (foutput, 0, &fit_out));
  SC3E (sc3_array_index (fpredef, 0, &fit_pre));
  for (i = 0; i < nface_out; ++i, ++fit_out, ++fit_pre) {
    SC3E_DEMAND (fit_out->tree_boundary == fit_pre->tree_boundary,
                 "Face: tree boundary property mismatches");
    SC3E_DEMAND (fit_out->orientation == fit_pre->orientation,
                 "Face: orientation property mismatches");
    SC3E (sc3_array_get_elem_count (fit_out->sides, &nsides_out));
    SC3E (sc3_array_get_elem_count (fit_pre->sides, &nsides_pre));
    SC3E_DEMAND (nsides_out == nsides_pre, "#Sides mismatches");
    SC3E (sc3_array_index (fit_out->sides, 0, &sit_out));
    SC3E (sc3_array_index (fit_pre->sides, 0, &sit_pre));
    for (side = 0; side < nsides_pre; ++side, ++sit_out, ++sit_pre) {
      SC3E_DEMAND (sit_out->is_ghost == sit_pre->is_ghost,
                   "Face: ghoust property mismatch");
      SC3E_DEMAND (sit_out->ntree == sit_pre->ntree, "Faces' ntree differs");
      SC3E_DEMAND (sit_out->nface == sit_pre->nface, "Faces' nface differs");
      SC3E_DEMAND (sit_out->nquad == sit_pre->nquad, "Faces' nquad differs");
      SC3E_DEMAND (sit_out->quadrant == sit_pre->quadrant,
                   "Faces poiner to various quadrants");
    }
  }

  return NULL;
}

static sc3_error_t *
free_arrays (sc3_array_t * voutput, sc3_array_t * vpredef,
             sc3_array_t * foutput, sc3_array_t * fpredef)
{
  p4est3_iterate_face_info_t *fit_pre;
  p4est3_iterate_face_info_t *fit_out;
  int                 i, nfaces_pre;
  int                 nfaces_out;

  SC3A_IS (sc3_array_is_setup, voutput);
  SC3A_IS (sc3_array_is_setup, vpredef);
  SC3A_IS (sc3_array_is_setup, fpredef);
  SC3A_IS (sc3_array_is_setup, foutput);

  SC3E (sc3_array_destroy (&voutput));
  SC3E (sc3_array_destroy (&vpredef));

  SC3E (sc3_array_get_elem_count (fpredef, &nfaces_pre));
  SC3E (sc3_array_get_elem_count (foutput, &nfaces_out));
  SC3E_DEMAND (nfaces_pre == nfaces_out, "#Sides mismatches");
  SC3E (sc3_array_index (fpredef, 0, &fit_pre));
  SC3E (sc3_array_index (foutput, 0, &fit_out));
  for (i = 0; i < nfaces_pre; ++i, ++fit_out, ++fit_pre) {
    SC3E (sc3_array_destroy (&fit_pre->sides));
    SC3E (sc3_array_destroy (&fit_out->sides));
  }
  SC3E (sc3_array_destroy (&fpredef));
  SC3E (sc3_array_destroy (&foutput));
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
  setup_t             st, *t = &st;
  p4est3_t           *p3;
  sc3_array_t        *vpredef, *fpredef;
  callback_data_t     sud, *user_data = &sud;
  p4est3_quadrant_vtable_t vtable, *qvt = &vtable;
  int                 l_face, r_face, ori;
  int                 nfaces, nori;

  /* v3 standard procedure to isolate memory allocation contexts */
  t->mainalloc = sc3_allocator_nothread ();
  SC3X (sc3_MPI_Init (&argc, &argv));
  t->mpicomm = SC3_MPI_COMM_WORLD;
  SC3X (sc3_MPI_Comm_rank (t->mpicomm, &t->mpirank));
  SC3X (make_allocator (t));
  SC3X (p4est3_quadrant_vtable_p4est (qvt, 0));
  nfaces = 2 * qvt->dim;
  nori = 1 << (qvt->dim - 1);

  t->level = 3;
  t->num_trees = 2;

#ifdef P4EST_ENABLE_DEBUG
  printf ("l = %d, t = %d\n", t->level, t->num_trees);
#endif /* P4EST_ENABLE_DEBUG */

  for (ori = 0; ori < nori; ++ori) {
    for (l_face = 0; l_face < nfaces; ++l_face) {
      for (r_face = 0; r_face < nfaces; ++r_face) {
        SC3X (make_connectivity (t, qvt->dim, l_face, r_face, ori));
        SC3X (make_new_p4est3 (&p3, t->alloc, t->conn, qvt, t->level));

        SC3X (make_result_arrays
              (t, p3, qvt, &user_data->volumes, &vpredef, &user_data->faces,
               &fpredef));
        SC3X (perform_test (t, p3, qvt, user_data));
        SC3X (compare_results
              (t, p3, qvt, user_data->volumes, vpredef, user_data->faces,
               fpredef));

        /*destroy forest, that was referenced for others */
        SC3X (free_arrays
              (user_data->volumes, vpredef, user_data->faces, fpredef));
        SC3X (p4est3_destroy (&p3));
        SC3X (p4est3_connectivity_destroy (&t->conn));
      }
    }
  }
  SC3X (free_allocator (&t->alloc));
  SC3X (sc3_MPI_Finalize ());
  return 0;
}
