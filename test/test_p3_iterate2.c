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
#ifndef P4_TO_P8
#include <p4est3_p4est.h>
#else
#include <p4est3_p8est.h>
#endif

/* *INDENT-OFF* */
const int           face_inner_nb_2d[4][4] =
{{ 0, 1 },
 { 0, 2 },
 { 1, 3 },
 { 1, 2 }};

const int           face_outer_nb_2d[8] =
{ 0, 2, 1, 3, 0, 1, 2, 3};

const int           face_inner_nface_2d[4][4] =
{{ 1, 0 },
 { 3, 2 },
 { 3, 2 },
 { 1, 0 }};

const int           face_outer_nface_2d[8] =
{ 0, 0, 1, 1, 2, 2, 3, 3 };
 /* *INDENT-ON* */

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
  sc3_array_t *volumes;
  sc3_array_t *faces;
}
callback_data_t;


static sc3_error_t *
make_allocator (setup_t *t)
{
  SC3A_IS (sc3_allocator_is_setup, t->mainalloc);
  SC3E (sc3_allocator_new (t->mainalloc, &t->alloc));
  SC3E (sc3_allocator_setup (t->alloc));

  return NULL;
}

static sc3_error_t *
make_connectivity (setup_t * t, int dim)
{
  SC3E (p4est3_connectivity_new (t->alloc, &t->conn));
  SC3E (p4est3_connectivity_set_dim (t->conn, dim));
  SC3E (p4est3_connectivity_set_num_trees (t->conn, t->num_trees));
  SC3E (p4est3_connectivity_setup (t->conn));

  return NULL;
}

static sc3_error_t *
volume_callback (p4est3_iterate_volume_info_t * vi)
{
  p4est3_iterate_volume_info_t * idx;
  sc3_array_t *out = ((callback_data_t *) vi->user_data)->volumes;
  SC3A_IS (sc3_array_is_setup, out);
  SC3E (sc3_array_push (out, &idx));
  *idx = *vi;

  return NULL;
}

static sc3_error_t *
face_callback (p4est3_iterate_face_info_t * fi)
{
  int nsides, side;
  p4est3_iterate_face_info_t * idx;
  p4est3_iterate_face_side_t * in_side, *out_side;
  sc3_array_t *out = ((callback_data_t *) fi->user_data)->faces;
  SC3A_IS (sc3_array_is_setup, out);
  SC3A_IS (sc3_array_is_setup, fi->sides);
  SC3E (sc3_array_push (out, &idx));
  SC3A_IS (sc3_array_is_setup, idx->sides);
  idx->p3 = fi->p3;
  idx->orientation = fi->orientation;
  idx->tree_boundary = fi->tree_boundary;

  SC3E (sc3_array_get_elem_count (idx->sides, &nsides));
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
make_result_arrays (setup_t *t, p4est3_t *p3, p4est3_quadrant_vtable_t *qvt,
                    sc3_array_t ** voutput, sc3_array_t ** vpredef,
                    sc3_array_t ** foutput, sc3_array_t ** fpredef)
{
  const int nquad = (1 << (qvt->dim * t->level));
  const int nface = qvt->dim * (nquad + (1 << t->level));
  p4est3_iterate_face_info_t *fit;
  p4est3_iterate_volume_info_t *vit;
  p4est3_iterate_face_side_t *sit;
  int i, side;
  char *qit;
  SC3E (array_new (t->alloc, sizeof (p4est3_iterate_volume_info_t), nquad,
        nquad, voutput));
  SC3E (array_new (t->alloc, sizeof (p4est3_iterate_volume_info_t), nquad,
        nquad, vpredef));
  SC3E (array_new (t->alloc, sizeof (p4est3_iterate_face_info_t), nface,
        nface, foutput));
  SC3E (array_new (t->alloc, sizeof (p4est3_iterate_face_info_t), nface,
        nface, fpredef));

  /* allocate memory for sides arrays at faces */
  /* side arrays for foutput have size 0, since we will extend them at
     corresponding callback */
  SC3E (sc3_array_index (*foutput, 0, &fit));
  for (i = 0; i < nface; ++i) {
    SC3E (array_new (t->alloc, sizeof (p4est3_iterate_face_side_t), 0, 0,
                     &fit->sides));
  }

  /* now we'll fill in the reference arrays */
  SC3E (sc3_array_index (*vpredef, 0, &vit));
  SC3E (p4est3_get_quadrants (p3, &qit));
  for (i = 0; i < nquad; ++i, ++vit, qit += qvt->quadrant_size) {
    vit->ntree = 0;
    vit->nquad = i;
    vit->quadrant = (void *) qit;
  }

  SC3E (p4est3_get_quadrants (p3, &qit));
  SC3E (sc3_array_index (*fpredef, 0, &fit));
  if (qvt->dim == 2) {
    /* fill in inner faces */
    for (i = 0; i < 4; ++i) {
      fit->tree_boundary = 0;
      fit->orientation = 0;
      SC3E (array_new (t->alloc, sizeof (p4est3_iterate_face_side_t), 2, 2,
                      &fit->sides));
      SC3E (sc3_array_index (fit->sides, 0, &sit));
      for (side = 0; side < 2; ++side, ++sit) {
        sit->is_ghost = 0;
        sit->ntree = 0;
        sit->nface = face_inner_nface_2d[i][side];
        sit->nquad = face_inner_nb_2d[i][side];
        sit->quadrant = (void *) (qit + sit->nquad);
      }
    }
      /* fill in outer faces */
    for (i = 0; i < 8; ++i) {
      fit->tree_boundary = 1;
      fit->orientation = 0;
      SC3E (array_new (t->alloc, sizeof (p4est3_iterate_face_side_t), 1, 1,
                      &fit->sides));
      SC3E (sc3_array_index (fit->sides, 0, &sit));
      sit->is_ghost = 0;
      sit->ntree = 0;
      sit->nface = face_outer_nface_2d[i];
      sit->nquad = face_outer_nb_2d[i];
      sit->quadrant = qit + sit->nquad;
    }
  }
  else {
    SC3E_UNREACH ("invalid dimension");
  }
  return NULL;
}

static sc3_error_t *
perform_test (setup_t * t, p4est3_t * p3, p4est3_quadrant_vtable_t * qvt,
              callback_data_t *user_data)
{
  SC3E (p4est3_iterate_codim
        (p3, P4EST3_ITERATE_VOLUME || P4EST3_ITERATE_FACE,
         volume_callback, face_callback, NULL, user_data));
  return NULL;
}

static sc3_error_t *
compare_results (setup_t *t, p4est3_t *p3, p4est3_quadrant_vtable_t *qvt,
                 sc3_array_t * voutput, sc3_array_t * vpredef,
                 sc3_array_t * foutput, sc3_array_t * fpredef)
{
  const int nquad = (1 << (qvt->dim * t->level));
  const int nface = qvt->dim * (nquad + (1 << t->level));
  p4est3_iterate_volume_info_t *vit_out, *vit_pre;
  p4est3_iterate_face_info_t *fit_out, *fit_pre;
  p4est3_iterate_face_side_t *sit_out, *sit_pre;
  int i, side, nsides_out, nsides_pre;

  SC3E (sc3_array_index (voutput, 0, &vit_out));
  SC3E (sc3_array_index (vpredef, 0, &vit_pre));
  for (i = 0; i < nquad; ++i, ++vit_out, ++vit_pre) {
    SC3E_DEMAND (vit_out->ntree == vit_pre->ntree, "Volume's ntree differs");
    SC3E_DEMAND (vit_out->nquad == vit_pre->nquad, "Volume's nquad differs");
    SC3E_DEMAND (vit_out->quadrant == vit_pre->quadrant,
                 "Volumes poiner to various quadrants");
  }

  SC3E (sc3_array_index (foutput, 0, &fit_out));
  SC3E (sc3_array_index (fpredef, 0, &fit_pre));
  for (i = 0; i < nface; ++i, ++fit_out, ++fit_pre) {
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
free_arrays (sc3_array_t *voutput, sc3_array_t *vpredef, sc3_array_t *foutput,
             sc3_array_t *fpredef)
{
  p4est3_iterate_face_info_t *fit_out, *fit_pre;
  int nfaces_out, nfaces_pre;
  int i;
  SC3A_IS (sc3_array_is_setup, voutput);
  SC3A_IS (sc3_array_is_setup, vpredef);
  SC3A_IS (sc3_array_is_setup, foutput);
  SC3A_IS (sc3_array_is_setup, fpredef);

  SC3E (sc3_array_destroy (&voutput));
  SC3E (sc3_array_destroy (&vpredef));
  SC3E (sc3_array_get_elem_count (voutput, &nfaces_out));
  SC3E (sc3_array_get_elem_count (vpredef, &nfaces_pre));
  SC3E_DEMAND (nfaces_pre == nfaces_out, "#Sides mismatches");
  SC3E (sc3_array_index (foutput, 0, &fit_out));
  SC3E (sc3_array_index (fpredef, 0, &fit_pre));
  for (i = 0; i < nfaces_out; ++i, ++fit_out, ++fit_pre) {
    SC3E (sc3_array_destroy (&fit_out->sides));
    SC3E (sc3_array_destroy (&fit_pre->sides));
  }
  SC3E (sc3_array_destroy (&foutput));
  SC3E (sc3_array_destroy (&fpredef));
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
  p4est3_quadrant_vtable_t vtable, *qvt = &vtable;
  p4est3_t           *p3;
  sc3_array_t        *vpredef, *fpredef;
  callback_data_t     sud, *user_data = &sud;

  /* v3 standard procedure to isolate memory allocation contexts */
  t->mainalloc = sc3_allocator_nothread ();
  SC3X (sc3_MPI_Init (&argc, &argv));
  t->mpicomm = SC3_MPI_COMM_WORLD;
  SC3X (sc3_MPI_Comm_rank (t->mpicomm, &t->mpirank));
  SC3X (make_allocator (t));
  SC3X (p4est3_quadrant_vtable_p4est (qvt, 0));

  t->level = 1;
  t->num_trees = 1;

#ifdef P4EST_ENABLE_DEBUG
  printf ("l = %d, t = %d\n", t->level, t->num_trees);
#endif /* P4EST_ENABLE_DEBUG */
  SC3X (make_connectivity (t, qvt->dim));
  SC3X (make_new_p4est3 (&p3, t->alloc, t->conn, qvt, t->level));

  SC3X (make_result_arrays
        (t, p3, qvt, &user_data->volumes, &vpredef, &user_data->faces,
         &fpredef));
  SC3X (perform_test (t, p3, qvt, user_data));
  SC3X (compare_results
        (t, p3, qvt, user_data->volumes, vpredef, user_data->faces, fpredef));

  /*destroy forest, that was referenced for others*/
  SC3X (p4est3_destroy (&p3));
  SC3X (p4est3_connectivity_destroy (&t->conn));
  SC3X (free_arrays (user_data->volumes, vpredef, user_data->faces, fpredef));

  SC3X (free_allocator (&t->alloc));
  SC3X (sc3_MPI_Finalize ());
  return 0;
}
