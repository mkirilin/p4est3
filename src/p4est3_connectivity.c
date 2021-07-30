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
#include <sc3_array.h>
#include <sc3_refcount.h>

/* *INDENT-OFF* */
static const int    face_permutation_refs[6][6] =
{{ 0, 1, 1, 0, 0, 1 },
 { 2, 0, 0, 1, 1, 0 },
 { 2, 0, 0, 1, 1, 0 },
 { 0, 2, 2, 0, 0, 1 },
 { 0, 2, 2, 0, 0, 1 },
 { 2, 0, 0, 2, 2, 0 }};
static const int    face_permutation_sets[3][4] =
{{ 1, 2, 5, 6 },
 { 0, 3, 4, 7 },
 { 0, 4, 3, 7 }};
static const int    face_permutations[8][4] =
{{ 0, 1, 2, 3 },
 { 0, 2, 1, 3 },
 { 1, 0, 3, 2 },
 { 1, 3, 0, 2 },
 { 2, 0, 3, 1 },
 { 2, 3, 0, 1 },
 { 3, 1, 2, 0 },
 { 3, 2, 1, 0 }};
/* *INDENT-ON* */

struct p4est3_connectivity
{
  sc3_refcount_t      rc;
  sc3_allocator_t    *alloc;
  int                 setup;

  /* this connectivity may be wrapping a virtual implementation */
  p4est3_connectivity_vtable_t scvt, *cvt;
  void               *slf;

  /* parameters fixed after setup call */
  int                 dim;
  int                 num_faces;
  int                 num_orient;
  p4est3_topidx       num_trees;
  int                 half_children;

  int                *face_corner;  /**< Corresdondance between face and
                                         quadrant's corners adjacent to it */
};

static sc3_error_t *
p4est3_connectivity_set_face_child_id (const p4est3_connectivity_t * c)
{
  int                *fc = c->face_corner;
  SC3A_IS (p4est3_connectivity_is_new, c);
  if (c->dim == 2) {
    SC3E (sc3_allocator_calloc (c->alloc, 8, sizeof (int), &fc));
    /* *INDENT-OFF* */
    fc[0] = 0; fc[1] = 2;
    fc[2] = 1; fc[3] = 3;
    fc[4] = 0; fc[5] = 1;
    fc[6] = 2; fc[7] = 3;
    /* *INDENT-ON* */
  }
  else if (c->dim == 3) {
    SC3E (sc3_allocator_calloc (c->alloc, 24, sizeof (int), &fc));
    /* *INDENT-OFF* */
    fc[0] = 0; fc[1] = 2; fc[2] = 4; fc[3] = 6;
    fc[4] = 1; fc[5] = 3; fc[6] = 5; fc[7] = 7;
    fc[8] = 0; fc[9] = 1; fc[10] = 4; fc[11] = 5;
    fc[12] = 2; fc[13] = 3; fc[14] = 6; fc[15] = 7;
    fc[16] = 0; fc[17] = 1; fc[18] = 2; fc[19] = 3;
    fc[20] = 4; fc[21] = 5; fc[22] = 6; fc[23] = 7;
    /* *INDENT-ON* */
  }
  else {
    SC3E_UNREACH ("wrong dimension");
  }
  return NULL;
}

int                 p4est3_connectivity_vtable_is_valid
  (const p4est3_connectivity_vtable_t * cvt, char *reason)
{
  SC3E_TEST (cvt != NULL, reason);
  SC3E_TEST (0 < cvt->dim && cvt->dim <= 3, reason);
  SC3E_TEST (0 < cvt->num_trees, reason);
  SC3E_YES (reason);
}

int
p4est3_connectivity_is_valid (const p4est3_connectivity_t * c, char *reason)
{
  SC3E_TEST (c != NULL, reason);
  SC3E_IS (sc3_refcount_is_valid, &c->rc, reason);
  SC3E_IS (sc3_allocator_is_setup, c->alloc, reason);
  if (c->cvt != NULL) {
    SC3E_IS (p4est3_connectivity_vtable_is_valid, c->cvt, reason);
    SC3E_TEST (c->cvt->dim == c->dim, reason);
    SC3E_TEST (c->cvt->num_trees == c->num_trees, reason);
  }
  SC3E_TEST (0 < c->dim && c->dim <= 3, reason);
  SC3E_TEST (0 < c->num_trees, reason);
  if (c->setup) {
    SC3E_TEST (c->num_faces == 2 * c->dim, reason);
    SC3E_TEST (c->num_orient == 2 * (c->dim - 1), reason);
    SC3E_TEST (c->half_children == 1 << (c->dim - 1), reason);
  }
  SC3E_YES (reason);
}

int
p4est3_connectivity_is_new (const p4est3_connectivity_t * c, char *reason)
{
  SC3E_IS (p4est3_connectivity_is_valid, c, reason);
  SC3E_TEST (!c->setup, reason);
  SC3E_YES (reason);
}

int
p4est3_connectivity_is_setup (const p4est3_connectivity_t * c, char *reason)
{
  SC3E_IS (p4est3_connectivity_is_valid, c, reason);
  SC3E_TEST (c->setup, reason);
  SC3E_YES (reason);
}

sc3_error_t        *
p4est3_connectivity_new (sc3_allocator_t * alloc, p4est3_connectivity_t ** pc)
{
  p4est3_connectivity_t *c;

  SC3E_RETVAL (pc, NULL);
  SC3A_IS (sc3_allocator_is_setup, alloc);

  SC3E (sc3_allocator_ref (alloc));
  SC3E (sc3_allocator_calloc_one (alloc, sizeof (p4est3_connectivity_t), &c));
  SC3E (sc3_refcount_init (&c->rc));
  c->alloc = alloc;
  c->dim = 2;
  c->num_trees = 1;
  SC3A_IS (p4est3_connectivity_is_new, c);

  *pc = c;
  return NULL;
}

sc3_error_t        *
p4est3_connectivity_set_vtable (p4est3_connectivity_t * c,
                                p4est3_connectivity_vtable_t * cvt, void *slf)
{
  SC3A_IS (p4est3_connectivity_is_new, c);
  SC3A_IS (p4est3_connectivity_vtable_is_valid, cvt);

  /* make shallow copy of virtual table */
  c->dim = cvt->dim;
  c->num_trees = cvt->num_trees;
  *(c->cvt = &c->scvt) = *cvt;
  c->slf = slf;
  return NULL;
}

sc3_error_t        *
p4est3_connectivity_set_dim (p4est3_connectivity_t * c, int dim)
{
  SC3A_IS (p4est3_connectivity_is_new, c);
  SC3A_CHECK (0 < dim && dim <= 3);
  c->dim = dim;
  return NULL;
}

sc3_error_t        *
p4est3_connectivity_set_num_trees (p4est3_connectivity_t * c,
                                   p4est3_topidx num_trees)
{
  SC3A_IS (p4est3_connectivity_is_new, c);
  SC3A_CHECK (num_trees > 0);
  c->num_trees = num_trees;
  return NULL;
}

sc3_error_t        *
p4est3_connectivity_setup (p4est3_connectivity_t * c)
{
  SC3A_IS (p4est3_connectivity_is_new, c);

  c->num_faces = 2 * c->dim;
  c->num_orient = 2 * (c->dim - 1);
  c->half_children = 1 << (c->dim - 1);

  SC3E (p4est3_connectivity_set_face_child_id (c));

  c->setup = 1;
  SC3A_IS (p4est3_connectivity_is_setup, c);
  return NULL;
}

sc3_error_t        *
p4est3_connectivity_ref (p4est3_connectivity_t * c)
{
  SC3A_IS (p4est3_connectivity_is_setup, c);
  SC3E (sc3_refcount_ref (&c->rc));
  return NULL;
}

sc3_error_t        *
p4est3_connectivity_unref (p4est3_connectivity_t * c)
{
  int                 waslast;

  SC3A_IS (p4est3_connectivity_is_setup, c);
  SC3E (sc3_refcount_unref (&c->rc, &waslast));

  /* This is a hard check that we do not unref below a count of one. */
  SC3E_DEMAND (!waslast, "Connectivity unrefd below a count of one");
  return NULL;
}

sc3_error_t        *
p4est3_connectivity_destroy (p4est3_connectivity_t ** pc)
{
  sc3_allocator_t    *alloc;
  p4est3_connectivity_t *c;

  SC3E_INULLP (pc, c);
  SC3A_IS (p4est3_connectivity_is_valid, c);

  /* This is a hard check for a reference count of exactly one. */
  SC3E_DEMIS (sc3_refcount_is_last, &c->rc);

  /* destruction callback if one was provided */
  if (c->cvt != NULL && c->cvt->destroy != NULL) {
    SC3E (c->cvt->destroy (c->slf));
  }

  /* remove allocation */
  alloc = c->alloc;
  SC3E (sc3_allocator_free (alloc, c));
  SC3E (sc3_allocator_unref (&alloc));

  /* nothing is left */
  return NULL;
}

sc3_error_t        *
p4est3_connectivity_get_dim (const p4est3_connectivity_t * c, int *pdim)
{
  SC3A_IS (p4est3_connectivity_is_setup, c);
  SC3A_CHECK (pdim != NULL);

  *pdim = c->dim;
  return NULL;
}

sc3_error_t        *
p4est3_connectivity_get_num_trees (const p4est3_connectivity_t * c,
                                   p4est3_topidx * pnum_trees)
{
  SC3A_IS (p4est3_connectivity_is_setup, c);
  SC3A_CHECK (pnum_trees != NULL);

  *pnum_trees = c->num_trees;
  return NULL;
}

sc3_error_t        *
p4est3_connectivity_get_face (const p4est3_connectivity_t * c,
                              p4est3_topidx * which_tree, int *nface,
                              int *orient)
{
  /* formally check arguments */
  SC3A_IS (p4est3_connectivity_is_setup, c);
  SC3A_CHECK (which_tree != NULL);
  SC3A_CHECK (nface != NULL);
  SC3E_RETVAL (orient, 0);

  /* verify input values */
  SC3A_CHECK (0 <= *which_tree && *which_tree < c->num_trees);
  SC3A_CHECK (0 <= *nface && *nface < c->num_faces);

  /* return physical boundary unless specified otherwise */
  if (c->cvt != NULL && c->cvt->get_face != NULL) {
    SC3E (c->cvt->get_face (c->slf, which_tree, nface, orient));
  }

  return NULL;
}

sc3_error_t        *
p4est3_connectivity_get_face_child_id (const p4est3_connectivity_t * c,
                                       int nface, int *i)
{
  SC3A_IS (p4est3_connectivity_is_setup, c);
  SC3A_CHECK (i != NULL);

  SC3A_CHECK (0 <= nface && nface < c->num_faces);
  SC3A_CHECK (0 <= *i && *i < c->num_orient);

  *i = c->face_corner[nface * c->num_orient + *i];
  SC3A_CHECK (0 <= *i && *i < (1 << c->dim));
  return NULL;
}

sc3_error_t        *
p4est3_connectivity_face_neighbor_face_corner (const p4est3_connectivity_t *
                                               c, int *fc, int f, int nf,
                                               int o)
{
  int                 pref, pset;

  SC3A_CHECK (fc != NULL);
  SC3A_CHECK (0 <= *fc && *fc < c->half_children);
  SC3A_CHECK (0 <= f && f < c->num_faces);
  SC3A_CHECK (0 <= nf && nf < c->num_faces);
  SC3A_CHECK (0 <= o && o < c->half_children);

  if (c->dim == 2) {
    *fc = *fc ^ o;
  }
  else if (c->dim == 3) {
    pref = face_permutation_refs[f][nf];
    pset = face_permutation_sets[pref][o];
    *fc = face_permutations[pset][*fc];
  }
  else {
    SC3E_UNREACH ("wrong dimension");
  }
  SC3A_CHECK (0 <= *fc && *fc < c->half_children);
  return NULL;
}

sc3_error_t        *
p4est3_connectivity_new_unitsquare (sc3_allocator_t * alloc,
                                    p4est3_connectivity_t ** pc)
{
  p4est3_connectivity_t *c;

  SC3E_RETVAL (pc, NULL);
  SC3E (p4est3_connectivity_new (alloc, &c));
  SC3E (p4est3_connectivity_set_dim (c, 2));
  SC3E (p4est3_connectivity_setup (c));
  SC3A_IS (p4est3_connectivity_is_setup, c);

  *pc = c;
  return NULL;
}

sc3_error_t        *
p4est3_connectivity_new_unitcube (sc3_allocator_t * alloc,
                                  p4est3_connectivity_t ** pc)
{
  p4est3_connectivity_t *c;

  SC3E_RETVAL (pc, NULL);
  SC3E (p4est3_connectivity_new (alloc, &c));
  SC3E (p4est3_connectivity_set_dim (c, 3));
  SC3E (p4est3_connectivity_setup (c));
  SC3A_IS (p4est3_connectivity_is_setup, c);

  *pc = c;
  return NULL;
}
