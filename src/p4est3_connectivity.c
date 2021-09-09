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
};

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
p4est3_connectivity_get_face_transform (const p4est3_connectivity_t * c,
                                        int32_t iface, int *itree,
                                        sc3_array_t * transform)
{
  SC3A_CHECK (itree != NULL);

  int                 itree_neighbor = *itree;
  int                 iface_neighbor = iface;
  int                 orient, *iter;
  SC3A_IS (p4est3_connectivity_is_setup, c);
  SC3A_IS (sc3_array_is_setup, transform);
  SC3A_CHECK (0 <= iface && iface < c->num_faces);
  SC3A_CHECK (0 <= *itree && *itree < c->num_trees);

  SC3E (p4est3_connectivity_get_face
        (c, &itree_neighbor, &iface_neighbor, &orient));

  *itree = *itree == itree_neighbor ? -1 : itree_neighbor;
  if (*itree == itree_neighbor && iface == iface_neighbor) {
    SC3A_CHECK (orient == 0);
    *itree = -1;
    return NULL;
  }

  SC3A_CHECK (0 <= iface_neighbor && iface_neighbor < c->num_faces);
  SC3A_CHECK (0 <= orient && orient <= c->half_children);

#ifdef P4EST_ENABLE_DEBUG
  int                 ecount;
  SC3E (sc3_array_get_elem_count (transform, &ecount));
  /* check is transform has a desired length */
  SC3A_CHECK (ecount == 9);
#endif

  SC3E (sc3_array_index (transform, 0, &iter));
  if (c->dim == 2) {
    iter[2] = iface / 2;
    iter[1] = 0;
    iter[0] = 1 - iter[2];
    iter[5] = iface_neighbor / 2;
    iter[4] = 0;
    iter[3] = 1 - iter[5];
    iter[6] = orient;
    iter[7] = 0;
    iter[8] = 2 * (iface & 1) + (iface_neighbor & 1);
  }
  else if (c->dim == 3) {
    int                 reverse;

#ifdef P4EST_ENABLE_DEBUG
    int                 i;
    int                *my_axis;
    int                *target_axis;
#endif

    iter[0] = iface < 2 ? 1 : 0;
    iter[1] = iface < 4 ? 2 : 1;
    iter[2] = iface / 2;
    reverse =
      face_permutation_refs[0][iface] ^
      face_permutation_refs[0][iface_neighbor] ^ (orient == 0 || orient == 3);
    iter[3 + reverse] = iface_neighbor < 2 ? 1 : 0;
    iter[3 + !reverse] = itree_neighbor < 4 ? 2 : 1;
    iter[5] = iface_neighbor / 2;
    reverse = (face_permutation_refs[iface][iface_neighbor] == 1);
    iter[6 + reverse] = (orient & 1);
    iter[6 + !reverse] = (orient >> 1);
    iter[8] = 2 * (iface & 1) + (iface_neighbor & 1);

#ifdef P4EST_ENABLE_DEBUG
    SC3E (sc3_array_index (transform, 0, &my_axis));
    SC3E (sc3_array_index (transform, 3, &target_axis));
    for (i = 0; i < 3; ++i) {
      SC3A_CHECK (0 <= my_axis[i] && my_axis[i] < 3);
      SC3A_CHECK (0 <= target_axis[i] && target_axis[i] < 3);
    }
    SC3A_CHECK (my_axis[0] != my_axis[1] &&
                my_axis[0] != my_axis[2] && my_axis[1] != my_axis[2]);
    SC3A_CHECK (target_axis[0] != target_axis[1] &&
                target_axis[0] != target_axis[2] &&
                target_axis[1] != target_axis[2]);
#endif

  }

  return NULL;
}

sc3_error_t        *
p4est3_connectivity_get_face_child_id (const p4est3_connectivity_t * c,
                                       int nface, int fcorner, int *childid)
{
  SC3A_IS (p4est3_connectivity_is_setup, c);
  SC3A_CHECK (childid != NULL);

  SC3A_CHECK (0 <= nface && nface < c->num_faces);
  SC3A_CHECK (0 <= fcorner && fcorner < c->num_orient);

  if (c->dim == 2) {
    *childid = face_corners_2d[nface][fcorner];
  }
  else {
    *childid = face_corners_3d[nface][fcorner];
  }

  SC3A_CHECK (0 <= *childid && *childid < (1 << c->dim));
  return NULL;
}

sc3_error_t        *
p4est3_connectivity_get_neighbor_face_corner (const p4est3_connectivity_t * c,
                                              int f, int nf, int o, int *fc)
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
