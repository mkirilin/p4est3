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
  int                 enable_faces;
  p4est3_topidx       num_trees;
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
    SC3E_TEST (!c->cvt->enable_faces == !c->enable_faces, reason);
    SC3E_TEST (c->cvt->num_trees == c->num_trees, reason);
  }
  SC3E_TEST (0 < c->dim && c->dim <= 3, reason);
  SC3E_TEST (0 < c->num_trees, reason);
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
  c->enable_faces = 1;
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
  c->enable_faces = cvt->enable_faces;
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
p4est3_connectivity_set_enable_faces (p4est3_connectivity_t * c,
                                      int enable_faces)
{
  SC3A_IS (p4est3_connectivity_is_new, c);
  c->enable_faces = enable_faces;
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
p4est3_connectivity_new_unitsquare (sc3_allocator_t * alloc,
                                    p4est3_connectivity_t ** pc)
{
  p4est3_connectivity_t *c;

  SC3E_RETVAL (pc, NULL);
  SC3E (p4est3_connectivity_new (alloc, &c));
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
  SC3E (p4est3_connectivity_setup (c));
  SC3A_IS (p4est3_connectivity_is_setup, c);

  *pc = c;
  return NULL;
}
