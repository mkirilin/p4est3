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
  p4est3_topidx       num_trees;
};

int
p4est3_connectivity_is_valid (const p4est3_connectivity_t * c, char *reason)
{
  SC3E_TEST (c != NULL, reason);
  SC3E_IS (sc3_refcount_is_valid, &c->rc, reason);
  SC3E_IS (sc3_allocator_is_setup, c->alloc, reason);
  if (c->cvt != NULL) {
    /* is there anything we should check? */
  }
  else {
    SC3E_TEST (c->num_trees > 0, reason);
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
  SC3E_ALLOCATOR_CALLOC (alloc, p4est3_connectivity_t, 1, c);
  SC3E (sc3_refcount_init (&c->rc));
  c->alloc = alloc;
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
  SC3A_CHECK (cvt != NULL);

  /* make deep copy of virtual table */
  *(c->cvt = &c->scvt) = *cvt;
  c->slf = slf;
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
p4est3_connectivity_unref (p4est3_connectivity_t ** pc)
{
  int                 waslast;
  sc3_allocator_t    *alloc;
  p4est3_connectivity_t *c;

  SC3E_INOUTP (pc, c);
  SC3A_IS (p4est3_connectivity_is_valid, c);

  SC3E (sc3_refcount_unref (&c->rc, &waslast));
  if (waslast) {
    *pc = NULL;

    /* destruction callback if one was provided */
    if (c->cvt != NULL && c->cvt->destroy != NULL) {
      SC3E (c->cvt->destroy (c->slf));
    }

    /* remove allocation */
    alloc = c->alloc;
    SC3E_ALLOCATOR_FREE (alloc, p4est3_connectivity_t, c);
    SC3E (sc3_allocator_unref (&alloc));
  }
  return NULL;
}

sc3_error_t        *
p4est3_connectivity_destroy (p4est3_connectivity_t ** pc)
{
  p4est3_connectivity_t *c;

  SC3E_INULLP (pc, c);
  SC3E_DEMIS (sc3_refcount_is_last, &c->rc);
  SC3E (p4est3_connectivity_unref (&c));

  SC3A_CHECK (c == NULL);
  return NULL;
}

sc3_error_t        *
p4est3_connectivity_get_num_trees (const p4est3_connectivity_t * c,
                                   p4est3_topidx * pnum_trees)
{
  SC3A_IS (p4est3_connectivity_is_setup, c);
  SC3A_CHECK (pnum_trees != NULL);

  if (c->cvt != NULL) {
    SC3A_CHECK (c->cvt->get_num_trees != NULL);
    SC3E (c->cvt->get_num_trees (c->slf, pnum_trees));
  }
  else {
    *pnum_trees = c->num_trees;
  }
  return NULL;
}

/* TODO this is demo/convenience code; move away */

typedef struct p4est3_connectivity_ntslf
{
  sc3_allocator_t    *alloc;
  p4est3_topidx       num_trees;

  /* no need really to allocate this here since it is deep copied */
  p4est3_connectivity_vtable_t scvt;
}
p4est3_connectivity_ntslf_t;

static sc3_error_t *
p4est3_connectivity_gnt (void *vslf, p4est3_topidx * pnt)
{
  p4est3_connectivity_ntslf_t *slf = (p4est3_connectivity_ntslf_t *) vslf;
  SC3A_CHECK (slf != NULL);
  SC3A_CHECK (pnt != NULL);

  *pnt = slf->num_trees;
  return NULL;
}

static sc3_error_t *
p4est3_connectivity_dstr (void *vslf)
{
  p4est3_connectivity_ntslf_t *slf = (p4est3_connectivity_ntslf_t *) vslf;
  SC3A_CHECK (slf != NULL);

  SC3E_ALLOCATOR_FREE (slf->alloc, p4est3_connectivity_ntslf_t, slf);
  return NULL;
}

sc3_error_t        *
p4est3_connectivity_new_num_trees (sc3_allocator_t * alloc,
                                   p4est3_topidx num_trees,
                                   p4est3_connectivity_t ** pc)
{
  p4est3_connectivity_ntslf_t *slf;
  p4est3_connectivity_t *c;

  SC3E_RETVAL (pc, NULL);
  SC3A_CHECK (num_trees > 0);

  /* create virtual structure */
  SC3E_ALLOCATOR_CALLOC (alloc, p4est3_connectivity_ntslf_t, 1, slf);
  slf->alloc = alloc;
  slf->num_trees = num_trees;
  slf->scvt.get_num_trees = p4est3_connectivity_gnt;
  slf->scvt.destroy = p4est3_connectivity_dstr;

  /* create connectivity */
  SC3E (p4est3_connectivity_new (alloc, &c));
  SC3E (p4est3_connectivity_set_vtable (c, &slf->scvt, slf));
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
