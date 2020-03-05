/*
  This file is part of p4est, version 3
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

#ifndef P4_TO_P8
#include <p4est_bits.h>
#include <p4est_p4est3.h>
#else
#include <p8est_bits.h>
#include <p8est_p4est3.h>
#endif

static sc3_error_t *
p4est3_connectivity_p4est_get_num_trees (void *vslf,
                                         p4est3_topidx * pnum_trees)
{
  p4est_connectivity_t *c4 = (p4est_connectivity_t *) vslf;

  SC3A_CHECK (c4 != NULL);
  SC3A_CHECK (pnum_trees != NULL);
  *pnum_trees = c4->num_trees;
  return NULL;
}

static sc3_error_t *
p4est3_connectivity_p4est_destroy (void *vslf)
{
  p4est_connectivity_t *c4 = (p4est_connectivity_t *) vslf;

  SC3A_CHECK (c4 != NULL);
  p4est_connectivity_destroy (c4);
  return NULL;
}

sc3_error_t        *
p4est3_connectivity_new_p4est (sc3_allocator_t * alloc,
                               p4est_connectivity_t * c4, int autodestroy,
                               p4est3_connectivity_t ** pc)
{
  p4est3_connectivity_t *c;
  p4est3_connectivity_vtable_t scvt, *cvt = &scvt;

  SC3E_RETVAL (pc, NULL);

  /* create virtual structure */
  memset (cvt, 0, sizeof (*cvt));
  cvt->get_num_trees = p4est3_connectivity_p4est_get_num_trees;
  if (autodestroy) {
    cvt->destroy = p4est3_connectivity_p4est_destroy;
  }

  /* create connectivity */
  SC3E (p4est3_connectivity_new (alloc, &c));
  SC3E (p4est3_connectivity_set_vtable (c, cvt, c4));
  SC3E (p4est3_connectivity_setup (c));
  SC3A_IS (p4est3_connectivity_is_setup, c);

  *pc = c;
  return NULL;
}

static int
p4est_vtable_max_level (void)
{
  return P4EST_QMAXLEVEL;
}

static int
p4est_vtable_num_children (void)
{
  return P4EST_CHILDREN;
}

static              size_t
p4est_quadrant_vtable_size (void)
{
  return sizeof (p4est_quadrant_t);
}

static int
p4est_quadrant_vtable_is_valid (const void *q, char *reason)
{
  SC3E_TEST (p4est_quadrant_is_valid ((const p4est_quadrant_t *) q), reason);
  SC3E_YES (reason);
}

static sc3_error_t *
p4est_quadrant_vtable_level (const void *q, int *l)
{
  SC3A_CHECK (q != NULL);
  SC3A_CHECK (l != NULL);

  *l = ((const p4est_quadrant_t *) q)->level;
  return NULL;
}

static sc3_error_t *
p4est_quadrant_vtable_child_id (const void *q, int *j)
{
  SC3A_CHECK (j != NULL);

  *j = p4est_quadrant_child_id ((const p4est_quadrant_t *) q);
  return NULL;
}

static sc3_error_t *
p4est_quadrant_vtable_ancestor_id (const void *q, int i, int *j)
{
  SC3A_CHECK (j != NULL);

  *j = p4est_quadrant_ancestor_id ((const p4est_quadrant_t *) q, i);
  return NULL;
}

static sc3_error_t *
p4est_quadrant_vtable_root (void *r)
{
  p4est_quadrant_set_morton ((p4est_quadrant_t *) r, 0, 0);
  return NULL;
}

static sc3_error_t *
p4est_quadrant_vtable_parent (const void *q, void *r)
{
  p4est_quadrant_parent
    ((const p4est_quadrant_t *) q, (p4est_quadrant_t *) r);
  return NULL;
}

static sc3_error_t *
p4est_quadrant_vtable_copy (const void *q, void *r)
{
  p4est_quadrant_copy ((const p4est_quadrant_t *) q, (p4est_quadrant_t *) r);
  return NULL;
}

static sc3_error_t *
p4est_quadrant_vtable_predecessor (const void *q, void *r)
{
  p4est_quadrant_predecessor
    ((const p4est_quadrant_t *) q, (p4est_quadrant_t *) r);
  return NULL;
}

static sc3_error_t *
p4est_quadrant_vtable_successor (const void *q, void *r)
{
  p4est_quadrant_successor
    ((const p4est_quadrant_t *) q, (p4est_quadrant_t *) r);
  return NULL;
}

static sc3_error_t *
p4est_quadrant_vtable_child (const void *q, int i, void *r)
{
  p4est_quadrant_child
    ((const p4est_quadrant_t *) q, (p4est_quadrant_t *) r, i);
  return NULL;
}

static sc3_error_t *
p4est_quadrant_vtable_ancestor (const void *q, int l, void *r)
{
  p4est_quadrant_ancestor
    ((const p4est_quadrant_t *) q, l, (p4est_quadrant_t *) r);
  return NULL;
}

static sc3_error_t *
p4est_quadrant_vtable_first_descendant (const void *q, int l, void *r)
{
  p4est_quadrant_first_descendant
    ((const p4est_quadrant_t *) q, (p4est_quadrant_t *) r, l);
  return NULL;
}

static sc3_error_t *
p4est_quadrant_vtable_last_descendant (const void *q, int l, void *r)
{
  p4est_quadrant_last_descendant
    ((const p4est_quadrant_t *) q, (p4est_quadrant_t *) r, l);
  return NULL;
}

static sc3_error_t *
p4est_quadrant_vtable_morton (int level, p4est_gloidx_t id, void *r)
{
  p4est_quadrant_set_morton ((p4est_quadrant_t *) r, level, (uint64_t) id);
  return NULL;
}

void
p4est_quadrant_vtable (p4est3_quadrant_vtable_t * qvt, int id)
{
  if (qvt == NULL) {
    return;
  }
  memset (qvt, 0, sizeof (p4est3_quadrant_vtable_t));
  qvt->id = id;
  qvt->dim = P4EST_DIM;
  qvt->max_level = p4est_vtable_max_level;
  qvt->num_children = p4est_vtable_num_children;
  qvt->quadrant_size = p4est_quadrant_vtable_size;
  qvt->quadrant_is_valid = p4est_quadrant_vtable_is_valid;
  qvt->quadrant_level = p4est_quadrant_vtable_level;
  qvt->quadrant_child_id = p4est_quadrant_vtable_child_id;
  qvt->quadrant_ancestor_id = p4est_quadrant_vtable_ancestor_id;
  qvt->quadrant_root = p4est_quadrant_vtable_root;
  qvt->quadrant_copy = p4est_quadrant_vtable_copy;
  qvt->quadrant_parent = p4est_quadrant_vtable_parent;
  qvt->quadrant_predecessor = p4est_quadrant_vtable_predecessor;
  qvt->quadrant_successor = p4est_quadrant_vtable_successor;
  qvt->quadrant_child = p4est_quadrant_vtable_child;
  qvt->quadrant_ancestor = p4est_quadrant_vtable_ancestor;
  qvt->quadrant_first_descendant = p4est_quadrant_vtable_first_descendant;
  qvt->quadrant_last_descendant = p4est_quadrant_vtable_last_descendant;
  qvt->quadrant_morton = p4est_quadrant_vtable_morton;
}
