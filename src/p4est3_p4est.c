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
#include <p4est_algorithms.h>
#include <p4est3_p4est.h>
#else
#include <p8est_bits.h>
#include <p8est_algorithms.h>
#include <p4est3_p8est.h>
#endif

static sc3_error_t *
p4est3_connectivity_p4est_destroy (void *cslf)
{
  p4est_connectivity_t *c4 = (p4est_connectivity_t *) cslf;

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

  /* verify arguments */
  SC3E_RETVAL (pc, NULL);
  SC3A_IS (sc3_allocator_is_valid, alloc);
  SC3A_CHECK (c4 != NULL && p4est_connectivity_is_valid (c4));

  /* create virtual structure */
  memset (cvt, 0, sizeof (*cvt));
  cvt->dim = P4EST_DIM;
  cvt->num_trees = c4->num_trees;
  if (autodestroy) {
    cvt->destroy = p4est3_connectivity_p4est_destroy;
  }

  /* create connectivity */
  SC3E (p4est3_connectivity_new (alloc, &c));

  /* legal to pass stack variable because the virtual table is shallow copied */
  SC3E (p4est3_connectivity_set_vtable (c, cvt, c4));

  /* finalize connectivity */
  SC3E (p4est3_connectivity_setup (c));
  SC3A_IS (p4est3_connectivity_is_setup, c);

  *pc = c;
  return NULL;
}

typedef struct p4est3_p4est_self
{
  int                 autodestroy;
  sc3_allocator_t    *alloc;
  p4est3_connectivity_t *c3;
  p4est_t            *p4;
}
p4est3_p4est_self_t;

static sc3_error_t *
p4est3_p4est_destroy (void *pslf)
{
  p4est3_p4est_self_t *slf = (p4est3_p4est_self_t *) pslf;

  /* leave p4->connectivity alone */
  SC3A_CHECK (slf != NULL);
  SC3A_CHECK (slf->p4 != NULL);
  if (slf->autodestroy) {
    p4est_destroy (slf->p4);
  }
  SC3E (p4est3_connectivity_destroy (&slf->c3));
  SC3E (sc3_allocator_free (slf->alloc, slf));
  return NULL;
}

sc3_error_t        *
p4est3_new_p4est (sc3_allocator_t * alloc, p4est_t * p4,
                  int autodestroy, p4est3_t ** pp3)
{
  p4est3_p4est_self_t *slf;
  p4est3_t           *p3;
  p4est3_vtable_t     spvt, *pvt = &spvt;
  p4est3_quadrant_vtable_t sqvt;

  /* verify arguments */
  SC3E_RETVAL (pp3, NULL);
  SC3A_IS (sc3_allocator_is_valid, alloc);
  SC3A_CHECK (p4 != NULL && p4est_is_valid (p4));

  /* wrap connectivity into a p4est3_connectivity_t object and build context */
  SC3E (sc3_allocator_malloc (alloc, sizeof (p4est3_p4est_self_t), &slf));
  SC3E (p4est3_connectivity_new_p4est (alloc, p4->connectivity, 0, &slf->c3));
  slf->autodestroy = autodestroy;
  slf->alloc = alloc;
  slf->p4 = p4;

  /* create forest and quadrant virtual tables on the stack */
  memset (pvt, 0, sizeof (*pvt));
  pvt->dim = P4EST_DIM;
  pvt->c3 = slf->c3;
  pvt->qvt = &sqvt;
  SC3E (p4est3_quadrant_vtable_p4est (pvt->qvt, 0));
  pvt->destroy = p4est3_p4est_destroy;

  /* create forest */
  SC3E (p4est3_new (alloc, &p3));

  /* legal to pass stack variable because the virtual table is shallow copied */
  SC3E (p4est3_set_vtable (p3, pvt, slf));

  /* finalize forest */
  SC3E (p4est3_setup (p3));
  SC3A_IS (p4est3_is_setup, p3);

  *pp3 = p3;
  return NULL;
}

static              p4est3_gloidx
p4est_quadrant_vtable_num_uniform (int level)
{
  if (level < 0) {
    return -1;
  }
  return p4est3_glopow (P4EST_CHILDREN, level);
}

static int
p4est_quadrant_vtable_is_valid (const void *q, char *reason)
{
  SC3E_TEST (p4est_quadrant_is_valid ((const p4est_quadrant_t *) q), reason);
  SC3E_YES (reason);
}

static int
p4est_quadrant_vtable_is_equal (const void *q1, const void *q2, char *reason)
{
  SC3E_TEST (p4est_quadrant_is_equal ((const p4est_quadrant_t *) q1,
                                      (const p4est_quadrant_t *) q2), reason);
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
p4est_quadrant_vtable_coordinates (const void *q, int n, int *j)
{
  const p4est_quadrant_t *quad = (const p4est_quadrant_t *) q;
  SC3A_CHECK (n == P4EST_DIM);
  SC3A_CHECK (j != NULL);
  j[0] = (int) quad->x;
  j[1] = (int) quad->y;
#ifdef P4_TO_P8
  j[2] = (int) quad->z;
#endif
  return NULL;
}

static sc3_error_t *
p4est_quadrant_vtable_compare (const void *q1, const void *q2, int *j)
{
  SC3A_CHECK (j != NULL);

  *j = p4est_quadrant_compare ((const p4est_quadrant_t *) q1,
                               (const p4est_quadrant_t *) q2);
  return NULL;
}

static sc3_error_t *
p4est_quadrant_vtable_root (void *r)
{
  p4est_quadrant_set_morton ((p4est_quadrant_t *) r, 0, 0);
  return NULL;
}

static sc3_error_t *
p4est_quadrant_vtable_copy (const void *q, void *r)
{
  p4est_quadrant_copy ((const p4est_quadrant_t *) q, (p4est_quadrant_t *) r);
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

sc3_error_t        *
p4est3_quadrant_vtable_p4est (p4est3_quadrant_vtable_t * qvt, int id)
{
  /* check arguments */
  SC3A_CHECK (qvt != NULL);
  SC3A_CHECK (id >= 0);
  memset (qvt, 0, sizeof (p4est3_quadrant_vtable_t));

  /* populate scalar members */
  qvt->id = id;
  qvt->dim = P4EST_DIM;
  qvt->max_level = P4EST_QMAXLEVEL;
  qvt->max_children = P4EST_CHILDREN;
  qvt->quadrant_size = sizeof (p4est_quadrant_t);

  /* populate member functions */
  qvt->quadrant_is_valid = p4est_quadrant_vtable_is_valid;
  qvt->quadrant_is_equal = p4est_quadrant_vtable_is_equal;
  qvt->quadrant_num_uniform = p4est_quadrant_vtable_num_uniform;
  qvt->quadrant_level = p4est_quadrant_vtable_level;
  qvt->quadrant_child_id = p4est_quadrant_vtable_child_id;
  qvt->quadrant_ancestor_id = p4est_quadrant_vtable_ancestor_id;
  qvt->quadrant_coordinates = p4est_quadrant_vtable_coordinates;
  /* quadrant_num_children is not necessary */
  qvt->quadrant_compare = p4est_quadrant_vtable_compare;
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

  /* verify correctness */
  SC3A_IS (p4est3_quadrant_vtable_is_valid, qvt);
  return NULL;
}
