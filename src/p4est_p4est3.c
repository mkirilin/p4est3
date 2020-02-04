/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2010 The University of Texas System
  Additional copyright (C) 2011 individual authors
  Written by Carsten Burstedde, Lucas C. Wilcox, and Tobin Isaac

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

#include <p4est_bits.h>
#include <p4est_p4est3.h>

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
  SC3E_TEST (p4est_quadrant_is_valid ((p4est_quadrant_t *) q), reason);
  SC3E_YES (reason);
}

static sc3_error_t *
p4est_quadrant_vtable_root (void *q)
{
  p4est_quadrant_set_morton ((p4est_quadrant_t *) q, 0, 0);
  return NULL;
}

static sc3_error_t *
p4est_quadrant_vtable_child (const void *q, void *r, int i)
{
  p4est_quadrant_child
    ((const p4est_quadrant_t *) q, (p4est_quadrant_t *) r, i);
  return NULL;
}

static sc3_error_t *
p4est_quadrant_vtable_parent (const void *q, void *r)
{
  p4est_quadrant_parent
    ((const p4est_quadrant_t *) q, (p4est_quadrant_t *) r);
  return NULL;
}

void
p4est_quadrant_vtable (p4est3_quadrant_vtable_t * qvt)
{
  if (qvt == NULL) {
    return;
  }
  qvt->max_level = p4est_vtable_max_level;
  qvt->num_children = p4est_vtable_num_children;
  qvt->quadrant_size = p4est_quadrant_vtable_size;
  qvt->quadrant_is_valid = p4est_quadrant_vtable_is_valid;
  qvt->quadrant_root = p4est_quadrant_vtable_root;
  qvt->quadrant_child = p4est_quadrant_vtable_child;
  qvt->quadrant_parent = p4est_quadrant_vtable_parent;
  qvt->quadrant_successor = (p4est3_quadrant_successor_t) NULL;
  qvt->quadrant_predecessor = (p4est3_quadrant_predecessor_t) NULL;
}
