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
#include <p4est3_quadrant_vtable.h>

static              size_t
p4est_quadrant_size (void)
{
  return sizeof (p4est_quadrant_t);
}

static void
p4est_quadrant_root (p4est_quadrant_t * q)
{
  P4EST_ASSERT (q != NULL);
  p4est_quadrant_set_morton (q, 0, 0);
}

void
p4est3_quadrant_vtable_p4est (p4est3_quadrant_vtable_t * qvt)
{
  P4EST_ASSERT (qvt != NULL);

  qvt->quadrant_size = p4est_quadrant_size;
  qvt->quadrant_root = (p4est3_quadrant_root_t) p4est_quadrant_root;
  qvt->quadrant_child = (p4est3_quadrant_child_t) p4est_quadrant_child;
  qvt->quadrant_parent = (p4est3_quadrant_parent_t) p4est_quadrant_parent;
  qvt->quadrant_successor = (p4est3_quadrant_successor_t) NULL;
  qvt->quadrant_predecessor = (p4est3_quadrant_predecessor_t) NULL;
}

size_t
p4est3_quadrant_size (p4est3_quadrant_vtable_t * qvt)
{
  P4EST_ASSERT (qvt != NULL && qvt->quadrant_size != NULL);

  return qvt->quadrant_size ();
}

void
p4est3_quadrant_root (p4est3_quadrant_vtable_t * qvt, void *q)
{
  P4EST_ASSERT (qvt != NULL && qvt->quadrant_root != NULL);

  qvt->quadrant_root (q);
}

void
p4est3_quadrant_child (p4est3_quadrant_vtable_t * qvt,
                       const void *q, void *r, int i)
{
  P4EST_ASSERT (qvt != NULL && qvt->quadrant_child != NULL);

  qvt->quadrant_child (q, r, i);
}

void
p4est3_quadrant_parent (p4est3_quadrant_vtable_t * qvt,
                        const void *q, void *r)
{
  P4EST_ASSERT (qvt != NULL && qvt->quadrant_parent != NULL);

  qvt->quadrant_parent (q, r);
}
