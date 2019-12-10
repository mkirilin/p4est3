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

#include <p4est3_quadrant_vtable.h>

size_t
p4est3_quadrant_size (p4est3_quadrant_vtable_t * qvt)
{
  if (qvt == NULL || qvt->quadrant_size == NULL) {
    return 0;
  }
  return qvt->quadrant_size ();
}

int
p4est3_quadrant_is_valid (p4est3_quadrant_vtable_t * qvt, const void *q)
{
  if (qvt == NULL) {
    return 0;
  }
  return qvt->quadrant_is_valid == NULL || qvt->quadrant_is_valid (q);
}

sc3_error_t        *
p4est3_quadrant_root (p4est3_quadrant_vtable_t * qvt, void *r)
{
  SC3A_CHECK (qvt != NULL && qvt->quadrant_root != NULL);

  return qvt->quadrant_root (r);
}

sc3_error_t        *
p4est3_quadrant_child (p4est3_quadrant_vtable_t * qvt,
                       const void *q, void *r, int i)
{
  SC3A_CHECK (qvt != NULL && qvt->quadrant_child != NULL);

  return qvt->quadrant_child (q, r, i);
}

sc3_error_t        *
p4est3_quadrant_parent (p4est3_quadrant_vtable_t * qvt,
                        const void *q, void *r)
{
  SC3A_CHECK (qvt != NULL && qvt->quadrant_parent != NULL);

  return qvt->quadrant_parent (q, r);
}
