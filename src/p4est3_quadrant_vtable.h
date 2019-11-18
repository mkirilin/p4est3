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

#ifndef P4EST3_QUADRANT_VTABLE
#define P4EST3_QUADRANT_VTABLE

#include <p4est_base.h>

SC_EXTERN_C_BEGIN;

/* *INDENT-OFF* */
typedef size_t      (*p4est3_quadrant_size_t) (void);
typedef void        (*p4est3_quadrant_out_t) (void *r);
typedef void        (*p4est3_quadrant_in_out_t) (const void *q, void *r);

typedef p4est3_quadrant_out_t p4est3_quadrant_root_t;
typedef void        (*p4est3_quadrant_child_t) (const void *q, void *r,
                                                int i);
typedef p4est3_quadrant_in_out_t p4est3_quadrant_parent_t;
typedef p4est3_quadrant_in_out_t p4est3_quadrant_successor_t;
typedef p4est3_quadrant_in_out_t p4est3_quadrant_predecessor_t;
/* *INDENT-ON* */

typedef struct p4est3_quadrant_vtable
{
  p4est3_quadrant_size_t quadrant_size;
  p4est3_quadrant_root_t quadrant_root;
  p4est3_quadrant_child_t quadrant_child;
  p4est3_quadrant_parent_t quadrant_parent;
  p4est3_quadrant_successor_t quadrant_successor;
  p4est3_quadrant_predecessor_t quadrant_predecessor;
}
p4est3_quadrant_vtable_t;

void                p4est3_quadrant_vtable_p4est (p4est3_quadrant_vtable_t *
                                                  qvt);

size_t              p4est3_quadrant_size (p4est3_quadrant_vtable_t * qvt);
void                p4est3_quadrant_root (p4est3_quadrant_vtable_t * qvt,
                                          void *q);
void                p4est3_quadrant_child (p4est3_quadrant_vtable_t * qvt,
                                           const void *q, void *r, int i);
void                p4est3_quadrant_parent (p4est3_quadrant_vtable_t * qvt,
                                            const void *q, void *r);

SC_EXTERN_C_END;

#endif /* !P4EST3_QUADRANT_VTABLE */
