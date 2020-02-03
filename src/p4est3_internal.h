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

#ifndef P4EST3_INTERNAL_H
#define P4EST3_INTERNAL_H

#include <p4est3.h>
#include <sc3_refcount.h>

struct p4est3
{
  sc3_refcount_t      rc;
  sc3_allocator_t    *alloc;
  int                 setup;

  sc3_MPI_Comm_t      mpicomm;
  int                 commdup;
  int                 mpisize;
  int                 mpirank;

  p4est3_quadrant_vtable_t *qvt;
  p4est3_topidx       num_trees;
  int                 level;
};

#ifdef __cplusplus
extern              "C"
{
#if 0
}
#endif
#endif

#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif /* !P4EST3_INTERNAL_H */
