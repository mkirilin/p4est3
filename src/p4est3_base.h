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

#ifndef P4EST3_BASE_H
#define P4EST3_BASE_H

#include <sc3_alloc.h>
#include <sc3_error.h>
#include <stdint.h>

#ifdef __cplusplus
extern              "C"
{
#if 0
}
#endif
#endif

#define P3A_CHECK(x) SC3A_CHECK (x)
#define P3A_STACK(x) SC3A_STACK (x)
#define P3E(f) SC3E (f)
#define P3E_DEMAND(f) SC3E_DEMAND (f)

typedef int         p4est3_topidx;
typedef int         p4est3_locidx;
typedef long        p4est3_gloidx;

#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif /* !P4EST3_BASE_H */
