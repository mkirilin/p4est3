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

#ifndef P4EST3_BITS_MORT_H
#define P4EST3_BITS_MORT_H

#include <p4est.h>
#include <p4est3_quadrant_vtable.h>

#define P4EST3_QUADRANT_MORT_LEN(n, l) ((uint64_t) (n) << P4EST_DIM * (P4EST_MAXLEVEL - (l)))

typedef struct p4est3_quadrant_mort
{
  /*@{ */
  uint64_t            coords;  /**< coordinates */
  /*@} */
  int8_t              level,    /**< level of refinement */
                      pad8;     /**< padding */
  int16_t             pad16;    /**< padding */
}
p4est3_quadrant_mort_t;

void                p4est3_quadrant_mort_vtable (p4est3_quadrant_vtable_t
                                                 * qvt);
#endif
