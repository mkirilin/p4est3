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

#ifndef P4EST3_H
#define P4EST3_H

#include <p4est3_base.h>

#ifdef __cplusplus
extern              "C"
{
#if 0
}
#endif
#endif

#define P4EST3_CHILDREN 4

/*------------------------- the connectivity -------------------------*/

typedef struct p4est3_connectivity_attr p4est3_connectivity_attr_t;
typedef struct p4est3_connectivity p4est3_connectivity_t;

/* attributes: set comm, maybe more */
/* ... */

p4est3_connectivity_t *p4est3_connectivity_new (p4est3_connectivity_attr_t *
                                                ca);
p4est3_connectivity_t *p4est3_connectivity_new_unitcube (void);

/* function-based construction of connectivity */
/* ... */

void                p4est3_connectivity_setup (p4est3_connectivity_t * conn);
void                p4est3_connectivity_ref (p4est3_connectivity_t * c3);
void                p4est3_connectivity_unref (p4est3_connectivity_t ** c3);
void                p4est3_connectivity_destroy (p4est3_connectivity_t ** c3);

/*------------------------- the p4est object -------------------------*/

typedef struct p4est3_attr p4est3_attr_t;
typedef struct p4est3 p4est3_t;

/* attributes: set minlevel, user data size, etc. */
/* ... */

p4est3_t           *p4est3_new (p4est3_connectivity_t * c3,
                                p4est3_attr_t * pa);

/*----------------------- accessing quadrants ------------------------*/

#if 0

p4est3_quadrant_t  *p4est3_quadrant_range (p4est3_t * p3,
                                           p4est3_topidx_t tbegin,
                                           p4est3_topidx_t tend,
                                           p4est3_locidx_t qbegin,
                                           p4est3_locidx_t qend);
void                p4est3_quadrant_restore (p4est3_t * p3,
                                             p4est3_quadrant_t * q3,
                                             p4est3_topidx_t tbegin,
                                             p4est3_topidx_t tend,
                                             p4est3_locidx_t qbegin,
                                             p4est3_locidx_t qend);

typedef struct p4est3_access_attr p4est3_access_attr_t;
typedef struct p4est3_access p4est3_access_t;

p4est3_access_t    *p4est3_iter_new (p4est3_t * p3,
                                     p4est3_access_attr_t * ia);
p4est3_quadrant_t  *p4est3_iter_end (p4est3_access_t * pi);
void                p4est3_iter_inc (p4est3_access_t * pi);

p4est3_access_t    *p4est3_access_new (p4est3_t * p3,
                                       p4est3_access_attr_t * ia);
void                p4est3_access_ref (p4est3_access_t * a3);
void                p4est3_access_unref (p4est3_access_t ** a3);
void                p4est3_access_destroy (p4est3_access_t ** a3);

p4est3_locidx_t     p4est3_access_get_length (p4est3_access_t * a3);
p4est3_quadrant_t  *p4est3_access_get_begin (p4est3_access_t * a3);
p4est3_quadrant_t  *p4est3_access_index (p4est3_access_t * a3,
                                         p4est3_locidx_t li);

#endif /* 0 */

#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif /* !P4EST3_H */
