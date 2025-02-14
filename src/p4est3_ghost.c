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

#include <p4est3_ghost.h>

#ifdef __cplusplus
extern              "C"
{
#if 0
}
#endif
#endif

typedef struct
{
  p4est_ghost_t     *ghosts;
}
p4est3_ghost_fill_data_t;

static sc3_error_t *
p4est3_ghost_fill_callback (p4est3_iterate_face_info_t *fi)
{
  return NULL;
}

sc3_error_t        *
p4est3_ghost_fill_p4est (p4est3_t * p3, p4est_ghost_t ** ghost)
{
  /*TODO: Allocate memory for ghosts */
  p4est3_ghost_fill_data_t data, *d = &data;
  /* ... */




  SC3E (p4est3_iterate_face (p3, NULL, p4est3_ghost_fill_callback, d));

  return NULL;
}

#ifdef __cplusplus
#if 0
{
#endif
}
#endif