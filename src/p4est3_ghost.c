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
#include <p4est3_internal.h>

#ifdef __cplusplus
extern              "C"
{
#if 0
}
#endif
#endif

typedef struct
{
  p4est_ghost_t     *ghost;
}
p4est3_ghost_fill_data_t;

static sc3_error_t *
p4est3_ghost_fill_callback (p4est3_iterate_face_info_t *fi)
{
  p4est3_ghost_fill_data_t *d = (p4est3_ghost_fill_data_t *) fi->user_data;
  p4est_ghost_t      *ghost = d->ghost;
  p4est3_iterate_face_side_t *fside[2], *gside;
  size_t               nsides;

  ghost = (void *) ghost; /* temporarily to avoid unused variable */

  SC3E(sc3_array_get_elem_count (fi->sides, &nsides));
  SC3A_CHECK(nsides == 2 || nsides == 1);

  if (nsides == 1) {
    /* Nothing to do here. There are no ghosts on a boundary. */
    return NULL;
  }

  /** Check if exactly one side is a ghost */
  sc3_array_index (fi->sides, 0, &fside[0]);
  sc3_array_index (fi->sides, 1, &fside[1]);
  SC3A_CHECK(fside[0]->is_ghost != -1 || fside[1]->is_ghost != -1);
  SC3A_CHECK(fside[0]->is_ghost != 1 && fside[1]->is_ghost != 1);

  if (fside[0]->is_ghost == 0 && fside[1] == 0) {
    /* It's not a ghost. Nothing to do here. */
    return NULL;
  }

  gside = fside[0]->is_ghost == 1 ? fside[0] : fside[1]; /*< ghost side*/
  SC3A_CHECK(gside->is_ghost == 1);

  /** Add ghost to \c ghost->ghosts array */
  /** How to know at what location of the array to place the quadrant?
   * 1. Simple and memory efficient solution is add any ghost to the array and
   *    sort it at the end. It costs additional O(n log n) operations in the
   *    worst case. Which is alright because the Iterator is O(n log n) anyway.
   * TODO: Double check the exact complexity of the Iterator.
   * 2. Iterate in the order of increasing ghosts id. Is it even possible?
  */

  /** Fill its \c piggy3 field */

  /** Contribute to a structure tracking \c tree_offsets */
  /** Contribute to a structure tracking \c proc_offsets */

  /** Do the same for mirrors */

  return NULL;
}

sc3_error_t        *
p4est3_ghost_fill_p4est (p4est3_t * p3, p4est_ghost_t * ghost)
{
  /*TODO: Allocate memory for ghosts outside and before this function call */
  p4est3_ghost_fill_data_t data, *d = &data;
  /* ... */

  d->ghost = ghost;

  ghost->mpisize = p3->mpisize;
  ghost->num_trees = p3->num_trees;
  ghost->btype = P4EST_CONNECT_FACE;

  /* Might be NULL for integration with Dune */
  ghost->mirror_proc_offsets = NULL;
  ghost->mirror_proc_fronts = NULL;
  ghost->mirror_proc_front_offsets = NULL;


  SC3E (p4est3_iterate_face (p3, NULL, p4est3_ghost_fill_callback, d));

  /** Fill \c tree_offsets */
  /** Fill \c proc_offsets */

  return NULL;
}

#ifdef __cplusplus
#if 0
{
#endif
}
#endif