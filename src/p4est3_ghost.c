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
#include <p4est3_search.h>

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
  p4est_quadrant_t    q;
  size_t              nsides;
  p4est_topidx_t     *ghosts_in_tree;
  p4est_gloidx_t     *ghosts_in_proc, proc_owner, global_qid;
  int                 coords[P4EST_DIM], level;

  /* TODO: Check if the ghost is already in the layer via hash table.
           If so, skip it. */

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

  /* Convert p3 quad to p2 quad throughout coordinates because of qvt */
  SC3E (p4est3_quadrant_coordinates (fi->p3->qvt, gside->quadrant, coords));
  SC3E (p4est3_quadrant_level (fi->p3->qvt, gside->quadrant, &q.level));
  q.x = coords[0];
  q.y = coords[1];
#ifdef P4_TO_P8
  q.z = coords[2];
#endif
  /* Find ghost proc owner */
  proc_owner = fi->p3->mpirank;
  global_qid = (p4est3_gloidx) gside->nquad + fi->p3->gtroffset[gside->ntree];
  /** TODO: We just did it in Iterator to fill callback data.
   *        Now we do it again here. Think on a way to optimize it. */
  SC3E (p4est3_search_lower_bound64
         (global_qid, fi->p3->goffset, fi->p3->mpisize + 1, &proc_owner));
  if (fi->p3->goffset[proc_owner] > global_qid) {
    SC3A_CHECK (proc_owner > 0);
    proc_owner--;
  }

  /** Fill its \c piggy3 field */
  q.p.piggy3.which_tree = gside->ntree;
  q.p.piggy3.local_num = global_qid - fi->p3->goffset[proc_owner];

  /* Push back to ghosts array */
  *(p4est_quadrant_t *) sc_array_push(&(ghost->ghosts)) = q;

  /** Contribute to a structure tracking \c tree_offsets */

  ghosts_in_tree =
    (p4est_locidx_t *) sc_array_index (ghost->tree_offsets, gside->ntree);
  (*ghosts_in_tree)++;

  /** Contribute to a structure tracking \c proc_offsets */
  ghosts_in_proc =
    (p4est_locidx_t *) sc_array_index (ghost->proc_offsets, proc_owner);
  (*ghosts_in_proc)++;

  /** Do the same for mirrors */
  /* Check 1st todo */

  return NULL;
}

sc3_error_t        *
p4est3_ghost_fill_p4est (p4est3_t * p3, p4est_ghost_t * ghost)
{
  /*TODO: Allocate memory for ghosts outside and before this function call */
  /* Ensure tree_ and proc_offsets are pre-initialized by 0 */
  p4est3_ghost_fill_data_t data, *d = &data;
  int                i;
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

  /** Accumulate \c tree_offsets */
  for (i = 1; i < ghost->num_trees; i++) {
    (*(p4est_locidx_t *) sc_array_index (ghost->tree_offsets, i)) +=
      (*(p4est_locidx_t *) sc_array_index (ghost->tree_offsets, i - 1));
  }

  /** Accumulate \c proc_offsets */
  for (i = 1; i < ghost->mpisize; i++) {
    (*(p4est_locidx_t *) sc_array_index (ghost->proc_offsets, i)) +=
      (*(p4est_locidx_t *) sc_array_index (ghost->proc_offsets, i - 1));
  }

  return NULL;
}

#ifdef __cplusplus
#if 0
{
#endif
}
#endif