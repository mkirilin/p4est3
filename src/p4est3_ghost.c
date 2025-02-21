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

/* A ghost quadrant index is uniquefied by adding the owning process. */
typedef struct ghost_hash_key
{
  p4est_locidx_t      qid;    /* local index within owner process */
  p4est_locidx_t      owner;
}
ghost_hash_key_t;

/* Provide an allocator for the key data as well as the hash map itself */
typedef struct ghost_hash_data
{
  sc_mempool_t       *ckeys;    /* memory pool for allocating the hash keys */
  sc_hash_t          *chash;    /* the hash map links keys without copying */
  p4est_locidx_t      added;    /* count each coordinate point just once */
  p4est_locidx_t      duped;    /* count attempts to add more than once */
}
ghost_hash_data_t;

/* Calculate a hash function for a ghost index */
static unsigned
ghost_hash_fn (const void *v, const void *u)
{
  uint32_t           q, o, z;
  const ghost_hash_key_t *k = (ghost_hash_key_t *) v;

  P4EST_ASSERT (k != NULL);
  q = (uint32_t) k->qid;
  o = (uint32_t) k->owner;
  z = (uint32_t) 0;

  sc_hash_final(q, o, z);

  return (unsigned) z;
}

/* Determine whether two ghosts are equal */
static int
ghost_equal_fn (const void *v1, const void *v2, const void *u)
{
  const ghost_hash_key_t *k1 = (ghost_hash_key_t *) v1;
  const ghost_hash_key_t *k2 = (ghost_hash_key_t *) v2;

  P4EST_ASSERT (k1 != NULL);
  P4EST_ASSERT (k2 != NULL);

  return (k1->qid == k2->qid && k1->owner == k2->owner);
}



typedef struct p4est3_ghost_fill_data
{
  p4est_ghost_t     *ghost;
  ghost_hash_data_t *hdata;
}
p4est3_ghost_fill_data_t;

static sc3_error_t *
p4est3_ghost_fill_callback (p4est3_iterate_face_info_t *fi)
{
  p4est3_ghost_fill_data_t *d = (p4est3_ghost_fill_data_t *) fi->user_data;
  p4est_ghost_t      *ghost = d->ghost;
  p4est3_iterate_face_side_t *fside[2], *gside, *mside;
  p4est_quadrant_t    q;
  size_t              nsides;
  p4est_gloidx_t      proc_owner, global_qid;
  ghost_hash_key_t   *k;
  int                 coords[P4EST_DIM], level;
  void              **found;

  /* TODO: 1. Check if the ghost and mirrors are already in the layer via hash table.
           If so, skip the nesessary entity. */

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
  mside = fside[0]->is_ghost == 0 ? fside[0] : fside[1]; /*< mirror side*/
  SC3A_CHECK(gside->is_ghost == 1 && mside->is_ghost == 0);

  /************* GHOST **************/

  /** Add ghost to \c ghost->ghosts array */
  /** How to know at what location of the array to place the quadrant?
   * - Simple and memory efficient solution is add any ghost to the array and
   *   sort it at the end. It costs additional O(n log n) operations in the
   *   worst case. Which is alright because the Iterator is O(n log n) anyway.
   * TODO: Double check the exact complexity of the Iterator.
  */

  /* Convert p3 quad to p2 quad */
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

  /* Check if considered ghost is unique */
  k = (ghost_hash_key_t *) sc_mempool_alloc (d->hdata->ckeys);
  k->qid = global_qid - fi->p3->goffset[proc_owner];
  k->owner = proc_owner;
  if (sc_hash_insert_unique (d->hdata->chash, k, &found)) {
    /* The key is newly linked into the hash table: count it */
    P4EST_ASSERT (*found == k);
    P4EST_INFOF ("First time adding ghost %ld, proc %ld\n",
                 (long) k->qid, (long) k->owner);
    d->hdata->added++;

    /** Fill its \c piggy3 field */
    q.p.piggy3.which_tree = gside->ntree;
    q.p.piggy3.local_num = k->qid;

    /* Push back to ghosts array */
    *(p4est_quadrant_t *) sc_array_push(&(ghost->ghosts)) = q;

    /** Contribute to a structure tracking \c tree_offsets */
    (ghost->tree_offsets[gside->ntree + 1])++;

    /** Contribute to a structure tracking \c proc_offsets */
    (ghost->proc_offsets[proc_owner + 1])++;
  }
  else {
    /* The key for this ghost had already been stored earlier */
    P4EST_ASSERT (*found != k);
    sc_mempool_free (d->hdata->ckeys, k);
    d->hdata->duped++;
  }


  /************* MIRROR **************/

  /* Convert p3 quad to p2 quad */
  SC3E (p4est3_quadrant_coordinates (fi->p3->qvt, mside->quadrant, coords));
  SC3E (p4est3_quadrant_level (fi->p3->qvt, mside->quadrant, &q.level));
  q.x = coords[0];
  q.y = coords[1];
#ifdef P4_TO_P8
  q.z = coords[2];
#endif

  global_qid = (p4est3_gloidx) mside->nquad + fi->p3->gtroffset[mside->ntree];
  /** Fill its \c piggy3 field */
  q.p.piggy3.which_tree = mside->ntree;
  q.p.piggy3.local_num = global_qid - fi->p3->goffset[fi->p3->mpirank];

  /* Push back to mirrors array */
  *(p4est_quadrant_t *) sc_array_push(&(ghost->mirrors)) = q;

  /** Contribute to a structure tracking \c mirror_tree_offsets */
  (ghost->mirror_tree_offsets[mside->ntree + 1])++;

  /** Contribute to a structure tracking \c mirror_proc_offsets */
  (ghost->mirror_proc_offsets[proc_owner + 1])++;

  /**  */
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
  ghost_hash_data_t  shdata, *hdata = &shdata;

  /* hash table for ghosts checking */
  hdata->ckeys = sc_mempool_new (sizeof (ghost_hash_key_t));
  hdata->chash = sc_hash_new (ghost_hash_fn, ghost_equal_fn, hdata, NULL);
  hdata->added = hdata->duped = 0;

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
  for (i = 1; i < ghost->num_trees + 1; i++) {
    (*(p4est_locidx_t *) sc_array_index (ghost->tree_offsets, i)) +=
      (*(p4est_locidx_t *) sc_array_index (ghost->tree_offsets, i - 1));
  }

  /** Accumulate \c proc_offsets */
  for (i = 1; i < ghost->mpisize + 1; i++) {
    (*(p4est_locidx_t *) sc_array_index (ghost->proc_offsets, i)) +=
      (*(p4est_locidx_t *) sc_array_index (ghost->proc_offsets, i - 1));
  }

  /** Sort \c ghosts */
  /** Sort \c mirrors. */

  /* clean up memory */
  sc_hash_destroy (hdata->chash);
  sc_mempool_destroy (hdata->ckeys);
  P4EST_PRODUCTINF ("Added %ld ghosts, duplicates %ld\n",
                    (long) hdata->added, (long) hdata->duped);

  return NULL;
}

#ifdef __cplusplus
#if 0
{
#endif
}
#endif