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

#ifndef P4_TO_P8
#include <p4est3_ghost_p4est.h>
#include <p4est3_p4est.h>
#include <p4est_bits.h>

#else
#include <p4est3_ghost_p8est.h>
#include <p4est3_p8est.h>
#include <p8est_bits.h>

#endif

#include <p4est3_internal.h>
#include <p4est3_search.h>
#include <stdlib.h>             /* for qsort */

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
  p4est_locidx_t      qid;      /* for ghost -- local index within owner process.
                                   for mirror -- local index and a process to
                                   what it is mirrored */
  int                 proc;
  p4est_locidx_t      i;        /* for mirror -- index within the mirror */
}
ghost_hash_key_t;

/* Provide an allocator for the key data as well as the hash map itself */
typedef struct ghost_hash_data
{
  sc_mempool_t       *ckeys;    /* memory pool for allocating the hash keys */
  sc_hash_t          *chash;    /* the hash map links keys without copying */
  p4est_locidx_t      added;    /* count each quadrant just once */
  p4est_locidx_t      duped;    /* count attempts to add more than once */
}
ghost_hash_data_t;

/* Context structure for mirror index comparison */
typedef struct
{
  p4est_quadrant_t   *mirrors;
}
mirror_compare_context_t;

/* Calculate a hash function for a ghost index */
static unsigned
ghost_hash_fn (const void *v, const void *u)
{
  uint32_t            q, o, z;
  const ghost_hash_key_t *k = (ghost_hash_key_t *) v;

  P4EST_ASSERT (k != NULL);
  q = (uint32_t) k->qid;
  o = (uint32_t) k->proc;
  z = (uint32_t) 0;

  sc_hash_final (q, o, z);

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

  return (k1->qid == k2->qid && k1->proc == k2->proc);
}

static int
compare_locidx (const void *a, const void *b)
{
  p4est_locidx_t      arg1 = *(const p4est_locidx_t *) a;
  p4est_locidx_t      arg2 = *(const p4est_locidx_t *) b;

  return (arg1 > arg2) - (arg1 < arg2);
}

typedef struct p4est3_ghost_fill_data
{
  p4est_ghost_t      *ghost;
  ghost_hash_data_t  *ghost_hdata;
  ghost_hash_data_t  *mirror_hdata;
  sc_array_t        **p2m;
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
  p4est_gloidx_t      p_own, global_qid;
  ghost_hash_key_t   *k, *k_unique_p;
  p4est3_quadrant_vtable_t *qvt_standard;
  void              **found, **found_unique_p;
  int                 coords[P4EST_DIM], level;

#ifdef P4EST_ENABLE_DEBUG
  int                 is_found = 0;
#endif

  SC3E (sc3_array_get_elem_count (fi->sides, &nsides));
  SC3A_CHECK (nsides == 2 || nsides == 1);

  if (nsides == 1) {
    /* Nothing to do here. There are no ghosts on a boundary. */
    return NULL;
  }

  sc3_array_index (fi->sides, 0, &fside[0]);
  sc3_array_index (fi->sides, 1, &fside[1]);

  /* Check both sides are initialized */
  SC3A_CHECK (fside[0]->is_ghost != -1 && fside[1]->is_ghost != -1);

  /* Check at least one side is not ghost */
  SC3A_CHECK (fside[0]->is_ghost == 0 || fside[1]->is_ghost == 0);

  if (fside[0]->is_ghost == 0 && fside[1]->is_ghost == 0) {
    /* No ghosts here. Nothing to do here. */
    return NULL;
  }

  gside = fside[0]->is_ghost == 1 ? fside[0] : fside[1];        /*< ghost side */
  mside = fside[0]->is_ghost == 0 ? fside[0] : fside[1];        /*< mirror side */
  SC3A_CHECK (gside->is_ghost == 1 && mside->is_ghost == 0);

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
  SC3E (p4est3_quadrant_level (fi->p3->qvt, gside->quadrant, &level));

  SC3E (p4est3_quadrant_vtable_p4est (&qvt_standard));
  SC3E (p4est3_quadrant_quadrant (qvt_standard, coords, level, &q));

  /* Find ghost proc owner */
  p_own = fi->p3->mpirank;
  global_qid = (p4est3_gloidx) gside->nquad + fi->p3->gtroffset[gside->ntree];
  /** TODO: We just did it in Iterator to fill callback data.
   *        Now we do it again here. Think on a way to optimize it. */
  SC3E (p4est3_search_lower_bound64
        (global_qid, fi->p3->goffset, fi->p3->mpisize + 1, &p_own));
  if (fi->p3->goffset[p_own] > global_qid) {
    SC3A_CHECK (p_own > 0);
    p_own--;
  }

  /* Check if considered ghost is unique */
  k = (ghost_hash_key_t *) sc_mempool_alloc (d->ghost_hdata->ckeys);
  k->qid = (p4est_locidx_t) (global_qid - fi->p3->goffset[p_own]);
  k->proc = (int) p_own;
  if (sc_hash_insert_unique (d->ghost_hdata->chash, k, &found)) {
    /* The key is newly linked into the hash table: count it */
    P4EST_ASSERT (*found == k);
    //P4EST_INFOF ("First time adding ghost %ld, proc %d\n",
    //             (long) k->qid, k->proc);
    d->ghost_hdata->added++;

    /** Fill its \c piggy3 field */
    q.p.piggy3.which_tree = gside->ntree;
    q.p.piggy3.local_num = k->qid;

    /* Push back to ghosts array */
    *(p4est_quadrant_t *) sc_array_push (&(ghost->ghosts)) = q;

    /** Contribute to a structure tracking \c tree_offsets */
    (ghost->tree_offsets[gside->ntree + 1])++;

    /** Contribute to a structure tracking \c proc_offsets */
    (ghost->proc_offsets[p_own + 1])++;
  }
  else {
    /* The key for this ghost had already been stored earlier */
    P4EST_ASSERT (*found != k);
    sc_mempool_free (d->ghost_hdata->ckeys, k);
    d->ghost_hdata->duped++;
  }

  /************* MIRROR **************/

  /* Convert p3 quad to p2 quad */
  SC3E (p4est3_quadrant_coordinates (fi->p3->qvt, mside->quadrant, coords));
  SC3E (p4est3_quadrant_level (fi->p3->qvt, mside->quadrant, &level));

  SC3E (p4est3_quadrant_quadrant (qvt_standard, coords, level, &q));

  global_qid = (p4est3_gloidx) mside->nquad + fi->p3->gtroffset[mside->ntree];

  /* Check if considered mirror is unique */
  k = (ghost_hash_key_t *) sc_mempool_alloc (d->mirror_hdata->ckeys);
  k->qid = (p4est_locidx_t) (global_qid - fi->p3->goffset[fi->p3->mpirank]);
  k->proc = fi->p3->mpirank;

  k_unique_p = (ghost_hash_key_t *) sc_mempool_alloc (d->mirror_hdata->ckeys);
  *k_unique_p = *k;
  k_unique_p->proc = (int) p_own;

  if (sc_hash_insert_unique (d->mirror_hdata->chash, k, &found)) {
    /* The key is newly linked into the hash table: count it */
    P4EST_ASSERT (*found == k);
    //P4EST_INFOF ("First time adding mirror %ld, proc %d\n",
    //             (long) k->qid, k->proc);
    (*(ghost_hash_key_t **) found)->i = d->mirror_hdata->added++;
#ifdef P4EST_ENABLE_DEBUG
    is_found = 1;
#endif

    /** Fill its \c piggy3 field */
    q.p.piggy3.which_tree = mside->ntree;
    q.p.piggy3.local_num = k->qid;

    /* Push back to mirrors array */
    *(p4est_quadrant_t *) sc_array_push (&(ghost->mirrors)) = q;

    /** Contribute to a structure tracking \c mirror_tree_offsets */
    (ghost->mirror_tree_offsets[mside->ntree + 1])++;
  }
  else {
    /* The key for this mirror had already been stored earlier */
    P4EST_ASSERT (*found != k);
    sc_mempool_free (d->mirror_hdata->ckeys, k);
    d->mirror_hdata->duped++;
  }

  if (sc_hash_insert_unique
      (d->mirror_hdata->chash, k_unique_p, &found_unique_p)) {
    /* The key is newly linked into the hash table: count it */
    P4EST_ASSERT (*found_unique_p == k_unique_p);
    //P4EST_INFOF ("First time adding mirror %ld, proc %ld\n",
    //             (long) k_unique_p->qid, (long) k_unique_p->proc);
    *(p4est_locidx_t *) sc_array_push (d->p2m[p_own]) =
      (*(ghost_hash_key_t **) found)->i;
  }
  else {
    /* if we uniquely inserted k before, this case is not possible */
    P4EST_ASSERT (is_found == 0);

    /* The key for this mirror had already been stored earlier */
    P4EST_ASSERT (*found_unique_p != k_unique_p);
    sc_mempool_free (d->mirror_hdata->ckeys, k_unique_p);
    d->mirror_hdata->duped++;
  }

  return NULL;
}

/* Sort ghost quadrants contained in ghost->ghosts */
static void
sort_ghost_quadrants (sc_array_t *ghosts)
{
  qsort (ghosts->array, ghosts->elem_count, ghosts->elem_size,
         p4est_quadrant_compare_piggy);
}

/* Compare function for mirror indices using context */
static int
compare_mirror_indices_context (const void *a, const void *b, void *ctx)
{
  p4est_locidx_t      idx1 = *(const p4est_locidx_t *) a;
  p4est_locidx_t      idx2 = *(const p4est_locidx_t *) b;
  mirror_compare_context_t *context = (mirror_compare_context_t *) ctx;
  return p4est_quadrant_compare_piggy (&context->mirrors[idx1],
                                       &context->mirrors[idx2]);
}

/* Quicksort implementation that allows context to be passed to comparison function */
static void
qsort_with_context (void *base, size_t nmemb, size_t size,
                    int (*compar) (const void *, const void *, void *),
                    void *context)
{
  char               *i, *j;
  char               *left, *right, *pivot_ptr;
  char               *pivot_val, *tmp, *inner_tmp;
  size_t              right_elements, left_elements;

  if (nmemb <= 1) {
    return;
  }

  /* Allocate memory for pivot and temp buffer */
  pivot_val = P4EST_ALLOC(char, size);
  tmp = P4EST_ALLOC(char, size);

  /* Simple insertion sort for small arrays */
  if (nmemb <= 16) {
    inner_tmp = P4EST_ALLOC(char, size);
    for (i = (char *) base + size; i < (char *) base + nmemb * size;
         i += size) {
      memcpy (tmp, i, size);
      for (j = i - size; j >= (char *) base &&
           compar (j, tmp, context) > 0; j -= size) {
        memcpy (j + size, j, size);
      }
      memcpy (j + size, tmp, size);
    }
    P4EST_FREE(inner_tmp);
    P4EST_FREE(pivot_val);
    P4EST_FREE(tmp);
    return;
  }

  /* Quicksort for larger arrays */
  left = (char *) base;
  right = (char *) base + (nmemb - 1) * size;

  /* Choose pivot (middle element) and copy its value */
  pivot_ptr = (char *) base + (nmemb / 2) * size;
  memcpy (pivot_val, pivot_ptr, size);

  while (left <= right) {
    /* Use pivot_val instead of pivot pointer */
    while (left < (char *) base + nmemb * size
           && compar (left, pivot_val, context) < 0)
      left += size;
    while (right >= (char *) base && compar (right, pivot_val, context) > 0)
      right -= size;

    if (left <= right) {
      /* Swap left and right */
      if (left != right) {
        memcpy (tmp, left, size);
        memcpy (left, right, size);
        memcpy (right, tmp, size);
      }
      left += size;
      right -= size;
    }
  }

  /* Recursively sort sub-arrays */
  right_elements = (right - (char *) base) / size + 1;
  if (right_elements > 0) {
    qsort_with_context (base, right_elements, size, compar, context);
  }

  left_elements = nmemb - (left - (char *) base) / size;
  if (left_elements > 0) {
    qsort_with_context (left, left_elements, size, compar, context);
  }

  /* Free allocated memory */
  P4EST_FREE(pivot_val);
  P4EST_FREE(tmp);
}

/* Sort the mirrors in a ghost layer and update the p2m arrays accordingly */
static void
sort_mirror_quadrants (p4est3_ghost_fill_data_t *d)
{
  int                 ii;
  void               *sorted;
  size_t              i, j, count;
  sc_array_t         *arr, *mirrors = &d->ghost->mirrors;
  size_t              nmirrors = mirrors->elem_count;
  p4est_locidx_t     *perm = P4EST_ALLOC (p4est_locidx_t, nmirrors);
  p4est_locidx_t     *inv = P4EST_ALLOC (p4est_locidx_t, nmirrors);
  p4est_locidx_t     *index_ptr, old_idx;
  mirror_compare_context_t context;

  if (nmirrors > 0) {

    /* Initialize permutation: identity */
    for (i = 0; i < nmirrors; i++) {
      perm[i] = (p4est_locidx_t) i;
    }

    /* Set context with pointer to mirrors array */
    context.mirrors = (p4est_quadrant_t *) mirrors->array;

    /* Sort permutation array using our context-based comparison */
    qsort_with_context (perm, nmirrors, sizeof (p4est_locidx_t),
                        compare_mirror_indices_context, &context);

    /* Build inverse permutation: inverse[old_index] = new_index */
    for (i = 0; i < nmirrors; i++) {
      inv[perm[i]] = (p4est_locidx_t) i;
    }

    /* Create a new array for sorted mirrors */
    sorted = P4EST_ALLOC (p4est_quadrant_t, nmirrors);
    for (i = 0; i < nmirrors; i++) {
      /* Copy mirror at old index perm[i] into sorted[i] */
      memcpy ((char *) sorted + i * mirrors->elem_size,
              (char *) mirrors->array +
              perm[i] * mirrors->elem_size, mirrors->elem_size);
    }

    /* Replace mirrors array with sorted order */
    P4EST_FREE (mirrors->array);
    mirrors->array = (char *) sorted;

    /* Update each sc_array in d->p2m */
    for (ii = 0; ii < d->ghost->mpisize; ii++) {
      arr = d->p2m[ii];
      count = arr->elem_count;
      for (j = 0; j < count; j++) {
        index_ptr = (p4est_locidx_t *) sc_array_index (arr, j);
        old_idx = *index_ptr;
        /* Map old index to new */
        *index_ptr = inv[old_idx];
      }
      /* sort newly mapped elements with qsort */
      qsort (arr->array, count, arr->elem_size, compare_locidx);
    }
  }
  P4EST_FREE (perm);
  P4EST_FREE (inv);
}

/* Merge the per-processor arrays of mirror indices into a single array */
static void
merge_mirror_proc_arrays (p4est3_t *p3, p4est_ghost_t *ghost,
                          sc_array_t **p2m)
{
  int                 i;
  p4est_locidx_t      total_mirrors = 0;
  p4est_locidx_t      offset = 0;
  size_t              j, count;

  /* Calculate total size needed for the merged array */
  for (i = 0; i < p3->mpisize; i++) {
    total_mirrors += (p4est_locidx_t) p2m[i]->elem_count;
    ghost->mirror_proc_offsets[i + 1] =
      ghost->mirror_proc_offsets[i] + (p4est_locidx_t) p2m[i]->elem_count;
  }

  /* Allocate memory for the merged array */
  ghost->mirror_proc_mirrors = P4EST_ALLOC (p4est_locidx_t, total_mirrors);

  /* Copy data from p2m arrays to mirror_proc_mirrors */
  offset = 0;
  for (i = 0; i < p3->mpisize; i++) {
    count = p2m[i]->elem_count;

    /* Copy indices from this p2m array */
    for (j = 0; j < count; j++) {
      ghost->mirror_proc_mirrors[offset + j] =
        *((p4est_locidx_t *) sc_array_index (p2m[i], j));
    }

    offset += (p4est_locidx_t) count;
  }
}

sc3_error_t        *
p4est3_ghost_fill_p4est (p4est3_t *p3, p4est_ghost_t **ptr_ghost)
{
  p4est3_ghost_fill_data_t data, *d = &data;
  int                 i;
  ghost_hash_data_t   sghost_hdata, *ghost_hdata = &sghost_hdata;
  ghost_hash_data_t   smirror_hdata, *mirror_hdata = &smirror_hdata;
  sc_array_t        **p2m;
  p4est_ghost_t      *ghost;

  /*--------------------------------------------------------------*/
  /************************ ALLOCATIONS ***************************/
  /*--------------------------------------------------------------*/
  ghost = P4EST_ALLOC (p4est_ghost_t, 1);
  sc_array_init (&(ghost)->ghosts, sizeof (p4est_quadrant_t));
  (ghost)->tree_offsets = P4EST_ALLOC (p4est_locidx_t, p3->num_trees + 1);
  (ghost)->proc_offsets = P4EST_ALLOC (p4est_locidx_t, p3->mpisize + 1);

  sc_array_init (&(ghost)->mirrors, sizeof (p4est_quadrant_t));
  (ghost)->mirror_tree_offsets =
    P4EST_ALLOC (p4est_locidx_t, p3->num_trees + 1);
  (ghost)->mirror_proc_mirrors = NULL;
  (ghost)->mirror_proc_offsets =
    P4EST_ALLOC (p4est_locidx_t, p3->mpisize + 1);

  /* hash table for ghosts checking */
  ghost_hdata->ckeys = sc_mempool_new (sizeof (ghost_hash_key_t));
  ghost_hdata->chash =
    sc_hash_new (ghost_hash_fn, ghost_equal_fn, ghost_hdata, NULL);
  ghost_hdata->added = ghost_hdata->duped = 0;

  /* hash table for mirrors checking */
  mirror_hdata->ckeys = sc_mempool_new (sizeof (ghost_hash_key_t));
  mirror_hdata->chash =
    sc_hash_new (ghost_hash_fn, ghost_equal_fn, mirror_hdata, NULL);
  mirror_hdata->added = mirror_hdata->duped = 0;

  /* c-array of sc_array_t * to store mirrors in a proc, that form
     mirrors_proc_mirrors later */
  p2m = P4EST_ALLOC (sc_array_t *, p3->mpisize);
  for (i = 0; i < p3->mpisize; i++) {
    p2m[i] = sc_array_new (sizeof (p4est_locidx_t));
  }

  d->ghost = ghost;
  d->ghost_hdata = ghost_hdata;
  d->mirror_hdata = mirror_hdata;
  d->p2m = p2m;

  ghost->mpisize = p3->mpisize;
  ghost->num_trees = p3->num_trees;
  ghost->btype = P4EST_CONNECT_FACE;

  /* Might be NULL for integration with Dune */
  ghost->mirror_proc_fronts = NULL;
  ghost->mirror_proc_front_offsets = NULL;

  /* Initialize tree_offsets and proc_offsets */
  memset (ghost->tree_offsets, 0,
          (ghost->num_trees + 1) * sizeof (p4est_locidx_t));

  memset (ghost->proc_offsets, 0,
          (ghost->mpisize + 1) * sizeof (p4est_locidx_t));

  memset (ghost->mirror_tree_offsets, 0,
          (ghost->num_trees + 1) * sizeof (p4est_locidx_t));

  memset (ghost->mirror_proc_offsets, 0,
          (ghost->mpisize + 1) * sizeof (p4est_locidx_t));

  /*--------------------------------------------------------------*/
  /************************** ITERATE *****************************/
  /*--------------------------------------------------------------*/

  SC3E (p4est3_iterate_face (p3, NULL, p4est3_ghost_fill_callback, d));

  /* Clean up hash tables */
  sc_hash_destroy (ghost_hdata->chash);
  sc_hash_destroy (mirror_hdata->chash);
  sc_mempool_destroy (ghost_hdata->ckeys);
  sc_mempool_destroy (mirror_hdata->ckeys);
  P4EST_PRODUCTIONF ("Added %ld ghosts, duplicates %ld\n",
                     (long) ghost_hdata->added, (long) ghost_hdata->duped);

  /*--------------------------------------------------------------*/
  /******************* POST-ITERATE PROCESSING ********************/
  /*--------------------------------------------------------------*/

  /** Accumulate \c tree_offsets */
  for (i = 1; i < ghost->num_trees + 1; i++) {
    ghost->tree_offsets[i] += ghost->tree_offsets[i - 1];
    ghost->mirror_tree_offsets[i] += ghost->mirror_tree_offsets[i - 1];
  }

  /** Accumulate \c proc_offsets */
  for (i = 1; i < ghost->mpisize + 1; i++) {
    ghost->proc_offsets[i] += ghost->proc_offsets[i - 1];
    ghost->mirror_proc_offsets[i] += ghost->mirror_proc_offsets[i - 1];
  }

  /** Sort \c ghosts */
  sort_ghost_quadrants (&(ghost->ghosts));

  /** Sort \c mirrors and update \c d->p2m arrays */
  sort_mirror_quadrants (d);

  /** Merge \c d->p2m arrays to \c mirror_proc_mirrors */
  merge_mirror_proc_arrays (p3, ghost, p2m);

  /* clean up memory temporary mirror_proc_mirrors sub-arrays */
  for (i = 0; i < p3->mpisize; i++) {
    /* Free this p2m array */
    sc_array_destroy (p2m[i]);
  }
  P4EST_FREE (p2m);

  SC3E_RETVAL (ptr_ghost, ghost);
  return NULL;
}

#ifdef __cplusplus
#if 0
{
#endif
}
#endif
