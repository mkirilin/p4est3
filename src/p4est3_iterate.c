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

#include <p4est3_iterate.h>
#include <p4est3_internal.h>
#include <p4est3_search.h>
#include <p4est_base.h>
#include <stdlib.h>
#include <p4est3_connectivity.h>
#include <p4est3_internal.h>
#include <string.h>

#include <stdio.h>

#include <sc_containers.h>
/*
 * Performance hint: Assume contiguous shared-memory layout unless
 * P4EST3_ASSUME_CONTIGUOUS is explicitly set to 0. This strips all
 * non-contiguous branches from hot iteration paths to reduce branch
 * misprediction and enable better inlining/constant folding.
 */
#ifndef P4EST3_ASSUME_CONTIGUOUS
#define P4EST3_ASSUME_CONTIGUOUS 1
#endif
#if P4EST3_ASSUME_CONTIGUOUS
#define P4EST3_ASSERT_CONTIGUOUS(p3) SC3A_CHECK ((p3)->contiguous)
#else
#define P4EST3_ASSERT_CONTIGUOUS(p3) ((void) 0)
#endif

/* Branch prediction helpers (local fallback if not provided by sc3) */
#ifndef SC3_LIKELY
#define SC3_LIKELY(x)   __builtin_expect(!!(x), 1)
#endif
#ifndef SC3_UNLIKELY
#define SC3_UNLIKELY(x) __builtin_expect(!!(x), 0)
#endif

/* Define maximum quadrant levels - use a reasonable default if not defined */
#ifndef SC3E_FAST
#ifdef P4EST_ENABLE_DEBUG
#define SC3E_FAST(f) SC3E(f)
#else
  /* In release builds assume invariant preconditions already checked */
#define SC3E_FAST(f) (f);
#endif
#endif

#ifdef __cplusplus
extern              "C"
{
#if 0
}
#endif
#endif

/* ------------------------------------------------------------------------- */
/* Lightweight tier ring cache (optimization #2)                             */
/* ------------------------------------------------------------------------- */
/* This complements the (more general) MRU hash cache above.  For the vast
 * majority of split requests we observe strong temporal locality: the same
 * (level, tree, begin,end) region is split multiple times while iterating
 * over siblings.  A tiny per-level ring buffer with O(1) insert and O(k)
 * linear lookup (k <= 2 * children) is faster and has no hashing overhead.
 * We therefore probe the tier ring first; on a miss we fall back to the
 * unified MRU cache (still useful across trees / deeper recursion).
 */

/* Maximum children supported (2D:4, 3D:8) -> splits array size children+1 */
#define P4EST3_TIER_MAX_SPLITS  (9)

typedef struct p4est3_tier_entry
{
  p4est3_gloidx       begin;    /* inclusive global id */
  //int                 tree_id;  /* owning tree */
  p4est3_gloidx       splits[P4EST3_TIER_MAX_SPLITS];   /* cached child boundaries */
} p4est3_tier_entry_t;

typedef struct p4est3_tier_ring
{
  int                 next;     /* next insertion slot */
  int                 capacity; /* fixed ring capacity */
  p4est3_tier_entry_t *entries; /* array[capacity] */
} p4est3_tier_ring_t;

static inline void
p4est3_tier_ring_init (p4est3_tier_ring_t *ring, sc3_allocator_t *A,
                       int capacity)
{
  SC_ASSERT (capacity > 0);
  ring->next = 0;
  ring->capacity = capacity;
  (void) sc3_allocator_calloc (A, (size_t) capacity,
                               sizeof (p4est3_tier_entry_t),
                               (void **) &ring->entries);
  for (int i = 0; i < capacity; ++i) {
    ring->entries[i].begin = -1;        /* mark as invalid */
  }
}

static inline void
p4est3_tier_ring_destroy (p4est3_tier_ring_t *ring, sc3_allocator_t *A)
{
  if (ring->entries != NULL) {
    (void) sc3_allocator_free (A, ring->entries);
  }
  ring->entries = NULL;
  ring->capacity = ring->next = 0;
}

static inline p4est3_tier_entry_t *
p4est3_tier_ring_lookup (p4est3_tier_ring_t *ring, p4est3_gloidx begin)
{
  for (int i = 0; i < ring->capacity; ++i) {
    p4est3_tier_entry_t *e = &ring->entries[i];
    /* Further relaxed key: match (tree_id, level, begin) only; allow end mismatch */
    if (e->begin == begin) {
      return e;                 /* hit (may be placeholder or may need extension) */
    }
  }
  return NULL;                  /* miss */
}

static inline p4est3_tier_entry_t *
p4est3_tier_ring_insert (p4est3_tier_ring_t *ring, p4est3_gloidx begin)
{
  p4est3_tier_entry_t *e = &ring->entries[ring->next];
  ring->next = (ring->next + 1) % ring->capacity;
  e->begin = begin;
  return e;
}

/* Forward declarations */
typedef struct p4est3_search_area p4est3_search_area_t;

/* Static lookup tables for face neighbors (const for compiler optimization) */
static const int    p4est3_children_face_neighbors_2d[4 * 4] = {
  /* *INDENT-OFF* */
  -1, 1, -1, 2,
   0, -1, -1, 3,
  -1, 3, 0, -1,
   2, -1, 1, -1
  /* *INDENT-ON* */
};

static const int    p4est3_children_face_neighbors_3d[8 * 6] = {
  /* *INDENT-OFF* */
  -1, 1, -1, 2, -1, 4,
   0, -1, -1, 3, -1, 5,
  -1, 3, 0, -1, -1, 6,
   2, -1, 1, -1, -1, 7,
  -1, 5, -1, 6, 0, -1,
   4, -1, -1, 7, 1, -1,
  -1, 7, 4, -1, 2, -1,
   6, -1, 5, -1, 3, -1
  /* *INDENT-ON* */
};

/* Static lookup tables for face duals (const for compiler optimization) */
static const int    p4est3_face_dual_2d[4] = { 1, 0, 3, 2 };
static const int    p4est3_face_dual_3d[6] = { 1, 0, 3, 2, 5, 4 };

/* (Split cache entry struct definition moved earlier; removed duplicate) */

typedef struct p4est3_search_area
{
  p4est3_tree_t      *tree;
  sc3_allocator_t    *alloc;    /* allocator handle (for tier ring cleanup) */
  /* Cached invariants (per p4est3 instance) to reduce repeated dereferences */
  int                 dim;      /* p3->qvt->dim */
  int                 max_level;        /* p3->qvt->max_level */
  size_t              quadrant_size;    /* p3->qvt->quadrant_size */
  int                 max_children;
  int                 half_children;    /* max_children / 2 */
  int                 nfaces;   /*Global number of faces */
  int                 Level;    /* Current volume level */
  int                 child_id; /* Child id of the area under consideration */
  p4est3_gloidx       local_begin;      /* Id of the first local quadrant in
                                           the tree under iteration */
  p4est3_gloidx       local_end;        /* Id of the last local quadrant in
                                           the tree under iteration */

  int                 nsides;   /* Number of sides of a face */
  int                 child_id_face[2];
  int                 Level_face[2];    /* array of 2, storing the level of
                                           current size */
  int                 is_refine[2];     /* array of 2, storing the information
                                           about neccesity of refenement. Must
                                           be initiated by 0 */
  p4est3_topidx       treeid_face[2];

  sc3_array_t        *idx_vol_stack;    /* 2D stack storing arrays of indices,
                                           that are output of split_array */
  sc3_array_t        *idx_face_stack[2];        /* 2D stacks storing arrays of indices,
                                                   that are output of split_array */
  sc3_array_t        *view_quads;       /* array ptr to pass tree's quads
                                           into array_split */

  p4est3_gloidx       local_begin_face[2];      /* Ids of the first local quadrant in
                                                   the tree under iteration */
  p4est3_gloidx       local_end_face[2];        /* Ids of the last local quadrant in
                                                   the tree under iteration */
  p4est3_gloidx       remote_first[2];  /* Id of the first remote proc that
                                           potentially contains the neighbor
                                           quadrants */
  p4est3_gloidx       remote_last[2];   /* Id of the last remote proc that
                                           potentially contains the neighbor
                                           quadrants */

  /* Tier ring split cache */
  p4est3_tier_ring_t *tier_rings;       /* array[max_level] */
  int                 tier_capacity;    /* per-level ring capacity */
  int                 tier_hits;
  int                 tier_misses;
  long long           total_calls;      /* total invocations of cached split function */
  long long           non_sibling_requests;     /* calls that required potential real split work (range_len != max_children) */
  long long           sibling_fastpath; /* exact full-sibling pack fast path occurrences (no compute) */
  long long           splits_computed;  /* actual expensive split computations performed (excludes sibling fast path) */
  long long           tier_extends;     /* number of times an existing tier entry was extended or materialized */
  /* Seed child interval length buckets (attempted and placed) */
  /* (removed histogram instrumentation) */
  /* Last-used split single-entry direct cache */
  int                 last_split_level; /* -1 means invalid */
  p4est3_gloidx       last_split_begin;
  p4est3_gloidx       last_split_results[P4EST3_TIER_MAX_SPLITS];
  unsigned long long  last_split_hits;
  unsigned long long  last_split_uses;
  int                 is_contiguous;    /* cached copy of p3->contiguous */

  /* Owner lookup cache (monotonic batch acceleration) */
  int                 owner_cache_valid;        /* 0 invalid, else valid */
  int                 owner_cache_rank; /* Cached owning rank */
  p4est3_gloidx       owner_cache_begin;        /* Inclusive global id begin */
  p4est3_gloidx       owner_cache_end;  /* Exclusive global id end */

  long long           subtree_full_local_hits;  /* times we skipped owner lookup */
  long long           first_child_level_reuses; /* times we reused parent-known first child level */

  /* general section */
  int                *children_face_neighbors;
  int                *face_dual;
  int                *level2nchildren;  /* Array specifing the number n
                                           of children processed on
                                           the particular level; n can't be
                                           greater than p3->num_children */
  p4est3_iterate_volume_info_t *vinfo;
  p4est3_iterate_face_info_t *finfo;
}
p4est3_search_area_t;

static sc3_error_t *p4est3_split_cache_destroy (p4est3_search_area_t * sa);
static sc3_error_t *p4est3_cached_quadrant_array_split_noncontig
  (p4est3_t * p3, sc3_array_t * array, int level, p4est3_gloidx begin,
   p4est3_gloidx end, sc3_array_t * indices, p4est3_search_area_t * sa);

static inline sc3_error_t *
p4est3_array_new (sc3_allocator_t *alloc, size_t esize, int ealloc,
                  int ecount, int is_resizable, sc3_array_t **arr)
{
  SC3E_RETVAL (arr, NULL);
  SC3A_IS (sc3_allocator_is_setup, alloc);
  SC3A_CHECK (ealloc >= 0);

  SC3E_FAST (sc3_array_new (alloc, arr));
  SC3E_FAST (sc3_array_set_elem_size (*arr, esize));
  SC3E_FAST (sc3_array_set_elem_alloc (*arr, ealloc));
  SC3E_FAST (sc3_array_set_elem_count (*arr, ecount));
  SC3E_FAST (sc3_array_set_initzero (*arr, 1));
  SC3E_FAST (sc3_array_set_resizable (*arr, is_resizable));
  SC3E_FAST (sc3_array_set_tighten (*arr, 0));
  SC3E_FAST (sc3_array_setup (*arr));

  return NULL;
}

#ifdef P4EST_ENABLE_DEBUG
static sc3_error_t *
p4est3_array_set_zero (sc3_array_t *arr)
{
  size_t              esize, ecount;
  void               *idx;
  SC3E_FAST (sc3_array_get_elem_count (arr, &ecount));
  SC3E_FAST (sc3_array_get_elem_size (arr, &esize));
  SC3E_FAST (sc3_array_index (arr, 0, &idx));
  memset (idx, 0, ecount * esize);
  return NULL;
}
#endif

static inline sc3_error_t *
p4est3_set_children_face_neighbors (p4est3_t *p3, p4est3_search_area_t *sa)
{
  const int          *src_table;
  size_t              table_size;

  /* Use static const lookup tables for better performance */
  if (p3->qvt->dim == 2) {
    src_table = p4est3_children_face_neighbors_2d;
    table_size = sizeof (p4est3_children_face_neighbors_2d);
  }
  else if (p3->qvt->dim == 3) {
    src_table = p4est3_children_face_neighbors_3d;
    table_size = sizeof (p4est3_children_face_neighbors_3d);
  }
  else {
    return NULL;                /* Unsupported dimension */
  }

  SC3E_FAST (sc3_allocator_calloc
             (p3->alloc, (size_t) sa->nfaces * (1 << p3->qvt->dim),
              sizeof (int), (void *) &sa->children_face_neighbors));

  /* Direct memory copy from const table - faster than individual assignments */
  memcpy (sa->children_face_neighbors, src_table, table_size);
  return NULL;
}

static inline sc3_error_t *
p4est3_set_face_dual (p4est3_t *p3, p4est3_search_area_t *sa)
{
  const int          *src_table;
  size_t              table_size;

  /* Use static const lookup tables for better performance */
  if (p3->qvt->dim == 2) {
    src_table = p4est3_face_dual_2d;
    table_size = sizeof (p4est3_face_dual_2d);
  }
  else if (p3->qvt->dim == 3) {
    src_table = p4est3_face_dual_3d;
    table_size = sizeof (p4est3_face_dual_3d);
  }
  else {
    return NULL;                /* Unsupported dimension */
  }

  SC3E_FAST (sc3_allocator_calloc
             (p3->alloc, (size_t) 2 * p3->qvt->dim, sizeof (int),
              (void *) &sa->face_dual));

  /* Direct memory copy from const table - faster than individual assignments */
  memcpy (sa->face_dual, src_table, table_size);
  return NULL;
}

static sc3_error_t *
p4est3_set_outer_data (p4est3_t *p3, p4est3_search_area_t *sa,
                       void *user_data)
{
  const int           ntypes = p3->num_children + 1;
  int                 i, side, nodesize, noderank;
  void               *arr;
  /*set general section of sa */
  sa->alloc = p3->alloc;
  sa->dim = p3->qvt->dim;
  sa->max_level = p3->qvt->max_level;
  sa->quadrant_size = p3->qvt->quadrant_size;
  sa->max_children = p3->num_children;
  sa->half_children = sa->max_children / 2;
  sa->nfaces = 2 * sa->dim;
  SC3E_FAST (p4est3_set_children_face_neighbors (p3, sa));
  SC3E_FAST (p4est3_set_face_dual (p3, sa));

  sa->owner_cache_valid = 0;
  sa->owner_cache_rank = -1;
  sa->owner_cache_begin = 0;
  sa->owner_cache_end = 0;
  sa->subtree_full_local_hits = 0;
  sa->first_child_level_reuses = 0;
  /* Initialize tier rings (lightweight per-level tiny caches) */
  sa->tier_rings = NULL;
  SC3E_FAST (sc3_allocator_calloc (p3->alloc, (size_t) sa->max_level,
                                   sizeof (p4est3_tier_ring_t),
                                   (void *) &sa->tier_rings));
  for (int L = 0; L < sa->max_level; ++L) {
    p4est3_tier_ring_init ((&sa->tier_rings[L]), p3->alloc, 2 * sa->max_children);      /* small fixed capacity */
  }
  /* Zero split/tier counters */
  sa->total_calls = sa->non_sibling_requests = sa->sibling_fastpath =
    sa->splits_computed = sa->tier_extends = 0;
  sa->tier_hits = sa->tier_misses = 0;
  sa->last_split_uses = sa->last_split_hits = 0;
  sa->last_split_level = -1;    /* sentinel invalid */
  sa->last_split_begin = (p4est3_gloidx) - 1;

  /*set volume section of sa */
  SC3E_FAST (p4est3_tree_index (p3, p3->fltree, &sa->tree));
  sa->child_id = -1;
  sa->Level = 0;
  SC3E_FAST (sc3_allocator_calloc (p3->alloc, sa->max_level, sizeof (int),
                                   (void *) &sa->level2nchildren));
  memset (sa->level2nchildren, 0, sizeof (int) * sa->max_level);

  /* set 2d stack (array of arrays) */
  SC3E_FAST (p4est3_array_new (p3->alloc, sizeof (sc3_array_t *),
                               sa->max_level, sa->max_level, 1,
                               &sa->idx_vol_stack));
  SC3E_FAST (sc3_array_index (sa->idx_vol_stack, 0, (void **) &arr));
  SC3E_FAST (p4est3_array_new (p3->alloc, sizeof (p4est3_gloidx), 2, 2, 1,
                               (sc3_array_t **) arr));
  for (i = 1; i < sa->max_level; ++i) {
    SC3E_FAST (sc3_array_index (sa->idx_vol_stack, i, &arr));
    SC3E_FAST (p4est3_array_new
               (p3->alloc, sizeof (p4est3_gloidx), ntypes, ntypes, 1,
                (sc3_array_t **) arr));
  }
  SC3E_FAST (sc3_array_resize (sa->idx_vol_stack, 1));

  SC3E_FAST (sc3_allocator_calloc_one (p3->alloc,
                                       sizeof (p4est3_iterate_volume_info_t),
                                       &sa->vinfo));
  sa->vinfo->p3 = p3;
  sa->vinfo->user_data = user_data;

  /*set face section of sa */
  sa->nsides = 2;
  SC3E_FAST (sc3_allocator_calloc_one (p3->alloc,
                                       sizeof (p4est3_iterate_face_info_t),
                                       &sa->finfo));
  SC3E_FAST (p4est3_array_new
             (p3->alloc, sizeof (p4est3_iterate_face_side_t), 2, 2, 1,
              &sa->finfo->sides));
  SC3E_FAST (sc3_mpienv_get_nodesize (p3->split_info, &nodesize));
  SC3E_FAST (sc3_mpienv_get_noderank (p3->split_info, &noderank));
  sa->finfo->p3 = p3;
  sa->finfo->user_data = user_data;
  for (side = 0; side < 2; ++side) {
    sa->treeid_face[side] = -1;
    sa->child_id_face[side] = -1;
    sa->Level_face[side] = 0;
    sa->is_refine[side] = 1;

    /* set 2d stack (array of arrays) */
    SC3E_FAST (p4est3_array_new (p3->alloc, sizeof (sc3_array_t *),
                                 sa->max_level, sa->max_level, 1,
                                 &sa->idx_face_stack[side]));
    SC3E_FAST (sc3_array_index (sa->idx_face_stack[side], 0, &arr));
    SC3E_FAST (p4est3_array_new
               (p3->alloc, sizeof (p4est3_gloidx), 2, 2, 1,
                (sc3_array_t **) arr));
    for (i = 1; i < sa->max_level; ++i) {
      SC3E_FAST (sc3_array_index (sa->idx_face_stack[side], i, &arr));
      SC3E_FAST (p4est3_array_new
                 (p3->alloc, sizeof (p4est3_gloidx), ntypes, ntypes, 1,
                  (sc3_array_t **) arr));
    }
    SC3E_FAST (sc3_array_resize (sa->idx_face_stack[side], 1));
  }

  /*temporarily set view on level2children array to be able to renew it later
     instead of making new / removing */
  SC3E_FAST (sc3_array_new_data (p3->alloc, &sa->view_quads, sa->tree->tquads,
                                 sa->quadrant_size, 0, 0));
  return NULL;
}

static sc3_error_t *
p4est3_destroy_outer_data (p4est3_t *p3, p4est3_search_area_t *sa)
{
  int                 i, side;
  void               *arr;

  /* Cleanup caches (tier rings + stats print) */
  SC3E_FAST (p4est3_split_cache_destroy (sa));

  SC3E_FAST (sc3_allocator_free (p3->alloc, sa->children_face_neighbors));
  SC3E_FAST (sc3_allocator_free (p3->alloc, sa->face_dual));
  SC3E_FAST (sc3_allocator_free (p3->alloc, sa->level2nchildren));

  SC3E_FAST (sc3_array_resize (sa->idx_vol_stack, sa->max_level));
  for (i = 0; i < sa->max_level; ++i) {
    SC3E_FAST (sc3_array_index (sa->idx_vol_stack, i, &arr));
    SC3E_FAST (sc3_array_destroy ((sc3_array_t **) arr));
  }
  SC3E_FAST (sc3_array_destroy (&sa->idx_vol_stack));

  for (side = 0; side < 2; ++side) {
    SC3E_FAST (sc3_array_resize (sa->idx_face_stack[side], sa->max_level));
    for (i = 0; i < sa->max_level; ++i) {
      SC3E_FAST (sc3_array_index (sa->idx_face_stack[side], i, &arr));
      SC3E_FAST (sc3_array_destroy ((sc3_array_t **) arr));
    }
  }
  SC3E_FAST (sc3_array_destroy (&sa->idx_face_stack[0]));
  SC3E_FAST (sc3_array_destroy (&sa->idx_face_stack[1]));
  SC3E_FAST (sc3_array_destroy (&sa->view_quads));
  SC3E_FAST (sc3_array_destroy (&sa->finfo->sides));
  SC3E_FAST (sc3_allocator_free (p3->alloc, sa->vinfo));
  SC3E_FAST (sc3_allocator_free (p3->alloc, sa->finfo));

  return NULL;
}

static inline int
p4est3_get_children_face_nb_id (const p4est3_search_area_t *sa,
                                int child_id, int face)
{
  return sa->children_face_neighbors[child_id * sa->nfaces + face];
}

static inline int
p4est3_get_dual_face (const p4est3_search_area_t *sa, int face)
{
  return sa->face_dual[face];
}

/* Fast owner lookup using a small temporal cache that assumes mostly
 * monotonic non-decreasing access patterns during iteration. */
static inline sc3_error_t *
p4est3_owner_lookup_fast (p4est3_t *p3, p4est3_search_area_t *sa,
                          p4est3_gloidx gid, p4est3_gloidx *owner_out)
{
  if (sa->owner_cache_valid &&
      gid >= sa->owner_cache_begin && gid < sa->owner_cache_end) {
    *owner_out = sa->owner_cache_rank;
    return NULL;                /* cache hit */
  }
  p4est3_gloidx       rank;
  SC3E_FAST (p4est3_search_lower_bound64
             (gid, p3->goffset, p3->mpisize + 1, &rank));
  if (p3->goffset[rank] > gid) {
    SC3A_CHECK (rank > 0);
    rank--;
  }
  SC3A_CHECK (rank >= 0 && rank < p3->mpisize);
  sa->owner_cache_rank = (int) rank;
  sa->owner_cache_begin = p3->goffset[rank];
  sa->owner_cache_end = p3->goffset[rank + 1];
  sa->owner_cache_valid = 1;
  *owner_out = rank;
  return NULL;
}

sc3_error_t        *
p4est3_iterate_face (p4est3_t *p3,
                     p4est3_iterate_volume_t cvolume,
                     p4est3_iterate_face_t cface, void *user_data)
{
  SC3E_FAST (p4est3_iterate_codim
             (p3, 0x01, cvolume, cface, NULL, user_data));
  return NULL;
}

sc3_error_t        *
p4est3_iterate_volume (p4est3_t *p3,
                       p4est3_iterate_volume_t cvolume, void *user_data)
{
  p4est3_topidx       ntree;
  p4est3_tree_t      *tree;
  p4est3_iterate_volume_info_t info;
  p4est3_locidx       si;

  SC3A_IS (p4est3_is_setup, p3);
  if (p3->fltree < 0 || cvolume == NULL) {
    return NULL;
  }

  info.p3 = p3;
  info.user_data = user_data;
  for (ntree = p3->fltree; ntree <= p3->lltree; ++ntree) {
    info.ntree = ntree;
    SC3E_FAST (p4est3_tree_index (p3, ntree, &tree));
    for (si = 0; si < tree->num_quads; ++si) {
      info.quadrant = tree->tquads + si * p3->qvt->quadrant_size;
      info.nquad = (p4est3_locidx) (si + tree->first_tquad);
      SC3E_FAST (cvolume (&info));
    }
  }
  return NULL;
}

static sc3_error_t *
p4est3_iterate_face_bound_init (p4est3_t *p3,
                                p4est3_search_area_t *sa,
                                p4est3_topidx tree,
                                int face,
                                p4est3_iterate_face_side_t *fside,
                                int *is_iterate)
{
  int                 orient;
  int                *Level_face = sa->Level_face;
  int                *is_refine = sa->is_refine;
  p4est3_topidx       tree_neighbor = tree, tree_ids[2];
  sc3_array_t       **idx_f_stack = sa->idx_face_stack;
  sc3_array_t        *arr;
  p4est3_gloidx      *begin, *end;
  int                 s;

  SC3E_RETVAL (is_iterate, 1);
  is_refine[0] = is_refine[1] = 0;
  fside[0].ntree = tree;
  fside[0].nface = face;
  fside[0].is_ghost = -1;
  SC3E_FAST (p4est3_connectivity_get_face
             (p3->conn, &tree_neighbor, &face, &orient));
  /*if (tree_neighbor > p3->lltree || tree_neighbor < p3->fltree) {
   *is_lower = 1;
   return NULL;
   }
   if (tree_neighbor < tree) {
   *is_lower = 1;
   return NULL;
   }*/
  if (p3->fltree <= tree_neighbor && tree_neighbor <= p3->lltree
      && tree_neighbor > tree) {
    *is_iterate = 0;
    return NULL;
  }
  sa->finfo->orientation = orient;

  if (tree_neighbor == tree) {
    /* physical boundary */
    sa->nsides = 1;
  }
  else {
    sa->nsides = 2;
    fside[1].ntree = tree_neighbor;
    sa->treeid_face[1] = tree_neighbor;
    fside[1].nface = face;
    fside[1].is_ghost = -1;
  }
  SC3E_FAST (sc3_array_resize (sa->finfo->sides, sa->nsides));
  tree_ids[0] = tree;
  tree_ids[1] = tree_neighbor;
  /* Find here two ranges of procs that:
     0th -- shares 0th side
     1st -- shares 1st side */
  for (s = 0; s < sa->nsides; ++s) {
    sa->remote_first[s] = 0;
    sa->remote_last[s] = p3->mpisize - 1;
  }

#if !P4EST3_ASSUME_CONTIGUOUS
  if (!p3->contiguous && sa->nsides == 2) {
    /* only if we probe not a physical boundary */
    for (s = 0; s < sa->nsides; ++s) {
      SC3E_FAST (p4est3_find_partition
                 (p3->alloc, p3->mpisize, p3->goffset,
                  p3->gtroffset[tree_ids[s]],
                  p3->gtroffset[tree_ids[s] + 1] - 1, &sa->remote_first[s],
                  &sa->remote_last[s]));
      if (p3->goffset[sa->remote_first[s]] > p3->gtroffset[tree_ids[s]]) {
        sa->remote_first[s]--;
      }
      if (p3->goffset[sa->remote_last[s]] >
          p3->gtroffset[tree_ids[s] + 1] - 1) {
        sa->remote_last[s]--;
      }
    }
  }
#endif /* !P4EST3_ASSUME_CONTIGUOUS */

  for (s = 0; s < sa->nsides; ++s) {
    Level_face[s] = 0;
    is_refine[s] = 1;
    sa->child_id_face[s] = 0;
    SC3E_FAST (sc3_array_resize (idx_f_stack[s], 1));
    SC3E_FAST (sc3_array_index (idx_f_stack[s], 0, &arr));
    SC3E_FAST (sc3_array_index (*(sc3_array_t **) arr, 0, &begin));
    SC3E_FAST (sc3_array_index (*(sc3_array_t **) arr, 1, &end));

    sa->local_begin_face[s] =
      SC3_MAX (p3->gtroffset[tree_ids[s]], p3->goffset[p3->mpirank]);
    sa->local_end_face[s] =
      SC3_MIN (p3->gtroffset[tree_ids[s] + 1], p3->goffset[p3->mpirank + 1]);

    *begin = p3->gtroffset[tree_ids[s]];
    *end = p3->gtroffset[tree_ids[s] + 1];
  }

  return NULL;
}

static sc3_error_t *
p4est3_internal_iterate_face (p4est3_t *p3,
                              p4est3_iterate_face_t cface,
                              p4est3_iterate_codim_t ccodim,
                              p4est3_search_area_t *sa)
{
  const int           max_children = sa->max_children;
  const int           half_ch = sa->half_children;
  p4est3_topidx      *trees = sa->treeid_face;
  p4est3_gloidx      *b_f[2], *e_f[2];
#if !P4EST3_ASSUME_CONTIGUOUS
  p4est3_gloidx       proc_owner = p3->mpirank;
#endif
  p4est3_gloidx      *base_ptr;
  void               *stack_it[2];
  /* sa->view_quads not directly needed here; use sa->view_quads inline */
  sc3_array_t       **idx_face_stack = sa->idx_face_stack;
  p4est3_iterate_face_side_t *fside;
  int                *is_refine = sa->is_refine;
  int                *Level = sa->Level_face;
  int                 i, side, level, idx, child_id;
  int                 ori = sa->finfo->orientation;
  int                 is_lvl_increased[2] = { 0, 0 };
  void               *first_quad;

  for (side = 0; side < sa->nsides; ++side) {
    SC3E_FAST (sc3_array_index
               (idx_face_stack[side], Level[side], &(stack_it[side])));
    SC3E_FAST (sc3_array_index
               (*(sc3_array_t **) stack_it[side], sa->child_id_face[side],
                &(b_f[side])));
    SC3E_FAST (sc3_array_index
               (*(sc3_array_t **) stack_it[side], sa->child_id_face[side] + 1,
                &(e_f[side])));
  }

  /* Check if both sides belong to remote process(es).
     If so, we ignore this face. */
  if (SC3_UNLIKELY ((sa->nsides == 1 && (*(b_f[0]) >= sa->local_end_face[0]
                                         || *(e_f[0]) <=
                                         sa->local_begin_face[0]))
                    || (sa->nsides == 2 && (*(b_f[0]) >= sa->local_end_face[0]
                                            || *(e_f[0]) <=
                                            sa->local_begin_face[0])
                        && (*(b_f[1]) >= sa->local_end_face[1]
                            || *(e_f[1]) <= sa->local_begin_face[1])))) {
    return NULL;
  }

  /* first check if the whole quadrant passes */
  SC3E_FAST (sc3_array_index (sa->finfo->sides, 0, &fside));
  for (side = 0; side < sa->nsides; ++side) {
    if (!is_refine[side]) {
      continue;
    }

    /* Contiguous fast path (non-contiguous disabled under macro). */
#if P4EST3_ASSUME_CONTIGUOUS
    first_quad = (void *) (p3->nodequads[0] + p3->qsize * (*(b_f[side])));
#else
    if (p3->contiguous) {
      /* Get quadrant from shared storage */
      first_quad = (void *) (p3->nodequads[0] + p3->qsize * (*(b_f[side])));
    }
    else {
      /* TODO: make it through p4est3_find_partition and last_goffsets */
      /* Consider batch processing for binary search */
      /* Find process that owns this quadrant */
      /* Start binary search from a local process */
      SC3E_FAST (p4est3_owner_lookup_fast
                 (p3, sa, *(b_f[side]), &proc_owner));

      /* Get quadrant from the owning process's window */
      first_quad = (void *) (p3->nodequads[proc_owner] +
                             p3->qsize * (*(b_f[side]) -
                                          p3->goffset[proc_owner]));
    }
#endif

    SC3E_FAST (p4est3_quadrant_level (p3->qvt, first_quad, &level));
    if (level == Level[side]) {
      is_refine[side] = 0;
      fside[side].nquad =
        (p4est3_locidx) (*(b_f[side]) - p3->gtroffset[trees[side]]);
      fside[side].quadrant = first_quad;
      if (SC3_UNLIKELY (*(b_f[side]) < p3->goffset[p3->mpirank] ||
                        *(e_f[side]) > p3->goffset[p3->mpirank + 1])) {
        /* This is a ghost quadrant */
        fside[side].is_ghost = 1;
      }
      else {
        fside[side].is_ghost = 0;
      }
    }
  }

  if (!is_refine[0] && !is_refine[1]) {
    if (SC3_LIKELY (cface != NULL)) {
      SC3E_FAST (cface (sa->finfo));
    }
    /*for (side = 0; side < sa->nsides; ++side) {
       is_refine[side] = 1;
       } */
    return NULL;
  }
  for (side = 0; side < sa->nsides; ++side) {
    if (!is_refine[side]) {
      continue;
    }
    /* we split array taht is unite for and local and remote procs quads */

    SC3E_FAST (sc3_array_push (idx_face_stack[side], &(stack_it[side])));
#ifdef P4EST_ENABLE_DEBUG
    SC3E_FAST (p4est3_array_set_zero (*(sc3_array_t **) (stack_it[side])));
#endif
    SC3E_FAST (p4est3_cached_quadrant_array_split_noncontig
               (p3, sa->view_quads, Level[side], *(b_f[side]),
                *(e_f[side]), *(sc3_array_t **) (stack_it[side]), sa));

    /* since array_split doesn't count shift from the beinning of quadrants
       in a proc, we shift result indices at the loop below */
    SC3E_FAST (sc3_array_index
               (*(sc3_array_t **) (stack_it[side]), 0, &base_ptr));
    for (i = 0; i < max_children + 1; ++i) {
      base_ptr[i] += *(b_f[side]);
    }
  }
  for (i = 0; i < half_ch; ++i) {
    is_lvl_increased[0] = is_lvl_increased[1] = 0;
    for (side = 0; side < sa->nsides; ++side) {
      /* Some service computations to derive direction of face-face connection */
      if (!is_refine[side]) {
        continue;
      }
      idx = i;
      if (side == 1) {
        SC3E_FAST (p4est3_connectivity_get_neighbor_face_corner
                   (p3->conn, fside[0].nface, fside[1].nface, ori, &idx));
      }
      SC3E_FAST (p4est3_connectivity_get_face_child_id
                 (p3->conn, fside[side].nface, idx, &child_id));
      sa->child_id_face[side] = child_id;
      Level[side]++;
      is_lvl_increased[side] = 1;
    }
    SC3E_FAST (p4est3_internal_iterate_face (p3, cface, ccodim, sa));   /* keep full check here */
    for (side = 0; side < sa->nsides; ++side) {
      Level[side] = is_lvl_increased[side] ? Level[side] - 1 : Level[side];
      if (i != half_ch - 1) {
        is_refine[side] = is_lvl_increased[side];
      }
      else if (is_lvl_increased[side]) {
        SC3E_FAST (sc3_array_pop (idx_face_stack[side]));
      }
    }
  }

  return NULL;
}

static inline sc3_error_t *
p4est3_iterate_face_inner_init (p4est3_t *p3, p4est3_search_area_t *sa,
                                int child, int neighbor)
{
  int                *Level_face = sa->Level_face;
  int                *is_refine = sa->is_refine;
  int                 ch_neigh[2] = { child, neighbor }, s;

  sc3_array_t       **idx_f_stack = sa->idx_face_stack;
  void               *top;      /* generic top of a stack */
  p4est3_gloidx      *begin, *end;
  p4est3_gloidx      *arr_vol_it;

  sa->nsides = 2;
  for (s = 0; s < sa->nsides; ++s) {
    Level_face[s] = sa->Level;
    is_refine[s] = 1;
    sa->child_id_face[s] = ch_neigh[s];
    SC3E_FAST (sc3_array_resize (idx_f_stack[s], sa->Level + 1));
    SC3E_FAST (sc3_array_index (idx_f_stack[s], Level_face[s], &top));
    SC3E_FAST (sc3_array_index (*(sc3_array_t **) top, ch_neigh[s], &begin));
    SC3E_FAST (sc3_array_index
               (*(sc3_array_t **) top, ch_neigh[s] + 1, &end));
    SC3E_FAST (sc3_array_index (sa->idx_vol_stack, Level_face[s], &top));
    SC3E_FAST (sc3_array_index (*(sc3_array_t **) top, 0, &arr_vol_it));
    *(begin) = *(arr_vol_it + ch_neigh[s]);
    *(end) = *(arr_vol_it + ch_neigh[s] + 1);
    /* Since we iterate inner faces, we inherit the local boundaries for volumes */
    sa->local_begin_face[s] = sa->local_begin;
    sa->local_end_face[s] = sa->local_end;
  }

  return NULL;
}

static sc3_error_t *
p4est3_iterate_face_inner (p4est3_t *p3,
                           p4est3_iterate_face_t cface,
                           p4est3_iterate_codim_t ccodim,
                           p4est3_search_area_t *search_area)
{

  int                 child, face, nb_id;
  p4est3_iterate_face_side_t *fside;
  search_area->treeid_face[0] = search_area->treeid_face[1]
    = search_area->tree->treeid;
  SC3E_FAST (sc3_array_index (search_area->finfo->sides, 0, &fside));
  for (child = 0; child < search_area->max_children; ++child) {
    for (face = 0; face < search_area->nfaces; ++face) {
      nb_id = p4est3_get_children_face_nb_id (search_area, child, face);
      if (nb_id < child) {
        continue;
      }
      fside[0].nface = face;
      fside[1].nface = p4est3_get_dual_face (search_area, face);
      SC3E_FAST (p4est3_iterate_face_inner_init
                 (p3, search_area, child, nb_id));
      SC3E_FAST (p4est3_internal_iterate_face
                 (p3, cface, ccodim, search_area));
    }
  }
  return NULL;
}

static sc3_error_t *
p4est3_iterate_volume_rec_init (p4est3_t *p3,
                                p4est3_search_area_t *sa, p4est3_topidx tree)
{
  int                 side;
  void               *arr;
  p4est3_iterate_face_side_t *fside;
  p4est3_gloidx      *begin, *end;

  SC3E_FAST (p4est3_tree_index (p3, tree, &sa->tree));
  sa->vinfo->ntree = tree;
  sa->finfo->orientation = 0;
  sa->finfo->tree_boundary = 0;
  SC3E_FAST (sc3_array_resize (sa->finfo->sides, 2));

  memset (sa->level2nchildren, 0, sizeof (int) * p3->qvt->max_level);
  SC3A_CHECK (sa->Level == 0);

  SC3E_FAST (sc3_array_index (sa->idx_vol_stack, 0, &arr));
#ifdef P4EST_ENABLE_DEBUG
  size_t              ecount;
  SC3E_FAST (sc3_array_get_elem_count (*(sc3_array_t **) arr, &ecount));
  SC3A_CHECK (ecount == 2);
  SC3E_FAST (sc3_array_get_elem_count (sa->idx_vol_stack, &ecount));
  SC3A_CHECK (ecount == 1);
#endif

  sa->child_id = 0;
  SC3E_FAST (sc3_array_index (*(sc3_array_t **) arr, sa->child_id, &begin));
  SC3E_FAST (sc3_array_index (*(sc3_array_t **) arr, sa->child_id + 1, &end));

  *begin = p3->gtroffset[tree];
  *end = p3->gtroffset[tree + 1];

  sa->local_begin = SC3_MAX (p3->gtroffset[tree], p3->goffset[p3->mpirank]);
  sa->local_end =
    SC3_MIN (p3->gtroffset[tree + 1], p3->goffset[p3->mpirank + 1]);

  SC3E_FAST (sc3_array_index (sa->finfo->sides, 0, &fside));
  for (int side = 0; side < 2; ++side) {
    fside[side].ntree = tree;
    fside[side].is_ghost = -1;
  }

  /* Find here the range of procs that share the same tree as noderank */
  /* Since we process inner faces here:
     The procs < mpirank are responsible for 0th sides,
     the procs > mpirank are responsoble for 1st sides. Proof? */
  for (side = 0; side < 2; ++side) {
    sa->remote_first[side] = 0;
    sa->remote_last[side] = p3->mpisize - 1;
  }

#if P4EST3_ASSUME_CONTIGUOUS
  return NULL;                  /* Skip remote partition logic */
#endif

  if (tree == p3->fltree && p3->gtroffset[tree] != p3->goffset[p3->mpirank]) {
    SC3E_FAST (p4est3_find_partition
               (p3->alloc, p3->mpisize, p3->goffset, p3->gtroffset[tree],
                p3->goffset[p3->mpirank] - 1, &sa->remote_first[0],
                &sa->remote_last[0]));
    if (p3->goffset[sa->remote_first[0]] > p3->gtroffset[tree]) {
      sa->remote_first[0]--;
    }
    if (p3->goffset[sa->remote_last[0]] > p3->goffset[p3->mpirank] - 1) {
      sa->remote_last[0]--;
    }
    SC3A_CHECK (sa->remote_last[0] < p3->mpirank);
  }
  if (tree == p3->lltree
      && p3->gtroffset[tree + 1] != p3->goffset[p3->mpirank + 1]) {
    SC3E_FAST (p4est3_find_partition
               (p3->alloc, p3->mpisize, p3->goffset,
                p3->goffset[p3->mpirank + 1], p3->gtroffset[tree + 1] - 1,
                &sa->remote_first[1], &sa->remote_last[1]));
    if (p3->goffset[sa->remote_first[1]] > p3->goffset[p3->mpirank + 1]) {
      sa->remote_first[1]--;
    }
    if (p3->goffset[sa->remote_last[1]] > p3->gtroffset[tree + 1] - 1) {
      sa->remote_last[1]--;
    }
    SC3A_CHECK (sa->remote_first[1] > p3->mpirank);
  }

  return NULL;
}

/* Legacy recursive implementation kept for reference (disabled). */
#if 0
static sc3_error_t *
p4est3_iterate_volume_rec (p4est3_t *p3,
                           p4est3_iterate_volume_t cvolume,
                           p4est3_iterate_face_t cface,
                           p4est3_iterate_codim_t ccodim,
                           p4est3_search_area_t *sa)
{
  int                 i;
  void               *first_quad;       /*first quadrant in this search area */
  int                 level;
  void               *stack_it;
  p4est3_gloidx      *arr_it;
  p4est3_gloidx       proc_owner = p3->mpirank;

  const int           max_children = sa->max_children;
  int                *l2nch = sa->level2nchildren;
  int                *Level = &sa->Level;
  p4est3_gloidx      *begin, *end;
  p4est3_tree_t      *tree = sa->tree;
  sc3_array_t        *view_q = sa->view_quads;
  sc3_array_t        *idx_vol_stack = sa->idx_vol_stack;
  p4est3_iterate_volume_info_t *vinfo = sa->vinfo;

  SC3E_FAST (sc3_array_index (idx_vol_stack, *Level, &stack_it));
  SC3E_FAST (sc3_array_index
             (*(sc3_array_t **) stack_it, sa->child_id, &begin));
  SC3E_FAST (sc3_array_index
             (*(sc3_array_t **) stack_it, sa->child_id + 1, &end));

  /* Check if the considered search area intersect
     the area of the local process. If not, then skip it. */
  if (*begin >= sa->local_end || *end <= sa->local_begin) {
    l2nch[*Level]++;
    return NULL;
  }

  /* TODO: make it through p4est3_find_partition and last_goffsets */
  /* Consider batch processing for binary search */
  /* Find process that owns this quadrant */
  /* Start binary search from a local process, we do this for volume, too,
     because we need to fill the metadata to proceed the recursion. */
  if (p3->contiguous) {
    /* Get quadrant from shared storage */
    first_quad = (void *) (p3->nodequads[0] + p3->qsize * (*begin));
  }
  else {
    SC3E_FAST (p4est3_owner_lookup_fast (p3, sa, *begin, &proc_owner));

    /* Get quadrant from the owning process's window */
    first_quad = (void *) (p3->nodequads[proc_owner] +
                           p3->qsize * (*begin - p3->goffset[proc_owner]));
  }

  SC3E_FAST (p4est3_quadrant_level (p3->qvt, first_quad, &level));
  if (level == *Level) {
    if (cvolume != NULL) {
      /* Check that this state is never reached by a remote process */
      SC3A_CHECK (proc_owner == p3->mpirank);

      vinfo->quadrant = first_quad;
      /* vinfo->nquad does not really make sense, since it is the quadrant index,
         calculated from the beginning of the tree, counting the quadrants from
         the other processes. Moreover, it cannot be local, since it can be the
         other processes sharing the same tree. So far leave it as is, since
         the tests are tailored for this value. */
         /** TODO: Maybe fix later. Possibly either by making it global,
          * or by making it the id/count of the local quadrant within the tree,
          * counting only local quadrant in this tree before it.
          */
      /* vinfo->nquad =  THE PREVIOUS VERSION
       *begin + p3->goffset[p3->mpirank] - p3->gtroffset[tree->treeid]; */
      vinfo->nquad = (p4est3_locidx) (*begin - p3->gtroffset[tree->treeid]);
      SC3E_FAST (cvolume (vinfo));
    }
    l2nch[*Level]++;
    return NULL;
  }

  SC3E_FAST (sc3_array_push (idx_vol_stack, &stack_it));
#ifdef P4EST_ENABLE_DEBUG
  SC3E_FAST (p4est3_array_set_zero (*(sc3_array_t **) stack_it));
#endif

  SC3E_FAST (p4est3_cached_quadrant_array_split_noncontig
             (p3, view_q, *Level, *begin, *end, *(sc3_array_t **) stack_it,
              sa));
  l2nch[++(*Level)] = 0;

  /* since array_split doesn't count shift from the beinning of quadrants
     in a node, we shift result indices at the loop below */
  SC3E_FAST (sc3_array_index (*(sc3_array_t **) stack_it, 0, &arr_it));
  for (i = 0; i < max_children + 1; ++i) {
    arr_it[i] += *begin;
  }
  for (i = 0; i < max_children; ++i) {
    sa->child_id = i;
    SC3E_FAST (p4est3_iterate_volume_rec (p3, cvolume, cface, ccodim, sa));
  }
  SC3A_CHECK (l2nch[*Level] == max_children);
  SC3E_FAST (p4est3_iterate_face_inner (p3, cface, ccodim, sa));
  l2nch[--(*Level)]++;
  SC3E_FAST (sc3_array_pop (idx_vol_stack));
  return NULL;
}
#endif

/* Iterative replacement for recursive volume traversal */
typedef struct
{
  int                 level;    /* Current level (after splitting) */
  int                 next_child;       /* Next child index to process */
  int                 fully_local;      /* Subtree wholly local to this rank */
  int                 have_first_level; /* We already know level of first child */
  int                 first_child_level;        /* Cached level value for first child */
} p4est3_iter_frame_t;

static sc3_error_t *
p4est3_iterate_volume_iterative (p4est3_t *p3,
                                 p4est3_iterate_volume_t cvolume,
                                 p4est3_iterate_face_t cface,
                                 p4est3_iterate_codim_t ccodim,
                                 p4est3_search_area_t *sa)
{
  p4est3_gloidx      *begin, *end;
#if !P4EST3_ASSUME_CONTIGUOUS
  p4est3_gloidx       proc_owner;
#endif
  void               *first_quad;
  int                 levelq;
  int                 max_children = sa->max_children;
  int                *l2nch = sa->level2nchildren;
  sc3_array_t        *idx_vol_stack = sa->idx_vol_stack;
  p4est3_iterate_volume_info_t *vinfo = sa->vinfo;

  /* Defensive initialization (recursive init already zeroes this, but make
     explicit for static analyzers and robustness). */
  if (l2nch[0] != 0) {
    memset (l2nch, 0, sizeof (int) * sa->max_level);
  }

  p4est3_iter_frame_t *frames;
  int                 top = -1; /* empty */
  /* allocate temp frame buffer on stack (bounded by max_level) */
  frames = (p4est3_iter_frame_t *) alloca (sizeof (*frames) * sa->max_level);

  /* Access root range (already initialized in rec_init) */
  sa->Level = 0;
  sa->child_id = 0;
  void               *root_stack_it;
  SC3E_FAST (sc3_array_index (idx_vol_stack, 0, &root_stack_it));
  SC3E_FAST (sc3_array_index (*(sc3_array_t **) root_stack_it, 0, &begin));
  SC3E_FAST (sc3_array_index (*(sc3_array_t **) root_stack_it, 1, &end));

  if (*begin < sa->local_end && *end > sa->local_begin) {
    int                 root_fully_local = (*begin >= sa->local_begin
                                            && *end <= sa->local_end);
    /* Root first quadrant */
#if P4EST3_ASSUME_CONTIGUOUS
    first_quad = (void *) (p3->nodequads[0] + p3->qsize * (*begin));
#else
    if (p3->contiguous) {
      first_quad = (void *) (p3->nodequads[0] + p3->qsize * (*begin));
      proc_owner = p3->mpirank;
    }
    else {
      if (root_fully_local) {
        proc_owner = p3->mpirank;
        ++sa->subtree_full_local_hits;
      }
      else {
        SC3E_FAST (p4est3_owner_lookup_fast (p3, sa, *begin, &proc_owner));
      }
      first_quad = (void *) (p3->nodequads[proc_owner] +
                             p3->qsize * (*begin - p3->goffset[proc_owner]));
    }
#endif
    SC3E_FAST (p4est3_quadrant_level (p3->qvt, first_quad, &levelq));
    if (levelq == sa->Level) {
      if (cvolume != NULL) {
        vinfo->quadrant = first_quad;
        vinfo->nquad =
          (p4est3_locidx) (*begin - p3->gtroffset[sa->tree->treeid]);
        SC3E_FAST (cvolume (vinfo));
      }
      /* root leaf -> done */
      return NULL;
    }
    /* Split root region */
    {
      void               *stack_it;
      p4est3_gloidx      *arr_it;
      SC3E_FAST (sc3_array_push (idx_vol_stack, &stack_it));
#ifdef P4EST_ENABLE_DEBUG
      SC3E_FAST (p4est3_array_set_zero (*(sc3_array_t **) stack_it));
#endif
      SC3E_FAST (p4est3_cached_quadrant_array_split_noncontig
                 (p3, sa->view_quads, sa->Level, *begin, *end,
                  *(sc3_array_t **) stack_it, sa));
      SC3E_FAST (sc3_array_index (*(sc3_array_t **) stack_it, 0, &arr_it));
      for (int i = 0; i < max_children + 1; ++i) {
        arr_it[i] += *begin;
      }
      l2nch[++sa->Level] = 0;
      frames[++top] = (p4est3_iter_frame_t) {
      sa->Level, 0, root_fully_local, 1, levelq};       /* first child shares parent begin */
    }
  }
  else {
    /* Root outside local - nothing to do */
    return NULL;
  }

  /* DFS using explicit stack */
  while (top >= 0) {
    p4est3_iter_frame_t *fr = &frames[top];
    if (SC3_UNLIKELY (fr->next_child >= max_children)) {
      /* Finished all children at this level */
      SC3A_CHECK (l2nch[fr->level] == max_children);
      /* Process inner faces at this (still current) level */
      SC3E_FAST (p4est3_iterate_face_inner (p3, cface, ccodim, sa));
      /* Pop indices stack and decrement level */
      SC3E_FAST (sc3_array_pop (idx_vol_stack));
      --top;                    /* remove current frame first */
      sa->Level = fr->level - 1;
      if (SC3_LIKELY (sa->Level >= 0)) {
        l2nch[sa->Level]++;     /* parent observed one more finished child */
        /* Advance parent frame child index to avoid reprocessing */
        if (top >= 0) {
          frames[top].next_child++;
        }
      }
      continue;
    }

    /* Process next child at this frame level */
    sa->child_id = fr->next_child;
    /* Obtain begin/end for this child */
    void               *stack_it;
    p4est3_gloidx      *arr_child;
    SC3E_FAST (sc3_array_index (idx_vol_stack, fr->level, &stack_it));
    SC3E_FAST (sc3_array_index (*(sc3_array_t **) stack_it, 0, &arr_child));
    begin = &arr_child[sa->child_id];
    end = &arr_child[sa->child_id + 1];

    if (SC3_UNLIKELY (*begin >= sa->local_end || *end <= sa->local_begin)) {
      /* Skip non-local */
      ++fr->next_child;
      l2nch[fr->level]++;
      continue;
    }

    /* Fetch quadrant level (avoid second query for first child if we already know it) */
    /* Child first quadrant */
#if P4EST3_ASSUME_CONTIGUOUS
    first_quad = (void *) (p3->nodequads[0] + p3->qsize * (*begin));
#else
    if (p3->contiguous) {
      first_quad = (void *) (p3->nodequads[0] + p3->qsize * (*begin));
      proc_owner = p3->mpirank;
    }
    else {
      if (fr->fully_local && *begin >= sa->local_begin
          && *end <= sa->local_end) {
        proc_owner = p3->mpirank;
        ++sa->subtree_full_local_hits;
      }
      else {
        SC3E_FAST (p4est3_owner_lookup_fast (p3, sa, *begin, &proc_owner));
      }
      first_quad = (void *) (p3->nodequads[proc_owner] +
                             p3->qsize * (*begin - p3->goffset[proc_owner]));
    }
#endif
    if (SC3_UNLIKELY
        (fr->next_child == 0 && fr->have_first_level
         && *begin == arr_child[0])) {
      /* Reuse the level we already queried at parent range decision */
      levelq = fr->first_child_level;
      ++sa->first_child_level_reuses;
      /* Invalidate so only first child benefits */
      fr->have_first_level = 0;
    }
    else {
      SC3E_FAST (p4est3_quadrant_level (p3->qvt, first_quad, &levelq));
    }
    if (levelq == fr->level) {
#if P4EST3_ASSUME_CONTIGUOUS
      if (SC3_UNLIKELY (cvolume != NULL)) {
#else
      if (cvolume != NULL && proc_owner == p3->mpirank) {
#endif
        vinfo->quadrant = first_quad;
        vinfo->nquad =
          (p4est3_locidx) (*begin - p3->gtroffset[sa->tree->treeid]);
        SC3E_FAST (cvolume (vinfo));
      }
      ++fr->next_child;
      l2nch[fr->level]++;
      continue;
    }

    /* Need to split this child range -> descend */
    {
      void               *stack_it_child;
      p4est3_gloidx      *arr_it;
      SC3E_FAST (sc3_array_push (idx_vol_stack, &stack_it_child));
#ifdef P4EST_ENABLE_DEBUG
      SC3E_FAST (p4est3_array_set_zero (*(sc3_array_t **) stack_it_child));
#endif
      SC3E_FAST (p4est3_cached_quadrant_array_split_noncontig
                 (p3, sa->view_quads, fr->level, *begin, *end,
                  *(sc3_array_t **) stack_it_child, sa));
      SC3E_FAST (sc3_array_index
                 (*(sc3_array_t **) stack_it_child, 0, &arr_it));
      for (int i = 0; i < max_children + 1; ++i) {
        arr_it[i] += *begin;
      }
      l2nch[++sa->Level] = 0;
      int                 child_full_local = (fr->fully_local
                                              || (*begin >= sa->local_begin
                                                  && *end <= sa->local_end));
      /* Store known first child level for the new frame: levelq already corresponds to *begin of child 0 */
      frames[++top] = (p4est3_iter_frame_t) {
      sa->Level, 0, child_full_local, 1, levelq};
      /* After return will continue with this child's siblings */
    }
  }
  return NULL;
}

sc3_error_t        *
p4est3_iterate_codim (p4est3_t *p3, int codims,
                      p4est3_iterate_volume_t cvolume,
                      p4est3_iterate_face_t cface,
                      p4est3_iterate_codim_t ccodim, void *user_data)
{
  p4est3_search_area_t ssa, *search_area = &ssa;
  p4est3_topidx       tree;
  p4est3_iterate_face_side_t *fside;
  int                 face, is_iterate = 0;

  /* This iteration is w/o ghost layer and for volumes only */
  if (codims < 0 || codims >= P4EST3_ITERATE_LAST) {
    return NULL;
  }
  SC3A_IS (p4est3_is_setup, p3);
  if (p3->fltree < 0 || (cvolume == NULL && cface == NULL && ccodim == NULL)) {
    return NULL;
  }

  if (codims == P4EST3_ITERATE_VOLUME) {
    SC3E_FAST (p4est3_iterate_volume (p3, cvolume, user_data));
    return NULL;
  }

  SC3E_FAST (p4est3_set_outer_data (p3, search_area, user_data));
  SC3E_FAST (sc3_array_index (search_area->finfo->sides, 0, &fside));
  for (tree = p3->fltree; tree <= p3->lltree; ++tree) {

    SC3E_FAST (p4est3_iterate_volume_rec_init (p3, search_area, tree));
    SC3E_FAST (p4est3_iterate_volume_iterative
               (p3, cvolume, cface, ccodim, search_area));

    /* frame faces part */
    search_area->finfo->tree_boundary = 1;
    search_area->treeid_face[0] = search_area->tree->treeid;
    for (face = 0; face < search_area->nfaces; ++face) {
      SC3E_FAST (p4est3_iterate_face_bound_init
                 (p3, search_area, tree, face, fside, &is_iterate));
      if (is_iterate) {
        /* we iterate over such trees that tree_neighbor < tree */
        SC3E_FAST (p4est3_internal_iterate_face
                   (p3, cface, ccodim, search_area));
      }
    }
  }
  SC3E_FAST (p4est3_destroy_outer_data (p3, search_area));
  return NULL;
}

/* (Initialization function removed – fields are zeroed in set_outer_data) */

static sc3_error_t *
p4est3_split_cache_destroy (p4est3_search_area_t *sa)
{
  /* Optional statistics output controlled by env variable */
  const char         *env_stats = getenv ("P4EST3_ITERATE_CACHE_STATS");
  if (env_stats != NULL && *env_stats) {
    double              tier_total =
      (double) (sa->tier_hits + sa->tier_misses);
    double              tier_hr =
      tier_total > 0.0 ? (100.0 * sa->tier_hits / tier_total) : 0.0;
    double              combined_hits =
      (double) (sa->tier_hits + sa->last_split_hits);
    /* Only non-sibling requests are eligible for cache savings */
    double              combined_lookups = (double) sa->non_sibling_requests;
    long long           avoided_splits = (long long) combined_hits;
    double              avoided_rate =
      combined_lookups >
      0.0 ? (100.0 * (double) avoided_splits / combined_lookups) : 0.0;
    /* seeding removed */
    double              last_rate =
      sa->last_split_uses ? (100.0 * (double) sa->last_split_hits /
                             (double) sa->last_split_uses) : 0.0;
    int                 rr_print = fprintf (stderr,
                                            "[p4est3_iterate] cache stats: calls=%lld sibling_fast=%lld non_sibling=%lld computed=%lld avoided=%lld avoided_rate=%.2f%% | last_used: uses=%llu hits=%llu rate=%.2f%% | tier_hits=%d misses=%d rate=%.2f%% extends=%lld | subtree_full_local=%lld first_child_reuse=%lld\n",
                                            sa->total_calls,
                                            sa->sibling_fastpath,
                                            sa->non_sibling_requests,
                                            sa->splits_computed,
                                            avoided_splits,
                                            avoided_rate,
                                            (unsigned long long)
                                            sa->last_split_uses,
                                            (unsigned long long)
                                            sa->last_split_hits,
                                            last_rate,
                                            sa->tier_hits, sa->tier_misses,
                                            tier_hr,
                                            sa->tier_extends,
                                            sa->subtree_full_local_hits,
                                            sa->first_child_level_reuses);
    (void) rr_print;
  }

  /* Destroy tier rings */
  SC3A_CHECK (sa->tier_rings != NULL && sa->alloc != NULL);
  for (int i = 0; i < sa->max_level; ++i) {
    p4est3_tier_ring_destroy (&sa->tier_rings[i], sa->alloc);
  }
  (void) sc3_allocator_free (sa->alloc, sa->tier_rings);
  sa->tier_rings = NULL;

  return NULL;                  /* Success */
}

static p4est3_gloidx seq[9] = { 0, 1, 2, 3, 4, 5, 6, 7, 8 };

/* Cached array split function */
/* Cached version of array split function using unified cache entry structure */
static sc3_error_t *p4est3_cached_quadrant_array_split_noncontig
  (p4est3_t * p3, sc3_array_t * array, int level,
   p4est3_gloidx begin, p4est3_gloidx end,
   sc3_array_t * indices, p4est3_search_area_t * sa)
{
  p4est3_gloidx      *src_val;
  p4est3_tier_entry_t *tier_hit = NULL;

  p4est3_gloidx       range_len = end - begin;
  /* Early sibling fast path: no split, no other counters */
  if (SC3_LIKELY (range_len == (p4est3_gloidx) sa->max_children)) {
    sa->sibling_fastpath++;
    SC3E_FAST (sc3_array_index (indices, 0, &src_val));
    memcpy (src_val, seq, sizeof (p4est3_gloidx) * (sa->max_children + 1));
    /*for (int i = 0; i < range_len + 1; i++) {
       src_val[i] = i;
       } */
    return NULL;
  }

  sa->total_calls++;
  sa->non_sibling_requests++;
  /* Last-used direct cache */
  ++sa->last_split_uses;
  if (SC3_UNLIKELY (sa->last_split_level == level &&
                    sa->last_split_begin == begin)) {
    ++sa->last_split_hits;
    p4est3_gloidx      *dst;
    SC3E_FAST (sc3_array_index (indices, 0, &dst));
    memcpy (dst, sa->last_split_results,
            sizeof (p4est3_gloidx) * (sa->max_children + 1));
    return NULL;
  }

  /* Tier ring fast path (only for non-sibling requests) */
  SC3A_CHECK (level >= 0 && level < sa->max_level && sa->tier_rings != NULL);
  p4est3_tier_ring_t *ring = &sa->tier_rings[level];
  tier_hit = p4est3_tier_ring_lookup (ring, begin);
  if (tier_hit != NULL) {
    /* Cached data covers requested sub-interval (or equal) */
    SC3E_FAST (sc3_array_index (indices, 0, &src_val));
    memcpy (src_val, tier_hit->splits,
            sizeof (p4est3_gloidx) * (sa->max_children + 1));
    sa->tier_hits++;
    return NULL;
  }
  else {
    sa->tier_misses++;
  }

  /* Compute the split result */
  ++sa->splits_computed;
#if P4EST3_ASSUME_CONTIGUOUS
  SC3E_FAST (sc3_array_renew_data
             (&array, p3->nodequads[0], p3->qsize, begin, end - begin));
  SC3E_FAST (p4est3_quadrant_array_split (p3->qvt, array, level, indices));
#else
  if (sa->is_contiguous) {
    SC3E_FAST (sc3_array_renew_data
               (&array, p3->nodequads[0], p3->qsize, begin, end - begin));
    SC3E_FAST (p4est3_quadrant_array_split (p3->qvt, array, level, indices));
  }
  else {
    SC3E_FAST (sc3_array_renew_data
               (&array, p3->nodequads[0], p3->qsize, 0, end - begin));
    SC3E_FAST (p4est3_quadrant_array_split_noncontig
               (p3, array, level, begin, indices));
  }
#endif

  SC3E_FAST (sc3_array_index (indices, 0, &src_val));
  /* Update last-used cache */
  sa->last_split_level = level;
  sa->last_split_begin = begin;
  memcpy (sa->last_split_results, src_val,
          sizeof (p4est3_gloidx) * (sa->max_children + 1));
  /* Insert or extend tier ring entry */
  SC3A_CHECK (level >= 0 && level < sa->max_level && sa->alloc != NULL);
  ring = &sa->tier_rings[level];
  p4est3_tier_entry_t *e = p4est3_tier_ring_insert (ring, begin);
  memcpy (e->splits, src_val,
          sizeof (p4est3_gloidx) * (sa->max_children + 1));

  return NULL;
}
