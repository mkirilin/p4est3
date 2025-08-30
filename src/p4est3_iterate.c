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
#include <stdio.h>

/* Include MRU cache functionality for array split optimization */
#include <sc_containers.h>

/* Branch prediction macros (fallback) */
#ifndef SC_LIKELY
#if defined(__GNUC__) || defined(__clang__)
#define SC_LIKELY(x)   __builtin_expect(!!(x), 1)
#define SC_UNLIKELY(x) __builtin_expect(!!(x), 0)
#else
#define SC_LIKELY(x)   (x)
#define SC_UNLIKELY(x) (x)
#endif
#endif

/* Define maximum quadrant levels - use a reasonable default if not defined */
#ifndef SC3E_FAST
#ifdef P4EST_ENABLE_DEBUG
#define SC3E_FAST(f) SC3E(f)
#else
  /* In release builds assume invariant preconditions already checked */
#define SC3E_FAST(f) do { (void) (f); } while (0)
#endif
#endif

#ifndef P4EST3_ITER_CACHE_LVL
#define P4EST3_ITER_CACHE_LVL 32
#endif

#ifdef __cplusplus
extern              "C"
{
#if 0
}
#endif
#endif

/* Define drop function type since it's not available in sc */
typedef void        (*sc_drop_function_t) (void *data, const void *user);

/* Split cache entry structure (embedded in MRU nodes) */
typedef struct p4est3_split_cache_entry
{
  p4est3_gloidx       first_quad_id;
  p4est3_gloidx       last_quad_id;
  int                 level;
  int                 tree_id;
  unsigned long long  composite_key;
  union
  {
    p4est3_gloidx       split_results2d[5];
    p4est3_gloidx       split_results3d[9];
  };
  int                 reuse_count;
} p4est3_split_cache_entry_t;

/* Doubly-linked list node for MRU cache (embeds entry) */
typedef struct sc_dlink
{
  void               *data;     /* Points to &entry */
  struct sc_dlink    *next;
  struct sc_dlink    *prev;
  p4est3_split_cache_entry_t entry;     /* Embedded cache entry */
} sc_dlink_t;

/* MRU cache structure definition (adapted from p4est_dune.c) */
typedef struct sc_hash_mru
{
  /* functions provided by the user */
  sc_hash_function_t  hash_fn;
  sc_equal_function_t equal_fn;
  sc_drop_function_t  drop_fn;
  void               *user;

  /* internal container objects */
  sc_hash_t          *hash;
  sc_mempool_t       *pool;
  sc_dlink_t         *first, *last;

  /* counters and statistics */
  size_t              maxcount;
  size_t              count;
  size_t              num_inserted;
  size_t              insert_missed;
  size_t              num_removed;
  size_t              remove_missed;
}
sc_hash_mru_t;

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
  p4est3_gloidx       end;      /* exclusive global id last computed (may differ from request) */
  int                 tree_id;  /* owning tree */
  int                 level;    /* level at which split was done */
  int                 valid;    /* entry has computed splits */
  int                 seeded;   /* entry was pre-seeded (placeholder) */
  p4est3_gloidx       splits[P4EST3_TIER_MAX_SPLITS];   /* cached child boundaries */
} p4est3_tier_entry_t;

typedef struct p4est3_tier_ring
{
  int                 next;     /* next insertion slot */
  int                 size;     /* number of valid entries */
  int                 capacity; /* fixed ring capacity */
  p4est3_tier_entry_t *entries; /* array[capacity] */
} p4est3_tier_ring_t;

static inline void
p4est3_tier_ring_init (p4est3_tier_ring_t *ring, sc3_allocator_t *A,
                       int capacity)
{
  ring->next = 0;
  ring->size = 0;
  ring->capacity = capacity;
  if (capacity <= 0) {
    ring->entries = NULL;
    return;
  }
  (void) sc3_allocator_calloc (A, (size_t) capacity,
                               sizeof (p4est3_tier_entry_t),
                               (void **) &ring->entries);
}

static inline void
p4est3_tier_ring_destroy (p4est3_tier_ring_t *ring, sc3_allocator_t *A)
{
  if (ring->entries != NULL) {
    (void) sc3_allocator_free (A, ring->entries);
  }
  ring->entries = NULL;
  ring->capacity = ring->size = ring->next = 0;
}

static inline p4est3_tier_entry_t *
p4est3_tier_ring_lookup (p4est3_tier_ring_t *ring, int tree_id,
                         int level, p4est3_gloidx begin, p4est3_gloidx end)
{
  if (SC_UNLIKELY (ring->size == 0)) {
    return NULL;
  }
  /* Linear scan – ring->size is tiny (<= 16) */
  for (int i = 0; i < ring->size; ++i) {
    p4est3_tier_entry_t *e = &ring->entries[i];
    /* Further relaxed key: match (tree_id, level, begin) only; allow end mismatch */
    if (e->tree_id == tree_id && e->level == level && e->begin == begin) {
      /* If stored covers a superset (end >= requested end) treat as hit outright */
      /* If requested extends further (end > e->end) we will later expand e->end after computing */
      return e;                 /* hit (may be placeholder or may need extension) */
    }
  }
  return NULL;                  /* miss */
}

static inline p4est3_tier_entry_t *
p4est3_tier_ring_insert (p4est3_tier_ring_t *ring, int tree_id,
                         int level, p4est3_gloidx begin, p4est3_gloidx end)
{
  if (SC_UNLIKELY (ring->capacity == 0)) {
    return NULL;                /* disabled */
  }
  p4est3_tier_entry_t *e = &ring->entries[ring->next];
  ring->next = (ring->next + 1) % ring->capacity;
  if (ring->size < ring->capacity) {
    ring->size++;
  }
  e->begin = begin;
  e->end = end;
  e->tree_id = tree_id;
  e->level = level;
  e->valid = 1;
  e->seeded = 0;
  return e;
}

/* Forward declarations for MRU cache functions */
static sc_hash_mru_t *sc_hash_mru_new (sc_hash_function_t hash_fn,
                                       sc_equal_function_t equal_fn,
                                       sc_drop_function_t drop_fn,
                                       void *user, size_t maxcount);
static void         sc_hash_mru_destroy (sc_hash_mru_t * mru);
static int          sc_hash_mru_insert_unique (sc_hash_mru_t * mru,
                                               void *v, void ***found);

/* MRU cache helper functions */
static unsigned int sc_hash_mru_hash (const void *v, const void *u);
static int          sc_hash_mru_is_equal (const void *v1, const void *v2,
                                          const void *u);
static void         sc_hash_mru_consolidate (sc_hash_mru_t * mru);

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

/* Hash and comparison functions for MRU cache */
static unsigned int
p4est3_split_cache_hash (const void *entry, const void *user)
{
  const p4est3_split_cache_entry_t *cache_entry =
    (const p4est3_split_cache_entry_t *) entry;

  if (cache_entry->composite_key) {
    /* Fast hash from composite key */
    unsigned long long  x = cache_entry->composite_key;
    /* 64->32 finalizer (splitmix32 style) */
    x ^= x >> 33;
    x *= 0xff51afd7ed558ccdULL;
    x ^= x >> 33;
    x *= 0xc4ceb9fe1a85ec53ULL;
    x ^= x >> 33;
    return (unsigned int) x;
  }

  unsigned int        a, b, c;
  a = b = c = 0x9e3779b9;       // Initialize with golden ratio

  // Mix in the key fields
  a += (unsigned int) cache_entry->first_quad_id;
  b += (unsigned int) (cache_entry->first_quad_id >> 32);
  c += (unsigned int) cache_entry->last_quad_id;

  // First mixing round
  sc_hash_mix (a, b, c);

  // Add more fields
  a += (unsigned int) (cache_entry->last_quad_id >> 32);
  b += (unsigned int) cache_entry->level;
  c += (unsigned int) cache_entry->tree_id;

  // Final mixing
  sc_hash_final (a, b, c);

  return c;
}

static int
p4est3_split_cache_equal (const void *entry1, const void *entry2,
                          const void *user)
{
  const p4est3_split_cache_entry_t *e1 =
    (const p4est3_split_cache_entry_t *) entry1;
  const p4est3_split_cache_entry_t *e2 =
    (const p4est3_split_cache_entry_t *) entry2;

  if (e1->composite_key || e2->composite_key) {
    return e1->composite_key == e2->composite_key &&
      e1->first_quad_id == e2->first_quad_id &&
      e1->last_quad_id == e2->last_quad_id;
  }

  return (e1->first_quad_id == e2->first_quad_id &&
          e1->last_quad_id == e2->last_quad_id &&
          e1->level == e2->level && e1->tree_id == e2->tree_id);
}

/* MRU cache implementation functions */
static inline unsigned int
sc_hash_mru_hash (const void *v, const void *u)
{
  const sc_hash_mru_t *mru = (const sc_hash_mru_t *) u;
  const sc_dlink_t   *lynk = (const sc_dlink_t *) v;

  SC_ASSERT (mru != NULL);
  SC_ASSERT (mru->hash_fn != NULL);

  return mru->hash_fn (lynk->data, mru->user);
}

static inline int
sc_hash_mru_is_equal (const void *v1, const void *v2, const void *u)
{
  const sc_hash_mru_t *mru = (const sc_hash_mru_t *) u;
  const sc_dlink_t   *lynk1 = (const sc_dlink_t *) v1;
  const sc_dlink_t   *lynk2 = (const sc_dlink_t *) v2;

  SC_ASSERT (mru != NULL);
  SC_ASSERT (mru->equal_fn != NULL);

  return mru->equal_fn (lynk1->data, lynk2->data, mru->user);
}

static void
sc_hash_mru_consolidate (sc_hash_mru_t *mru)
{
  sc_dlink_t         *drop;

  /* verify preconditions */
  SC_ASSERT (mru != NULL);
  SC_ASSERT (mru->pool->elem_count == mru->count);
  SC_ASSERT (mru->hash->elem_count == mru->count);

  /* drop superfluous objects */
  while (mru->count > mru->maxcount) {
    drop = mru->first;
    SC_ASSERT (drop != NULL);
    SC_ASSERT (drop->prev == NULL);

    /* first remove element from hash */
    P4EST_EXECUTE_ASSERT_TRUE (sc_hash_remove (mru->hash, drop, NULL));

    /* call the user's drop handler */
    if (mru->drop_fn != NULL) {
      mru->drop_fn (drop->data, mru->user);
    }

    /* drop oldest list entry */
    mru->first = drop->next;
    if (mru->first == NULL) {
      SC_ASSERT (mru->count == 1);
      mru->last = NULL;
    }
    else {
      SC_ASSERT (drop->next->prev == drop);
      mru->first->prev = NULL;
    }

    /* update memory and count */
    sc_mempool_free (mru->pool, drop);
    --mru->count;
  }
}

static inline sc_hash_mru_t *
sc_hash_mru_new (sc_hash_function_t hash_fn, sc_equal_function_t equal_fn,
                 sc_drop_function_t drop_fn, void *user, size_t maxcount)
{
  sc_hash_mru_t      *mru;

  SC_ASSERT (hash_fn != NULL);
  SC_ASSERT (equal_fn != NULL);

  mru = SC_ALLOC_ZERO (sc_hash_mru_t, 1);
  mru->hash_fn = hash_fn;
  mru->equal_fn = equal_fn;
  mru->drop_fn = drop_fn;
  mru->user = user;

  mru->hash = sc_hash_new (sc_hash_mru_hash, sc_hash_mru_is_equal, mru, NULL);
  mru->pool = sc_mempool_new (sizeof (sc_dlink_t));

  mru->maxcount = maxcount;

  return mru;
}

static void
sc_hash_mru_destroy (sc_hash_mru_t *mru)
{
  /* verify preconditions */
  SC_ASSERT (mru != NULL);
  SC_ASSERT (mru->pool->elem_count == mru->count);
  SC_ASSERT (mru->hash->elem_count == mru->count);

  /* call drop handler on remaining items */
  if (mru->drop_fn != NULL) {
    sc_dlink_t         *head = mru->first;

    /* walk through the list from oldest to newest */
    while (head != NULL) {
      mru->drop_fn (head->data, mru->user);
      head = head->next;

      /* returning to mempool would be redundant here */
    }
  }

  /* free all stored list elements */
  sc_hash_destroy (mru->hash);

  /* free the hash structure itself */
  sc_mempool_destroy (mru->pool);

  /* free this object */
  SC_FREE (mru);
}

static int
sc_hash_mru_insert_unique (sc_hash_mru_t *mru, void *v, void ***found)
{
  int                 inserted;
  void              **lfound;
  sc_dlink_t          key, *lkey = &key;
  sc_dlink_t         *add;

  /* verify preconditions */
  SC_ASSERT (mru != NULL);
  SC_ASSERT (mru->pool->elem_count == mru->count);
  SC_ASSERT (mru->hash->elem_count == mru->count);

  /* construct hash key */
  lkey->data = v;
  inserted = sc_hash_insert_unique (mru->hash, lkey, &lfound);
  if (inserted) {

    /* this object is newly added */
    add = (sc_dlink_t *) sc_mempool_alloc (mru->pool);
    /* Copy provided key/value (currently only key fields valid) into embedded entry */
    add->entry = *(p4est3_split_cache_entry_t *) v;
    add->data = &add->entry;
    add->next = NULL;
    if (mru->last == NULL) {

      /* the list was empty before */
      SC_ASSERT (mru->first == NULL && mru->count == 0);
      (mru->first = add)->prev = NULL;
    }
    else {

      /* append to the list */
      SC_ASSERT (mru->last->next == NULL && mru->count > 0);
      (mru->last->next = add)->prev = mru->last;
    }

    /* update memory and counters */
    *(sc_dlink_t **) lfound = mru->last = add;
    ++mru->count;
    ++mru->insert_missed;
  }
  else {

    /* this object exists already */
    add = *(sc_dlink_t **) lfound;
    /* add->data already points to embedded entry */
    if (add != mru->last) {

      /* remove it from its place */
      SC_ASSERT (add->next != NULL);
      add->next->prev = add->prev;
      if (add->next->prev == NULL) {

        /* we are removing the first element */
        SC_ASSERT (add == mru->first);
        mru->first = add->next;
      }
      else {

        /* we keep the first element */
        add->prev->next = add->next;
      }

      /* and append it to the end */
      (add->prev = mru->last)->next = add;
      (mru->last = add)->next = NULL;
    }
  }
  ++mru->num_inserted;

  /* return data location if so desired */
  if (found != NULL) {
    *found = &add->data;
  }

  /* indicate pre-existing object and return */
  sc_hash_mru_consolidate (mru);
  return inserted;
}

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

  /* MRU cache for array split optimization */
  sc_hash_mru_t      *split_cache[P4EST3_ITER_CACHE_LVL];       /* Legacy per-level (unused in unified mode) */
  sc_hash_mru_t      *split_cache_unified;      /* Unified cache across all levels */
  /* cache_entry_pool removed: entries embedded in MRU nodes */
  int                 cache_max_size;   /* Maximum cache size per level */
  int                 cache_base_size;  /* Baseline size for dynamic resizing */
  int                 cache_hits;       /* Statistics: cache hits */
  int                 cache_misses;     /* Statistics: cache misses */
  int                 cache_op_counter; /* Operations since last resize check */
  int                 cache_resize_interval;    /* How often to re-evaluate size */

  /* Tier ring split cache (fast path before MRU) */
  p4est3_tier_ring_t *tier_rings;       /* array[max_level] */
  int                 tier_capacity;    /* per-level ring capacity */
  int                 tier_hits;
  int                 tier_misses;
  int                 tier_seed_total;  /* number of tier entries pre-seeded */
  int                 tier_seed_hits;   /* how many seeded entries later filled & hit */
  int                 tier_seed_attempted;      /* attempted seeds before filtering */
  int                 tier_seed_skipped;        /* children skipped due to filter */
  int                 tier_seed_evicted;        /* placeholder overwritten before materialization */
  int                 tier_placeholder_lookups; /* found placeholder during lookup */
  long long           split_total;      /* total split requests */
  long long           sibling_fastpath; /* exact full-sibling pack fast path count */
  long long           splits_computed;  /* actual expensive split computations performed */
  int                 tier_extensions;  /* times an existing tier entry was extended (end grew) */
  int                 tier_subinterval_hits;    /* hits where cached end > requested end */
  /* (removed histogram instrumentation) */
  /* Last-used split single-entry direct cache */
  int                 last_split_level; /* -1 means invalid */
  p4est3_topidx       last_split_tree;
  p4est3_gloidx       last_split_begin;
  p4est3_gloidx       last_split_end;
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

/* Cache drop function */
static inline void
p4est3_split_cache_drop (void *entry, const void *user)
{
  (void) entry;
  (void) user;                  /* no-op for embedded entries */
}

/* Forward declarations for MRU cache functions */
static inline void  p4est3_split_cache_drop (void *entry, const void *user);
static inline sc3_error_t *p4est3_split_cache_init (p4est3_t * p3,
                                                    p4est3_search_area_t * sa,
                                                    int max_cache_size);
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

  SC3E (sc3_array_new (alloc, arr));
  SC3E (sc3_array_set_elem_size (*arr, esize));
  SC3E (sc3_array_set_elem_alloc (*arr, ealloc));
  SC3E (sc3_array_set_elem_count (*arr, ecount));
  SC3E (sc3_array_set_initzero (*arr, 1));
  SC3E (sc3_array_set_resizable (*arr, is_resizable));
  SC3E (sc3_array_set_tighten (*arr, 0));
  SC3E (sc3_array_setup (*arr));

  return NULL;
}

#ifdef P4EST_ENABLE_DEBUG
static sc3_error_t *
p4est3_array_set_zero (sc3_array_t *arr)
{
  size_t              esize, ecount;
  void               *idx;
  SC3E (sc3_array_get_elem_count (arr, &ecount));
  SC3E (sc3_array_get_elem_size (arr, &esize));
  SC3E (sc3_array_index (arr, 0, &idx));
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

  SC3E (sc3_allocator_calloc
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

  SC3E (sc3_allocator_calloc
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
  SC3E (p4est3_set_children_face_neighbors (p3, sa));
  SC3E (p4est3_set_face_dual (p3, sa));

  /* Initialize MRU cache for array split optimization */
  /* Allow cache size to be configured via environment variable */
  int                 cache_size = 16;  /* Default (reduced) unified cache size */
  SC3E (p4est3_split_cache_init (p3, sa, cache_size));
  sa->cache_hits = 0;
  sa->cache_misses = 0;
  sa->owner_cache_valid = 0;
  sa->owner_cache_rank = -1;
  sa->owner_cache_begin = 0;
  sa->owner_cache_end = 0;
  sa->subtree_full_local_hits = 0;
  sa->first_child_level_reuses = 0;

  /*set volume section of sa */
  SC3E (p4est3_tree_index (p3, p3->fltree, &sa->tree));
  sa->child_id = -1;
  sa->Level = 0;
  SC3E (sc3_allocator_calloc (p3->alloc, sa->max_level, sizeof (int),
                              (void *) &sa->level2nchildren));
  memset (sa->level2nchildren, 0, sizeof (int) * sa->max_level);

  /* set 2d stack (array of arrays) */
  SC3E (p4est3_array_new (p3->alloc, sizeof (sc3_array_t *),
                          sa->max_level, sa->max_level, 1,
                          &sa->idx_vol_stack));
  SC3E (sc3_array_index (sa->idx_vol_stack, 0, (void **) &arr));
  SC3E (p4est3_array_new (p3->alloc, sizeof (p4est3_gloidx), 2, 2, 1,
                          (sc3_array_t **) arr));
  for (i = 1; i < sa->max_level; ++i) {
    SC3E (sc3_array_index (sa->idx_vol_stack, i, &arr));
    SC3E (p4est3_array_new (p3->alloc, sizeof (p4est3_gloidx), ntypes, ntypes,
                            1, (sc3_array_t **) arr));
  }
  SC3E (sc3_array_resize (sa->idx_vol_stack, 1));

  SC3E (sc3_allocator_calloc_one (p3->alloc,
                                  sizeof (p4est3_iterate_volume_info_t),
                                  &sa->vinfo));
  sa->vinfo->p3 = p3;
  sa->vinfo->user_data = user_data;

  /*set face section of sa */
  sa->nsides = 2;
  SC3E (sc3_allocator_calloc_one (p3->alloc,
                                  sizeof (p4est3_iterate_face_info_t),
                                  &sa->finfo));
  SC3E (p4est3_array_new (p3->alloc, sizeof (p4est3_iterate_face_side_t), 2,
                          2, 1, &sa->finfo->sides));
  SC3E (sc3_mpienv_get_nodesize (p3->split_info, &nodesize));
  SC3E (sc3_mpienv_get_noderank (p3->split_info, &noderank));
  sa->finfo->p3 = p3;
  sa->finfo->user_data = user_data;
  for (side = 0; side < 2; ++side) {
    sa->treeid_face[side] = -1;
    sa->child_id_face[side] = -1;
    sa->Level_face[side] = 0;
    sa->is_refine[side] = 1;

    /* set 2d stack (array of arrays) */
    SC3E (p4est3_array_new (p3->alloc, sizeof (sc3_array_t *),
                            sa->max_level, sa->max_level, 1,
                            &sa->idx_face_stack[side]));
    SC3E (sc3_array_index (sa->idx_face_stack[side], 0, &arr));
    SC3E (p4est3_array_new
          (p3->alloc, sizeof (p4est3_gloidx), 2, 2, 1, (sc3_array_t **) arr));
    for (i = 1; i < sa->max_level; ++i) {
      SC3E (sc3_array_index (sa->idx_face_stack[side], i, &arr));
      SC3E (p4est3_array_new
            (p3->alloc, sizeof (p4est3_gloidx), ntypes, ntypes, 1,
             (sc3_array_t **) arr));
    }
    SC3E (sc3_array_resize (sa->idx_face_stack[side], 1));
  }

  /*temporarily set view on level2children array to be able to renew it later
     instead of making new / removing */
  SC3E (sc3_array_new_data (p3->alloc, &sa->view_quads, sa->tree->tquads,
                            sa->quadrant_size, 0, 0));
  return NULL;
}

static sc3_error_t *
p4est3_destroy_outer_data (p4est3_t *p3, p4est3_search_area_t *sa)
{
  int                 i, side;
  void               *arr;

  /* Destroy MRU cache */
  SC3E (p4est3_split_cache_destroy (sa));

  SC3E (sc3_allocator_free (p3->alloc, sa->children_face_neighbors));
  SC3E (sc3_allocator_free (p3->alloc, sa->face_dual));
  SC3E (sc3_allocator_free (p3->alloc, sa->level2nchildren));

  SC3E (sc3_array_resize (sa->idx_vol_stack, sa->max_level));
  for (i = 0; i < sa->max_level; ++i) {
    SC3E (sc3_array_index (sa->idx_vol_stack, i, &arr));
    SC3E (sc3_array_destroy ((sc3_array_t **) arr));
  }
  SC3E (sc3_array_destroy (&sa->idx_vol_stack));

  for (side = 0; side < 2; ++side) {
    SC3E (sc3_array_resize (sa->idx_face_stack[side], sa->max_level));
    for (i = 0; i < sa->max_level; ++i) {
      SC3E (sc3_array_index (sa->idx_face_stack[side], i, &arr));
      SC3E (sc3_array_destroy ((sc3_array_t **) arr));
    }
  }
  SC3E (sc3_array_destroy (&sa->idx_face_stack[0]));
  SC3E (sc3_array_destroy (&sa->idx_face_stack[1]));
  SC3E (sc3_array_destroy (&sa->view_quads));
  SC3E (sc3_array_destroy (&sa->finfo->sides));
  SC3E (sc3_allocator_free (p3->alloc, sa->vinfo));
  SC3E (sc3_allocator_free (p3->alloc, sa->finfo));

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
  SC3E (p4est3_search_lower_bound64
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
  SC3E (p4est3_iterate_codim (p3, 0x01, cvolume, cface, NULL, user_data));
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
    SC3E (p4est3_tree_index (p3, ntree, &tree));
    for (si = 0; si < tree->num_quads; ++si) {
      info.quadrant = tree->tquads + si * p3->qvt->quadrant_size;
      info.nquad = (p4est3_locidx) (si + tree->first_tquad);
      SC3E (cvolume (&info));
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
  SC3E (p4est3_connectivity_get_face
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
  SC3E (sc3_array_resize (sa->finfo->sides, sa->nsides));
  tree_ids[0] = tree;
  tree_ids[1] = tree_neighbor;
  /* Find here two ranges of procs that:
     0th -- shares 0th side
     1st -- shares 1st side */
  for (s = 0; s < sa->nsides; ++s) {
    sa->remote_first[s] = 0;
    sa->remote_last[s] = p3->mpisize - 1;
  }

  if (!p3->contiguous && sa->nsides == 2) {
    /* only if we probe not a physical boundary */
    for (s = 0; s < sa->nsides; ++s) {
      SC3E (p4est3_find_partition
            (p3->alloc, p3->mpisize, p3->goffset, p3->gtroffset[tree_ids[s]],
             p3->gtroffset[tree_ids[s] + 1] - 1,
             &sa->remote_first[s], &sa->remote_last[s]));
      if (p3->goffset[sa->remote_first[s]] > p3->gtroffset[tree_ids[s]]) {
        sa->remote_first[s]--;
      }
      if (p3->goffset[sa->remote_last[s]] >
          p3->gtroffset[tree_ids[s] + 1] - 1) {
        sa->remote_last[s]--;
      }
    }
  }

  for (s = 0; s < sa->nsides; ++s) {
    Level_face[s] = 0;
    is_refine[s] = 1;
    sa->child_id_face[s] = 0;
    SC3E (sc3_array_resize (idx_f_stack[s], 1));
    SC3E (sc3_array_index (idx_f_stack[s], 0, &arr));
    SC3E (sc3_array_index (*(sc3_array_t **) arr, 0, &begin));
    SC3E (sc3_array_index (*(sc3_array_t **) arr, 1, &end));

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
  p4est3_gloidx       proc_owner = p3->mpirank;
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
               (*(sc3_array_t **) stack_it[side],
                sa->child_id_face[side], &(b_f[side])));
    SC3E_FAST (sc3_array_index
               (*(sc3_array_t **) stack_it[side],
                sa->child_id_face[side] + 1, &(e_f[side])));
  }

  /* Check if both sides belong to remote process(es).
     If so, we ignore this face. */
  if ((sa->nsides == 1 && (*(b_f[0]) >= sa->local_end_face[0]
                           || *(e_f[0]) <= sa->local_begin_face[0]))
      || (sa->nsides == 2 && (*(b_f[0]) >= sa->local_end_face[0]
                              || *(e_f[0]) <= sa->local_begin_face[0])
          && (*(b_f[1]) >= sa->local_end_face[1]
              || *(e_f[1]) <= sa->local_begin_face[1]))) {
    return NULL;
  }

  /* first check if the whole quadrant passes */
  SC3E (sc3_array_index (sa->finfo->sides, 0, &fside));
  for (side = 0; side < sa->nsides; ++side) {
    if (!is_refine[side]) {
      continue;
    }

    if (p3->contiguous) {
      /* Get quadrant from shared storage */
      first_quad = (void *) (p3->nodequads[0] + p3->qsize * (*(b_f[side])));
    }
    else {
      /* TODO: make it through p4est3_find_partition and last_goffsets */
      /* Consider batch processing for binary search */
      /* Find process that owns this quadrant */
      /* Start binary search from a local process */
      SC3E (p4est3_owner_lookup_fast (p3, sa, *(b_f[side]), &proc_owner));

      /* Get quadrant from the owning process's window */
      first_quad = (void *) (p3->nodequads[proc_owner] +
                             p3->qsize * (*(b_f[side]) -
                                          p3->goffset[proc_owner]));
    }

    SC3E (p4est3_quadrant_level (p3->qvt, first_quad, &level));
    if (level == Level[side]) {
      is_refine[side] = 0;
      fside[side].nquad =
        (p4est3_locidx) (*(b_f[side]) - p3->gtroffset[trees[side]]);
      fside[side].quadrant = first_quad;
      if (*(b_f[side]) < p3->goffset[p3->mpirank] ||
          *(e_f[side]) > p3->goffset[p3->mpirank + 1]) {
        /* This is a ghost quadrant */
        fside[side].is_ghost = 1;
      }
      else {
        fside[side].is_ghost = 0;
      }
    }
  }

  if (!is_refine[0] && !is_refine[1]) {
    if (cface != NULL) {
      SC3E (cface (sa->finfo));
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
    SC3E (p4est3_array_set_zero (*(sc3_array_t **) (stack_it[side])));
#endif
    SC3E (p4est3_cached_quadrant_array_split_noncontig
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
        SC3E (p4est3_connectivity_get_neighbor_face_corner
              (p3->conn, fside[0].nface, fside[1].nface, ori, &idx));
      }
      SC3E (p4est3_connectivity_get_face_child_id
            (p3->conn, fside[side].nface, idx, &child_id));
      sa->child_id_face[side] = child_id;
      Level[side]++;
      is_lvl_increased[side] = 1;
    }
    SC3E (p4est3_internal_iterate_face (p3, cface, ccodim, sa));        /* keep full check here */
    for (side = 0; side < sa->nsides; ++side) {
      Level[side] = is_lvl_increased[side] ? Level[side] - 1 : Level[side];
      if (i != half_ch - 1) {
        is_refine[side] = is_lvl_increased[side];
      }
      else if (is_lvl_increased[side]) {
        SC3E (sc3_array_pop (idx_face_stack[side]));
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
      SC3E (p4est3_iterate_face_inner_init (p3, search_area, child, nb_id));
      SC3E (p4est3_internal_iterate_face (p3, cface, ccodim, search_area));
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

  SC3E (p4est3_tree_index (p3, tree, &sa->tree));
  sa->vinfo->ntree = tree;
  sa->finfo->orientation = 0;
  sa->finfo->tree_boundary = 0;
  SC3E (sc3_array_resize (sa->finfo->sides, 2));

  memset (sa->level2nchildren, 0, sizeof (int) * p3->qvt->max_level);
  SC3A_CHECK (sa->Level == 0);

  SC3E (sc3_array_index (sa->idx_vol_stack, 0, &arr));
#ifdef P4EST_ENABLE_DEBUG
  size_t              ecount;
  SC3E (sc3_array_get_elem_count (*(sc3_array_t **) arr, &ecount));
  SC3A_CHECK (ecount == 2);
  SC3E (sc3_array_get_elem_count (sa->idx_vol_stack, &ecount));
  SC3A_CHECK (ecount == 1);
#endif

  sa->child_id = 0;
  SC3E (sc3_array_index (*(sc3_array_t **) arr, sa->child_id, &begin));
  SC3E (sc3_array_index (*(sc3_array_t **) arr, sa->child_id + 1, &end));

  *begin = p3->gtroffset[tree];
  *end = p3->gtroffset[tree + 1];

  sa->local_begin = SC3_MAX (p3->gtroffset[tree], p3->goffset[p3->mpirank]);
  sa->local_end =
    SC3_MIN (p3->gtroffset[tree + 1], p3->goffset[p3->mpirank + 1]);

  SC3E (sc3_array_index (sa->finfo->sides, 0, &fside));
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

  if (p3->contiguous) {
    /* The rest functional is valid only for non-contiguous shared memory */
    return NULL;
  }

  if (tree == p3->fltree && p3->gtroffset[tree] != p3->goffset[p3->mpirank]) {
    SC3E (p4est3_find_partition
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
    SC3E (p4est3_find_partition
          (p3->alloc, p3->mpisize, p3->goffset, p3->goffset[p3->mpirank + 1],
           p3->gtroffset[tree + 1] - 1,
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
    SC3E (p4est3_owner_lookup_fast (p3, sa, *begin, &proc_owner));

    /* Get quadrant from the owning process's window */
    first_quad = (void *) (p3->nodequads[proc_owner] +
                           p3->qsize * (*begin - p3->goffset[proc_owner]));
  }

  SC3E (p4est3_quadrant_level (p3->qvt, first_quad, &level));
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
      SC3E (cvolume (vinfo));
    }
    l2nch[*Level]++;
    return NULL;
  }

  SC3E_FAST (sc3_array_push (idx_vol_stack, &stack_it));
#ifdef P4EST_ENABLE_DEBUG
  SC3E (p4est3_array_set_zero (*(sc3_array_t **) stack_it));
#endif

  SC3E (p4est3_cached_quadrant_array_split_noncontig
        (p3, view_q, *Level, *begin, *end, *(sc3_array_t **) stack_it, sa));
  l2nch[++(*Level)] = 0;

  /* since array_split doesn't count shift from the beinning of quadrants
     in a node, we shift result indices at the loop below */
  SC3E_FAST (sc3_array_index (*(sc3_array_t **) stack_it, 0, &arr_it));
  for (i = 0; i < max_children + 1; ++i) {
    arr_it[i] += *begin;
  }
  for (i = 0; i < max_children; ++i) {
    sa->child_id = i;
    SC3E (p4est3_iterate_volume_rec (p3, cvolume, cface, ccodim, sa));
  }
  SC3A_CHECK (l2nch[*Level] == max_children);
  SC3E (p4est3_iterate_face_inner (p3, cface, ccodim, sa));
  l2nch[--(*Level)]++;
  SC3E (sc3_array_pop (idx_vol_stack));
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

/* Small 2-slot level cache query */

static sc3_error_t *
p4est3_iterate_volume_iterative (p4est3_t *p3,
                                 p4est3_iterate_volume_t cvolume,
                                 p4est3_iterate_face_t cface,
                                 p4est3_iterate_codim_t ccodim,
                                 p4est3_search_area_t *sa)
{
  p4est3_gloidx      *begin, *end;
  p4est3_gloidx       proc_owner;
  void               *first_quad;
  int                 levelq;
  int                 max_children = sa->max_children;
  int                *l2nch = sa->level2nchildren;
  sc3_array_t        *idx_vol_stack = sa->idx_vol_stack;
  p4est3_iterate_volume_info_t *vinfo = sa->vinfo;

  /* Defensive initialization (recursive init already zeroes this, but make
     explicit for static analyzers and robustness). */
  if (SC_UNLIKELY (l2nch[0] != 0)) {
    memset (l2nch, 0, sizeof (int) * sa->max_level);
  }

  p4est3_iter_frame_t *frames;
  int                 top = -1; /* empty */
  /* allocate temp frame buffer on stack (bounded by max_level) */
  frames = (p4est3_iter_frame_t *) alloca (sizeof (*frames) * sa->max_level);

  /* Initialize owner cache */
  proc_owner = p3->mpirank;

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
        SC3E (p4est3_owner_lookup_fast (p3, sa, *begin, &proc_owner));
      }
      first_quad = (void *) (p3->nodequads[proc_owner] +
                             p3->qsize * (*begin - p3->goffset[proc_owner]));
    }
    SC3E (p4est3_quadrant_level (p3->qvt, first_quad, &levelq));
    if (levelq == sa->Level) {
      if (cvolume != NULL) {
        vinfo->quadrant = first_quad;
        vinfo->nquad =
          (p4est3_locidx) (*begin - p3->gtroffset[sa->tree->treeid]);
        SC3E (cvolume (vinfo));
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
      SC3E (p4est3_array_set_zero (*(sc3_array_t **) stack_it));
#endif
      SC3E (p4est3_cached_quadrant_array_split_noncontig
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
    if (fr->next_child >= max_children) {
      /* Finished all children at this level */
      SC3A_CHECK (l2nch[fr->level] == max_children);
      /* Process inner faces at this (still current) level */
      SC3E (p4est3_iterate_face_inner (p3, cface, ccodim, sa));
      /* Pop indices stack and decrement level */
      SC3E_FAST (sc3_array_pop (idx_vol_stack));
      --top;                    /* remove current frame first */
      sa->Level = fr->level - 1;
      if (sa->Level >= 0) {
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

    if (*begin >= sa->local_end || *end <= sa->local_begin) {
      /* Skip non-local */
      ++fr->next_child;
      l2nch[fr->level]++;
      continue;
    }

    /* Fetch quadrant level (avoid second query for first child if we already know it) */
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
        SC3E (p4est3_owner_lookup_fast (p3, sa, *begin, &proc_owner));
      }
      first_quad = (void *) (p3->nodequads[proc_owner] +
                             p3->qsize * (*begin - p3->goffset[proc_owner]));
    }
    if (fr->next_child == 0 && fr->have_first_level && *begin == arr_child[0]) {
      /* Reuse the level we already queried at parent range decision */
      levelq = fr->first_child_level;
      ++sa->first_child_level_reuses;
      /* Invalidate so only first child benefits */
      fr->have_first_level = 0;
    }
    else {
      SC3E (p4est3_quadrant_level (p3->qvt, first_quad, &levelq));
    }
    if (levelq == fr->level) {
      if (cvolume != NULL && proc_owner == p3->mpirank) {
        vinfo->quadrant = first_quad;
        vinfo->nquad =
          (p4est3_locidx) (*begin - p3->gtroffset[sa->tree->treeid]);
        SC3E (cvolume (vinfo));
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
      SC3E (p4est3_array_set_zero (*(sc3_array_t **) stack_it_child));
#endif
      SC3E (p4est3_cached_quadrant_array_split_noncontig
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
    SC3E (p4est3_iterate_volume (p3, cvolume, user_data));
    return NULL;
  }

  SC3E (p4est3_set_outer_data (p3, search_area, user_data));
  SC3E (sc3_array_index (search_area->finfo->sides, 0, &fside));
  for (tree = p3->fltree; tree <= p3->lltree; ++tree) {

    SC3E (p4est3_iterate_volume_rec_init (p3, search_area, tree));
    SC3E (p4est3_iterate_volume_iterative
          (p3, cvolume, cface, ccodim, search_area));

    /* frame faces part */
    search_area->finfo->tree_boundary = 1;
    search_area->treeid_face[0] = search_area->tree->treeid;
    for (face = 0; face < search_area->nfaces; ++face) {
      SC3E (p4est3_iterate_face_bound_init
            (p3, search_area, tree, face, fside, &is_iterate));
      if (is_iterate) {
        /* we iterate over such trees that tree_neighbor < tree */
        SC3E (p4est3_internal_iterate_face (p3, cface, ccodim, search_area));
      }
    }
  }
  SC3E (p4est3_destroy_outer_data (p3, search_area));
  return NULL;
}

/* MRU cache initialization function */
static inline sc3_error_t *
p4est3_split_cache_init (p4est3_t *p3, p4est3_search_area_t *sa,
                         int max_cache_size)
{
  int                 i;
  sa->cache_base_size = max_cache_size;
  sa->cache_max_size = max_cache_size;
  sa->cache_hits = 0;
  sa->cache_misses = 0;
  sa->cache_op_counter = 0;
  sa->cache_resize_interval = 512;      /* heuristic */

  /* Initialize unified cache (per-level disabled) */
  for (i = 0; i < P4EST3_ITER_CACHE_LVL; ++i) {
    sa->split_cache[i] = NULL;
  }
  sa->split_cache_unified =
    sc_hash_mru_new (p4est3_split_cache_hash, p4est3_split_cache_equal,
                     p4est3_split_cache_drop, sa, sa->cache_max_size);

  /* Initialize tier rings: capacity heuristic mirrors legacy implementation */
  sa->tier_capacity =
    (p3->mpisize == 1 ? p3->num_children : 2 * p3->num_children);
  (void) sc3_allocator_calloc (p3->alloc, (size_t) sa->max_level,
                               sizeof (p4est3_tier_ring_t),
                               (void **) &sa->tier_rings);
  for (i = 0; i < sa->max_level; ++i) {
    p4est3_tier_ring_init (&sa->tier_rings[i], p3->alloc, sa->tier_capacity);
  }
  sa->tier_hits = sa->tier_misses = 0;
  sa->tier_seed_total = 0;
  sa->tier_seed_hits = 0;
  sa->tier_seed_attempted = 0;
  sa->tier_seed_skipped = 0;
  sa->tier_seed_evicted = 0;
  sa->tier_placeholder_lookups = 0;
  sa->split_total = 0;
  sa->sibling_fastpath = 0;
  sa->splits_computed = 0;
  sa->tier_extensions = 0;
  sa->tier_subinterval_hits = 0;
  /* initialize last-used split cache */
  sa->last_split_level = -1;
  sa->last_split_tree = -1;
  sa->last_split_begin = -1;
  sa->last_split_end = -2;
  sa->last_split_hits = 0ULL;
  sa->last_split_uses = 0ULL;
  sa->is_contiguous = p3->contiguous;
  /* histogram removed */

  return NULL;
}

/* MRU cache cleanup function */
static sc3_error_t *
p4est3_split_cache_destroy (p4est3_search_area_t *sa)
{
  /* Optional statistics output controlled by env variable */
  const char         *env_stats = getenv ("P4EST3_ITERATE_CACHE_STATS");
  if (env_stats != NULL && *env_stats) {
    double              tier_total =
      (double) (sa->tier_hits + sa->tier_misses);
    double              mru_total =
      (double) (sa->cache_hits + sa->cache_misses);
    double              tier_hr =
      tier_total > 0.0 ? (100.0 * sa->tier_hits / tier_total) : 0.0;
    double              mru_hr =
      mru_total > 0.0 ? (100.0 * sa->cache_hits / mru_total) : 0.0;
    double              combined_hits =
      (double) (sa->tier_hits + sa->cache_hits);
    double              combined_lookups =
      (double) (sa->split_total - sa->sibling_fastpath);
    double              combined_rate = combined_lookups > 0.0 ?
      (100.0 * combined_hits / combined_lookups) : 0.0;
    long long           avoided_splits = (long long) combined_hits;
    double              avoided_rate = combined_lookups > 0.0 ?
      (100.0 * (double) avoided_splits / combined_lookups) : 0.0;
    double              seed_fill_rate = sa->tier_seed_total ?
      (100.0 * (double) sa->tier_seed_hits / (double) sa->tier_seed_total) :
      0.0;
    double              last_rate = sa->last_split_uses ?
      (100.0 * (double) sa->last_split_hits /
       (double) sa->last_split_uses) : 0.0;
    int                 rr_print = fprintf (stderr,
                                            "[p4est3_iterate] cache stats: total=%lld fastpath=%lld computed=%lld avoided=%lld avoided_rate=%.2f%% combined_rate=%.2f%% | last_used: uses=%llu hits=%llu rate=%.2f%% | tier_hits=%d misses=%d rate=%.2f%% ext=%d subhits=%d | seeds: placed=%d hits=%d fill=%.2f%% attempted=%d skipped=%d evicted=%d ph_lookups=%d | mru_hits=%d misses=%d rate=%.2f%% cap=%zu | subtree_full_local=%lld first_child_reuse=%lld\n",
                                            sa->split_total,
                                            sa->sibling_fastpath,
                                            sa->splits_computed,
                                            avoided_splits,
                                            avoided_rate,
                                            combined_rate,
                                            (unsigned long long)
                                            sa->last_split_uses,
                                            (unsigned long long)
                                            sa->last_split_hits,
                                            last_rate,
                                            sa->tier_hits, sa->tier_misses,
                                            tier_hr, sa->tier_extensions,
                                            sa->tier_subinterval_hits,
                                            sa->tier_seed_total,
                                            sa->tier_seed_hits,
                                            seed_fill_rate,
                                            sa->tier_seed_attempted,
                                            sa->tier_seed_skipped,
                                            sa->tier_seed_evicted,
                                            sa->tier_placeholder_lookups,
                                            sa->cache_hits, sa->cache_misses,
                                            mru_hr,
                                            sa->split_cache_unified ?
                                            sa->split_cache_unified->
                                            maxcount : 0UL,
                                            sa->subtree_full_local_hits,
                                            sa->first_child_level_reuses);
    /* histogram print removed */
    (void) rr_print;
  }
  /* Print cache statistics before cleanup */
//  if (sa->cache_hits > 0 || sa->cache_misses > 0) {
//    printf ("MRU Cache Statistics: Hits=%d, Misses=%d, Hit Rate=%.2f%%\n",
//            sa->cache_hits, sa->cache_misses,
//            100.0 * sa->cache_hits / (sa->cache_hits + sa->cache_misses));
//  }
//
//  /* Print detailed per-level cache statistics */
//  for (i = 0; i < P4EST3_ITER_CACHE_LVL; i++) {
//    if (sa->split_cache[i] != NULL && sa->split_cache[i]->count > 0) {
//      printf ("  Level %d: Cache entries=%lu, Insertions=%lu\n",
//              i, sa->split_cache[i]->count, sa->split_cache[i]->num_inserted);
//    }
//  }

  /* Destroy MRU cache for each level */
  if (sa->split_cache_unified != NULL) {
    sc_hash_mru_destroy (sa->split_cache_unified);
    sa->split_cache_unified = NULL;
  }

  /* Destroy tier rings */
  if (sa->tier_rings != NULL) {
    if (sa->alloc != NULL) {
      for (int i = 0; i < sa->max_level; ++i) {
        p4est3_tier_ring_destroy (&sa->tier_rings[i], sa->alloc);
      }
      (void) sc3_allocator_free (sa->alloc, sa->tier_rings);
    }
    sa->tier_rings = NULL;
  }

  /* Embedded entries freed with MRU nodes */

  return NULL;                  /* Success */
}

/* Cached array split function */
/* Cached version of array split function using unified cache entry structure */
static sc3_error_t *p4est3_cached_quadrant_array_split_noncontig
  (p4est3_t * p3, sc3_array_t * array, int level,
   p4est3_gloidx begin, p4est3_gloidx end,
   sc3_array_t * indices, p4est3_search_area_t * sa)
{
  p4est3_split_cache_entry_t search_entry;
  p4est3_split_cache_entry_t *cache_entry = NULL;
  void              **found;
  int                 inserted;
  sc_hash_mru_t      *cache;
  p4est3_gloidx      *src_val;
  p4est3_tier_entry_t *tier_hit = NULL;
  int                 tier_extend = 0;  /* need to extend stored end after computing */
  /* (tier_extend used later to decide extension; suppress unused-value warning in some analyzers) */
  if (0) {
    tier_extend = tier_extend;
  }

  sa->split_total++;
  /* Range length (only needed for sibling fast path check) */
  p4est3_gloidx       range_len = end - begin;
  /* Last-used direct cache: exact reuse only */
  if (sa->last_split_level == level &&
      sa->last_split_tree == sa->tree->treeid &&
      sa->last_split_begin == begin) {
    ++sa->last_split_uses;
    if (sa->last_split_end == end) {
      ++sa->last_split_hits;    /* exact match */
      p4est3_gloidx      *dst;
      SC3E (sc3_array_index (indices, 0, &dst));
      memcpy (dst, sa->last_split_results,
              sizeof (p4est3_gloidx) * (sa->max_children + 1));
      return NULL;
    }
  }
  else {
    ++sa->last_split_uses;      /* count attempt with different begin/tree/level */
  }

  /* Tier ring fast path (per-level, extremely small & hot) */
  if (SC_LIKELY
      (level >= 0 && level < sa->max_level && sa->tier_rings != NULL)) {
    p4est3_tier_ring_t *ring = &sa->tier_rings[level];
    tier_hit =
      p4est3_tier_ring_lookup (ring, sa->tree->treeid, level, begin, end);
    if (SC_LIKELY (tier_hit != NULL)) {
      if (tier_hit->valid && tier_hit->end >= end) {
        if (tier_hit->end > end) {
          ++sa->tier_subinterval_hits;
        }
        /* Cached data covers requested sub-interval (or equal) */
        SC3E (sc3_array_index (indices, 0, &src_val));
        memcpy (src_val, tier_hit->splits,
                sizeof (p4est3_gloidx) * (sa->max_children + 1));
        sa->tier_hits++;
        return NULL;
      }
      if (!tier_hit->valid) {
        ++sa->tier_placeholder_lookups; /* will compute below */
      }
      else if (tier_hit->end < end) {
        /* Have a prefix; we will recompute and then extend */
        tier_extend = 1;
      }
    }
    if (tier_hit == NULL) {
      sa->tier_misses++;
    }
  }

  /* Sibling pack fast path: exactly full set of children; direct split */
  if (SC_UNLIKELY (range_len == (p4est3_gloidx) sa->max_children)) {
    sa->sibling_fastpath++;
    if (sa->is_contiguous) {
      SC3E (sc3_array_renew_data (&array, p3->nodequads[0], p3->qsize,
                                  begin, range_len));
      return p4est3_quadrant_array_split (p3->qvt, array, level, indices);
    }
    else {
      SC3E (sc3_array_renew_data (&array, p3->nodequads[0], p3->qsize,
                                  0, range_len));
      return p4est3_quadrant_array_split_noncontig
        (p3, array, level, begin, indices);
    }
  }

  /* Check if caching is available for this level */
  /* Choose unified cache */
  if (sa->split_cache_unified == NULL) {
    sa->cache_misses++;
    if (sa->is_contiguous) {
      SC3E (sc3_array_renew_data
            (&array, p3->nodequads[0], p3->qsize, begin, end - begin));
      return p4est3_quadrant_array_split (p3->qvt, array, level, indices);
    }
    else {
      /* here we start with the very beginning of not necessary local node
         quadrants, because of our specialized array_split_noncontig function */
      SC3E (sc3_array_renew_data
            (&array, p3->nodequads[0], p3->qsize, 0, end - begin));
      return p4est3_quadrant_array_split_noncontig
        (p3, array, level, begin, indices);
    }
  }

  cache = sa->split_cache_unified;

  SC3A_CHECK (sa->tree != NULL);
  /* Create search entry with key information */
  search_entry.first_quad_id = begin;
  search_entry.last_quad_id = end;
  search_entry.level = level;
  search_entry.tree_id = sa->tree->treeid;
  search_entry.composite_key =
    (((unsigned long long) (unsigned int) level) << 56) ^
    (((unsigned long long) (unsigned int) sa->tree->treeid) << 40) ^
    ((unsigned long long) begin << 3) ^ (unsigned long long) (end - begin);
  search_entry.reuse_count = 0;

  /* Try to find in cache */
  inserted = sc_hash_mru_insert_unique (cache, &search_entry, &found);

  if (!inserted) {
    /* Cache hit - get the cached entry and copy its split result */
    cache_entry = (p4est3_split_cache_entry_t *) * found;
    cache_entry->reuse_count++;
    sa->cache_hits++;

    /* Copy cached split indices to the output array */
    /* Copy the cached data directly */
    SC3E (sc3_array_index (indices, 0, &src_val));
    if (sa->dim == 2) {
      memcpy (src_val, cache_entry->split_results2d,
              sizeof (p4est3_gloidx) * (sa->max_children + 1));
    }
    else {
      memcpy (src_val, cache_entry->split_results3d,
              sizeof (p4est3_gloidx) * (sa->max_children + 1));
    }
    /* Success - used cached result */
    return NULL;
  }

  /* Cache miss - allocate persistent cache entry and compute result */
  sa->cache_misses++;
  /* This path will result in a real split computation below */
  if (++sa->cache_op_counter == sa->cache_resize_interval) {
    /* Dynamic resize heuristic: expand if hit rate high, shrink if low */
    double              hr = (sa->cache_hits + sa->cache_misses) ?
      (double) sa->cache_hits / (double) (sa->cache_hits +
                                          sa->cache_misses) : 0.0;
    size_t              target = sa->cache_max_size;
    if (hr > 0.75 && sa->cache_max_size < sa->cache_base_size * 8) {
      target = (size_t) (sa->cache_max_size * 1.5) + 1;
    }
    else if (hr < 0.30 && sa->cache_max_size > sa->cache_base_size) {
      target = (size_t) (sa->cache_max_size / 1.5) + 1;
      if (target < (size_t) sa->cache_base_size)
        target = sa->cache_base_size;
    }
    if (target != (size_t) sa->cache_max_size) {
      sa->split_cache_unified->maxcount = target;       /* simple adjust */
      sa->cache_max_size = (int) target;
      sc_hash_mru_consolidate (sa->split_cache_unified);
    }
    sa->cache_op_counter = 0;
  }

  /* Retrieve embedded entry pointer from MRU node (inserted case) */
  cache_entry = (p4est3_split_cache_entry_t *) * found;
  *cache_entry = search_entry;
  cache_entry->reuse_count = 1;

  /* Compute the split result */
  ++sa->splits_computed;
  if (sa->is_contiguous) {
    SC3E (sc3_array_renew_data
          (&array, p3->nodequads[0], p3->qsize, begin, end - begin));
    SC3E (p4est3_quadrant_array_split (p3->qvt, array, level, indices));
  }
  else {
    SC3E (sc3_array_renew_data
          (&array, p3->nodequads[0], p3->qsize, 0, end - begin));
    SC3E (p4est3_quadrant_array_split_noncontig
          (p3, array, level, begin, indices));
  }

  /* Copy the computed split indices to the cache */
  SC3E (sc3_array_index (indices, 0, &src_val));
  if (sa->dim == 2) {
    memcpy (cache_entry->split_results2d, src_val,
            sizeof (p4est3_gloidx) * (sa->max_children + 1));
  }
  else {
    memcpy (cache_entry->split_results3d, src_val,
            sizeof (p4est3_gloidx) * (sa->max_children + 1));
  }
  /* Update last-used cache */
  sa->last_split_level = level;
  sa->last_split_tree = sa->tree->treeid;
  sa->last_split_begin = begin;
  sa->last_split_end = end;
  memcpy (sa->last_split_results, src_val,
          sizeof (p4est3_gloidx) * (sa->max_children + 1));

  /* Insert or extend tier ring entry */
  if (level >= 0 && level < sa->max_level && sa->tier_rings != NULL) {
    p4est3_tier_ring_t *ring = &sa->tier_rings[level];
    p4est3_tier_entry_t *e = NULL;
    if (tier_hit != NULL && (tier_extend || !tier_hit->valid)) {
      e = tier_hit;             /* extend or materialize existing */
    }
    else {
      /* search again for overwrite (cheap small ring) */
      for (int i = 0; i < ring->size; ++i) {
        p4est3_tier_entry_t *e2 = &ring->entries[i];
        if (e2->tree_id == sa->tree->treeid && e2->level == level
            && e2->begin == begin) {
          e = e2;
          break;
        }
      }
      if (e == NULL) {
        e =
          p4est3_tier_ring_insert (ring, sa->tree->treeid, level, begin, end);
      }
    }
    if (e != NULL) {
      if (e->seeded) {
        ++sa->tier_seed_hits;   /* seeded prediction realized (even if extension) */
      }
      if (e->end < end) {
        ++sa->tier_extensions;
        e->end = end;
      }
      memcpy (e->splits, src_val,
              sizeof (p4est3_gloidx) * (sa->max_children + 1));
      e->valid = 1;
      e->seeded = 0;
    }

    /* Pre-seed child intervals for next level (conservative: first large child only) */
    if (level + 1 < sa->max_level) {
      p4est3_tier_ring_t *next_ring = &sa->tier_rings[level + 1];
      int                 nchildren = sa->max_children;
      int                 seeded_one = 0;
      for (int cid = 0; cid < nchildren; ++cid) {
        p4est3_gloidx       cbeg = src_val[cid] + begin;
        p4est3_gloidx       cend = src_val[cid + 1] + begin;
        if (cend <= cbeg) {
          continue;
        }
        ++sa->tier_seed_attempted;      /* attempted consideration */
        /* Filter: require interval bigger than fastpath ( > max_children ) */
        if ((cend - cbeg) <= (p4est3_gloidx) sa->max_children) {
          ++sa->tier_seed_skipped;
          continue;
        }
        if (seeded_one) {       /* only first qualifying child */
          ++sa->tier_seed_skipped;
          continue;
        }
        /* Check if already present */
        int                 have = 0;
        for (int ti = 0; ti < next_ring->size; ++ti) {
          p4est3_tier_entry_t *te = &next_ring->entries[ti];
          if (te->tree_id == sa->tree->treeid && te->level == level + 1
              && te->begin == cbeg) {
            have = 1;
            break;
          }
        }
        if (!have) {
          /* Track eviction if overwriting a seeded placeholder */
          int                 pos = next_ring->next;
          p4est3_tier_entry_t *victim = next_ring->entries + pos;
          if (next_ring->size == next_ring->capacity && victim->seeded
              && !victim->valid) {
            ++sa->tier_seed_evicted;
          }
          p4est3_tier_entry_t *se =
            p4est3_tier_ring_insert (next_ring, sa->tree->treeid, level + 1,
                                     cbeg, cend);
          if (se != NULL) {
            se->end = cend;
            se->valid = 0;      /* placeholder */
            se->seeded = 1;
            ++sa->tier_seed_total;
            seeded_one = 1;     /* stop after first */
          }
        }
        else {
          ++sa->tier_seed_skipped;      /* already present counts as skip */
        }
        if (seeded_one) {
          break;                /* conservative strategy: only first large child */
        }
      }
    }
  }

  return NULL;
}
