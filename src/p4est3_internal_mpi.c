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

#include <p4est3_internal.h>
#include <sc3_omp.h>

sc3_error_t        *
p4est3_internal_setup_comm (p4est3_t * p3)
{
  int                 headsize, headrank;
  int                 p, next, *ofs;
  int                 dispunit;
  int                *nodesizemem;
  sc3_MPI_Aint_t      nodeabytes;

  /* this is a special-purpose function to simplify p4est3_setup */
  SC3A_CHECK (p3 != NULL);

  /* query input communicator */
  SC3E (sc3_MPI_Comm_size (p3->mpicomm, &p3->mpisize));
  SC3E (sc3_MPI_Comm_rank (p3->mpicomm, &p3->mpirank));

  /* create one communicator on each shared-memory node */
  SC3E (sc3_MPI_Comm_split_type (p3->mpicomm, SC3_MPI_COMM_TYPE_SHARED,
                                 0, SC3_MPI_INFO_NULL, &p3->nodecomm));
  SC3E (sc3_MPI_Comm_size (p3->nodecomm, &p3->nodesize));
  SC3E (sc3_MPI_Comm_rank (p3->nodecomm, &p3->noderank));

  /* create communicator that contains the first rank on each node */
  SC3E (sc3_MPI_Comm_split (p3->mpicomm, p3->noderank == 0 ? 0 :
                            SC3_MPI_UNDEFINED, 0, &p3->headcomm));
  SC3A_CHECK ((p3->noderank != 0) == (p3->headcomm == SC3_MPI_COMM_NULL));
  if (p3->noderank == 0) {
    SC3E (sc3_MPI_Comm_size (p3->headcomm, &headsize));
    SC3E (sc3_MPI_Comm_rank (p3->headcomm, &headrank));
    nodeabytes = (2 + 2 * headsize + 1) * sizeof (int);
  }
  else {
    headsize = headrank = 0;
    nodeabytes = 0;
  }

  /* create info structure to allow for per-rank allocation */
  SC3E (sc3_MPI_Info_create (&p3->info_noncontig));
  SC3E (sc3_MPI_Info_set
        (p3->info_noncontig, "alloc_shared_noncontig", "true"));

  /* allocate shared memory for information on node and head communicators */
  SC3E (sc3_MPI_Win_allocate_shared
        (nodeabytes, sizeof (int),
         p3->info_noncontig, p3->nodecomm, &nodesizemem, &p3->nodesizewin));
  if (p3->noderank == 0) {
    SC3E (sc3_MPI_Win_lock (SC3_MPI_LOCK_EXCLUSIVE, 0, SC3_MPI_MODE_NOCHECK,
                            p3->nodesizewin));
    nodesizemem[0] = p3->num_nodes = headsize;
    nodesizemem[1] = p3->node_num = headrank;
    p3->node_sizes = &nodesizemem[2];
    p3->node_frank = p3->mpirank;

    /* allgather information about all nodes and compute offsets */
    SC3E (sc3_MPI_Allgather (&p3->nodesize, 1, SC3_MPI_INT,
                             p3->node_sizes, 1, SC3_MPI_INT, p3->headcomm));
    *(ofs = p3->node_offsets = &nodesizemem[2 + headsize]) = 0;
    for (p = 0; p < headsize; ++p) {
      next = *ofs + p3->node_sizes[p];
      *++ofs = next;
    }
    SC3A_CHECK (p3->node_offsets[headrank] == p3->mpirank);
    SC3A_CHECK (p3->node_offsets[headsize] == p3->mpisize);

    /* make sure shared memory contents are consistent */
    SC3E (sc3_MPI_Win_unlock (0, p3->nodesizewin));
    SC3E (sc3_MPI_Barrier (p3->nodecomm));
  }
  else {
    SC3E (sc3_MPI_Win_shared_query (p3->nodesizewin, 0,
                                    &nodeabytes, &dispunit, &nodesizemem));
    SC3A_CHECK (nodeabytes >= (sc3_MPI_Aint_t) sizeof (int));
    SC3A_CHECK (dispunit == (int) sizeof (int));
    SC3A_CHECK (nodesizemem != NULL);

    /* access shared memory written by other process */
    SC3E (sc3_MPI_Barrier (p3->nodecomm));
    p3->num_nodes = nodesizemem[0];
    SC3A_CHECK (nodeabytes ==
                (sc3_MPI_Aint_t) ((2 + 2 * p3->num_nodes + 1) *
                                  sizeof (int)));
    p3->node_num = nodesizemem[1];
    p3->node_sizes = &nodesizemem[2];
    p3->node_offsets = &nodesizemem[2 + p3->num_nodes];
  }
  SC3A_CHECK (p3->node_frank == p3->node_offsets[p3->node_num]);

  return NULL;
}

sc3_error_t        *
p4est3_internal_setup_cut (p4est3_t * p3,
                           p4est3_gloidx num_uniform, int qsize)
{
#ifdef P4EST_ENABLE_DEBUG
  int                 dp;
#endif
  int                 beginr, endr;
  int                 dispunit;
  char               *gfposmem;
  p4est3_topidx      *gftreemem;
  p4est3_gloidx      *goffsetmem, num_global;
  sc3_MPI_Aint_t      gftreebytes, gfposbytes, goffsetbytes, tempbytes;
  sc3_omp_esync_t     esync, *s = &esync;

  /* this is a special-purpose function to simplify p4est3_setup */
  SC3A_CHECK (p3 != NULL);
  SC3A_CHECK (0 <= p3->mpirank && p3->mpirank < p3->mpisize);
  SC3A_CHECK (0 <= p3->noderank && p3->noderank < p3->nodesize);
  SC3A_CHECK (num_uniform > 0);
  SC3A_CHECK (qsize > 0);

  /* create shared partition arrays */
  gftreebytes = (p3->mpisize + 1) * sizeof (p4est3_topidx);
  SC3E (sc3_MPI_Win_allocate_shared
        (p3->noderank == 0 ? gftreebytes : 0, sizeof (p4est3_topidx),
         p3->info_noncontig, p3->nodecomm, &gftreemem, &p3->gftreewin));
  gfposbytes = (p3->mpisize + 1) * qsize;
  SC3E (sc3_MPI_Win_allocate_shared
        (p3->noderank == 0 ? gfposbytes : 0, qsize,
         p3->info_noncontig, p3->nodecomm, &gfposmem, &p3->gfposwin));
  goffsetbytes = (p3->mpisize + 1) * sizeof (p4est3_gloidx);
  SC3E (sc3_MPI_Win_allocate_shared
        (p3->noderank == 0 ? goffsetbytes : 0, sizeof (p4est3_gloidx),
         p3->info_noncontig, p3->nodecomm, &goffsetmem, &p3->goffsetwin));
  if (p3->noderank > 0) {
    SC3E (sc3_MPI_Win_shared_query (p3->gftreewin, 0,
                                    &tempbytes, &dispunit, &gftreemem));
    SC3A_CHECK (gftreebytes == tempbytes);
    SC3A_CHECK (dispunit == sizeof (p4est3_topidx));
    SC3A_CHECK (gftreemem != NULL);
    SC3E (sc3_MPI_Win_shared_query (p3->gfposwin, 0,
                                    &tempbytes, &dispunit, &gfposmem));
    SC3A_CHECK (gfposbytes == tempbytes);
    SC3A_CHECK (dispunit == qsize);
    SC3A_CHECK (gfposmem != NULL);
    SC3E (sc3_MPI_Win_shared_query (p3->goffsetwin, 0,
                                    &tempbytes, &dispunit, &goffsetmem));
    SC3A_CHECK (goffsetbytes == tempbytes);
    SC3A_CHECK (dispunit == (int) sizeof (p4est3_gloidx));
    SC3A_CHECK (goffsetmem != NULL);
  }

  /* compute global partition information fairly across node ranks */
  SC3E (sc3_MPI_Win_lock (SC3_MPI_LOCK_SHARED, 0, SC3_MPI_MODE_NOCHECK,
                          p3->gftreewin));
  SC3E (sc3_MPI_Win_lock (SC3_MPI_LOCK_SHARED, 0, SC3_MPI_MODE_NOCHECK,
                          p3->gfposwin));
  SC3E (sc3_MPI_Win_lock (SC3_MPI_LOCK_SHARED, 0, SC3_MPI_MODE_NOCHECK,
                          p3->goffsetwin));
  num_global = p3->num_trees * num_uniform;
  beginr = sc3_intcut (p3->mpisize + 1, p3->nodesize, p3->noderank);
  endr = sc3_intcut (p3->mpisize + 1, p3->nodesize, p3->noderank + 1);
  SC3E (sc3_omp_esync_init (s));
#pragma omp parallel
  {
    int                 beginrt = beginr;
    int                 endrt = endr;
    int                 pt;
    char               *qptr;
    char               *temp = p3->temp_quad[sc3_omp_thread_num ()];
    sc3_error_t        *e = NULL;

    /* parallelize process loop across threads */
    sc3_omp_thread_intrange (&beginrt, &endrt);
    qptr = gfposmem + beginrt * qsize;
    for (pt = beginrt; pt < endrt; ++pt) {
      gftreemem[pt] =
        (goffsetmem[pt] =
         p4est3_glocut (num_global, p3->mpisize, pt)) / num_uniform;
      SC3E_SET (e, p4est3_quadrant_morton
                (p3->qvt, p3->level,
                 goffsetmem[pt] - gftreemem[pt] * num_uniform, temp));
      SC3E_NULL_SET (e, p4est3_quadrant_first_descendant
                     (p3->qvt, temp, p3->qmaxlevel, qptr));
      SC3E_NULL_BREAK (e);
      qptr += qsize;
    }
    sc3_omp_esync (s, &e);
  }
  SC3E (sc3_omp_esync_summary (s));
  SC3E (sc3_MPI_Win_unlock (0, p3->gftreewin));
  SC3E (sc3_MPI_Win_unlock (0, p3->gfposwin));
  SC3E (sc3_MPI_Win_unlock (0, p3->goffsetwin));
  SC3E (sc3_MPI_Barrier (p3->nodecomm));
#ifdef P4EST_ENABLE_DEBUG
  for (dp = 0; dp <= p3->mpisize; ++dp) {
    SC3A_CHECK (gftreemem[dp] == goffsetmem[dp] / num_uniform);
    SC3A_CHECK (goffsetmem[dp] ==
                p4est3_glocut (num_global, p3->mpisize, dp));
  }
  SC3A_CHECK (gftreemem[p3->mpisize] == p3->num_trees);
  SC3A_CHECK (goffsetmem[p3->mpisize] == num_global);
#endif

  /* assign further object members */
  p3->qsize = qsize;
  p3->global_num_quads = num_global;
  p3->goffset = goffsetmem;
  p3->gftree = gftreemem;
  p3->gfpos = gfposmem;
  return NULL;
}

sc3_error_t        *
p4est3_tree_index (p4est3_t * p3, p4est3_topidx tt, p4est3_tree_t ** tree)
{
  void               *vt;

  SC3A_CHECK (p3 != NULL);
  SC3A_CHECK (p3->fltree <= tt && tt <= p3->lltree);
  SC3A_CHECK (tree != NULL);

  SC3E (sc3_array_index (p3->trees, tt - p3->fltree, &vt));
  *tree = (p4est3_tree_t *) vt;

  return NULL;
}

sc3_error_t        *
p4est3_internal_setup_tree (p4est3_t * p3, p4est3_gloidx num_uniform)
{
  int                 n;
  int                 dispunit;
  char               *quadmem, *nqmem;
  p4est3_topidx       tt;
  p4est3_gloidx       first_quad, end_quad, tt_offset, next_offset;
  p4est3_tree_t      *tree;
  sc3_MPI_Aint_t      quadbytes, tempbytes;

  /* this is a special-purpose function to simplify p4est3_setup */
  SC3A_CHECK (p3 != NULL);
  SC3A_CHECK (0 <= p3->mpirank && p3->mpirank < p3->mpisize);
  SC3A_CHECK (0 <= p3->noderank && p3->noderank < p3->nodesize);

  /* determine local trees */
  first_quad = p3->goffset[p3->mpirank];
  end_quad = p3->goffset[p3->mpirank + 1];
  SC3A_CHECK (end_quad - first_quad <= P4EST3_LOCIDX_MAX);
  if ((p3->local_num_quads = (p4est3_locidx) (end_quad - first_quad)) == 0) {
    p3->fltree = -1;
    p3->lltree = -2;
    p3->nltrees = 0;
  }
  else {
    p3->fltree = (p4est3_topidx) (first_quad / num_uniform);
    p3->lltree = (p4est3_topidx) ((end_quad - 1) / num_uniform);
    SC3A_CHECK (p3->fltree == p3->gftree[p3->mpirank]);
    p3->nltrees = p3->lltree - p3->fltree + 1;
  }

  /* create shared quadrant storage */
  SC3E_ALLOCATOR_MALLOC (p3->alloc, char *, p3->nodesize, p3->nodequads);
  quadbytes = p3->local_num_quads * p3->qsize;
  SC3E (sc3_MPI_Win_allocate_shared
        (quadbytes, p3->qsize,
         p3->info_noncontig, p3->nodecomm, &quadmem, &p3->quadwin));
  for (n = 0; n < p3->nodesize; ++n) {
    SC3E (sc3_MPI_Win_shared_query (p3->quadwin, n,
                                    &tempbytes, &dispunit, &nqmem));
    SC3A_CHECK (tempbytes == (sc3_MPI_Aint_t)
                ((p3->goffset[p3->node_frank + n + 1] -
                  p3->goffset[p3->node_frank + n]) * p3->qsize));
    SC3A_CHECK (dispunit == p3->qsize);
    SC3A_CHECK (nqmem != NULL || tempbytes == 0);
    p3->nodequads[n] = nqmem;
  }
  p3->quads = quadmem;
  SC3A_CHECK (p3->nodequads[p3->noderank] == p3->quads);

  /* populate tree metadata */
  SC3E (sc3_array_new (p3->alloc, &p3->trees));
  SC3E (sc3_array_set_elem_size (p3->trees, sizeof (p4est3_tree_t)));
  SC3E (sc3_array_set_elem_count (p3->trees, p3->nltrees));
  SC3E (sc3_array_setup (p3->trees));
  next_offset = 0;
  tt_offset = p3->fltree * num_uniform;
  for (tt = p3->fltree; tt <= p3->lltree; ++tt) {
    SC3E (p4est3_tree_index (p3, tt, &tree));
    tree->treeid = tt;
    tree->quad_offset = next_offset;
    tree->first_tquad = (tt == p3->fltree) ? first_quad - tt_offset : 0;
    if (tt == p3->lltree) {
      /* this is the last iteration: no need to update tt_offset */
      tree->end_tquad = end_quad - tt_offset;
    }
    else {
      /* TODO: double check whether tt_offset is needed here or not. */
      /* tree->end_tquad = tt_offset += num_uniform; */
      /* Temporary solution */
      tree->end_tquad = num_uniform;
      tt_offset += num_uniform;
    }
    /* by construction each local tree contains at least one element */
    SC3A_CHECK (0 <= tree->first_tquad &&
                tree->first_tquad < tree->end_tquad);
    tree->num_quads = tree->end_tquad - tree->first_tquad;
    SC3A_CHECK (0 < tree->num_quads && tree->num_quads <= num_uniform);
    next_offset = tree->quad_offset + tree->num_quads;
    tree->tquads = p3->quads + tree->quad_offset * p3->qsize;
  }
  SC3A_CHECK (tt_offset == p3->lltree * num_uniform);
  SC3A_CHECK (next_offset == p3->local_num_quads);
  return NULL;
}

static sc3_error_t *
p4est3_lowest_children (p4est3_t * p3, const void *q,
                        sc3_array_t * levelq, char **threadq)
{
  int                 i;
  void               *child;
  int                 level;
  SC3E (p4est3_quadrant_level (p3->qvt, q, &level));

  if (level < p3->level - 1) {
    SC3E (sc3_array_index (levelq, level + 1, &child));
    for (i = 0; i < p3->num_children; ++i) {
      SC3E (p4est3_quadrant_child (p3->qvt, q, i, child));
      SC3E (p4est3_lowest_children (p3, child, levelq, threadq));
    }
  }
  else if (level == p3->level - 1) {
    for (i = 0; i < p3->num_children; ++i, *threadq += p3->qsize) {
      SC3E (p4est3_quadrant_child (p3->qvt, q, i, *threadq));
      SC3A_IS (p3->qvt->quadrant_is_valid, *threadq);
    }
  }
  else {
    SC3E (p4est3_quadrant_copy (p3->qvt, q, *threadq));
    *threadq += p3->qsize;
  }

  return NULL;
}

static sc3_error_t *
p4est3_recursive_partition (p4est3_t * p3, int level,
                            p4est3_locidx rf,
                            p4est3_locidx rl,
                            p4est3_locidx mf,
                            p4est3_locidx ml,
                            sc3_array_t * levelq, char **threadq)
{
  int                 i;
  p4est3_locidx       n_lowerq;
  void               *q;

  if (rf >= mf && rl <= ml) {
    SC3E (sc3_array_index (levelq, level, &q));
    SC3E (p4est3_quadrant_morton (p3->qvt, level,
                                  rf >> (p3->qvt->dim * (p3->level - level)),
                                  q));
    SC3E (p4est3_lowest_children (p3, q, levelq, threadq));
  }
  else if (rf > ml || rl < mf) {
    return NULL;
  }
  else {
    n_lowerq = (rl - rf + 1) / p3->num_children;
    rl = rf + n_lowerq - 1;
    for (i = 0; i < p3->num_children; ++i, rf += n_lowerq, rl += n_lowerq) {
      SC3E (p4est3_recursive_partition (p3, level + 1, rf, rl, mf, ml,
                                        levelq, threadq));
    }
  }

  return NULL;
}

static sc3_error_t *
p4est3_recursive_partition_child (p4est3_t * p3, int level,
                                  p4est3_locidx rf,
                                  p4est3_locidx rl,
                                  p4est3_locidx mf,
                                  p4est3_locidx ml,
                                  sc3_array_t * levelq, char **threadq)
{
  int                 i;
  p4est3_locidx       n_lowerq;
  void               *q, *r;

  if (rf >= mf && rl <= ml) {
    SC3E (sc3_array_index (levelq, level, &q));
    SC3E (p4est3_lowest_children (p3, q, levelq, threadq));
  }
  else if (rf > ml || rl < mf) {
    return NULL;
  }
  else {
    n_lowerq = (rl - rf + 1) / p3->num_children;
    rl = rf + n_lowerq - 1;
    for (i = 0; i < p3->num_children; ++i, rf += n_lowerq, rl += n_lowerq) {
      SC3A_CHECK (level < p3->level);
      SC3E (sc3_array_index (levelq, level, &q));
      SC3E (sc3_array_index (levelq, level + 1, &r));
      SC3E (p4est3_quadrant_child (p3->qvt, q, i, r));
      SC3E (p4est3_recursive_partition_child (p3, level + 1, rf, rl, mf, ml,
                                              levelq, threadq));
    }
  }

  return NULL;
}

static sc3_error_t *
p4est3_region (p4est3_t * p3, void *a, void *b, sc3_array_t * region)
{
  int                 la, lb, lc;
  int                 i, j, j1, j2, ecount;
  p4est3_gloidx       aid, bid, cid;
  void               *c, *q;
  sc3_array_t        *testq, *buff;
  sc3_allocator_t    *alloc;

  p4est3_quadrant_level (p3->qvt, a, &la);
  p4est3_quadrant_level (p3->qvt, b, &lb);
  SC3A_CHECK (la == lb);

  SC3E (p4est3_quadrant_linear_id (p3->qvt, a, p3->level, &aid));
  SC3E (p4est3_quadrant_linear_id (p3->qvt, b, p3->level, &bid));

  SC3E (sc3_array_index
        (p3->talloc, sc3_omp_thread_num (), (void **) &(alloc)));
  SC3E (sc3_array_new (*(sc3_allocator_t **) alloc, &testq));
  SC3E (sc3_array_set_elem_size (testq, p3->qsize));
  SC3E (sc3_array_set_elem_count (testq, p3->num_children + 1));
  SC3E (sc3_array_set_resizable (testq, 1));
  SC3E (sc3_array_setup (testq));

  SC3E (sc3_array_new (*(sc3_allocator_t **) alloc, &buff));
  SC3E (sc3_array_set_elem_size (buff, p3->qsize));
  SC3E (sc3_array_set_elem_count (buff, 2));
  SC3E (sc3_array_setup (buff));

  SC3E (p4est3_quadrant_compare (p3->qvt, a, b, &j1));
  SC3A_CHECK (j1 < 0);
  SC3E (sc3_array_index (testq, 0, &c));
  SC3E (p4est3_nearest_common_ancestor (p3->qvt, a, b, c));
  for (i = 0; i < p3->num_children; ++i) {
    SC3E (sc3_array_index (testq, i + 1, &q));
    SC3E (p4est3_quadrant_child (p3->qvt, c, i, q));
  }
  SC3E (sc3_array_get_elem_count (testq, &ecount));
  SC3A_CHECK (ecount == p3->num_children + 1);
  for (i = 1; i < ecount; ++i) {
    SC3E (sc3_array_index (testq, i, &c));
    SC3E (p4est3_quadrant_linear_id (p3->qvt, c, p3->level, &cid));
    SC3E (p4est3_quadrant_level (p3->qvt, c, &lc));
    SC3E (p4est3_quadrant_is_ancestor (p3->qvt, c, b, &j1));
    if ((aid < cid || (aid == cid && la <= lc))
        && (cid < bid || (cid == bid && lc < lb))
        && !j1) {
      SC3E (sc3_array_push (region, c));
    }
    else {
      SC3E (p4est3_quadrant_is_ancestor (p3->qvt, c, a, &j1));
      SC3E (p4est3_quadrant_is_ancestor (p3->qvt, c, b, &j2));
      if (j1 || j2) {
        SC3E (sc3_array_index (buff, 0, &q));
        SC3E (p4est3_quadrant_copy (p3->qvt, c, q));
        c = q;
        for (j = 0; j < p3->num_children; ++j) {
          SC3E (sc3_array_index (buff, 1, &q));
          SC3E (p4est3_quadrant_child (p3->qvt, c, j, q));
          SC3E (sc3_array_push (testq, q));
        }
        ecount += p3->num_children;
      }
    }
  }

  SC3E (sc3_array_destroy (&testq));
  SC3E (sc3_array_destroy (&buff));
  return NULL;
}

static sc3_error_t *
p4est3_region_end (p4est3_t * p3, void *a, void *b, sc3_array_t * region)
{
  int                 la, lc;
  int                 i, j, j1, ecount;
  p4est3_gloidx       aid, cid;
  void               *c, *q;
  sc3_array_t        *testq, *buff;
  sc3_allocator_t    *alloc;

  p4est3_quadrant_level (p3->qvt, a, &la);

  SC3E (p4est3_quadrant_linear_id (p3->qvt, a, p3->level, &aid));

  SC3E (sc3_array_index
        (p3->talloc, sc3_omp_thread_num (), (void **) &(alloc)));
  SC3E (sc3_array_new (*(sc3_allocator_t **) alloc, &testq));
  SC3E (sc3_array_set_elem_size (testq, p3->qsize));
  SC3E (sc3_array_set_elem_count (testq, 2));
  SC3E (sc3_array_set_resizable (testq, 1));
  SC3E (sc3_array_setup (testq));

  SC3E (sc3_array_new (*(sc3_allocator_t **) alloc, &buff));
  SC3E (sc3_array_set_elem_size (buff, p3->qsize));
  SC3E (sc3_array_set_elem_count (buff, 2));
  SC3E (sc3_array_setup (buff));

  SC3E (sc3_array_index (testq, 1, &c));
  SC3E (p4est3_nearest_common_ancestor (p3->qvt, a, b, c));

  SC3E (sc3_array_get_elem_count (testq, &ecount));
  SC3A_CHECK (ecount == 2);
  for (i = 1; i < ecount; ++i) {
    SC3E (sc3_array_index (testq, i, &c));
    SC3E (p4est3_quadrant_linear_id (p3->qvt, c, p3->level, &cid));
    SC3E (p4est3_quadrant_level (p3->qvt, c, &lc));
    if (aid < cid || (aid == cid && la <= lc)) {
      SC3E (sc3_array_push (region, c));
    }
    else {
      SC3E (p4est3_quadrant_is_ancestor (p3->qvt, c, a, &j1));
      if (j1) {
        SC3E (sc3_array_index (buff, 0, &q));
        SC3E (p4est3_quadrant_copy (p3->qvt, c, q));
        c = q;
        for (j = 0; j < p3->num_children; ++j) {
          SC3E (sc3_array_index (buff, 1, &q));
          SC3E (p4est3_quadrant_child (p3->qvt, c, j, q));
          SC3E (sc3_array_push (testq, q));
        }
        ecount += p3->num_children;
      }
    }
  }

  SC3E (sc3_array_destroy (&testq));
  SC3E (sc3_array_destroy (&buff));
  return NULL;
}

static sc3_error_t *
p4est3_recursive_partition_region (p4est3_t * p3, int is_region_end,
                                   p4est3_locidx mf, p4est3_locidx ml,
                                   sc3_array_t * levelq, char **threadq)
{
  sc3_allocator_t    *alloc;
  sc3_array_t        *region;
  void               *a, *b;
  char               *threadq_ptr;
  int                 rcount, i;
  p4est3_gloidx       id;

  SC3E (sc3_array_index (levelq, 0, &a));
  SC3E (sc3_array_index (levelq, 1, &b));
  SC3E (p4est3_quadrant_morton (p3->qvt, p3->level, mf, a));

  if (mf == ml) {
    SC3E (p4est3_quadrant_copy (p3->qvt, a, *threadq));
    *threadq += p3->qsize;
    return NULL;
  }

  SC3E (sc3_array_index
        (p3->talloc, sc3_omp_thread_num (), (void **) &(alloc)));
  SC3E (sc3_array_new (*(sc3_allocator_t **) alloc, &region));
  SC3E (sc3_array_set_elem_size (region, p3->qsize));
  SC3E (sc3_array_set_resizable (region, 1));
  SC3E (sc3_array_setup (region));

  if (is_region_end) {
    SC3E (p4est3_quadrant_morton (p3->qvt, p3->level, ml, b));
    SC3E (p4est3_region_end (p3, a, b, region));
  }
  else {
    SC3E (p4est3_quadrant_morton (p3->qvt, p3->level, ml + 1, b));
    SC3E (p4est3_region (p3, a, b, region));
  }

  SC3E (sc3_array_get_elem_count (region, &rcount));
  for (i = 0; i < rcount; ++i) {
    SC3E (sc3_array_index (region, i, &a));
    SC3E (p4est3_quadrant_linear_id (p3->qvt, a, p3->level, &id));
    threadq_ptr = *threadq + p3->qsize * (id - mf);
    SC3E (p4est3_lowest_children (p3, a, levelq, &threadq_ptr));
  }
  *threadq += p3->qsize * (ml - mf + 1);

  SC3E (sc3_array_destroy (&region));
  return NULL;
}

/** Binary search a local quad number in the local trees */
static sc3_error_t *
p4est3_local_quad_tree (p4est3_t * p3,
                        p4est3_locidx local_num, p4est3_tree_t ** ptree)
{
  p4est3_topidx       mint, maxt, guess;

  /* sanity checks */
  SC3A_CHECK (p3 != NULL);
  SC3A_CHECK (0 <= local_num && local_num < p3->local_num_quads);
  SC3A_CHECK (ptree != NULL);

  /* begin search with local range of trees (inclusive) */
  mint = p3->fltree;
  maxt = p3->lltree;
  for (;;) {
    SC3A_CHECK (mint <= maxt);

    /* have we found our result? */
    if (mint == maxt) {
      SC3E (p4est3_tree_index (p3, mint, ptree));
      return NULL;
    }

    /* if not, it is important to look ahead of previous minimum */
    guess = (mint + maxt + 1) / 2;
    SC3E (p4est3_tree_index (p3, guess, ptree));
    if (local_num < (*ptree)->quad_offset) {
      /* the quadrant is on a lower tree */
      maxt = guess - 1;
    }
    else {
      /* the quadrant is on this or a higher tree */
      mint = guess;
    }
  }
}

/* TODO: char * is a good convention for type? */
static sc3_error_t *
p4est3_internal_populate_morton (p4est3_locidx tmine, p4est3_t * p3,
                                 p4est3_locidx * tq, p4est3_gloidx * gq,
                                 char **charq)
{
  for (; *tq < tmine; ++(*tq), ++(*gq), *charq += p3->qsize) {
    SC3E (p4est3_quadrant_morton (p3->qvt, p3->level, *gq, *charq));
  }
  return NULL;
}

static sc3_error_t *
p4est3_internal_populate_successor (p4est3_locidx tmine, p4est3_t * p3,
                                    p4est3_locidx * tq, p4est3_gloidx * gq,
                                    char **charq)
{
  if (*tq < tmine) {
    char               *cq_prev = *charq;
    SC3E (p4est3_quadrant_morton (p3->qvt, p3->level, *gq, *charq));
    ++(*tq);
    ++(*gq);
    *charq += p3->qsize;
    for (; *tq < tmine;
         ++(*tq), ++(*gq), cq_prev = *charq, *charq += p3->qsize) {
      SC3E (p4est3_quadrant_successor (p3->qvt, cq_prev, *charq));
    }
  }
  return NULL;
}

static sc3_error_t *
p4est3_internal_populate_recursive (p4est3_locidx tmine, p4est3_t * p3,
                                    p4est3_locidx * tq, p4est3_gloidx * gq,
                                    char **charq)
{
  if (*tq < tmine) {
    sc3_array_t        *levelq;
    sc3_allocator_t    *alloc;
    const p4est3_locidx rl = (1 << (p3->qvt->dim * p3->level)) - 1;
    const p4est3_locidx ml = *gq + (tmine - *tq) - 1;

    /* TODO: use per-thread allocotor here */
    SC3E (sc3_array_index
          (p3->talloc, sc3_omp_thread_num (), (void **) &(alloc)));
    SC3E (sc3_array_new (*(sc3_allocator_t **) alloc, &levelq));
    SC3E (sc3_array_set_elem_size (levelq, p3->qsize));
    SC3E (sc3_array_set_elem_alloc (levelq, p3->level + 1));
    SC3E (sc3_array_set_elem_count (levelq, p3->level + 1));
    SC3E (sc3_array_set_initzero (levelq, 1));
    SC3E (sc3_array_setup (levelq));
    SC3E (p4est3_recursive_partition (p3, 0, 0, rl, *gq, ml, levelq, charq));
    //SC3E (p4est3_recursive_partition_child (p3, 0, 0, rl, *gq, ml, levelq, charq));
    //SC3E (p4est3_recursive_partition_region (p3, ml == rl ? 1 : 0, *gq, ml, levelq, charq));
    SC3E (sc3_array_destroy (&levelq));
    *tq = tmine;
    *gq = ml;
  }
  return NULL;
}

static sc3_error_t *
p4est3_internal_populate (p4est3_locidx tmine, p4est3_t * p3,
                          p4est3_locidx * tq, p4est3_gloidx * gq,
                          char **charq)
{
  switch (p3->setup_mode) {
  case P4EST3_NEW_MORTON:
    SC3E (p4est3_internal_populate_morton (tmine, p3, tq, gq, charq));
    break;
  case P4EST3_NEW_SUCCESSOR:
    SC3E (p4est3_internal_populate_successor (tmine, p3, tq, gq, charq));
    break;
  case P4EST3_NEW_RECURSIVE:
    SC3E (p4est3_internal_populate_recursive (tmine, p3, tq, gq, charq));
    break;
  default:
    SC3E_UNREACH ("wrong setup mode");
  }
  return NULL;
}

sc3_error_t        *
p4est3_internal_setup_quadrants (p4est3_t * p3)
{
  sc3_omp_esync_t     esync, *s = &esync;
  int                 tcount;
  sc3_allocator_t    *malloc;

  /* this is a special-purpose function to simplify p4est3_setup */
  SC3A_CHECK (p3 != NULL && p3->quads != NULL);
  SC3A_CHECK (0 <= p3->mpirank && p3->mpirank < p3->mpisize);
  SC3A_CHECK (0 <= p3->noderank && p3->noderank < p3->nodesize);

  /* we work on the process-local window onte the quadrants */
  SC3E (sc3_MPI_Win_lock (SC3_MPI_LOCK_EXCLUSIVE, p3->noderank,
                          SC3_MPI_MODE_NOCHECK, p3->quadwin));
  SC3E (sc3_omp_esync_init (s));
#pragma omp parallel
  {
    const int           tnum = sc3_omp_num_threads ();
    const int           tid = sc3_omp_thread_num ();
    char               *charq;
    sc3_error_t        *e;
    p4est3_locidx       first_quad_num, end_quad_num, tmine, tq;
    p4est3_gloidx       gq;
    p4est3_tree_t      *tree;
    sc3_allocator_t    *alloc;

    /* TODO: if recursive mode is selected, create one allocator per thread
       derived from p3->alloc.
       Please see sc/example/v3basics/basics.c
       Would it make sense to allocate the per-thread allocators
       persistent through the lifetime of the p4est3 object.
     */

    /* find tree sub-range for each thread separately */
    first_quad_num = p4est3_loccut (p3->local_num_quads, tnum, tid);
    end_quad_num = p4est3_loccut (p3->local_num_quads, tnum, tid + 1);
    SC3E_SET (e, p4est3_local_quad_tree (p3, first_quad_num, &tree));

    if (p3->setup_mode == P4EST3_NEW_RECURSIVE) {
      if (tid == 0) {
        SC3E_NULL_SET (e, sc3_array_new (p3->alloc, &p3->talloc));
        SC3E_NULL_SET (e,
                       sc3_array_set_elem_size (p3->talloc,
                                                sizeof (sc3_allocator_t *)));
        SC3E_NULL_SET (e, sc3_array_set_elem_count (p3->talloc, tnum));
        SC3E_NULL_SET (e, sc3_array_setup (p3->talloc));
      }
#pragma omp barrier
      SC3E_NULL_SET (e,
                     sc3_array_index (p3->talloc, tid, (void **) &(alloc)));
#pragma omp critical
      {
        SC3E_NULL_SET (e,
                       sc3_allocator_new (p3->alloc,
                                          (sc3_allocator_t **) alloc));
        SC3E_NULL_SET (e, sc3_allocator_setup (*(sc3_allocator_t **) alloc));
        sc3_omp_esync_in_critical (s, &e);
      }
    }
    if (e == NULL) {
      tq = first_quad_num;
      gq = tree->first_tquad + (first_quad_num - tree->quad_offset);
      charq = p3->quads + first_quad_num * p3->qsize;

      /* loop over subset of local trees */
      for (;;) {
        SC3E_NULL_REQ (e, tree->quad_offset <= tq);
        SC3E_NULL_REQ (e, tq < tree->quad_offset + tree->num_quads);
        tmine = SC3_MIN (end_quad_num, tree->quad_offset + tree->num_quads);

        /* loop over quadrants in local tree with creating of quadrants
           by selected method */
        SC3E_SET (e, p4est3_internal_populate (tmine, p3, &tq, &gq, &charq));

        SC3E_NULL_REQ (e, tq <= end_quad_num);
        if (tq == end_quad_num) {
          break;
        }

        /* move forward to next tree */
        SC3E_NULL_SET (e, p4est3_tree_index (p3, tree->treeid + 1, &tree));
        SC3E_NULL_BREAK (e);
        gq = 0;
      }
    }
    sc3_omp_esync (s, &e);
  }
  SC3E (sc3_omp_esync_summary (s));
  if (p3->setup_mode == P4EST3_NEW_RECURSIVE) {
    SC3E (sc3_array_get_elem_count (p3->talloc, &tcount));
    for (int i = 0; i < tcount; ++i) {
      SC3E (sc3_array_index (p3->talloc, i, (void **) &(malloc)));
      SC3E (sc3_allocator_destroy ((sc3_allocator_t **) malloc));
    }
    SC3E (sc3_array_destroy (&p3->talloc));
  }
  SC3E (sc3_MPI_Win_unlock (p3->noderank, p3->quadwin));
  SC3E (sc3_MPI_Barrier (p3->nodecomm));
  return NULL;
}
