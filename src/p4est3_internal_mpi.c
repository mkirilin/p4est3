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
#include <p4est3_refine.h>
#include <sc3_omp.h>

#ifndef P4EST_ENABLE_OPENMP
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#endif

sc3_error_t        *
p4est3_internal_setup_comm (p4est3_t * p3)
{
  /* this is a special-purpose function to simplify p4est3_setup */
  SC3A_CHECK (p3 != NULL);

  /* make a new default mpi environment object */
  SC3E (sc3_mpienv_new (p3->alloc, &p3->split_info));

  /* set the same mpi environment's members as the current p4est has */
  SC3E (sc3_mpienv_set_comm (p3->split_info, p3->mpicomm, 0));
  SC3E (sc3_mpienv_set_shared (p3->split_info, p3->shared));

  /* query input communicator */
  SC3E (sc3_MPI_Comm_size (p3->mpicomm, &p3->mpisize));
  SC3E (sc3_MPI_Comm_rank (p3->mpicomm, &p3->mpirank));

  /* make calculations and setup mpi environment */
  SC3E (sc3_mpienv_setup (p3->split_info));
  return NULL;
}

sc3_error_t        *
p4est3_internal_setup_cut (p4est3_t * p3,
                           p4est3_gloidx num_uniform, int qsize)
{
#ifdef P4EST_ENABLE_DEBUG
  int                 dp;
#endif
  int                 noderank, nodesize;
  int                 beginr, endr;
  int                 dispunit;
  char               *gfposmem;
  p4est3_topidx      *gftreemem;
  p4est3_gloidx      *goffsetmem, num_global;
  sc3_MPI_Aint_t      gftreebytes, gfposbytes, goffsetbytes, tempbytes;
  sc3_MPI_Comm_t      nodecomm;
  sc3_MPI_Info_t      info_noncontig;
  sc3_omp_esync_t     esync, *s = &esync;

  SC3E (sc3_mpienv_get_noderank (p3->split_info, &noderank));
  SC3E (sc3_mpienv_get_nodesize (p3->split_info, &nodesize));

  /* this is a special-purpose function to simplify p4est3_setup */
  SC3A_CHECK (p3 != NULL);
  SC3A_CHECK (0 <= p3->mpirank && p3->mpirank < p3->mpisize);
  SC3A_CHECK (0 <= noderank && noderank < nodesize);
  SC3A_CHECK (num_uniform > 0);
  SC3A_CHECK (qsize > 0);

  SC3E (sc3_mpienv_get_nodecomm (p3->split_info, &nodecomm));
  SC3E (sc3_mpienv_get_info_noncont (p3->split_info, &info_noncontig));

  /* create shared partition arrays */
  gftreebytes = (p3->mpisize + 1) * sizeof (p4est3_topidx);
  SC3E (sc3_MPI_Win_allocate_shared
        (noderank == 0 ? gftreebytes : 0, sizeof (p4est3_topidx),
         info_noncontig, nodecomm, &gftreemem, &p3->gftreewin));
  gfposbytes = (p3->mpisize + 1) * qsize;
  SC3E (sc3_MPI_Win_allocate_shared
        (noderank == 0 ? gfposbytes : 0, qsize,
         info_noncontig, nodecomm, &gfposmem, &p3->gfposwin));
  goffsetbytes = (p3->mpisize + 1) * sizeof (p4est3_gloidx);
  SC3E (sc3_MPI_Win_allocate_shared
        (noderank == 0 ? goffsetbytes : 0, sizeof (p4est3_gloidx),
         info_noncontig, nodecomm, &goffsetmem, &p3->goffsetwin));
  if (noderank > 0) {
    SC3E (sc3_MPI_Win_shared_query (p3->gftreewin, 0,
                                    &tempbytes, &dispunit, &gftreemem));
    SC3A_CHECK (tempbytes >= gftreebytes);
    SC3A_CHECK (dispunit == sizeof (p4est3_topidx));
    SC3A_CHECK (gftreemem != NULL);
    SC3E (sc3_MPI_Win_shared_query (p3->gfposwin, 0,
                                    &tempbytes, &dispunit, &gfposmem));
    SC3A_CHECK (tempbytes >= gfposbytes);
    SC3A_CHECK (dispunit == qsize);
    SC3A_CHECK (gfposmem != NULL);
    SC3E (sc3_MPI_Win_shared_query (p3->goffsetwin, 0,
                                    &tempbytes, &dispunit, &goffsetmem));
    SC3A_CHECK (tempbytes >= goffsetbytes);
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
  beginr = sc3_intcut (p3->mpisize + 1, nodesize, noderank);
  endr = sc3_intcut (p3->mpisize + 1, nodesize, noderank + 1);
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
  SC3E (sc3_MPI_Barrier (nodecomm));
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
  int                 nodesize;
#ifdef P4EST_ENABLE_DEBUG
  int                 noderank, node_frank;
#endif
  char               *quadmem, *nqmem;
  p4est3_topidx       tt;
  p4est3_gloidx       first_quad, end_quad, tt_offset, next_offset;
  p4est3_tree_t      *tree;
  sc3_MPI_Aint_t      quadbytes, tempbytes;
  sc3_MPI_Comm_t      nodecomm;
  sc3_MPI_Info_t      info_noncontig;

#ifdef P4EST_ENABLE_DEBUG
  SC3E (sc3_mpienv_get_noderank (p3->split_info, &noderank));
  SC3E (sc3_mpienv_get_node_frank (p3->split_info, &node_frank));
#endif
  SC3E (sc3_mpienv_get_nodesize (p3->split_info, &nodesize));
  /* this is a special-purpose function to simplify p4est3_setup */
  SC3A_CHECK (p3 != NULL);
  SC3A_CHECK (0 <= p3->mpirank && p3->mpirank < p3->mpisize);
  SC3A_CHECK (0 <= noderank && noderank < nodesize);

  SC3E (sc3_mpienv_get_nodecomm (p3->split_info, &nodecomm));
  SC3E (sc3_mpienv_get_info_noncont (p3->split_info, &info_noncontig));

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
  SC3E (sc3_allocator_malloc (p3->alloc, nodesize * sizeof (char *),
                              &p3->nodequads));
  quadbytes = (sc3_MPI_Aint_t) p3->local_num_quads * p3->qsize;
  SC3E (sc3_MPI_Win_allocate_shared
        (quadbytes, p3->qsize,
         info_noncontig, nodecomm, &quadmem, &p3->quadwin));
  for (n = 0; n < nodesize; ++n) {
    SC3E (sc3_MPI_Win_shared_query (p3->quadwin, n,
                                    &tempbytes, &dispunit, &nqmem));
    SC3A_CHECK (tempbytes >= (sc3_MPI_Aint_t)
                ((p3->goffset[node_frank + n + 1] -
                  p3->goffset[node_frank + n]) * p3->qsize));
    SC3A_CHECK (dispunit == p3->qsize);
    SC3A_CHECK (nqmem != NULL || tempbytes == 0);
    p3->nodequads[n] = nqmem;
  }
  p3->quads = quadmem;
  SC3A_CHECK (p3->nodequads[noderank] == p3->quads);

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
      /* if it's not the last tree, we always end at top right corner */
      tt_offset += tree->end_tquad = num_uniform;
    }
    SC3A_CHECK (tree->end_tquad > 0);
    tree->last_tquad = tree->end_tquad - 1;
    /* by construction each local tree contains at least one element */
    SC3A_CHECK (0 <= tree->first_tquad &&
                tree->first_tquad <= tree->last_tquad &&
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
  SC3E (sc3_array_set_resizable (testq, 1));
  SC3E (sc3_array_setup (testq));

  SC3E (sc3_array_new (*(sc3_allocator_t **) alloc, &buff));
  SC3E (sc3_array_set_elem_size (buff, p3->qsize));
  SC3E (sc3_array_set_elem_count (buff, 2));
  SC3E (sc3_array_setup (buff));

  SC3E (p4est3_quadrant_compare (p3->qvt, a, b, &j1));
  SC3A_CHECK (j1 < 0);
  SC3E (sc3_array_push (testq, &c));
  SC3E (p4est3_nearest_common_ancestor (p3->qvt, a, b, c));
  for (i = 0; i < p3->num_children; ++i) {
    SC3E (sc3_array_push (testq, &q));
    SC3E (p4est3_quadrant_child (p3->qvt, c, i, q));
  }
  SC3E (sc3_array_get_elem_count (testq, &ecount));
  SC3A_CHECK (ecount == p3->num_children + 1);
  for (i = 1; i < ecount; ++i) {
    SC3E (sc3_array_index (testq, i, &c));
    SC3E (p4est3_quadrant_linear_id (p3->qvt, c, p3->level, &cid));
    SC3E (p4est3_quadrant_level (p3->qvt, c, &lc));
    SC3E (p4est3_quadrant_is_ancestor (p3->qvt, c, b, &j2));
    if ((aid < cid || (aid == cid && la <= lc))
        && (cid < bid || (cid == bid && lc < lb))
        && !j2) {
      SC3E (sc3_array_push (region, &q));
      SC3E (p4est3_quadrant_copy (p3->qvt, c, q));
    }
    else {
      SC3E (p4est3_quadrant_is_ancestor (p3->qvt, c, a, &j1));
      if (j1 || j2) {
        for (j = 0; j < p3->num_children; ++j) {
          SC3E (sc3_array_push (testq, &q));
          SC3E (p4est3_quadrant_child (p3->qvt, c, j, q));
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
  SC3E (sc3_array_set_resizable (testq, 1));
  SC3E (sc3_array_setup (testq));

  SC3E (sc3_array_new (*(sc3_allocator_t **) alloc, &buff));
  SC3E (sc3_array_set_elem_size (buff, p3->qsize));
  SC3E (sc3_array_set_elem_count (buff, 2));
  SC3E (sc3_array_setup (buff));

  SC3E (sc3_array_push (testq, &c));
  SC3E (p4est3_nearest_common_ancestor (p3->qvt, a, b, c));

  SC3E (sc3_array_get_elem_count (testq, &ecount));
  SC3A_CHECK (ecount == 1);
  for (i = 0; i < ecount; ++i) {
    SC3E (sc3_array_index (testq, i, &c));
    SC3E (p4est3_quadrant_linear_id (p3->qvt, c, p3->level, &cid));
    SC3E (p4est3_quadrant_level (p3->qvt, c, &lc));
    if (aid < cid || (aid == cid && la <= lc)) {
      SC3E (sc3_array_push (region, &q));
      SC3E (p4est3_quadrant_copy (p3->qvt, c, q));
    }
    else {
      SC3E (p4est3_quadrant_is_ancestor (p3->qvt, c, a, &j1));
      if (j1) {
        for (j = 0; j < p3->num_children; ++j) {
          SC3E (sc3_array_push (testq, &q));
          SC3E (p4est3_quadrant_child (p3->qvt, c, j, q));
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
    switch (p3->setup_mode) {
    case P4EST3_NEW_RECURSIVE:
      SC3E (p4est3_recursive_partition
            (p3, 0, 0, rl, *gq, ml, levelq, charq));
      break;
    case P4EST3_NEW_RECURSIVE_CHILD:
      SC3E (p4est3_recursive_partition_child
            (p3, 0, 0, rl, *gq, ml, levelq, charq));
      break;
    case P4EST3_NEW_RECURSIVE_REGION:
      SC3E (p4est3_recursive_partition_region
            (p3, ml == rl ? 1 : 0, *gq, ml, levelq, charq));
      break;
    default:
      SC3E_UNREACH ("wrong setup mode");
    }
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
  case P4EST3_NEW_RECURSIVE_CHILD:
    SC3E (p4est3_internal_populate_recursive (tmine, p3, tq, gq, charq));
    break;
  case P4EST3_NEW_RECURSIVE_REGION:
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
  int                 noderank, nodesize;
  sc3_allocator_t    *malloc;
  sc3_MPI_Comm_t      nodecomm;

  SC3E (sc3_mpienv_get_noderank (p3->split_info, &noderank));
  SC3E (sc3_mpienv_get_nodesize (p3->split_info, &nodesize));

  /* this is a special-purpose function to simplify p4est3_setup */
  SC3A_CHECK (p3 != NULL && p3->quads != NULL);
  SC3A_CHECK (0 <= p3->mpirank && p3->mpirank < p3->mpisize);
  SC3A_CHECK (0 <= noderank && noderank < nodesize);

  /* we work on the process-local window onte the quadrants */
  SC3E (sc3_MPI_Win_lock (SC3_MPI_LOCK_EXCLUSIVE, noderank,
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

    if (p3->setup_mode == P4EST3_NEW_RECURSIVE
        || p3->setup_mode == P4EST3_NEW_RECURSIVE_CHILD
        || p3->setup_mode == P4EST3_NEW_RECURSIVE_REGION) {
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
      charq = p3->quads + (p4est3_gloidx) first_quad_num *p3->qsize;

      /* loop over subset of local trees */
      for (;;) {
        SC3E_NULL_REQ (e, tree->quad_offset <= tq);
        SC3E_NULL_REQ (e, tq < tree->quad_offset + tree->num_quads);
        tmine =
          SC3_MIN (end_quad_num,
                   (p4est3_gloidx) tree->quad_offset + tree->num_quads);

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
  if (p3->setup_mode == P4EST3_NEW_RECURSIVE
      || p3->setup_mode == P4EST3_NEW_RECURSIVE_CHILD
      || p3->setup_mode == P4EST3_NEW_RECURSIVE_REGION) {
    SC3E (sc3_array_get_elem_count (p3->talloc, &tcount));
    for (int i = 0; i < tcount; ++i) {
      SC3E (sc3_array_index (p3->talloc, i, (void **) &(malloc)));
      SC3E (sc3_allocator_destroy ((sc3_allocator_t **) malloc));
    }
    SC3E (sc3_array_destroy (&p3->talloc));
  }
  SC3E (sc3_MPI_Win_unlock (noderank, p3->quadwin));
  SC3E (sc3_mpienv_get_nodecomm (p3->split_info, &nodecomm));
  SC3E (sc3_MPI_Barrier (nodecomm));
  return NULL;
}

/** TODO: add a coordinate-based translation to a quadrant virtual table. */
/* Where to put this function for the best? */
static sc3_error_t *
p4est3_internal_translate_quadrant (p4est3_quadrant_vtable_t * qvt_old,
                                    p4est3_quadrant_vtable_t * qvt_new,
                                    void *qin, void *qout)
{
  int                 level;
  p4est3_gloidx       id;

  SC3A_IS (p4est3_quadrant_vtable_is_valid, qvt_old);
  SC3A_IS (p4est3_quadrant_vtable_is_valid, qvt_new);
  SC3A_IS2 (p4est3_quadrant_is2_valid, qvt_old, qin);
  SC3A_IS2 (p4est3_quadrant_is2_valid, qvt_new, qout);

  if (qvt_old == qvt_new) {
    /* just hardcopy the quadrant */
    SC3E (p4est3_quadrant_copy (qvt_old, qin, qout));
  }
  else {
    SC3E (p4est3_quadrant_level (qvt_old, qin, &level));
    SC3A_CHECK (level <= qvt_new->max_level);
    SC3E (p4est3_quadrant_linear_id (qvt_old, qin, level, &id));
    SC3E (p4est3_quadrant_morton (qvt_new, level, id, qout));
  }

  return NULL;
}

sc3_error_t        *
p4est3_internal_setup_from_source (p4est3_t * p3)
{
  int                 i, nodesize, noderank;
  int                 dispunit, beginr, endr;
  p4est3_t           *old = p3->old;
  sc3_MPI_Aint_t      gftreebytes, tempbytes, gfposbytes;
  sc3_MPI_Info_t      info_noncontig;
  sc3_MPI_Comm_t      nodecomm;

  SC3A_IS (p4est3_is_new, p3);
  SC3A_CHECK (p3->old != NULL);
  SC3A_IS (p4est3_is_setup, p3->old);

  /* variables of internal state used during the whole lifetime */
  /* or is it == 0? */
  p3->accessed_conn = old->accessed_conn;

  /* virtual implementation variables */
  if (old->slf != NULL) {
    /** TODO: decide what to do if a virtual implementation exists */
  }
  /* variables set before p4est3_setup */
  /* this call also sets p4est3_t::commdup */
  SC3E (p4est3_set_comm (p3, old->mpicomm, 1));

  /* this call also sets p4est3_t::num_trees */
  SC3E (p4est3_set_connectivity (p3, old->conn));

  /* inherit qvt only if we didn't set it up visibly */
  if (p3->qvt == NULL) {
    /* this call also sets p4est3_t::sqvt */
    SC3E (p4est3_set_quadrant_vtable (p3, old->qvt));
  }
  /* the next two calls are not necessary,
     but leave them here for completeness */
  SC3E (p4est3_set_level (p3, old->level));
  SC3E (p4est3_set_setup_mode (p3, old->setup_mode));

  /* variables populated during p4est3_setup: communicator related */
  p3->mpisize = old->mpisize;
  p3->mpirank = old->mpirank;
  SC3E (p4est3_set_shared (p3, old->shared));
  /* p3->split_info == NULL at the current code state,
     but we check it just in case of future development */
  if (p3->split_info != NULL) {
    SC3E (sc3_mpienv_unref (&p3->split_info));
  }
  p3->split_info = old->split_info;
  SC3E (sc3_mpienv_ref (p3->split_info));

  /* variables populated during p4est3_setup: partition related */
 /** TODO: why is it int type while p4est3_quadrant_size returs size_t? */
  p3->qsize = (int) p4est3_quadrant_size (p3->qvt);
  p3->qmaxlevel = old->qmaxlevel;
  p3->num_children = old->num_children;
  p3->max_threads = old->max_threads;

  SC3E (sc3_mpienv_get_noderank (p3->split_info, &noderank));
  SC3E (sc3_mpienv_get_info_noncont (p3->split_info, &info_noncontig));
  SC3E (sc3_mpienv_get_nodecomm (p3->split_info, &nodecomm));

  gftreebytes = (p3->mpisize + 1) * sizeof (p4est3_topidx);
  SC3E (sc3_MPI_Win_allocate_shared
        (noderank == 0 ? gftreebytes : 0, sizeof (p4est3_topidx),
         info_noncontig, nodecomm, &p3->gftree, &p3->gftreewin));
  gfposbytes = (p3->mpisize + 1) * p3->qsize;
  SC3E (sc3_MPI_Win_allocate_shared
        (noderank == 0 ? gfposbytes : 0, p3->qsize,
         info_noncontig, nodecomm, &p3->gfpos, &p3->gfposwin));

  if (noderank > 0) {
    SC3E (sc3_MPI_Win_shared_query (p3->gftreewin, 0,
                                    &tempbytes, &dispunit, &p3->gftree));
    SC3A_CHECK (tempbytes >= gftreebytes);
    SC3A_CHECK (dispunit == sizeof (p4est3_topidx));
    SC3A_CHECK (p3->gftree != NULL);
    SC3E (sc3_MPI_Win_shared_query (p3->gfposwin, 0,
                                    &tempbytes, &dispunit, &p3->gfpos));
    SC3A_CHECK (tempbytes >= gfposbytes);
    SC3A_CHECK (dispunit == p3->qsize);
    SC3A_CHECK (p3->gfpos != NULL);
  }
  SC3E (sc3_MPI_Win_lock (SC3_MPI_LOCK_SHARED, 0, SC3_MPI_MODE_NOCHECK,
                          p3->gftreewin));
  SC3E (sc3_MPI_Win_lock (SC3_MPI_LOCK_SHARED, 0, SC3_MPI_MODE_NOCHECK,
                          p3->gfposwin));

  /* We make a hardcopy (with allocating new shared memory)
     of p4est3_t::gftreewin and p4est3_t::gfposwin so far.
     In the future we will introduce reference interface for them since they
     stay the same in a new tree. */
  SC3E (sc3_mpienv_get_nodesize (p3->split_info, &nodesize));
  beginr = sc3_intcut (p3->mpisize + 1, nodesize, noderank);
  endr = sc3_intcut (p3->mpisize + 1, nodesize, noderank + 1);
  for (i = 0; i < endr - beginr; ++i) {
    p3->gftree[i] = old->gftree[i];
    SC3E (p4est3_internal_translate_quadrant
          (old->qvt, p3->qvt, (void *) (old->gfpos + i * old->qsize),
           (void *) (p3->gfpos + i * p3->qsize)));
  }
  SC3E (sc3_MPI_Win_unlock (0, p3->gftreewin));
  SC3E (sc3_MPI_Win_unlock (0, p3->gfposwin));

  SC3E (sc3_allocator_malloc (p3->alloc, p3->max_threads * sizeof (char *),
                              &p3->temp_quad));
  for (i = 0; i < p3->max_threads; ++i) {
    SC3E (sc3_allocator_malloc (p3->alloc, p3->qsize, &p3->temp_quad[i]));
  }
  /* variables populated during p4est3_setup: tree and quadrant storage */
  p3->fltree = old->fltree;
  p3->lltree = old->lltree;
  p3->nltrees = old->nltrees;

  /* functions set before p4est3_setup */
  if (p3->crefine == NULL) {
    p3->crefine = old->crefine;
  }
  /** We set inside all the values left, namely:
   * p4est3_t::nodequads, local_num_quads, quadwin,
   * quads, trees, goffsetwin, goffset and global_num_quads.
  */
  SC3E (p4est3_refine (p3));

  p3->setup = 1;
  SC3A_IS (p4est3_is_setup, p3);
  return NULL;
}
