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

#include <p4est3_refine.h>
#include <p4est3_iterate.h>
#include <p4est3_internal.h>

#ifdef __cplusplus
extern              "C"
{
#if 0
}
#endif
#endif

typedef struct refine_callback_data
{
  int                 counter;
  int                 n_new;
  sc3_array_t        *pattern;
  p4est3_refine_callback_t crefine;
}
refine_callback_data_t;

typedef struct coarsen_callback_data
{
  int                 counter;
  sc3_array_t        *family; /**< Array of pointers to quadrants*/
  sc3_array_t        *pattern;
  int                 nsiblings;
  p4est3_coarsen_callback_t ccoarse;
}
coarsen_callback_data_t;

static sc3_error_t *
p4est3_refine_array_new (sc3_allocator_t * alloc, size_t esize, int ealloc,
                         int ecount, sc3_array_t ** arr)
{
  SC3E_RETVAL (arr, NULL);
  SC3A_IS (sc3_allocator_is_setup, alloc);
  SC3A_CHECK (ealloc >= 0);

  SC3E (sc3_array_new (alloc, arr));
  SC3E (sc3_array_set_elem_size (*arr, esize));
  SC3E (sc3_array_set_elem_alloc (*arr, ealloc));
  SC3E (sc3_array_set_elem_count (*arr, ecount));
  SC3E (sc3_array_set_initzero (*arr, 1));
  SC3E (sc3_array_setup (*arr));

  return NULL;
}

/* Specific volume iterator callback that processes volume info by refinement
   callback and returns (with user_data) population pattern
   (see refinement documentation).
   The pattern consists of possible values: 0 and num_children.
   0 means a simple copying (probably, with translation);
   num_children means creating all the children of the corresponding quadrant
   and inserting them into a new forest.*/
static sc3_error_t *
p4est3_refine_volume_callback (p4est3_iterate_volume_info_t * vi)
{
  int                 is_refine, level;
  refine_callback_data_t *cdata = (refine_callback_data_t *) vi->user_data;
  char               *pattern_it;

  /* pack data for refinement callback input */
  p4est3_refine_callback_info_t ri;
  ri.p3 = vi->p3;
  ri.ntree = vi->ntree;
  ri.quadrant = vi->quadrant;
  ri.qvt = vi->p3->qvt;
  ri.user_data = vi->p3->user_data;

  SC3A_CHECK (cdata->crefine != NULL);
  SC3E (cdata->crefine (&ri, &is_refine));
  SC3E (sc3_array_index (cdata->pattern, cdata->n_new, &pattern_it));
  SC3E (p4est3_quadrant_level (vi->p3->qvt, vi->quadrant, &level));
  is_refine = level == vi->p3->qmaxlevel ? 0 : is_refine;
  if (!is_refine) {
    *pattern_it = 0;
    cdata->counter++;
  }
  else {
    *pattern_it = vi->p3->num_children;
    cdata->counter += vi->p3->num_children;
  }
  cdata->n_new++;
  return NULL;
}

static sc3_error_t *
p4est3_coarsen_volume_zero (int child_id, coarsen_callback_data_t *cdata,
                            p4est3_iterate_volume_info_t * vi, int *is_ret)
{
  int i;
  char *pattern_it;
  void              **quad;
  p4est3_tree_t *tree;

  SC3E_RETVAL (is_ret, 0);
  if (child_id == 0) {
    /* the beginning of a new family,
       process previous (part of the) family */
    for (i = 0; i < cdata->nsiblings; ++i) {
      SC3E (sc3_array_index (cdata->pattern, cdata->counter, &pattern_it));
      *pattern_it = 0;
      cdata->counter++;
    }
    SC3E (sc3_array_index (cdata->family, 0, &quad));
    *quad = vi->quadrant;
    cdata->nsiblings = 1;
    /* Check if it is the last quadrant in a forest.
       If so, process quadrants in current family. */
    SC3E (p4est3_tree_index (vi->p3, vi->ntree, &tree));
    if (vi->nquad - tree->first_tquad + tree->quad_offset + 1 ==
        vi->p3->local_num_quads) {
      for (i = 0; i < cdata->nsiblings; ++i) {
        SC3E (sc3_array_index (cdata->pattern, cdata->counter, &pattern_it));
        *pattern_it = 0;
        cdata->counter++;
      }
    }
    *is_ret = 1;
  }
  return NULL;
}

static sc3_error_t *
p4est3_coarsen_volume_complete (coarsen_callback_data_t *cdata,
                                p4est3_iterate_volume_info_t * vi,
                                int *is_coarsen, int *is_ret)
{
  int i;
  char *pattern_it;
  p4est3_coarsen_callback_info_t ci;
  p4est3_tree_t *tree;
  SC3E_RETVAL (is_coarsen, 0);
  SC3E_RETVAL (is_ret, 0);

  if (cdata->nsiblings == vi->p3->num_children) {
    /* we have complete family, pack data and call coarsening */
    ci.p3 = vi->p3;
    ci.ntree = vi->ntree;
    ci.family = cdata->family;
    ci.qvt = vi->p3->qvt;
    ci.user_data = vi->p3->user_data;
    SC3A_CHECK (cdata->ccoarse != NULL);
    SC3E (cdata->ccoarse (&ci, is_coarsen));
    cdata->nsiblings = 0;
  }
  else {
    /* Check if it is the last quadrant in a forest.
       If so, process quadrants in current family. */
    SC3E (p4est3_tree_index (vi->p3, vi->ntree, &tree));
    if (vi->nquad - tree->first_tquad + tree->quad_offset + 1 ==
        vi->p3->local_num_quads) {
      for (i = 0; i < cdata->nsiblings; ++i) {
        SC3E (sc3_array_index (cdata->pattern, cdata->counter, &pattern_it));
        *pattern_it = 0;
        cdata->counter++;
      }
    }
    /* nothing left to do here, go to the next volume */
    *is_ret = 1;
  }
  return NULL;
}

static sc3_error_t *
p4est3_coarsen_volume_fill (coarsen_callback_data_t *cdata, int is_coarsen, int num_children)
{
  int i;
  char *pattern_it;
  /* fill in the level information */
  if (is_coarsen) {
    SC3E (sc3_array_index (cdata->pattern, cdata->counter, &pattern_it));
    *pattern_it = 1;
    cdata->counter++;
  }
  else {
    for (i = 0; i < num_children; ++i) {
      SC3E (sc3_array_index (cdata->pattern, cdata->counter, &pattern_it));
      *pattern_it = 0;
      cdata->counter++;
    }
  }
  return NULL;
}

/* Specific volume iterator callback that processes volume info by coarsining
   callback and returns (with user_data) population pattern
   (see coarsining documentation).
   The pattern consists of possible values: 0 and 1.
   0 means a simple copying (probably, with translation);
   1 means creating the parent for the next num_children
   corresponding quadrants and inserting them into a new forest.
   We set 1 only for the first quadrant in a family. */
static sc3_error_t *
p4est3_coarsen_volume_callback (p4est3_iterate_volume_info_t * vi)
{
  int                 is_coarsen, child_id, is_ret;
  char               *pattern_it;
  void              **quad;
  coarsen_callback_data_t *cdata = (coarsen_callback_data_t *) vi->user_data;

  /* Decide if we call coarse callback.
     We do this only if we find a whole family. */
  SC3E (p4est3_quadrant_child_id (vi->p3->qvt, vi->quadrant, &child_id));
  SC3E (p4est3_coarsen_volume_zero (child_id, cdata, vi, &is_ret));
  if (is_ret) {
    return NULL;
  }
  if (cdata->nsiblings != child_id) {
    /* cannot be a part of a complete family */
    SC3E (sc3_array_index (cdata->pattern, cdata->counter, &pattern_it));
    *pattern_it = 0;
    cdata->counter++;

    return NULL;
  }

  /* cdata->nsiblings == child_id: the volume is a part of the family */
  SC3E (sc3_array_index (cdata->family, cdata->nsiblings, &quad));
  *quad = vi->quadrant;
  cdata->nsiblings++;

  SC3E (p4est3_coarsen_volume_complete (cdata, vi, &is_coarsen, &is_ret));
  if (is_ret) {
    return NULL;
  }

  SC3E (p4est3_coarsen_volume_fill (cdata, is_coarsen, vi->p3->num_children));
  return NULL;
}

static sc3_error_t *
p4est3_translate_quadrant (const p4est3_quadrant_vtable_t * qvt_old,
                           const p4est3_quadrant_vtable_t * qvt_new,
                           const void *qin, void *qout, int32_t * c)
{
  int                 level;

  SC3A_IS (p4est3_quadrant_vtable_is_valid, qvt_old);
  SC3A_IS (p4est3_quadrant_vtable_is_valid, qvt_new);
  SC3A_IS2 (p4est3_quadrant_vtable_is2_valid, qvt_old, qin);
  SC3A_CHECK (qvt_old->dim == qvt_new->dim);

  if (qvt_old == qvt_new) {
    /* just hardcopy the quadrant */
    SC3E (p4est3_quadrant_copy (qvt_old, qin, qout));
  }
  else {
    SC3E (p4est3_quadrant_level (qvt_old, qin, &level));
    SC3A_CHECK (level <= qvt_new->max_level);
    SC3E (p4est3_quadrant_coordinates (qvt_old, qin, qvt_old->dim, c));
    SC3E (p4est3_quadrant_quadrant (qvt_new, c, level, qout));
  }

  SC3A_IS2 (p4est3_quadrant_vtable_is2_valid, qvt_new, qout);
  return NULL;
}

static sc3_error_t *
p4est3_pattern_populate_tree_coarse (p4est3_t * p3, p4est3_tree_t * tree,
                                     sc3_array_t * pattern, int *lt_offset,
                                     int32_t * c)
{
  int                 i = 0, n_new_quads = 0;
  char               *n_insert;
  void               *quad_old, *quad_new;
  p4est3_tree_t      *oldtree;

  SC3E (p4est3_tree_index (p3->old, tree->treeid, &oldtree));
  while (i < oldtree->num_quads) {
    SC3E (sc3_array_index (pattern, *lt_offset + n_new_quads, &n_insert));
    SC3A_CHECK ((int) *n_insert == 0 || (int) *n_insert == 1);
    quad_old = (void *) (oldtree->tquads + i * p3->old->qsize);
    if ((int) *n_insert == 0) {
      quad_new = (void *) (tree->tquads + n_new_quads * p3->qsize);
      SC3E (p4est3_translate_quadrant
            (p3->old->qvt, p3->qvt, quad_old, quad_new, c));
      n_new_quads++;
      i++;
    }
    else {
      SC3E (p4est3_translate_quadrant
            (p3->old->qvt, p3->qvt, quad_old, p3->temp_quad[0], c));
      quad_new = (void *) (tree->tquads + n_new_quads * p3->qsize);
      SC3E (p4est3_quadrant_parent (p3->qvt, p3->temp_quad[0], quad_new));
      n_new_quads++;
      i += p3->num_children;
    }
  }
  *lt_offset += n_new_quads;
  return NULL;
}

static sc3_error_t *
p4est3_pattern_populate_tree_ref (p4est3_t * p3, p4est3_tree_t * tree,
                                  sc3_array_t * pattern, int *lt_offset,
                                  int32_t * c)
{
  int                 i, nch, n_new_quads = 0;
  char               *n_insert;
  void               *quad_old, *quad_new;
  p4est3_tree_t      *oldtree;

  SC3E (p4est3_tree_index (p3->old, tree->treeid, &oldtree));
  for (i = 0; i < oldtree->num_quads; ++i) {
    SC3E (sc3_array_index (pattern, oldtree->quad_offset + i, &n_insert));
    SC3A_CHECK ((int) *n_insert == 0
                || (int) *n_insert == p3->old->num_children);
    quad_old = (void *) (oldtree->tquads + i * p3->old->qsize);
    if ((int) *n_insert == 0) {
      quad_new = (void *) (tree->tquads + n_new_quads * p3->qsize);
      SC3E (p4est3_translate_quadrant
            (p3->old->qvt, p3->qvt, quad_old, quad_new, c));
      n_new_quads++;
    }
    else {
      SC3E (p4est3_translate_quadrant
            (p3->old->qvt, p3->qvt, quad_old, p3->temp_quad[0], c));
      for (nch = 0; nch < p3->num_children; ++nch) {
        quad_new = (void *) (tree->tquads + n_new_quads * p3->qsize);
        SC3E (p4est3_quadrant_child
              (p3->qvt, p3->temp_quad[0], nch, quad_new));
        n_new_quads++;
      }
    }
  }
  *lt_offset += n_new_quads;
  return NULL;
}

static sc3_error_t *
p4est3_populate_tree_cpy (p4est3_t *p3, p4est3_tree_t * tree,
                          int *lt_offset, int32_t *c)
{
  int i;
  void *quad_old, *quad_new;
  p4est3_tree_t *oldtree;

  SC3E (p4est3_tree_index (p3->old, tree->treeid, &oldtree));
  for (i = 0; i < oldtree->num_quads; ++i) {
    quad_old = (void *) (oldtree->tquads + i * p3->old->qsize);
    quad_new = (void *) (tree->tquads + i * p3->qsize);
    SC3E (p4est3_translate_quadrant
          (p3->old->qvt, p3->qvt, quad_old, quad_new, c));
  }
  *lt_offset += oldtree->num_quads;
  return NULL;
}

sc3_error_t        *
p4est3_fill_from_source (p4est3_t * p3)
{
  int                 i, nodesize;
  int                 dispunit;
  int                 noderank;
  char               *quadmem, *nqmem;
  int32_t            *coords;
  p4est3_locidx       lt_offset, num_quads = 0;
  p4est3_locidx      *local_num_quads;  /**< Array of the numbers of quadrants at every rank */
  p4est3_locidx      *first_tree_quads; /**< Array of the numbers of quadrants at the first local tree */
  sc3_MPI_Info_t      info_noncontig;
  sc3_MPI_Comm_t      nodecomm;
  sc3_MPI_Aint_t      tempbytes, goffsetbytes;
  p4est3_tree_t      *tree;
  sc3_array_t        *pattern;

  coarsen_callback_data_t scdata, *cdata = &scdata;
  refine_callback_data_t srdata, *rdata = &srdata;

  /* We suppose to call this function after setting up routine */
  SC3A_CHECK (p3->old != NULL);
  SC3A_IS (p4est3_is_setup, p3->old);

  SC3E (sc3_mpienv_get_nodesize (p3->old->split_info, &nodesize));
  SC3E (sc3_mpienv_get_noderank (p3->old->split_info, &noderank));
  SC3E (sc3_mpienv_get_nodecomm (p3->old->split_info, &nodecomm));
  SC3E (sc3_mpienv_get_info_noncont (p3->old->split_info, &info_noncontig));
  /* Do some preliminary allocations */
  SC3E (sc3_allocator_calloc
        (p3->alloc, p3->mpisize, sizeof (p4est3_locidx), &local_num_quads));
  SC3E (sc3_allocator_calloc
        (p3->alloc, p3->mpisize, sizeof (p4est3_locidx), &first_tree_quads));

  SC3E (p4est3_refine_array_new
        (p3->alloc, sizeof (char), 0, p3->old->local_num_quads, &pattern));

  SC3E (p4est3_refine_array_new (p3->alloc, sizeof (p4est3_tree_t),
                                 p3->nltrees, p3->nltrees, &p3->trees));

  /* All the preparations are done. Begin with the algorithm. */
  /* A call to fill in level information, that is necessary for
     generation of a new forest's mesh */
  if (p3->crefine != NULL) {
    rdata->counter = 0;
    rdata->n_new = 0;
    rdata->pattern = pattern;
    rdata->crefine = p3->crefine;
    SC3E (p4est3_iterate_volume
          (p3->old, p4est3_refine_volume_callback, rdata));
    p3->local_num_quads = rdata->counter;
  }
  else if (p3->ccoarse != NULL) {
    cdata->counter = 0;
    cdata->pattern = pattern;
    SC3E (p4est3_refine_array_new
          (p3->alloc, sizeof (void *),
           p3->num_children, p3->num_children, &cdata->family));
    cdata->nsiblings = 0;
    cdata->ccoarse = p3->ccoarse;
    SC3E (p4est3_iterate_volume
          (p3->old, p4est3_coarsen_volume_callback, cdata));
    p3->local_num_quads = cdata->counter;
  }
  else {
    /*in this case we will perform a simple quadrant copying with translation*/
    p3->local_num_quads = p3->old->local_num_quads;
  }

  /* Here we allocate shared p4est3_t::quadwin and p4est3_t::nodequads.
     We will fill in the latter later. */
  SC3E (sc3_allocator_malloc (p3->alloc, nodesize * sizeof (char *),
                              &p3->nodequads));
  SC3E (sc3_MPI_Win_allocate_shared
        ((sc3_MPI_Aint_t) (p3->local_num_quads * p3->qsize), p3->qsize,
         info_noncontig, nodecomm, &quadmem, &p3->quadwin));
  for (i = 0; i < nodesize; ++i) {
    SC3E (sc3_MPI_Win_shared_query (p3->quadwin, i,
                                    &tempbytes, &dispunit, &nqmem));
    SC3A_CHECK (dispunit == p3->qsize);
    SC3A_CHECK (nqmem != NULL || tempbytes == 0);
    p3->nodequads[i] = nqmem;
  }
  p3->quads = quadmem;
  SC3A_CHECK (p3->nodequads[noderank] == p3->quads);

  /* Allocate shared memory for global offsets */
  goffsetbytes = (p3->mpisize + 1) * sizeof (p4est3_gloidx);
  SC3E (sc3_MPI_Win_allocate_shared
        (noderank == 0 ? goffsetbytes : 0, sizeof (p4est3_gloidx),
         info_noncontig, nodecomm, &p3->goffset, &p3->goffsetwin));
  if (noderank > 0) {
    SC3E (sc3_MPI_Win_shared_query (p3->goffsetwin, 0,
                                    &tempbytes, &dispunit, &p3->goffset));
    SC3A_CHECK (tempbytes >= goffsetbytes);
    SC3A_CHECK (dispunit == sizeof (p4est3_gloidx));
    SC3A_CHECK (p3->goffset != NULL);
  }
  /* We got a pattern of population, and now we populate it
     and set p4est3_tree_t:: treeid, quad_offset and num_quads.
     We also initialize p4est3_tree_t::first_tquad by 0.
     We work on the process-local window onto the quadrants */
  SC3E (sc3_allocator_calloc
        (p3->alloc, p3->qvt->dim, sizeof (int32_t), &coords));
  lt_offset = 0;
  for (i = p3->fltree; i <= p3->lltree; ++i) {
    SC3E (p4est3_tree_index (p3, i, &tree));
    tree->treeid = i;
    tree->quad_offset = lt_offset;
    tree->tquads = p3->quads + p3->qsize * tree->quad_offset;
    if (p3->crefine != NULL) {
      SC3E (p4est3_pattern_populate_tree_ref
            (p3, tree, rdata->pattern, &lt_offset, coords));
    }
    else if (p3->ccoarse != NULL) {
      SC3E (p4est3_pattern_populate_tree_coarse
            (p3, tree, cdata->pattern, &lt_offset, coords));
    }
    else {
      /*simply copying*/
      SC3E (p4est3_populate_tree_cpy (p3, tree, &lt_offset, coords));
    }
    tree->num_quads = lt_offset - tree->quad_offset;
    tree->first_tquad = 0;
  }

  /* This check here is only to avoid creating a new mpi datatype. */
  SC3A_CHECK (sizeof (p4est3_locidx) == sizeof (int));
  if (p3->fltree != -1) {
    SC3E (p4est3_tree_index (p3, p3->fltree, &tree));
    num_quads = tree->num_quads;
  }
  SC3E (sc3_MPI_Allgather
        (&p3->local_num_quads, 1, SC3_MPI_INT,
         local_num_quads, 1, SC3_MPI_INT, p3->mpicomm));
  SC3E (sc3_MPI_Allgather
        (&num_quads, 1, SC3_MPI_INT,
         first_tree_quads, 1, SC3_MPI_INT, p3->mpicomm));
  SC3E (sc3_MPI_Win_lock (SC3_MPI_LOCK_SHARED, 0, SC3_MPI_MODE_NOCHECK,
                          p3->goffsetwin));
  if (noderank == 0) {
    p3->goffset[0] = 0;
    for (i = 1; i < p3->mpisize + 1; ++i) {
      p3->goffset[i] = p3->goffset[i - 1]
        + (p4est3_gloidx) local_num_quads[i - 1];
    }
  }
  SC3E (sc3_MPI_Win_unlock (0, p3->goffsetwin));

  if (p3->mpirank != 0) {
    if (p3->fltree != -1) {
      for (i = p3->gftree[p3->mpirank - 1]; i == p3->fltree; --i) {
        tree->first_tquad += first_tree_quads[i];
      }
    }
  }
  for (i = p3->fltree; i <= p3->lltree; ++i) {
    SC3E (p4est3_tree_index (p3, i, &tree));
    tree->end_tquad = tree->first_tquad + tree->num_quads;
    tree->last_tquad =
      (tree->end_tquad = tree->first_tquad + tree->num_quads) - 1;
  }

  SC3E (sc3_array_destroy (&pattern));
  SC3E (sc3_allocator_free (p3->alloc, local_num_quads));
  SC3E (sc3_allocator_free (p3->alloc, first_tree_quads));
  SC3E (sc3_allocator_free (p3->alloc, coords));
  if (p3->crefine == NULL && p3->ccoarse != NULL) {
    SC3E (sc3_array_destroy (&cdata->family));
  }
  sc3_MPI_Barrier (p3->mpicomm);
  p3->global_num_quads = p3->goffset[p3->mpisize];
  return NULL;
}

#ifdef __cplusplus
#if 0
{
#endif
}
#endif
