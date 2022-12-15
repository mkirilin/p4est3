/*
  This file is part of p4est, version 3
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
#include <sc3_refcount.h>

int
p4est3_vtable_is_valid (const p4est3_vtable_t * pvt, char *reason)
{
  int                 cdim;

  SC3E_TEST (pvt != NULL, reason);
  SC3E_TEST (0 < pvt->dim && pvt->dim <= 3, reason);
  SC3E_TEST (pvt->mpicomm != SC3_MPI_COMM_NULL, reason);

  /* check connectivity and retrieve dimension */
  SC3E_IS (p4est3_connectivity_is_valid, pvt->c3, reason);
  SC3E_DO (p4est3_connectivity_get_dim (pvt->c3, &cdim), reason);

  /* internal consistency */
  SC3E_TEST (pvt->get_local_num_trees != NULL, reason);

  /* check quadrant table and compare dimension */
  SC3E_IS (p4est3_quadrant_vtable_is_valid, pvt->qvt, reason);
  SC3E_TEST (cdim == pvt->qvt->dim, reason);

  /* successfully done */
  SC3E_YES (reason);
}

int
p4est3_is_valid (const p4est3_t * p3, char *reason)
{
  SC3E_TEST (p3 != NULL, reason);
  SC3E_IS (sc3_refcount_is_valid, &p3->rc, reason);
  SC3E_IS (sc3_allocator_is_setup, p3->alloc, reason);

  if (!p3->setup) {
    SC3E_TEST (p3->accessed_conn == 0, reason);
    SC3E_TEST (p3->split_info == NULL, reason);
  }
  else {
    SC3E_TEST (p3->accessed_conn >= 0, reason);
    SC3E_IS (sc3_mpienv_is_valid, p3->split_info, reason);
    SC3E_TEST (p3->old == NULL, reason);
    SC3E_TEST (p3->crefine == NULL && p3->ccoarse == NULL, reason);
  }

  /* TODO check communicator and connectivity members */

  if (p3->pvt != NULL) {
    SC3E_IS (p4est3_vtable_is_valid, p3->pvt, reason);

    /* TODO check whatever else happens with a forest virtual table */
  }
  else {
    SC3E_TEST (p3->mpicomm != SC3_MPI_COMM_NULL, reason);
    SC3E_TEST (p3->level >= 0, reason);
    SC3E_TEST (p3->setup_mode < P4EST3_NEW_MODE_LAST, reason);

    if (!p3->setup) {
      SC3E_TEST (p3->mpisize == 0 && p3->mpirank == 0, reason);
    }
    else {
      SC3E_IS (p4est3_connectivity_is_setup, p3->conn, reason);
      SC3E_TEST (p3->num_trees > 0, reason);
      SC3E_TEST (p3->qvt != NULL, reason);

      /* TODO thoroughly test all member variables */
    }
  }

  SC3E_YES (reason);
}

int
p4est3_is_new (const p4est3_t * p3, char *reason)
{
  SC3E_IS (p4est3_is_valid, p3, reason);
  SC3E_TEST (!p3->setup, reason);
  SC3E_YES (reason);
}

int
p4est3_is_setup (const p4est3_t * p3, char *reason)
{
  SC3E_IS (p4est3_is_valid, p3, reason);
  SC3E_TEST (p3->setup, reason);
  SC3E_YES (reason);
}

sc3_error_t        *
p4est3_new (sc3_allocator_t * alloc, p4est3_t ** pp3)
{
  p4est3_t           *p3;

  SC3E_RETVAL (pp3, NULL);
  SC3A_IS (sc3_allocator_is_setup, alloc);

  SC3E (sc3_allocator_ref (alloc));
  SC3E (sc3_allocator_calloc_one (alloc, sizeof (p4est3_t), &p3));
  SC3E (sc3_refcount_init (&p3->rc));
  p3->alloc = alloc;
  p3->mpicomm = SC3_MPI_COMM_WORLD;
  p3->setup_mode = P4EST3_NEW_MORTON;
  p3->shared = 0;
  p3->family = 0;
  p3->partition = 0;
  SC3A_IS (p4est3_is_new, p3);

  *pp3 = p3;
  return NULL;
}

sc3_error_t        *
p4est3_set_vtable (p4est3_t * p3, p4est3_vtable_t * pvt, void *slf)
{
  SC3A_IS (p4est3_is_new, p3);
  SC3A_IS (p4est3_vtable_is_valid, pvt);

  /* make shallow copy of virtual table */
  *(p3->pvt = &p3->spvt) = *pvt;

  /* we make a dupliacte of the communicator */
  SC3E (p4est3_set_comm (p3, p3->pvt->mpicomm, 1));

  /* use objects passed in table now and overwrite them with our own */
  SC3E (p4est3_set_connectivity (p3, p3->pvt->c3));
  SC3E (p4est3_set_quadrant_vtable (p3, p3->pvt->qvt));
  p3->pvt->c3 = p3->conn;
  p3->pvt->qvt = p3->qvt;

  /* assign virtual context */
  p3->slf = slf;
  return NULL;
}

sc3_error_t        *
p4est3_set_comm (p4est3_t * p3, sc3_MPI_Comm_t comm, int dup)
{
  SC3A_IS (p4est3_is_new, p3);
  SC3A_CHECK (comm != SC3_MPI_COMM_NULL);

  /* remove previous communicator */
  if (p3->commdup) {
    SC3E (sc3_MPI_Comm_free (&p3->mpicomm));
  }

  /* register new communicator */
  if (dup) {
    SC3E (sc3_MPI_Comm_dup (comm, &p3->mpicomm));
    SC3E (sc3_MPI_Comm_set_errhandler (p3->mpicomm, SC3_MPI_ERRORS_RETURN));
  }
  else {
    p3->mpicomm = comm;
  }
  p3->commdup = dup;
  return NULL;
}

sc3_error_t        *
p4est3_set_connectivity (p4est3_t * p3, p4est3_connectivity_t * conn)
{
  SC3A_IS (p4est3_is_new, p3);
  SC3A_IS (p4est3_connectivity_is_setup, conn);

  if (p3->conn != NULL) {
    SC3E (p4est3_connectivity_unref (p3->conn));
  }
  p3->conn = conn;
  SC3E (p4est3_connectivity_ref (p3->conn));

  /* query connectivity for number of trees */
  SC3E (p4est3_connectivity_get_num_trees (p3->conn, &p3->num_trees));
  return NULL;
}

sc3_error_t        *
p4est3_set_quadrant_vtable (p4est3_t * p3, const p4est3_quadrant_vtable_t * qvt)
{
  SC3A_IS (p4est3_is_new, p3);
  SC3A_CHECK (qvt != NULL);

  /* make shallow copy of virtual table */
  p3->qvt = qvt;
  return NULL;
}

sc3_error_t        *
p4est3_set_level (p4est3_t * p3, int level)
{
  SC3A_IS (p4est3_is_new, p3);
  SC3A_CHECK (level >= 0);

  p3->level = level;
  return NULL;
}

static sc3_error_t *
p4est3_setup_vtable (p4est3_t * p3)
{
  /* TODO: make MPI communicator wrappers of sc and sc3 compatible */
  /* TODO: set as many p3 member variables as makes sense */

  return NULL;
}

sc3_error_t        *
p4est3_set_setup_mode (p4est3_t * p3, p4est3_setup_mode_t mode)
{
  SC3A_IS (p4est3_is_new, p3);
  SC3A_CHECK (0 <= mode && mode < P4EST3_NEW_MODE_LAST);

  p3->setup_mode = mode;
  return NULL;
}

sc3_error_t        *
p4est3_set_source (p4est3_t * p3, p4est3_t * old)
{
  SC3A_IS (p4est3_is_new, p3);
  SC3A_IS (p4est3_is_setup, old);
  SC3A_CHECK (p3 != old);

  if (p3->old != NULL) {
    SC3E (p4est3_unref (p3->old));
  }
  p3->old = old;
  SC3E (p4est3_ref (p3->old));

  return NULL;
}

sc3_error_t        *
p4est3_set_refine (p4est3_t * p3, p4est3_refine_callback_t crefine)
{
  SC3A_IS (p4est3_is_new, p3);
  p3->crefine = crefine;

  return NULL;
}

sc3_error_t        *
p4est3_set_coarsen (p4est3_t * p3, p4est3_coarsen_callback_t ccoarse)
{
  SC3A_IS (p4est3_is_new, p3);
  p3->ccoarse = ccoarse;

  return NULL;
}

sc3_error_t        *
p4est3_set_shared (p4est3_t * p3, int shared)
{
  SC3A_IS (p4est3_is_new, p3);
  p3->shared = shared;
  return NULL;
}

sc3_error_t        *
p4est3_set_family (p4est3_t * p3, int is_family)
{
  SC3A_IS (p4est3_is_new, p3);
  p3->family = is_family;

  return NULL;
}

sc3_error_t        *
p4est3_set_partition (p4est3_t * p3, int partition)
{
  SC3A_IS (p4est3_is_new, p3);
  p3->partition = partition;
  return NULL;
}

sc3_error_t        *
p4est3_setup (p4est3_t * p3)
{
  int                 cdim;
  int                 lev;
  int                 qsize;
  int                 ti;
  p4est3_gloidx       num_uniform, high_uniform;

  /*
   * We will extend the functionality of p4est3_set_* and _setup in the future.
   * They will be used to create a new forest by refinement, partitioning, etc.
   * The forest created is generally immutable after being setup.
   * Currently, we are just creating a uniformly refined forest.
   */

  SC3A_IS (p4est3_is_new, p3);

  /* Check conditions that arise due to omitting mandatory _set_ functions.
     Note that p4est3_set_vtable sets connectivity and quadrant vtable. */
  if (p3->old == NULL) {
    SC3E_DEMAND (p3->conn != NULL, "Connectivity must be set");
    SC3E_DEMAND (p3->qvt != NULL, "Quadrant virtual table must be set");
    SC3E (p4est3_connectivity_get_dim (p3->conn, &cdim));
    SC3E_DEMAND (cdim == p3->qvt->dim,
                 "Dimensions of connectivity and quadrant vtable must match");
    /* further pre-setup consistency checks */
    SC3A_CHECK (p3->num_trees > 0);
  }

  if (p3->pvt != NULL) {
    SC3E (p4est3_setup_vtable (p3));
  }
  if (p3->old != NULL) {
    SC3E (p4est3_internal_setup_from_source (p3));
    SC3E (p4est3_unref (p3->old));
    p3->old = NULL;
    p3->crefine = NULL;
    p3->ccoarse = NULL;
  }
  else {
    /* query input communicator and populate node and head communicators */
    SC3E (p4est3_internal_setup_comm (p3));

    /* determine a quadrant's size in memory */
    qsize = (int) p4est3_quadrant_size (p3->qvt);
    SC3A_CHECK (qsize > 0);

    /* determine principal and initial refinement level */
    p3->qmaxlevel = p4est3_quadrant_max_level (p3->qvt);
    SC3A_CHECK (p3->qmaxlevel >= 0);
    p3->level = SC3_MIN (p3->level, p3->qmaxlevel);

    /* TODO make sure that the number of local quadrants stays bounded */

    /* with number of children determine number of elements per tree */
    /* TODO use uniform_level function and consider variable num_children */
    p3->num_children = p4est3_quadrant_num_children (p3->qvt);
    SC3A_CHECK (p3->num_children > 0);
    high_uniform = P4EST3_GLOIDX_MAX / p3->num_children;
    for (num_uniform = 1, lev = 0;
         num_uniform <= high_uniform && lev < p3->level; ++lev) {
      /* we iterate so we do not roll over the gloidx limit */
      num_uniform *= p3->num_children;
    }
    p3->level = lev;
    SC3A_CHECK (p4est3_glopow (p3->num_children, p3->level) == num_uniform);

    /* allocate one temporary quadrant per thread */
    p3->max_threads = sc3_omp_max_threads ();
    SC3E (sc3_allocator_malloc (p3->alloc, p3->max_threads * sizeof (char *),
                                &p3->temp_quad));
    for (ti = 0; ti < p3->max_threads; ++ti) {
      SC3E (sc3_allocator_malloc (p3->alloc, qsize, &p3->temp_quad[ti]));
    }

    /* compute partition cuts and create shared partition arrays */
    SC3E (p4est3_internal_setup_cut (p3, num_uniform, qsize));

    /* create tree and quadrant metadata */
    SC3E (p4est3_internal_setup_tree (p3, num_uniform));

    /* create quadrants by the previously specified method */
    SC3E (p4est3_internal_setup_quadrants (p3));
  }

  /* we are done creating a valid forest */
  p3->setup = 1;
  SC3A_IS (p4est3_is_setup, p3);
  return NULL;
}

sc3_error_t        *
p4est3_ref (p4est3_t * p3)
{
  SC3A_IS (p4est3_is_setup, p3);
  SC3E (sc3_refcount_ref (&p3->rc));
  return NULL;
}

sc3_error_t        *
p4est3_unref (p4est3_t * p3)
{
  int                 waslast;

  SC3A_IS (p4est3_is_setup, p3);
  SC3E (sc3_refcount_unref (&p3->rc, &waslast));

  /* This is a hard check that we do not unref below a count of one. */
  SC3E_DEMAND (!waslast, "Forest unrefd below a count of one");
  return NULL;
}

sc3_error_t        *
p4est3_destroy (p4est3_t ** pp3)
{
  sc3_allocator_t    *alloc;
  p4est3_t           *p3;
  sc3_MPI_Win_t       gftreewin, gfposwin, goffsetwin;


  SC3E_INULLP (pp3, p3);
  SC3A_IS (p4est3_is_valid, p3);
  SC3A_CHECK (p3->accessed_conn == 0);

  /* This is a hard check for a reference count of exactly one. */
  SC3E_DEMIS (sc3_refcount_is_last, &p3->rc);

  if (p3->pvt != NULL) {
    /* the connectivity has been set from the virtual table */
    SC3E (p4est3_connectivity_unref (p3->conn));

    /* destruction callback if one was provided */
    if (p3->pvt->destroy != NULL) {
      SC3E (p3->pvt->destroy (p3->slf));
    }
  }
  else {
    /* free memory that has been populated during non-virtual setup */
    if (p3->setup) {
      int                 ti;

      /* free internal MPI objects */
      SC3E (p4est3_get_gftreewin (p3, &gftreewin));
      SC3E (p4est3_get_gfposwin (p3, &gfposwin));
      SC3E (p4est3_get_goffsetwin (p3, &goffsetwin));
      SC3E (sc3_MPI_Win_free (gftreewin));
      SC3E (sc3_MPI_Win_free (gfposwin));
      SC3E (sc3_MPI_Win_free (goffsetwin));
      SC3E (sc3_MPI_Win_free (&p3->quadwin));

      /* deallocate internal storage */
      for (ti = 0; ti < p3->max_threads; ++ti) {
        SC3E (sc3_allocator_free (p3->alloc, p3->temp_quad[ti]));
      }
      SC3E (sc3_allocator_free (p3->alloc, p3->temp_quad));

      SC3E (sc3_array_destroy (&p3->trees));
      SC3E (sc3_allocator_free (p3->alloc, p3->nodequads));
    }

    /* release data that has been referenced before setup */
    if (p3->conn != NULL) {
      SC3E (p4est3_connectivity_unref (p3->conn));
    }
    if (p3->commdup) {
      SC3E (sc3_MPI_Comm_free (&p3->mpicomm));
    }
    /* unref mpi environment related data */
    SC3E (sc3_mpienv_unref (&p3->split_info));
    if (p3->old != NULL) {
      SC3E (p4est3_unref (p3->old));
    }

    /* TODO: think about freeing all setup data also for virtual forest */
  }

  /* remove allocation */
  alloc = p3->alloc;
  SC3E (sc3_allocator_free (alloc, p3));
  SC3E (sc3_allocator_unref (&alloc));
  return NULL;
}

sc3_error_t        *
p4est3_access_connectivity (p4est3_t * p3, p4est3_connectivity_t ** pconn)
{
  SC3E_RETVAL (pconn, NULL);
  SC3A_IS (p4est3_is_setup, p3);

  SC3E (p4est3_connectivity_ref (p3->conn));
  *pconn = p3->conn;
  ++p3->accessed_conn;
  return NULL;
}

sc3_error_t        *
p4est3_restore_connectivity (p4est3_t * p3, p4est3_connectivity_t * conn)
{
  SC3A_IS (p4est3_is_setup, p3);
  SC3A_CHECK (conn == p3->conn);
  SC3A_CHECK (p3->accessed_conn > 0);

  SC3E (p4est3_connectivity_unref (p3->conn));
  --p3->accessed_conn;
  return NULL;
}

sc3_error_t        *
p4est3_get_local_num_trees (const p4est3_t * p3,
                            p4est3_topidx * first_local_tree,
                            p4est3_topidx * last_local_tree)
{
  SC3A_IS (p4est3_is_setup, p3);
  SC3A_CHECK (first_local_tree != NULL);
  SC3A_CHECK (last_local_tree != NULL);

  if (p3->pvt != NULL) {
    SC3A_CHECK (p3->pvt->get_local_num_trees != NULL);
    return p3->pvt->get_local_num_trees (p3->slf,
                                         first_local_tree, last_local_tree);
  }
  else {
    *first_local_tree = p3->fltree;
    *last_local_tree = p3->lltree;
    return NULL;
  }
}

sc3_error_t        *
p4est3_get_quadrants (const p4est3_t * p3, char **q)
{
  SC3A_IS (p4est3_is_setup, p3);

  *q = p3->quads;

  return NULL;
}

sc3_error_t        *
p4est3_get_global_num_quads (const p4est3_t * p3, p4est3_gloidx * n)
{
  if (n != NULL) {
    *n = 0L;
  }
  SC3A_IS (p4est3_is_setup, p3);

  if (n != NULL) {
    *n = p3->global_num_quads;
  }
  return NULL;
}

sc3_error_t        *
p4est3_get_local_num_quads (const p4est3_t * p3, p4est3_locidx * n)
{
  if (n != NULL) {
    *n = 0L;
  }
  SC3A_IS (p4est3_is_setup, p3);

  if (n != NULL) {
    *n = p3->local_num_quads;
  }
  return NULL;
}

sc3_error_t        *
p4est3_get_gftreewin (const p4est3_t *p3, sc3_MPI_Win_t *tw)
{
  SC3A_IS (p4est3_is_setup, p3);
  SC3A_CHECK (tw != NULL);

  *tw = p3->gtrees->gftreewin;
  return NULL;
}

sc3_error_t        *
p4est3_get_gfposwin (const p4est3_t *p3, sc3_MPI_Win_t *pw)
{
  SC3A_IS (p4est3_is_setup, p3);
  SC3A_CHECK (pw != NULL);

  *pw = p3->gposition->gfposwin;
  return NULL;
}

sc3_error_t        *
p4est3_get_goffsetwin (const p4est3_t *p3, sc3_MPI_Win_t *ow)
{
  SC3A_IS (p4est3_is_setup, p3);
  SC3A_CHECK (ow != NULL);

  *ow = p3->goffsets->goffsetwin;
  return NULL;
}