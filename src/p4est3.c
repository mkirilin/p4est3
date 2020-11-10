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
  SC3E_TEST (pvt != NULL, reason);
  SC3E_TEST (0 < pvt->dim && pvt->dim <= 3, reason);
  SC3E_YES (reason);
}

int
p4est3_is_valid (const p4est3_t * p3, char *reason)
{
  SC3E_TEST (p3 != NULL, reason);
  SC3E_IS (sc3_refcount_is_valid, &p3->rc, reason);
  SC3E_IS (sc3_allocator_is_setup, p3->alloc, reason);

  if (p3->pvt != NULL) {
    SC3E_IS (p4est3_vtable_is_valid, p3->pvt, reason);

    /* TODO check whatever else happens with a forest virtual table */
  }
  else {
    SC3E_TEST (p3->mpicomm != SC3_MPI_COMM_NULL, reason);
    SC3E_TEST (p3->level >= 0, reason);

    if (!p3->setup) {
      SC3E_TEST (p3->mpisize == 0 && p3->mpirank == 0, reason);
    }
    else {
      SC3E_IS (p4est3_connectivity_is_setup, p3->conn, reason);
      SC3E_TEST (p3->num_trees > 0, reason);
      SC3E_TEST (p3->qvt == &p3->sqvt, reason);

      SC3E_TEST (p3->nodesizewin != SC3_MPI_WIN_NULL, reason);
      SC3E_TEST (p3->headcomm != SC3_MPI_COMM_NULL
                 || p3->noderank > 0, reason);
      SC3E_TEST (p3->nodecomm != SC3_MPI_COMM_NULL, reason);

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
  p3->nodesizewin = SC3_MPI_WIN_NULL;
  p3->headcomm = SC3_MPI_COMM_NULL;
  p3->nodecomm = SC3_MPI_COMM_NULL;
  SC3A_IS (p4est3_is_new, p3);

  *pp3 = p3;
  return NULL;
}

sc3_error_t        *
p4est3_set_vtable (p4est3_t * p3, p4est3_vtable_t * pvt, void *slf)
{
  SC3A_IS (p4est3_is_new, p3);
  SC3A_IS (p4est3_vtable_is_valid, pvt);

  /* make deep copy of virtual table */
  *(p3->pvt = &p3->spvt) = *pvt;
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
p4est3_set_quadrant_vtable (p4est3_t * p3, p4est3_quadrant_vtable_t * qvt)
{
  SC3A_IS (p4est3_is_new, p3);
  SC3A_CHECK (qvt != NULL);

  /* make deep copy of virtual table */
  *(p3->qvt = &p3->sqvt) = *qvt;
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
p4est3_setup (p4est3_t * p3)
{
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

  if (p3->pvt != NULL) {
    SC3E (p4est3_setup_vtable (p3));
  }
  else {
    /* check conditions that arise due to omitting mandatory _set_ functions */
    SC3E_DEMAND (p3->conn != NULL, "Connectivity must be set");
    SC3E_DEMAND (p3->qvt == &p3->sqvt, "Quadrant virtual table must be set");

    /* further pre-setup consistency checks */
    SC3A_CHECK (p3->num_trees > 0);

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
    p3->num_children = p4est3_quadrant_max_children (p3->qvt);
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

    /* create quadrants by the morton method, which is presumably slowest */
    SC3E (p4est3_internal_setup_morton (p3));
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

  SC3E_INULLP (pp3, p3);
  SC3A_IS (p4est3_is_valid, p3);

  /* This is a hard check for a reference count of exactly one. */
  SC3E_DEMIS (sc3_refcount_is_last, &p3->rc);

  /* destruction callback if one was provided */
  if (p3->pvt != NULL && p3->pvt->destroy != NULL) {
    SC3E (p3->pvt->destroy (p3->slf));
  }

  /* free memory that has been populated during non-virtual setup */
  if (p3->pvt == NULL) {
    if (p3->setup) {
      int                 ti;

      /* free internal MPI objects */
      SC3E (sc3_MPI_Win_free (&p3->nodesizewin));
      SC3E (sc3_MPI_Win_free (&p3->gfposwin));
      SC3E (sc3_MPI_Win_free (&p3->gftreewin));
      SC3E (sc3_MPI_Win_free (&p3->goffsetwin));
      SC3E (sc3_MPI_Win_free (&p3->quadwin));
      if (p3->noderank == 0) {
        SC3E (sc3_MPI_Comm_free (&p3->headcomm));
      }
      SC3E (sc3_MPI_Comm_free (&p3->nodecomm));
      SC3E (sc3_MPI_Info_free (&p3->info_noncontig));

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
  }

  /* remove allocation */
  alloc = p3->alloc;
  SC3E (sc3_allocator_free (alloc, p3));
  SC3E (sc3_allocator_unref (&alloc));
  return NULL;
}
