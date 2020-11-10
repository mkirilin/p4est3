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
#include <p4est3_p4est.h>
#else
#include <p4est3_p8est.h>
#endif

typedef struct p4est3_connectivity_ntslf
{
  sc3_allocator_t    *alloc;

  /* no need really to allocate this here since it is deep copied */
  p4est3_connectivity_vtable_t scvt;
}
p4est3_connectivity_ntslf_t;

static sc3_error_t *
p4est3_connectivity_dstr (void *vslf)
{
  p4est3_connectivity_ntslf_t *slf = (p4est3_connectivity_ntslf_t *) vslf;
  SC3A_CHECK (slf != NULL);

  SC3E (sc3_allocator_free (slf->alloc, slf));
  return NULL;
}

static sc3_error_t *
basics_connectivity_new_virtual (sc3_allocator_t * alloc,
                                 p4est3_topidx num_trees,
                                 p4est3_connectivity_t ** pc)
{
  p4est3_connectivity_ntslf_t *slf;
  p4est3_connectivity_t *c;

  SC3E_RETVAL (pc, NULL);
  SC3A_IS (sc3_allocator_is_valid, alloc);
  SC3A_CHECK (num_trees > 0);

  /* create virtual structure */
  SC3E (sc3_allocator_calloc_one
        (alloc, sizeof (p4est3_connectivity_ntslf_t), &slf));
  slf->alloc = alloc;           /**< no need to ref alloc since conn_new does it */
  slf->scvt.dim = P4EST_DIM;
  slf->scvt.num_trees = num_trees;
  slf->scvt.destroy = p4est3_connectivity_dstr;

  /* create connectivity */
  SC3E (p4est3_connectivity_new (alloc, &c));
  SC3E (p4est3_connectivity_set_vtable (c, &slf->scvt, slf));
  SC3E (p4est3_connectivity_setup (c));
  SC3A_IS (p4est3_connectivity_is_setup, c);

  *pc = c;
  return NULL;
}

static sc3_error_t *
make_allocator (sc3_allocator_t * oa, sc3_allocator_t ** alloc)
{
  SC3A_IS (sc3_allocator_is_setup, oa);
  SC3E (sc3_allocator_new (oa, alloc));
  SC3E (sc3_allocator_setup (*alloc));
  return NULL;
}

static sc3_error_t *
test_p4est_new (sc3_allocator_t * alloc,
                sc3_MPI_Comm_t mpicomm, p4est3_quadrant_vtable_t * qvt,
                p4est3_topidx num_trees, int level)
{
  int                 i;
  p4est_connectivity_t *c4;
  p4est3_connectivity_t *conn;
  p4est3_t           *p3;

  SC3A_IS (sc3_allocator_is_setup, alloc);

  for (i = 0; i < 4; ++i) {
    /* create connectivity structure */
    fprintf (stderr, "Trying %d\n", i);
    switch (i) {
    case 0:
      /* default connectivity with one tree */
      SC3E (p4est3_connectivity_new (alloc, &conn));
      SC3E (p4est3_connectivity_set_dim (conn, P4EST_DIM));
      SC3E (p4est3_connectivity_set_num_trees (conn, num_trees));
      SC3E (p4est3_connectivity_setup (conn));
      break;
    case 1:
      /* virtual connectivity with one tree */
      SC3E (basics_connectivity_new_virtual (alloc, num_trees, &conn));
      break;
    case 2:
      /* use convenience constructors */
#ifndef P4_TO_P8
      SC3E (p4est3_connectivity_new_unitsquare (alloc, &conn));
#else
      SC3E (p4est3_connectivity_new_unitcube (alloc, &conn));
#endif
      break;
    case 3:
      /* wrapping a p4est connectivity */
#ifndef P4_TO_P8
      c4 = p4est_connectivity_new_unitsquare ();
#else
      c4 = p8est_connectivity_new_unitcube ();
#endif
      SC3E (p4est3_connectivity_new_p4est (alloc, c4, 1, &conn));
      break;
    default:
      SC3E_UNREACH ("Invalid example counter");
    }

    /* create p4est object with connectivity */
    SC3E (p4est3_new (alloc, &p3));
    SC3E (p4est3_set_comm (p3, mpicomm, 1));
    SC3E (p4est3_set_connectivity (p3, conn));
    SC3E (p4est3_set_quadrant_vtable (p3, qvt));
    SC3E (p4est3_set_level (p3, level));
    SC3E (p4est3_setup (p3));

    SC3E (p4est3_destroy (&p3));
    SC3E (p4est3_connectivity_destroy (&conn));
  }
  return NULL;
}

static sc3_error_t *
free_allocator (sc3_allocator_t ** alloc)
{
  SC3A_IS (sc3_allocator_is_setup, *alloc);
  SC3E (sc3_allocator_destroy (alloc));
  return NULL;
}

typedef struct v3basics
{
  sc3_MPI_Comm_t      mpicomm;
  int                 mpirank;
  sc3_allocator_t    *alloc;
  p4est3_quadrant_vtable_t sqvt_legacy, *qvt_legacy;

}
v3basics_t;

static sc3_error_t *
v3basics_prepare (v3basics_t * t)
{
  /* consistency checks */
  SC3A_CHECK (t != NULL);

  /* initialize global data */
  t->qvt_legacy = &t->sqvt_legacy;
  t->mpicomm = SC3_MPI_COMM_WORLD;
  SC3E (sc3_MPI_Comm_rank (t->mpicomm, &t->mpirank));
  SC3E (sc3_MPI_Comm_set_errhandler (t->mpicomm, SC3_MPI_ERRORS_RETURN));

  /* we don't need init calls for v3.  Just to check legacy wrapping */
  /* must not use SC3_MPI_COMM_WORLD due to incompatible non-mpi wrapping */
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  /* legacy wrapping for p4est quadrants */
  SC3E (p4est3_quadrant_vtable_p4est (t->qvt_legacy, 0));

  /* perspectively make one allocator for each thread */
  SC3E (make_allocator (sc3_allocator_nothread (), &t->alloc));
  return NULL;
}

static sc3_error_t *
v3basics_run (v3basics_t * t, p4est3_topidx num_trees, int level)
{
  SC3A_CHECK (t != NULL);
  SC3A_CHECK (num_trees > 0);

  SC3E (test_p4est_new
        (t->alloc, t->mpicomm, t->qvt_legacy, num_trees, level));
  return NULL;
}

static sc3_error_t *
v3basics_cleanup (v3basics_t * t)
{
  SC3A_CHECK (t != NULL);

  /* free resources allocated earlier */
  SC3E (free_allocator (&t->alloc));

  /* again, just to check legacy wrapping */
  SC3E_DEMAND (sc_finalize_noabort () == 0, "Legacy sc_finalize failed");
  return NULL;
}

/* It is generally a nice idea to make the error status collective.
   However, this will not work since MPI state prior to entering here
   may be inconsistent between ranks due to rank-specific error history. */
static int
v3basics_error_check (v3basics_t * t, sc3_error_t ** e)
{
  int                 retval;
  char                buffer[SC3_BUFSIZE];

  retval = sc3_error_check (e, buffer, SC3_BUFSIZE);
  if (retval) {
    fprintf (stderr, "Internal error to program on rank %d:\n"
             "%s\nThis rank %d will skip the rest.\n",
             t->mpirank, buffer, t->mpirank);
  }
  return retval < 0;
}

static void
v3basics_error_summary (v3basics_t * t, int arewedead)
{
  if (arewedead) {
    fprintf (stderr, "Ended rank %d on error.\n", t->mpirank);
    SC3X (sc3_MPI_Abort (t->mpicomm, SC3_MPI_ERR_OTHER));
  }
}

int
main (int argc, char **argv)
{
  int                 level;
  int                 arewedead;
  p4est3_topidx       num_trees;
  sc3_error_t        *e;
  v3basics_t          st, *t = &st;

  /* Generally needed for MPI.  No room for continuing on error. */
  SC3X (sc3_MPI_Init (&argc, &argv));

  /* wanna-be command line parameters */
  num_trees = 2;
  level = 3;

  /*** The way of using p4est3 in the following is one suggestion.
       Application may use shortcuts and crash on error, but here
       we try to report error conditions cleanly to calling code. ***/

  /* setup data structures to use */
  e = v3basics_prepare (t);
  arewedead = v3basics_error_check (t, &e);

  /* do something.  Supposing the use of p4est3 is part of a bigger program */
  if (!arewedead) {
    e = v3basics_run (t, num_trees, level);
    arewedead = v3basics_error_check (t, &e);
  }

  /* we will not try to cleanup if we must assume fatal inconsistencies */
  if (!arewedead) {
    e = v3basics_cleanup (t);
    arewedead = v3basics_error_check (t, &e);
  }

  /* print summary and abort if we have encountered a fatal inconsistency */
  v3basics_error_summary (t, arewedead);

  /* Generally needed for MPI.  No room for continuing on error. */
  SC3X (sc3_MPI_Finalize ());
  return 0;
}
