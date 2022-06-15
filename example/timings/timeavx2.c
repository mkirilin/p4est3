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
#include <p4est3_quadrant_yx.h>
#include <p4est3_p4est.h>
#else
#include <p4est3_quadrant_zyx.h>
#include <p4est3_p8est.h>
#endif

static sc3_error_t *
test_child (sc3_array_t * pull, p4est3_quadrant_vtable_t * qvt,
            p4est3_locidx n_quads, double *exec_time)
{
  p4est3_locidx       quad, put_ind;
  int                 child;
  void               *p, *q;
  double              t_b, t_e;

  SC3A_IS (sc3_array_is_valid, pull);
  SC3A_CHECK (qvt != NULL);
  SC3A_CHECK (exec_time != NULL);

  t_b = sc3_MPI_Wtime ();
  for (quad = 0, put_ind = 1; P4EST_CHILDREN * (quad + 1) < n_quads;
       ++quad, put_ind += P4EST_CHILDREN) {
    for (child = 0; child < P4EST_CHILDREN; ++child) {
      SC3E (sc3_array_index (pull, quad, &p));
      SC3E (sc3_array_index (pull, put_ind + child, &q));
      SC3E (p4est3_quadrant_child (qvt, p, child, q));
    }
  }

  for (child = 0; child < n_quads - put_ind; ++child) {
    SC3E (sc3_array_index (pull, quad, &p));
    SC3E (sc3_array_index (pull, put_ind + child, &q));
    SC3E (p4est3_quadrant_child (qvt, p, child, q));
  }
  t_e = sc3_MPI_Wtime ();

  *exec_time = t_e - t_b;
  return NULL;
}

static sc3_error_t *
test_parent (void *q, sc3_array_t * pull, p4est3_quadrant_vtable_t * qvt,
             p4est3_locidx n_quads, double *exec_time)
{
  p4est3_locidx       quad;
  void               *p;
  double              t_b, t_e;

  SC3A_CHECK (q != NULL);
  SC3A_IS (sc3_array_is_valid, pull);
  SC3A_CHECK (qvt != NULL);
  SC3A_CHECK (exec_time != NULL);

  t_b = sc3_MPI_Wtime ();
  for (quad = 1; quad < n_quads; ++quad) {
    SC3E (sc3_array_index (pull, quad, &p));
    SC3E (p4est3_quadrant_parent (qvt, p, q));
  }
  t_e = sc3_MPI_Wtime ();
  *exec_time = t_e - t_b;
  return NULL;
}

static sc3_error_t *
test_compare (sc3_array_t * pull, p4est3_quadrant_vtable_t * qvt,
              p4est3_locidx n_quads, double *exec_time)
{
  p4est3_locidx       quad;
  int                 j;
  void               *p, *q;
  double              t_b, t_e;

  SC3A_IS (sc3_array_is_valid, pull);
  SC3A_CHECK (qvt != NULL);
  SC3A_CHECK (exec_time != NULL);

  t_b = sc3_MPI_Wtime ();
  for (quad = 0; quad < n_quads; ++quad) {
    SC3E (sc3_array_index (pull, quad, &p));
    SC3E (sc3_array_index (pull, n_quads - quad - 1, &q));
    SC3E (p4est3_quadrant_compare (qvt, p, q, &j));
  }
  t_e = sc3_MPI_Wtime ();
  *exec_time = t_e - t_b;
  return NULL;
}

static sc3_error_t *
test_successor (void *q, sc3_array_t * pull, p4est3_quadrant_vtable_t * qvt,
                p4est3_locidx n_quads, double *exec_time)
{
  p4est3_locidx       quad;
  int                 i;
  void               *p;
  double              t_b, t_e;

  SC3A_CHECK (q != NULL);
  SC3A_IS (sc3_array_is_valid, pull);
  SC3A_CHECK (qvt != NULL);
  SC3A_CHECK (exec_time != NULL);

  t_b = sc3_MPI_Wtime ();
  for (quad = 1; quad < n_quads - P4EST_CHILDREN; quad += P4EST_CHILDREN) {
    for (i = 0; i < P4EST_CHILDREN - 1; ++i) {
      SC3E (sc3_array_index (pull, quad + i, &p));
      SC3E (p4est3_quadrant_successor (qvt, p, q));
    }
  }
  t_e = sc3_MPI_Wtime ();
  *exec_time = t_e - t_b;
  return NULL;
}

static void
print_time_info (double exec_avx, double exec_nonavx, const char *name)
{
  printf ("  %s: \n"
          "    Vectorized:        %g\n"
          "    Non-Vectorized:    %g\n"
          "    Vect/Non-Vect Ratio:  %g\n",
          name, exec_avx, exec_nonavx,
          exec_nonavx <= 0. ? 0. : exec_avx / exec_nonavx);
}

static sc3_error_t *
measure_child (sc3_array_t * v_pull2check, sc3_array_t * q_pull2check,
               p4est3_quadrant_vtable_t * qvt_avx,
               p4est3_quadrant_vtable_t * qvt, p4est3_locidx n_quads)
{
  double              exec_avx, exec_nonavx;

  SC3E (test_child (v_pull2check, qvt_avx, n_quads, &exec_avx));
  SC3E (test_child (q_pull2check, qvt, n_quads, &exec_nonavx));

  print_time_info (exec_avx, exec_nonavx, "Child");
  return NULL;
}

static sc3_error_t *
measure_parent (sc3_array_t * v_pull2check, sc3_array_t * q_pull2check,
                p4est3_quadrant_vtable_t * qvt_avx,
                p4est3_quadrant_vtable_t * qvt, p4est3_locidx n_quads)
{
  double              exec_avx, exec_nonavx;
  void               *v, *q;

  SC3E (sc3_array_index (v_pull2check, 0, &v));
  SC3E (sc3_array_index (q_pull2check, 0, &q));
  SC3E (test_parent (v, v_pull2check, qvt_avx, n_quads, &exec_avx));
  SC3E (test_parent (q, q_pull2check, qvt, n_quads, &exec_nonavx));

  print_time_info (exec_avx, exec_nonavx, "Parent");
  return NULL;
}

static sc3_error_t *
measure_compare (sc3_array_t * v_pull2check, sc3_array_t * q_pull2check,
                 p4est3_quadrant_vtable_t * qvt_avx,
                 p4est3_quadrant_vtable_t * qvt, p4est3_locidx n_quads)
{
  double              exec_avx, exec_nonavx;

  SC3E (test_compare (v_pull2check, qvt_avx, n_quads, &exec_avx));
  SC3E (test_compare (q_pull2check, qvt, n_quads, &exec_nonavx));

  print_time_info (exec_avx, exec_nonavx, "Compare");
  return NULL;
}

static sc3_error_t *
measure_successor (sc3_array_t * v_pull2check, sc3_array_t * q_pull2check,
                   p4est3_quadrant_vtable_t * qvt_avx,
                   p4est3_quadrant_vtable_t * qvt, p4est3_locidx n_quads)
{
  double              exec_avx, exec_nonavx;
  void               *v, *q;

  SC3E (sc3_array_index (v_pull2check, 0, &v));
  SC3E (sc3_array_index (q_pull2check, 0, &q));
  SC3E (test_successor (v, v_pull2check, qvt_avx, n_quads, &exec_avx));
  SC3E (test_successor (q, q_pull2check, qvt, n_quads, &exec_nonavx));

  print_time_info (exec_avx, exec_nonavx, "Successor");
  return NULL;
}

#if 0

static void
report_errors (sc3_error_t ** pe)
{
  char                eflat[SC3_BUFSIZE];

  if (pe != NULL && *pe != NULL) {
    sc3_error_destroy_noerr (pe, eflat);
    fprintf (stderr, "Error: %s\n", eflat);
  }
}

#endif

typedef struct timeavx2
{
  int                 mpirank;
  p4est3_locidx       n_quads;

  sc3_allocator_t    *alloc;
  sc3_array_t        *qarr;
  sc3_array_t        *qarr_avx;
  p4est3_quadrant_vtable_t sqvt, *qvt;
  const p4est3_quadrant_vtable_t *qvt_avx;
}
timeavx2_t;

static sc3_error_t *
timeavx2_prepare (timeavx2_t * t, int *retval)
{
  void               *p;

  SC3E_RETVAL (retval, -1);
  SC3A_CHECK (t != NULL);
  SC3A_CHECK (t->n_quads > 0);

  /* static initializers */
  t->qvt = &t->sqvt;

  /* the standard p4est2 virtual table always exists */
  p4est3_quadrant_vtable_p4est (t->qvt, 0);

  /* the AVX virtual table can only be set with hardware support */
  SC3E (p4est3_quadrant_yx_vtable (&t->qvt_avx));
  if (t->qvt_avx == NULL) {
    /* AVX is not supported by hardware
      or p4est is not build neither in 2D nor 3D*/
    if (t->mpirank == 0) {
      fprintf (stderr, "%s\nWill not proceed\n",
      "AVX is not supported by hardware "
      "or p4est is not build neither in 2D nor 3D\n");
    }
    /* return value has been initialized to failure above */
    return NULL;
  }

  /* create a toplevel allocator */
  SC3E (sc3_allocator_new (sc3_allocator_nocount (), &t->alloc));
  SC3E (sc3_allocator_setup (t->alloc));

  /* allocate quadrant arrays */
  SC3E (p4est3_quadrant_array_new (t->alloc, t->qvt, t->n_quads, &t->qarr));
  SC3E (p4est3_quadrant_array_new (t->alloc,
                                   t->qvt_avx, t->n_quads, &t->qarr_avx));

  /* initialize first element */
  SC3E (sc3_array_index (t->qarr, 0, &p));
  SC3E (p4est3_quadrant_root (t->qvt, p));
  SC3E (sc3_array_index (t->qarr_avx, 0, &p));
  SC3E (p4est3_quadrant_root (t->qvt_avx, p));

  /* clean and successful return */
  *retval = 0;
  return NULL;
}

static sc3_error_t *
timeavx2_measure (timeavx2_t * t)
{
  SC3A_CHECK (t != NULL);
  SC3A_CHECK (t->n_quads > 0);

  SC3E (measure_child (t->qarr_avx, t->qarr, t->qvt_avx, t->qvt, t->n_quads));
  SC3E (measure_parent
        (t->qarr_avx, t->qarr, t->qvt_avx, t->qvt, t->n_quads));
  SC3E (measure_compare
        (t->qarr_avx, t->qarr, t->qvt_avx, t->qvt, t->n_quads));
  SC3E (measure_successor
        (t->qarr_avx, t->qarr, t->qvt_avx, t->qvt, t->n_quads));
  return NULL;
}

static sc3_error_t *
timeavx2_cleanup (timeavx2_t * t)
{
  SC3A_CHECK (t != NULL);
  SC3A_CHECK (t->n_quads > 0);

  SC3E (sc3_array_destroy (&t->qarr_avx));
  SC3E (sc3_array_destroy (&t->qarr));
  SC3E (sc3_allocator_destroy (&t->alloc));
  return NULL;
}

int
main (int argc, char **argv)
{
  int                 retval;
  p4est3_locidx       n_quads;
#if 0
  p4est3_quadrant_vtable_t sqvt_avx, *qvt_avx = &sqvt_avx;
  p4est3_quadrant_vtable_t sqvt, *qvt = &sqvt;
  sc3_error_t        *e;
  sc3_array_t        *qarr_avx, *qarr;
  void               *p;
#endif
  timeavx2_t          st, *t = &st;

  /* MPI_Init comes first in a program.  We abort should this go wrong. */
  SC3X (sc3_MPI_Init (&argc, &argv));
  SC3X (sc3_MPI_Comm_rank (SC3_MPI_COMM_WORLD, &t->mpirank));

  /* interpret command line arguments */
  if (argc == 1) {
    n_quads = 0;
    for (int i = 0; i < 10; ++i) {
      n_quads += (1 << (i * P4EST_DIM));
    }
  }
  else {
    n_quads = atoll (argv[1]);
  }
  t->n_quads = n_quads = SC3_MAX (n_quads, 1);

  /* choose virtual tables and initialize resources */
  SC3X (timeavx2_prepare (t, &retval));

#if 0
  /* TODO create a dedicated allocator for this program */
  SC3E_NULL_SET (e, p4est3_quadrant_array_new (sc3_allocator_nocount (),
                                               qvt_avx, n_quads, &qarr_avx));
  SC3E_NULL_SET (e, p4est3_quadrant_array_new (sc3_allocator_nocount (),
                                               qvt, n_quads, &qarr));

  SC3E_NULL_SET (e, sc3_array_index (qarr_avx, 0, &p));
  SC3E_NULL_SET (e, p4est3_quadrant_root (qvt_avx, p));

  SC3E_NULL_SET (e, sc3_array_index (qarr, 0, &p));
  SC3E_NULL_SET (e, p4est3_quadrant_root (qvt, p));
#endif

  if (!retval) {
    if (t->mpirank == 0) {
      SC3X (timeavx2_measure (t));
    }

#if 0
    if (e == NULL) {
      printf ("Executing time: \n");
      SC3E_SET (e, measure_child (qarr_avx, qarr, qvt_avx, qvt, n_quads));
      report_errors (&e);

      SC3E_SET (e, measure_parent (qarr_avx, qarr, qvt_avx, qvt, n_quads));
      report_errors (&e);

      SC3E_SET (e, measure_compare (qarr_avx, qarr, qvt_avx, qvt, n_quads));
      report_errors (&e);

      SC3E_SET (e, measure_successor (qarr_avx, qarr, qvt_avx, qvt, n_quads));
      report_errors (&e);
    }
#endif

    SC3X (timeavx2_cleanup (t));
  }

#if 0
  SC3E_NULL_SET (e, sc3_array_destroy (&qarr_avx));
  SC3E_NULL_SET (e, sc3_array_destroy (&qarr));

  /* A program is not guaranteed to exist beyond MPI_Finalize. */
  report_errors (&e);
#endif

  /* MPI_Finalize comes last in a program.  We abort should this go wrong. */
  SC3X (sc3_MPI_Finalize ());
  return 0;
}
