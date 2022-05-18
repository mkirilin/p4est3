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

typedef struct time
{
  int                 mpirank;
  int max_level;
  int mpisize;
  int nquads_init;
  int max_cpu_ref;

  sc3_allocator_t    *alloc;
  sc3_array_t        *qarr;
  sc3_array_t        *qarr_avx;
  p4est3_quadrant_vtable_t sqvt, *qvt;
  p4est3_quadrant_vtable_t sqvt_avx, *qvt_avx;
  p4est3_quadrant_vtable_t sqvt_mort, *qvt_mort;
}
time_t;

static sc3_error_t *
time_prepare (time_t * t, int *retval)
{
  void               *p;

  SC3E_RETVAL (retval, -1);
  SC3A_CHECK (t != NULL);
  SC3A_CHECK (t->n_quads > 0);

  /* the standard p4est2 virtual table always exists */
  SC3E (p4est3_quadrant_vtable_p4est (&t->qvt));
  SC3E_DEMAND (t->qvt != NULL, "standard qvt: "
               "p4est is not build neither in 2D nor 3D");

  /* the AVX virtual table can only be set with hardware support */
  SC3E (p4est3_quadrant_yx_vtable (&t->qvt_avx));
  SC3E_DEMAND (t->qvt_avx != NULL, "AVX is not supported by hardware "
      "or p4est is not build neither in 2D nor 3D\n");

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
timeavx2_measure (time_t * t)
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
timeavx2_cleanup (time_t * t)
{
  SC3A_CHECK (t != NULL);
  SC3A_CHECK (t->n_quads > 0);

  SC3E (sc3_array_destroy (&t->qarr_avx));
  SC3E (sc3_array_destroy (&t->qarr));
  SC3E (sc3_allocator_destroy (&t->alloc));
  return NULL;
}

static sc3_error_t *
interpret_command_line (const int argc, const char **argv, time_t *t)
{
  int scale;
  SC3E (sc3_MPI_Comm_size (SC3_MPI_COMM_WORLD, &t->mpisize));
  if (argc == 1) {
    printf ("Execution without parameters. "
            "Default parameters are applied.\n"
            "Parameter's format: "
            "<MAX CPU REF> <SCALE> <MAX LEVEL>\n");
    t->max_cpu_ref = t->mpisize;
    t->max_level = 1;
    t->nquads_init = t->max_cpu_ref / t->mpisize;
  }
  if (argc >= 2){
    t->max_cpu_ref = atoi (argv[2]);
    SC3E_DEMAND (t->max_cpu_ref != 0, "MAX CPU REF is interpreted wrong");
  }
  if (argc >= 3){
    scale = atoi(argv[3]);
    SC3E_DEMAND (scale != 0, "SCALE is interpreted wrong");
    t->max_cpu_ref *= scale;
  }
  if (argc >= 4){
    t->max_level = atoi(argv[4]);
    SC3E_DEMAND (t->max_level != 0, "MAX LEVEL is interpreted wrong");
  }
  return NULL;
}

int
main (int argc, char **argv)
{
  int                 retval;
  p4est3_locidx       n_quads;
  time_t          st, *t = &st;

  /* MPI_Init comes first in a program.  We abort should this go wrong. */
  SC3X (sc3_MPI_Init (&argc, &argv));
  SC3X (sc3_MPI_Comm_rank (SC3_MPI_COMM_WORLD, &t->mpirank));

  SC3X (interpret_command_line (argc, argv, t));

  /* choose virtual tables and initialize resources */
  SC3X (time_prepare (t, &retval));

  if (!retval) {
    if (t->mpirank == 0) {
      SC3X (timeavx2_measure (t));
    }
    SC3X (timeavx2_cleanup (t));
  }
  /* MPI_Finalize comes last in a program.  We abort should this go wrong. */
  SC3X (sc3_MPI_Finalize ());
  return 0;
}
