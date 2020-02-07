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

#include <p4est_p4est3.h>
#include <p4est3.h>

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
  p4est3_t           *p3;

  SC3A_IS (sc3_allocator_is_setup, alloc);

  SC3E (p4est3_new (alloc, &p3));
  SC3E (p4est3_set_comm (p3, mpicomm, 1));
  SC3E (p4est3_set_vtable (p3, qvt));
  SC3E (p4est3_set_num_trees (p3, num_trees));
  SC3E (p4est3_set_level (p3, level));
  SC3E (p4est3_setup (p3));

  SC3E (p4est3_destroy (&p3));
  return NULL;
}

static sc3_error_t *
free_allocator (sc3_allocator_t ** alloc)
{
  SC3A_IS (sc3_allocator_is_setup, *alloc);
  SC3E (sc3_allocator_destroy (alloc));
  return NULL;
}

static void
report_errors (sc3_allocator_t * mainalloc, sc3_error_t ** pe)
{
  char                eflat[SC3_BUFSIZE];
  char                reason[SC3_BUFSIZE];

  if (pe != NULL && *pe != NULL) {
    /* TODO print error messages in a nicer way */
    sc3_error_destroy_noerr (pe, eflat);
    fprintf (stderr, "Error: %s\n", eflat);
  }

  if (!sc3_allocator_is_free (mainalloc, reason)) {
    fprintf (stderr, "Allocation error: %s\n", reason);
  }

#if 0
  e = sc3_error_destroy (pe);

  /* TODO synchronize e across MPI processes */

  if (e != NULL) {
    fprintf (stderr, "Errors remain\n");
    sc3_error_destroy (&e);
  }
#endif
}

int
main (int argc, char **argv)
{
  int                 level;
  p4est3_topidx       num_trees;
  sc3_allocator_t    *alloc, *mainalloc;
  sc3_error_t        *e;
  sc3_MPI_Comm_t      mpicomm;
  p4est3_quadrant_vtable_t vtable, *qvt = &vtable;

  mainalloc = sc3_allocator_nothread ();
  mpicomm = SC3_MPI_COMM_WORLD;
  p4est_quadrant_vtable (qvt);

  num_trees = 2;
  level = 3;

  SC3E_SET (e, sc3_MPI_Init (&argc, &argv));
  SC3E_NULL_SET (e, make_allocator (mainalloc, &alloc));
  SC3E_NULL_SET (e, test_p4est_new (alloc, mpicomm, qvt, num_trees, level));
  SC3E_NULL_SET (e, free_allocator (&alloc));

  SC3E_NULL_SET (e, sc3_MPI_Finalize ());
  report_errors (mainalloc, &e);

  return 0;
}
