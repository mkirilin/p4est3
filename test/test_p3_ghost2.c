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

#include <p4est_bits.h>
#include <p4est_ghost.h>

#include <p4est3.h>
#include <p4est3_ghost.h>
#include <p4est3_convert_p4est.h>

/*
static int refine_level = 5;

static int
refine_fn (p4est_t * p4est, p4est_topidx_t which_tree,
           p4est_quadrant_t * quadrant)
{
  int                 cid;

  if (which_tree == 2 || which_tree == 3) {
    return 0;
  }

  cid = p4est_quadrant_child_id (quadrant);

  if (cid == P4EST_CHILDREN - 1 ||
      (quadrant->x >= P4EST_LAST_OFFSET (P4EST_MAXLEVEL - 2) &&
       quadrant->y >= P4EST_LAST_OFFSET (P4EST_MAXLEVEL - 2)
#ifdef P4_TO_P8
       && quadrant->z >= P4EST_LAST_OFFSET (P4EST_MAXLEVEL - 2)
#endif
      )) {
    return 1;
  }
  if ((int) quadrant->level >= (refine_level - (int) (which_tree % 3))) {
    return 0;
  }
  if (quadrant->level == 1 && cid == 2) {
    return 1;
  }
  if (quadrant->x == P4EST_QUADRANT_LEN (2) &&
      quadrant->y == P4EST_LAST_OFFSET (2)) {
    return 1;
  }
  if (quadrant->y >= P4EST_QUADRANT_LEN (2)) {
    return 0;
  }

  return 1;
}*/

int main(int argc, char **argv) {
    int mpiret;
    sc_MPI_Comm mpicomm;
    //p4est_t *p4est;
    //p4est3_t *p4est3;
    //p4est_connectivity_t *conn;
    //p4est_ghost_t *ghost_p4est, *ghost_p4est3;
    //int num_cycles = 2;
    //int i;

    /* initialize MPI */
    mpiret = sc_MPI_Init(&argc, &argv);
    SC_CHECK_MPI(mpiret);
    mpicomm = sc_MPI_COMM_WORLD;

    sc_init(mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
    p4est_init(NULL, SC_LP_DEFAULT);

    //conn = p4est_connectivity_new_unitsquare();
    //p4est = p4est_new(mpicomm, conn, 0, NULL, NULL);

    /* refine to make the number of elements interesting */
    //p4est_refine(p4est, 1, refine_fn, NULL);

    /* balance the forest */
    //p4est_balance(p4est, P4EST_CONNECT_FULL, NULL);

    /* do a uniform partition */
    //p4est_partition(p4est, 0, NULL);

    /* create the ghost layer for p4est */
    //ghost_p4est = p4est_ghost_new(p4est, P4EST_CONNECT_FULL);

    /* convert p4est to p4est3 */
    //p4est3_new(mpicomm, &p4est3);
    //p4est3_convert_p4est(p4est, p4est3);

    /* create the ghost layer for p4est3 */
    //ghost_p4est3 = p4est_ghost_new(p4est, P4EST_CONNECT_FULL);
    //p4est3_ghost_fill_p4est(p4est3, ghost_p4est3);

    /* compare ghost layers */
    //SC_CHECK_ABORT(ghost_p4est->ghosts.elem_count == ghost_p4est3->ghosts.elem_count,
    //               "Ghost count mismatch between p4est and p4est3");

    //for (i = 0; i < ghost_p4est->ghosts.elem_count; ++i) {
    //    p4est_quadrant_t *q_p4est = p4est_quadrant_array_index(&ghost_p4est->ghosts, i);
    //    p4est_quadrant_t *q_p4est3 = p4est_quadrant_array_index(&ghost_p4est3->ghosts, i);
    //    SC_CHECK_ABORT(p4est_quadrant_is_equal(q_p4est, q_p4est3),
    //                   "Ghost quadrant mismatch between p4est and p4est3");
    //}

    /* clean up */
    //p4est_ghost_destroy(ghost_p4est);
    //p4est_ghost_destroy(ghost_p4est3);
    //p4est_destroy(p4est);
    //p4est3_destroy(&p4est3);
    //p4est_connectivity_destroy(conn);

    /* exit */
    sc_finalize();

    mpiret = sc_MPI_Finalize();
    SC_CHECK_MPI(mpiret);

    return 0;
}