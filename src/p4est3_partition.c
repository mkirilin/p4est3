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

#include <p4est3_partition.h>
#include <p4est3_internal.h>

#ifdef __cplusplus
extern              "C"
{
#if 0
}
#endif
#endif

static sc3_error_t *
p4est3_part_array_new (sc3_allocator_t * alloc, size_t esize, int ealloc,
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

sc3_error_t        *
p4est3_partition (p4est3_t * p3)
{
  /* at this stage we have a completely setup refined forest */
  int p;
  int nodesize, noderank, node_num;
  int *node_sizes, *node_offsets;
  char               *temp = p3->temp_quad[0];
  int32_t            *coords;
  p4est3_gloidx        prev_quadrant, next_quadrant, qcount_node;
  p4est3_gloidx        new_right_border, new_left_border;
  sc3_MPI_Comm_t nodecomm;

  /* We suppose to call this function after setting up routine */
  SC3A_CHECK (p3->old != NULL);
  SC3A_IS (p4est3_is_setup, p3->old);

  /* this function does nothing for processes without shared memory */
  if (nodesize == 1) {
    return;
  }
  if (p3->cweight == NULL) {
    /* Divide up the quadrants equally */
    /* Find a new right border for the local partition */
    SC3E (sc3_mpienv_get_nodesize (p3->split_info, &nodesize));
    SC3E (sc3_mpienv_get_node_num (p3->split_info, &node_num));
    SC3E (sc3_mpienv_get_node_offsets (p3->split_info, &node_offsets));
    SC3E (sc3_mpienv_get_noderank (p3->split_info, &noderank));
    SC3E (sc3_mpienv_get_nodecomm (p3->split_info, &nodecomm));

    qcount_node = p3->goffset[node_offsets[node_num + 1]]
                  - p3->goffset[node_offsets[node_num]];
    SC3E (sc3_MPI_Barrier (nodecomm));
    new_left_border = p4est3_glocut (qcount_node, noderank, nodesize);
    new_right_border = p4est3_glocut (qcount_node, noderank + 1, nodesize);

#ifdef P4EST_ENABLE_DEBUG
    if (noderank + 1 == nodesize) {
      SC3E (sc3_MPI_Win_lock (SC3_MPI_LOCK_SHARED, 0, SC3_MPI_MODE_NOCHECK,
                              p3->goffsets->goffsetwin));
      SC3A_CHECK
        (p3->goffset[p3->mpirank] == new_right_border);
      SC3E (sc3_MPI_Win_unlock (0, p3->goffsets->goffsetwin));
    }
    else if (noderank == 0) {
      SC3E (sc3_MPI_Win_lock (SC3_MPI_LOCK_SHARED, 0, SC3_MPI_MODE_NOCHECK,
                              p3->goffsets->goffsetwin));
      SC3A_CHECK (p3->goffset[p3->mpirank] == new_left_border);
      SC3E (sc3_MPI_Win_unlock (0, p3->goffsets->goffsetwin));
    }
#endif

    SC3E (sc3_MPI_Win_lock (SC3_MPI_LOCK_SHARED, 0, SC3_MPI_MODE_NOCHECK,
                            p3->goffsets->goffsetwin));
    p3->goffset[p3->mpirank] = new_left_border;
    SC3E (sc3_MPI_Win_unlock (0, p3->goffsets->goffsetwin));
    p3->local_num_quads = new_right_border - new_left_border;
    SC3E (sc3_MPI_Win_lock (SC3_MPI_LOCK_SHARED, 0, SC3_MPI_MODE_NOCHECK,
                            p3->gposition->gfposwin));
    SC3E(p4est3_quadrant_first_descendant
          (p3->qvt, p3->nodequads[0] + p3->qsize * new_left_border,
           p3->qmaxlevel, p3->gfpos[p3->mpirank]));
    if (noderank == 0) {
      SC3E (sc3_allocator_calloc
            (p3->alloc, p3->qvt->dim, sizeof (int32_t), &coords));
      SC3E (p4est3_translate_quadrant
            (p3->old->qvt, p3->qvt,
             p3->old->gfpos[node_offsets[node_num + 1]],
             p3->gfpos[node_offsets[node_num + 1]], coords));
      SC3E (sc3_allocator_free (p3->alloc, coords));
    }
    SC3E (sc3_MPI_Win_unlock (0, p3->gposition->gfposwin));
  }

}

#ifdef __cplusplus
#if 0
{
#endif
}
#endif