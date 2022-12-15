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

sc3_error_t        *
p4est3_partition (p4est3_t * p3)
{
  int p;
  int nodesize, noderank, node_num;
  int *node_sizes, *node_offsets;
  p4est3_locidx       *num_quadrants_in_proc;
  p4est3_gloidx        prev_quadrant, next_quadrant, qcount;
  p4est3_gloidx        new_right_border;

  /* We suppose to call this function after setting up routine */
  SC3A_CHECK (p3->old != NULL);
  SC3A_IS (p4est3_is_setup, p3->old);

  SC3E (sc3_mpienv_get_nodesize (p3->split_info, &nodesize));
  /* this function does nothing for processes without shared memory */
  if (nodesize == 1) {
    return;
  }

  if (p3->cweight == NULL) {
    /* Divide up the quadrants equally */
    /*SC3E (p4est3_part_array_new
          (p3->alloc, sizeof (p4est3_gloidx), nodesize + 1, nodesize + 1,
           num_quadrants_in_proc));
    for (p = 0, next_quadrant = 0; p < nodesize; ++p) {
      prev_quadrant = next_quadrant;
      next_quadrant = p4est3_glocut (p3->global_num_quads, p + 1, nodesize);
      qcount = next_quadrant - prev_quadrant;
      SC3A_CHECK (0 <= qcount && qcount <= (p4est3_gloidx) P4EST3_GLOIDX_MAX);
      num_quadrants_in_proc[p] = (p4est3_locidx) (qcount);
    }*/

    /* Find a new right border for the local partition */
    SC3E (sc3_mpienv_get_noderank (p3->split_info, &noderank));
    new_right_border
      = p4est3_glocut (p3->global_num_quads, noderank + 1, nodesize);
    /* Find to which process belongs the new right border */
    SC3E (sc3_mpienv_get_node_num (p3->split_info, &node_num));
    SC3E (sc3_mpienv_get_node_sizes (p3->split_info, &node_sizes));
    SC3E (sc3_mpienv_get_node_offsets (p3->split_info, &node_offsets));
    SC3E (p4est3_find_partition (p3->alloc, node_sizes[node_num], p3->goffset +,
                                 new_right_border, new_right_border,
                                 
                                 ))


  }

}

#ifdef __cplusplus
#if 0
{
#endif
}
#endif