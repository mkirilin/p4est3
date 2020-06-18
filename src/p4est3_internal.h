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

#ifndef P4EST3_INTERNAL_H
#define P4EST3_INTERNAL_H

#include <p4est3.h>
#include <sc3_array.h>
#include <sc3_refcount.h>

typedef struct p4est3_tree
{
  p4est3_topidx       treeid;
  p4est3_gloidx       first_tquad, end_tquad;
  p4est3_locidx       quad_offset;
  p4est3_locidx       num_quads;
  char               *tquads;
}
p4est3_tree_t;

struct p4est3
{
  sc3_refcount_t      rc;
  sc3_allocator_t    *alloc;
  int                 setup;

  sc3_MPI_Comm_t      mpicomm;
  int                 commdup;
  p4est3_connectivity_t *conn;

  sc3_MPI_Comm_t      nodecomm, headcomm;
  sc3_MPI_Info_t      info_noncontig;
  sc3_MPI_Win_t       nodesizewin;
  sc3_MPI_Win_t       gfposwin, gftreewin, goffsetwin;
  sc3_MPI_Win_t       quadwin;
  int                 mpisize, mpirank;
  int                 nodesize, noderank;
  int                 num_nodes;
  int                 node_num;
  int                 node_frank;
  int                *node_sizes;
  int                *node_offsets;

  p4est3_quadrant_vtable_t sqvt, *qvt;
  p4est3_topidx       num_trees;
  int                 level;
  int                 qmaxlevel;
  int                 num_children;
  int                 qsize;
  int                 setup_mode;

  int                 max_threads;
  char              **temp_quad;

  p4est3_locidx       local_num_quads;
  p4est3_gloidx       global_num_quads;
  p4est3_gloidx      *goffset;
  p4est3_topidx      *gftree;
  char               *gfpos;
  char              **nodequads, *quads;

  p4est3_topidx       fltree, lltree, nltrees;
  sc3_array_t        *trees;
};

typedef enum p4est3_setup_mode
{
  P4EST3_NEW_MORTON = 1,        /**< Set every quadrant by its Morton index */
  P4EST3_NEW_SUCCESSOR = 2,     /**< Set every quadrant by the previous one */
  P4EST3_NEW_RECURSIVE = 4,     /**< Recursive calling the child function */
  P4EST3_NEW_MAX_TYPE = 8       /**< Unused bounding value */
}
p4est3_setup_mode_t;

#ifdef __cplusplus
extern              "C"
{
#if 0
}
#endif
#endif

sc3_error_t        *p4est3_tree_index (p4est3_t * p3, p4est3_topidx tt,
                                       p4est3_tree_t ** tree);

sc3_error_t        *p4est3_internal_setup_comm (p4est3_t * p3);
sc3_error_t        *p4est3_internal_setup_cut (p4est3_t * p3,
                                               p4est3_gloidx num_uniform,
                                               int qsize);
sc3_error_t        *p4est3_internal_setup_tree (p4est3_t * p3,
                                                p4est3_gloidx num_uniform);
sc3_error_t        *p4est3_internal_setup_quadrants (p4est3_t * p3);

#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif /* !P4EST3_INTERNAL_H */
