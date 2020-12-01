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

/** \file p4est3_internal.h
 *
 * Private declarations for use within the library.
 * This file must never be included by public header files.
 *
 * Me make no provisions on stability or backwards compatibility.
 * There is usually no reason to include this file outside of the library.
 *
 * \ingroup p4est3
 */

#ifndef P4EST3_INTERNAL_H
#define P4EST3_INTERNAL_H

#include <sc3_array.h>
#include <sc3_refcount.h>
#include <p4est3.h>

/** Internal data for a process-local tree and the quadrants it contains. */
typedef struct p4est3_tree
{
  p4est3_topidx       treeid;   /**< Tree number is zero based. */
  p4est3_gloidx       first_tquad;
  p4est3_gloidx       end_tquad;
  p4est3_locidx       quad_offset;      /**< Local quadrants before this tree. */
  p4est3_locidx       num_quads;        /**< Local quadrants within this tree. */
  char               *tquads;   /**< Array of local tree quadrants. */
}
p4est3_tree_t;

/** This internal structure holds the members of a forest object.
 * Don't rely on its declaration in code outside the library. */
struct p4est3
{
  /* variables of internal state used during the whole lifetime */
  sc3_refcount_t      rc;       /**< Reference counter in use. */
  sc3_allocator_t    *alloc;    /**< Memory allocator in use. */
  int                 setup;    /**< Boolean: object is setup. */
  int                 accessed_conn;    /**< Number of currently active
                                             connectivity accesses. */

  /* this forest may be wrapping a virtual implementation */
  p4est3_vtable_t     spvt;     /**< Memory pointed to by \ref pvt. */
  p4est3_vtable_t     *pvt;     /**< If not NULL, forest virtual table. */
  void               *slf;      /**< Context to use with virtual forest. */

  /* variables set before \ref p4est3_setup */
  sc3_MPI_Comm_t      mpicomm;  /**< Valid MPI communicator. */
  int                 commdup;  /**< Boolean: communicator has been duped. */
  p4est3_connectivity_t *conn;  /**< Pointer to the relevant connectivity. */
  p4est3_topidx       num_trees;        /**< Number of trees in \ref conn. */
  p4est3_quadrant_vtable_t sqvt;        /**< Memory pointed to by \ref qvt.
                                             Stores virtual quadrant methods. */
  p4est3_quadrant_vtable_t *qvt;        /**< Always points to \ref sqvt. */
  int                 level;    /**< Configuration variable for initiel level. */

  /* variables populated during p4est3_setup: communicator related */
  sc3_MPI_Comm_t      nodecomm;         /**< All ranks of shared memory node. */
  sc3_MPI_Comm_t      headcomm;         /**< Contains first rank of each node. */
  sc3_MPI_Info_t      info_noncontig;   /**< Key "alloc_shared_noncontig" set. */
  sc3_MPI_Win_t       nodesizewin;      /**< Shared memory segment allocated
                                             on first rank of a node, available
                                             to all ranks on that node.  Its
                                             element count is (2 + 2 * \ref
                                             num_nodes + 1) integers.
                                             Its contents hold
 *                                  * number of nodes for this run
 *                                  * zero-based number of this node
 *                                  * for each node number of ranks on it
 *                                  * for each node and one beyond the
 *                                    number of ranks before it
 */

  /* variables populated during p4est3_setup: partition related */
  sc3_MPI_Win_t       gfposwin;
  sc3_MPI_Win_t       gftreewin;
  sc3_MPI_Win_t       goffsetwin;
  sc3_MPI_Win_t       quadwin;
  int                 mpisize;          /**< Size of forest communicator. */
  int                 mpirank;          /**< Rank in forest communicator. */
  int                 nodesize;         /**< Size of node communicator. */
  int                 noderank;         /**< Rank in node communicator. */
  int                 num_nodes;        /**< Number of shared memory nodes. */
  int                 node_num;         /**< Zero-based node number. */
  int                 node_frank;       /**< Rank within forest communicator
                                             of first rank on this node. */
  int                *node_sizes;       /**< For each node, number of its ranks. */
  int                *node_offsets;     /**< For each node and one beyond, the
                                             number of ranks before it. */

  int                 qmaxlevel;
  int                 num_children;
  int                 qsize;            /**< Store byte size of one quadrant. */

  int                 max_threads;      /**< Max threads from querying openmp. */
  char              **temp_quad;

  p4est3_locidx       local_num_quads;  /**< Count process-local quadrants. */
  p4est3_gloidx       global_num_quads; /**< Count all quadrants globally. */
  p4est3_gloidx      *goffset;
  p4est3_topidx      *gftree;
  char               *gfpos;
  char              **nodequads;
  char               *quads;

  p4est3_topidx       fltree;   /**< Number of first local tree.
                                     Relative to all trees in \ref conn. */
  p4est3_topidx       lltree;   /**< Number of last local tree inclusive. */
  p4est3_topidx       nltrees;  /**< Number of process-local trees. */
  sc3_array_t        *trees;    /**< Array of only the process-local trees. */
};

#ifdef __cplusplus
extern              "C"
{
#if 0
}
#endif
#endif

/** Index into array of local trees and return pointer to indexed tree.
 * \param [in] p3   Required to be non-NULL.
 *                  For speed, we do not check for validity each time.
 * \param [in] tt   Must be the index of a local tree, thus between \ref
 *                  p4est3_t::fltree and \ref p4est3_t::lltree inclusive.
 * \param [in,out] tree     Non-NULL on input.  On output, value is
 *                          populated with a pointer to local tree structure.
 * \return                  NULL on success, error object otherwise.
 */
sc3_error_t        *p4est3_tree_index (p4est3_t * p3, p4est3_topidx tt,
                                       p4est3_tree_t ** tree);

/** \cond P4EST_FALSE */
/* these functions are not documented on purpose */
sc3_error_t        *p4est3_internal_setup_comm (p4est3_t * p3);
sc3_error_t        *p4est3_internal_setup_cut (p4est3_t * p3,
                                               p4est3_gloidx num_uniform,
                                               int qsize);
sc3_error_t        *p4est3_internal_setup_tree (p4est3_t * p3,
                                                p4est3_gloidx num_uniform);
sc3_error_t        *p4est3_internal_setup_morton (p4est3_t * p3);
/** \endcond */

#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif /* !P4EST3_INTERNAL_H */
