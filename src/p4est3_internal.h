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
#include <sc3_mpienv.h>
#include <p4est3.h>

/* Use to choose a way of filling a tree with quadrants.*/
typedef enum p4est3_setup_mode
{
  P4EST3_NEW_MORTON,    /**< Set every quadrant by its Morton index */
  P4EST3_NEW_SUCCESSOR, /**< Set every quadrant by the previous one */
  P4EST3_NEW_RECURSIVE, /**< Recursive calling the child function */
  P4EST3_NEW_RECURSIVE_CHILD,
  P4EST3_NEW_RECURSIVE_REGION,
  P4EST3_NEW_MODE_LAST  /**< Unused bounding value */
}
p4est3_setup_mode_t;

/* Use to choose a way of filling a forest based on another one. */
typedef enum p4est3_source_setup
{
  P4EST3_SRC_COPY,     /**< Setup by a simple copying quadrants */
  P4EST3_SRC_REFINE,   /**< Setup with refinement where necessary */
  P4EST3_SRC_COARSE,   /**< Setup with coarsening where necessary */
  P4EST3_SRC_MODE_LAST  /**< Unused bounding value */
}
p4est3_source_setup_t;

/** Internal data for a process-local tree and the quadrants it contains. */
typedef struct p4est3_tree
{
  p4est3_topidx       treeid;   /**< Tree number between p4est3_t::fltree
                                     and p4est3_t::lltree inclusive. */
  p4est3_gloidx       first_tquad;      /**< First local quadrant in this tree
                                             counted from the very first
                                             (lower left) quadrant of this tree. */
  p4est3_gloidx       last_tquad;       /**< Last local quadrant in this tree
                                             (inclusive), counted from first
                                             (lower left) quadrant of this tree.
                                             Greater equal \ref first_tquad. */
  p4est3_gloidx       end_tquad;        /**< Last local quadrant in this tree
                                             (exclusive), counted from first
                                             (lower left) quadrant of this tree.
                                             Equals \ref last_tquad + 1. */
  p4est3_locidx       quad_offset;      /**< Local quadrants before this tree. */
  p4est3_locidx       num_quads;        /**< Local quadrants within this tree. */
  char               *tquads;   /**< Array of local quadrants in this tree.
                                     Subarray of \ref p4est3_t::quads. */
}
p4est3_tree_t;

/** This internal structure holds the members of a forest object.
 * Don't rely on its declaration in code outside the library. */
struct p4est3
{
  /* variables of internal state used during the whole lifetime */
  sc3_refcount_t      rc;       /**< Reference counter in use. */
  sc3_array_t        *talloc;   /**< Allocator for recoursive mode */
  sc3_allocator_t    *alloc;    /**< Memory allocator in use. */
  int                 setup;    /**< Boolean: object is setup. */
  int                 accessed_conn;    /**< Number of currently active
                                             connectivity accesses. */

  /* this forest may be wrapping a virtual implementation */
  p4est3_vtable_t     spvt;     /**< Memory pointed to by \ref pvt. */
  p4est3_vtable_t    *pvt;      /**< If not NULL, forest virtual table. */
  void               *slf;      /**< Context to use with virtual forest. */

  /* variables set before p4est3_setup */
  sc3_MPI_Comm_t      mpicomm;  /**< Valid MPI communicator. */
  int                 commdup;  /**< Boolean: communicator has been duped. */
  p4est3_connectivity_t *conn;  /**< Pointer to the relevant connectivity. */
  p4est3_topidx       num_trees;        /**< Number of trees in \ref conn. */
  const p4est3_quadrant_vtable_t *qvt;        /**< Always points to \ref sqvt. */
  int                 level;    /**< Configuration variable for initial level.
                                     Depending on the available memory and
                                     index space, may be reduced during
                                     \ref p4est3_setup. */
  p4est3_setup_mode_t setup_mode;       /**< Choose the method of quadrant creation*/
  p4est3_t           *old;      /**< Pointer to the setup forest */

  /* variables populated during p4est3_setup: communicator related */
  int                 mpisize;          /**< Size of forest communicator. */
  int                 mpirank;          /**< Rank in forest communicator. */
  int                 shared;           /**< MPI sharined memory enable/disable
                                             indicator. */
  sc3_mpienv_t       *split_info;       /**<  Pointer to a relevant MPI
                                              processes split related
                                              information. */

  /* variables populated during p4est3_setup: partition related */
  sc3_MPI_Win_t       gftreewin;        /**< Array of (\ref mpisize + 1) \ref
                                             p4est3_topidx integers for the
                                             global partition of trees. */
  sc3_MPI_Win_t       gfposwin;         /**< Array of (\ref mpisize + 1) times \ref
                                        qsize bytes for global first quadrant. */
  sc3_MPI_Win_t       goffsetwin;       /**< Array of (\ref mpisize + 1) \ref
                                        p4est3_gloidx for global quadrant offsets. */
  int                 qsize;            /**< Store byte size of one quadrant. */
  int                 qmaxlevel;        /**< Maximum allowed refinement level. */
  int                 num_children;     /**< Number of children for a quadrant. */
  int                 max_threads;      /**< Max threads from querying openmp. */
  char              **temp_quad;        /**< Quadrant work space, one per thread. */
  p4est3_locidx       local_num_quads;  /**< Count process-local quadrants. */
  p4est3_gloidx       global_num_quads; /**< Count all quadrants globally. */
  p4est3_gloidx      *goffset;          /**< Pointer to \ref goffsetwin's memory. */
  p4est3_topidx      *gftree;           /**< Pointer to \ref gftreewin's memory. */
  char               *gfpos;            /**< Pointer to \ref gfposwin's memory. */

  /* variables populated during p4est3_setup: tree and quadrant storage */
  sc3_MPI_Win_t       quadwin;          /**< Shared memory stores the quadrants
                                        for all ranks on this node in order.
                                        Each node rank's local subwindow is
                                        associated with its rank. */
  char              **nodequads;        /**< Array of \ref nodesize holds pointers
                                        to their respective first quadrants in
                                        the storage of \ref quadwin. */
  char               *quads;            /**< Pointer to first quadrant local
                                        to his process equals \ref
                                        nodequads[\ref noderank]. */
  sc3_array_t        *trees;    /**< Array of only the process-local trees. */
  p4est3_topidx       fltree;   /**< Number of first local tree, or -1
                                     if process holds no quadrants.
                                     Relative to all trees in \ref conn. */
  p4est3_topidx       lltree;   /**< Number of last local tree inclusive,
                                     or -2 if process holds no quadrants. */
  p4est3_topidx       nltrees;  /**< Number of trees with local quadrants. */

  /* functions set before p4est3_setup */
  p4est3_refine_callback_t crefine; /**< Refinemet callback function */
  p4est3_coarsen_callback_t ccoarse; /**< Coarsening  callback function */

  /* pointer to user data, p4est does not touch them */
  void               *user_data;
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
sc3_error_t        *p4est3_internal_setup_quadrants (p4est3_t * p3);
sc3_error_t        *p4est3_internal_setup_from_source (p4est3_t * p3);
/** \endcond */

/* TODO: document default value for all _set_ */
/** Set a way that creates quadrants in a tree in a setup p4est3 phase.
 * TODO: possibly rename function
 * \param [in,out] p3       The forest must not have been setup.
 * \param [in] mode         See \ref p4est3_setup_mode_t type for
 *                          available options. Default value is
 *                          P4EST3_NEW_MORTON.
 */
sc3_error_t        *p4est3_set_setup_mode (p4est3_t * p3,
                                           p4est3_setup_mode_t mode);

#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif /* !P4EST3_INTERNAL_H */
