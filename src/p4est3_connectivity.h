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

/** \file p4est3_connectivity.h
 *
 * The connectivity data structure describes neighborhood relations of trees.
 * Each tree root identifies a logical square (2D) or cube (3D) in the domain.
 * Multiple root cubes can be connected across faces, corners, and edges.
 * These connections include a relative orientation between neighbors.
 * They are purely topological: we do not specify the geometry here.
 * The forest is composed of all roots and their connections.
 *
 * To create a forest, one connectivity structure is required.
 *
 * File looks ok in general.  Issues remaining:
 *   - Doxygenate everything
 *   - Write test/test3_connectivity.c
 *   - Move connectivity_new_num_trees and implementation to the test
 *   - Treat all sc3 leak errors as fatal.
 *   - ref and unref merely count.  Deallocation happens in destroy.
 *
 * \ingroup p4est3
 */

#ifndef P4EST3_CONNECTIVITY_H
#define P4EST3_CONNECTIVITY_H

#include <p4est3_base.h>

#ifdef __cplusplus
extern              "C"
{
#if 0
}
#endif
#endif

#if 0
typedef int         (*p4est3_connectivity_is_t) (void *slf, char *reason);
#endif

/** General virtual function taking one in-out argument. */
typedef sc3_error_t *(*p4est3_connectivity_inout_t) (void *slf);

/** Virtual function with one \ref p4est3_topidx output argument. */
typedef sc3_error_t *(*p4est3_connectivity_get_topidx_t)
                    (void *slf, p4est3_topidx * ptopidx);

/** One way to create a connectivity is to provide a virtual table with state.
 * The members of this table must be set before passing it to \ref
 * p4est3_connectivity_set_vtable, where we make a deep copy.
 * This method is suited to wrap any third-party object into p4est.
 */
typedef struct p4est3_connectivity_vtable
{
  /** Necessary function to return a positive number of trees. */
  p4est3_connectivity_get_topidx_t get_num_trees;

  /** This function may be NULL, e.g. when no state requires destruction */
  p4est3_connectivity_inout_t destroy;
}
p4est3_connectivity_vtable_t;

/** Check whether a connectivity virtual table is valid, thus ready to use.
 * \param [in] cvt      Any pointer.  NULL is considered not valid.
 * \param [out] reason  May be NULL.  Otherwise, will be filled with the
 *                      empty string on validity or the issue found otherwise.
 * \return              Boolean value.
 */
int                 p4est3_connectivity_vtable_is_valid
  (const p4est3_connectivity_vtable_t * cvt, char *reason);

/** The connectivity object is an opaque structure. */
typedef struct p4est3_connectivity p4est3_connectivity_t;

/** Check whether a connectivity is valid.
 * \param [in] c        Any pointer.  NULL is considered not valid.
 * \param [out] reason  May be NULL.  Otherwise, will be filled with the
 *                      empty string on validity or the issue found otherwise.
 * \return              Boolean value.
 */
int                 p4est3_connectivity_is_valid (const p4est3_connectivity_t
                                                  * c, char *reason);

/** Check whether a connectivity is valid and not setup yet.
 * \param [in] c        Any pointer.  NULL is considered not new.
 * \param [out] reason  May be NULL.  Otherwise, will be filled with the
 *                      empty string on validity or the issue found otherwise.
 * \return              Boolean value.
 */
int                 p4est3_connectivity_is_new (const p4est3_connectivity_t *
                                                c, char *reason);

/** Check whether a connectivity is valid and setup.
 * \param [in] c        Any pointer.  NULL is considered not setup.
 * \param [out] reason  May be NULL.  Otherwise, will be filled with the
 *                      empty string on validity or the issue found otherwise.
 * \return              Boolean value.
 */
int                 p4est3_connectivity_is_setup (const p4est3_connectivity_t
                                                  * c, char *reason);

/** Create a connectivity under construction.
 * Its properties must further be set by \c p4est3_connectivit_set_*
 * and then \ref p4est3_connectivity_setup must be called to finalize construction.
 * \param [in,out] alloc   Allocator must be setup.  It is referenced
 *                         and kept around while connectivity is live.
 * \param [out] pc         Pointer to a pointer, the latter will be updated.
 * \return                 NULL on success, error object otherwise.
 */
sc3_error_t        *p4est3_connectivity_new (sc3_allocator_t * alloc,
                                             p4est3_connectivity_t ** pc);

/** Select the virtual table creation method for the connectivity.
 * \param [in,out] c    Connectivity under construction.
 * \param [in] cvt      Valid connectivity virtual table.  We make a deep copy.
 * \param [in] slf      Self (state) of virtual connectivity passed along.
 * \return              NULL on success, error object otherwise.
 */
sc3_error_t        *p4est3_connectivity_set_vtable
  (p4est3_connectivity_t * c, p4est3_connectivity_vtable_t * cvt, void *slf);

/** Set the number of trees that constitute this connectivity.
 * \param [in,out] c        Connectivity under construction.
 * \param [in] num_trees    Positive number of trees.
 * \return                  NULL on success, error object otherwise.
 */
sc3_error_t        *p4est3_connectivity_set_num_trees
  (p4est3_connectivity_t * c, p4est3_topidx num_trees);

/** Finalize a connectivity under construction for use with a forest.
 * \param [in,out] c        Connectivity must be valid but not yet setup.
 * \return                  NULL on success, error object otherwise.
 */
sc3_error_t        *p4est3_connectivity_setup (p4est3_connectivity_t * c);

sc3_error_t        *p4est3_connectivity_ref (p4est3_connectivity_t * c);
sc3_error_t        *p4est3_connectivity_unref (p4est3_connectivity_t ** c);

/** Destroy a connectivity that must be valid, but may or may not be setup.
 * It is required that all additional references have been dropped beforehand.
 * \param [in,out] c    Valid connectivity on input.  NULL on output.
 * \return              NULL on success, error object otherwise.
 */
sc3_error_t        *p4est3_connectivity_destroy (p4est3_connectivity_t ** c);

/** Query number of trees in a connectivity.
 * \param [in] c            Connectivity must be setup.
 * \param [out] pnum_trees  Not NULL.  Number of trees is placed here.
 * \return                  NULL on success, error object otherwise.
 */
sc3_error_t        *p4est3_connectivity_get_num_trees
  (const p4est3_connectivity_t * c, p4est3_topidx * pnum_trees);

/* TODO remove me */
sc3_error_t        *p4est3_connectivity_new_num_trees
  (sc3_allocator_t * alloc, p4est3_topidx num_trees,
   p4est3_connectivity_t ** pc);

/** Create a connectivity readily setup to represent the unit cube.
 * \param [in,out] alloc   Allocator must be setup.  It is referenced
 *                         and kept around while connectivity is live.
 * \param [out] pc         Connectivity is setup on output.
 * \return                 NULL on success, error object otherwise.
 */
sc3_error_t        *p4est3_connectivity_new_unitcube
  (sc3_allocator_t * alloc, p4est3_connectivity_t ** pc);

#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif /* !P4EST3_CONNECTIVITY_H */
