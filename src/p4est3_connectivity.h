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
 * We provide the convenience constructors \ref
 * p4est3_connectivity_new_unitsquare (2D) and \ref
 * p4est3_connectivity_new_unitcube (3D).
 *
 * The connectivity object is designed as a shell for a virtual implementation.
 * (The builtin construction is limited to setting the dimension and
 * a number of non-connected trees (all faces are boundary faces)).
 * This connectivity becomes more practically useful when setting the virtual
 * table to wrap an arbitrary externally maintained connectivity.
 * We provide wrappers for classic p4est and p8est connectivities
 * in the files \ref p4est3_p4est.h and \ref p4est3_p8est.h.
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

/** General virtual function taking one in-out argument. */
typedef sc3_error_t *(*p4est3_connectivity_inout_t) (void *slf);

/** In/out: tree number, face number; out: orientation (input 0).
 * Only for a tree connection the output arguments need to be updated.
 */
typedef             sc3_error_t
  * (*p4est3_connectivity_get_face_t) (void *slf, p4est3_topidx * which_tree,
                                       int *nface, int *orient);

/** One way to create a connectivity is to provide a virtual table with state.
 * The members of this table must be set before passing it to \ref
 * p4est3_connectivity_set_vtable, where we make a shallow copy.
 *
 * Whenever a non-NULL virtual table is set in a connectivity at the time of
 * \ref p4est3_connectivity_setup, it will override all other settings.
 *
 * This method is suited to wrap any compatible third-party object into p4est.
 */
typedef struct p4est3_connectivity_vtable
{
  int                 dim;      /**< Space dimension is 1, 2 or 3. */
  p4est3_topidx       num_trees;    /**< Number of trees is positive. */

  /** This function may be NULL, e.g.\ when no state requires destruction */
  p4est3_connectivity_inout_t destroy;
  /** Query face connection; if NULL report physical boundary always */
  p4est3_connectivity_get_face_t get_face;
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
 * The connectivity defaults to dimension 2 and one tree.
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
 * This overrides all other \c p4est3_connectivity_set_* calls made.
 * \param [in,out] c    Connectivity under construction.
 * \param [in] cvt      Valid connectivity virtual table.
 *                      We make a shallow copy, that is, copy all elements.
 *                      It is thus safe if \c *cvt lives on the stack.
 * \param [in] slf      Self (state) of virtual connectivity passed along.
 * \return              NULL on success, error object otherwise.
 */
sc3_error_t        *p4est3_connectivity_set_vtable
  (p4est3_connectivity_t * c, p4est3_connectivity_vtable_t * cvt, void *slf);

/** Set the spatial dimension of this connectivity.
 * \param [in,out] c        Connectivity under construction.
 * \param [in] dim          Dimension from 1 to 3.  Default is 2.
 * \return                  NULL on success, error object otherwise.
 */
sc3_error_t        *p4est3_connectivity_set_dim
  (p4est3_connectivity_t * c, int dim);

/** Set the number of trees that constitute this connectivity.
 * \param [in,out] c        Connectivity under construction.
 * \param [in] num_trees    Positive number of trees.  Default is 1.
 * \return                  NULL on success, error object otherwise.
 */
sc3_error_t        *p4est3_connectivity_set_num_trees
  (p4est3_connectivity_t * c, p4est3_topidx num_trees);

/** Finalize a connectivity under construction for use with a forest.
 * \param [in,out] c        Connectivity must be valid but not yet setup.
 *                          It is returned with a reference count of one.
 * \return                  NULL on success, error object otherwise.
 */
sc3_error_t        *p4est3_connectivity_setup (p4est3_connectivity_t * c);

/** Increase reference counter of a connectivity after setup.
 * \param [in,out] c        Must be setup.  Increase its reference counter.
 * \return                  NULL on success, error object otherwise.
 */
sc3_error_t        *p4est3_connectivity_ref (p4est3_connectivity_t * c);

/** Decrease reference counter of a connectivity after setup.
 * The lowest legal value for the reference counter is one, as after setup.
 * This function never destroys the object: it is not legal to unref below one.
 * The only way to deallocate a connectivity is \ref p4est3_connectivity_destroy.
 * \param [in,out] c        Must be setup and have reference counter greater one.
 *                          Decrease its reference counter.
 *                          When the count reaches one, nothing happens.
 * \return                  NULL on success, error object otherwise.
 */
sc3_error_t        *p4est3_connectivity_unref (p4est3_connectivity_t * c);

/** Destroy a connectivity that must be valid, but may or may not be setup.
 * It must have a reference count of exactly one.  Otherwise we return an error.
 * Thus, all additional references must have been dropped before calling.
 * \param [in,out] c    Valid connectivity with one reference.  NULL on output.
 * \return              NULL on success, error object otherwise.
 */
sc3_error_t        *p4est3_connectivity_destroy (p4est3_connectivity_t ** c);

/** Query spatial dimension of a connectivity.
 * \param [in] c            Connectivity must be setup.
 * \param [out] pdim        Not NULL.  Dimension is placed here.
 * \return                  NULL on success, error object otherwise.
 */
sc3_error_t        *p4est3_connectivity_get_dim
  (const p4est3_connectivity_t * c, int *pdim);

/** Query number of trees in a connectivity.
 * \param [in] c            Connectivity must be setup.
 * \param [out] pnum_trees  Not NULL.  Number of trees is placed here.
 * \return                  NULL on success, error object otherwise.
 */
sc3_error_t        *p4est3_connectivity_get_num_trees
  (const p4est3_connectivity_t * c, p4est3_topidx * pnum_trees);

/** Query a tree connection across a face.
 * For a physical boundary face we return same tree and same face.
 * \param [in] c            Connectivity must be setup.
 *                          If the connectivity is not face-enabled,
 *                          we always return a physical boundary.
 * \param [in,out] which_tree   On input, valid tree number.
 *                              On output, connecting tree's number.
 * \param [in,out] nface        On input, valid face number.
 *                              On output, connecting face's number.
 * \param [out] orient          On output, orientation of connection.
 * \return                  NULL on success, error object otherwise.
 */
sc3_error_t        *p4est3_connectivity_get_face
  (const p4est3_connectivity_t * c,
   p4est3_topidx * which_tree, int *nface, int *orient);

/** Query child id that touches a face \a nface at a face corner \a i.
 * \param [in] c            Connectivity must be setup.
 * \param [in] nface        Valid face number.
 * \param [in, out] i       On input, valid face corner number of the
 *                          face \a nface. 3D: 0..3, 2D: 0..1.
 *                          On output, child id of the quadrant touches both 
 *                          nface \a face and face corner \a i.
 * \return                  NULL on success, error object otherwise.
 */
sc3_error_t        *p4est3_connectivity_get_face_child_id
  (const p4est3_connectivity_t * c, int nface, int *i);

/** Transform a face corner across one of the adjacent faces into a neighbor tree.
 * This version expects the neighbor face and orientation separately.
 * \param [in] c            Connectivity must be setup.
 * \param [in,out] fc       On input, a face corner number in 0..3.
 *                          On output, the face corner number relative
 *                          to the neighbor's face.
 * \param [in] f     A face that the face corner \a fc is relative to.
 * \param [in] nf    A neighbor face that is on the other side of \f.
 * \param [in] o     The orientation between tree boundary faces \a f and \nf.
 * \return                  NULL on success, error object otherwise.
 */
sc3_error_t        *p4est3_connectivity_face_neighbor_face_corner
  (const p4est3_connectivity_t * c, int *fc, int f, int nf, int o);

/** Create a connectivity readily setup to represent the 2D unit square.
 * \param [in,out] alloc   Allocator must be setup.  It is referenced
 *                         and kept around while connectivity is live.
 * \param [out] pc         Connectivity is setup on output.
 * \return                 NULL on success, error object otherwise.
 */
sc3_error_t        *p4est3_connectivity_new_unitsquare
  (sc3_allocator_t * alloc, p4est3_connectivity_t ** pc);

/** Create a connectivity readily setup to represent the 3D unit cube.
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
