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

/** \file p4est3.h
 *
 * Main interface file to construct and interact with a version 3 forest.
 * The forest may be constructed piece by piece using the setter functions.
 * It may also be defined by populating and passing a forest virtual table.
 * The latter approach allows third-party objects to pass as a legal forest.
 *
 * \ingroup p4est3
 */

#ifndef P4EST3_H
#define P4EST3_H

#include <sc3_mpi.h>
#include <p4est3_connectivity.h>
#include <p4est3_quadrant_vtable.h>

#ifdef __cplusplus
extern              "C"
{
#if 0
}
#endif
#endif

/** General virtual function taking one in-out argument. */
typedef sc3_error_t *(*p4est3_inout_t) (void *slf);

/** General virtual function with one topidx out argument. */
typedef sc3_error_t *(*p4est3_out1t_t) (void *slf, p4est3_topidx * t1);

/** General virtual function with two topidx out argument. */
typedef sc3_error_t *(*p4est3_out2t_t) (const void *slf,
                                        p4est3_topidx * t1,
                                        p4est3_topidx * t2);

/** One way to create a forest is to provide a virtual table with state.
 * The members of this table must be set before passing it to \ref
 * p4est3_set_vtable, where we make a shallow copy.
 *
 * Whenever a non-NULL virtual table is set in a forest at the time of
 * \ref p4est3_setup, it will override all other settings.
 *
 * This method is suited to wrap any compatible third-party object into p4est.
 */

typedef struct p4est3_vtable
{
  int                 dim;      /**< Space dimension is 1, 2 or 3. */
  sc3_MPI_Comm_t      mpicomm;  /**< Valid MPI communicator. */
  p4est3_connectivity_t *c3;    /**< This connectivity must match the virtual
                                     forest to create.  It must be setup. */
  const p4est3_quadrant_vtable_t *qvt;  /**< Quadrant table must match forest.
                                             We make a shallow copy. */

  /** Query functions */
  p4est3_out2t_t      get_local_num_trees;      /**< Number of local trees. */

  /** This function may be NULL, e.g.\ when no state requires destruction */
  p4est3_inout_t      destroy;
}
p4est3_vtable_t;

/* p4est construction parameters: connectivity, uniform level, etc. */
/* While we're not ready defining the connectivity, use abstract trees. */

/** Check whether a forest virtual table is valid, thus ready to use.
 * \param [in] pvt      Any pointer.  NULL is considered not valid.
 * \param [out] reason  May be NULL.  Otherwise, will be filled with the
 *                      empty string on validity or the issue found otherwise.
 * \return              Boolean value.
 */
int                 p4est3_vtable_is_valid (const p4est3_vtable_t * pvt,
                                            char *reason);

/** The forest is an opaque structure. */
typedef struct p4est3 p4est3_t;

/** Pass context information about a local element to decide for refinement. */
typedef struct p4est3_refine_callback_info
{
  /* these variables stay the same during the refine process */
  p4est3_t           *p3;               /**< Pointer to the forest */
  p4est3_quadrant_vtable_t *qvt;        /**< Pointer to the quadrant virtual
                                             table of current implementation */
  void               *user_data;        /**< For convenience, the user data */

  /* these variables are specific to each quadrant asked for refinement */
  p4est3_topidx       ntree;            /**< Number of tree of quadrant */
  void               *quadrant;         /**< Pointer to the quadrant that
                                             may be refined */
}
p4est3_refine_callback_info_t;

/** Document this. */
typedef struct p4est3_coarsen_callback_info
{
  /* these variables stay the same during the refine process */
  p4est3_t           *p3;               /**< Pointer to the forest */
  p4est3_quadrant_vtable_t *qvt;        /**< Pointer to the quadrant virtual
                                             table of current implementation */
  void               *user_data;        /**< For convenience, the user data */

  /* these variables are specific to each quadrant asked for refinement */
  p4est3_topidx       ntree;            /**< Number of tree of quadrant family */
  sc3_array_t        *family;           /**< Array of quadrants that represent
                                             a family that may be coarsened */
}
p4est3_coarsen_callback_info_t;

/** Document this. */
typedef             sc3_error_t
  * (*p4est3_refine_callback_t) (p4est3_refine_callback_info_t * ci,
                                 int *is_refine);

/** Document this. */
typedef             sc3_error_t
  * (*p4est3_coarsen_callback_t) (p4est3_coarsen_callback_info_t * ci,
                                  int *is_coarsen);

/** Check whether a forest is valid (no matter if setup or not).
 * \param [in] p3       Forest pointer.  NULL is considered not valid.
 * \param [out] reason  May be NULL.  Otherwise, will be filled with the
 *                      empty string on validity or the issue found otherwise.
 * \return              Boolean value.
 */
int                 p4est3_is_valid (const p4est3_t * p3, char *reason);

/** Check whether a forest is valid and not yet setup.
 * \param [in] p3       Forest pointer.  NULL is considered not valid.
 * \param [out] reason  May be NULL.  Otherwise, will be filled with the
 *                      empty string on validity or the issue found otherwise.
 * \return              Boolean value.
 */
int                 p4est3_is_new (const p4est3_t * p3, char *reason);

/** Check whether a forest is valid and setup.
 * \param [in] p3      Forest pointer.  NULL is considered not valid.
 * \param [out] reason  May be NULL.  Otherwise, will be filled with the
 *                      empty string on validity or the issue found otherwise.
 * \return              Boolean value.
 */
int                 p4est3_is_setup (const p4est3_t * p3, char *reason);

/** Begin life cycle of a forest object.
 * \param [in,out] alloc    Allocator must be setup.  It is referenced
 *                          and kept around while forest is live.
 * \param [out] pp3         Pointer to a pointer, the latter will be updated.
 * \return                  NULL on success, error object otherwise.
 */
sc3_error_t        *p4est3_new (sc3_allocator_t * alloc, p4est3_t ** pp3);

/** Select the virtual table creation method for the forest.
 * This overrides all other \c p4est3_set_* calls made before or after.
 * \param [in,out] p3   Forest under construction.
 * \param [in] pvt      Valid forest virtual table.
 *                      We make a shallow copy, that is, copy all elements.
 *                      We call \ref p4est3_set_comm with duping true.
 *                      We call \ref p4est3_set_connectivity and \ref
 *                      p4est3_set_quadrant_vtable with the table contents.
 *                      It is thus safe if \c *pvt lives on the stack.
 * \param [in] slf      Self (state) of virtual forest passed along.
 * \return              NULL on success, error object otherwise.
 */
sc3_error_t        *p4est3_set_vtable (p4est3_t * p3,
                                       p4est3_vtable_t * pvt, void *slf);

/** Provide an MPI communicator to use.
 * The default after \ref p4est3_new is \c SC3_MPI_COMM_WORLD.
 * \param [in,out] p3       The forest must not have been setup.
 * \param [in] comm         This communicator replaces any previous one.
 *                          If it is dupd, we also set it to return errors.
 *                          We free a previous communicator when dupd.
 * \param [in] dup          If true, the input communicator is dupd
 *                          right here and set to return errors.
 * \return                  NULL on success, error object otherwise.
 */
sc3_error_t        *p4est3_set_comm (p4est3_t * p3,
                                     sc3_MPI_Comm_t comm, int dup);

/** Enable/disable use of MPI shared memory.
 * \param [in,out] p3       The forest must not have been setup.
 * \param [in] shared       True value indicates enabling, while false
 *                          is for disabling of MPI-3 shared memory.
 *                          Default value is false.
 * \return                  NULL on success, error object otherwise.
*/
sc3_error_t        *p4est3_set_shared (p4est3_t * p3, int shared);

/** Provide a connectivity to be used in creating the forest.
 * TODO: set 2D unit square as default.
 * This function is mandatory to call at least once before \ref p4est3_setup.
 * \param [in,out] p3       Forest object under construction.
 * \param [in] conn         Connectivity structure must be setup.
 * \return                  NULL on success, error object otherwise.
 */
sc3_error_t        *p4est3_set_connectivity (p4est3_t * p3,
                                             p4est3_connectivity_t * conn);

/** Set a virtual quadrant implementation to use in the forest.
 * TODO: set standard virtual table as default.
 * \param [in,out] p3       Forest under construction.
 * \param [in] qvt          Valid virtual quadrant table.
 *                          We make a shallow copy, that is, copy all elements.
 *                          It is thus safe if \c *qvt lives on the stack.
 * \return                  NULL on success, error object otherwise.
 */
sc3_error_t        *p4est3_set_quadrant_vtable
  (p4est3_t * p3, const p4est3_quadrant_vtable_t * qvt);

/** Set minimum refinement level on creation of the forest.
 * The default after \ref p4est3_new is 0.
 * \param [in,out] p3       Forest under construction.
 * \param [in] level        Range must match quadrant virtual table.
 * \return                  NULL on success, error object otherwise.
 */
sc3_error_t        *p4est3_set_level (p4est3_t * p3, int level);

/** Provide a forest to be used in setting up a new one
 * and a mode to get this forest.
 * \param [in,out] p3       New forest object under construction.
 * \param [in] old          Source forest object that data will be used on
 *                          the setting up stage.
 * \param [in] mode         See \ref p4est3_source_setup_t type for
 *                          available options. Default value is
 *                          P4EST3_COPY_MODE.
 * \return                  NULL on success, error object otherwise.
 */
sc3_error_t        *p4est3_set_source (p4est3_t * p3, p4est3_t * old);

/** Provide a function and data to be used as refinement contition.
 * Must be assign to the source forest.
 * \param [in,out] p3       New forest object under construction.
 * \param [in] crefine      Callback function prototype to decide
 *                          for refinement.  NULL value is possible,
 *                          in this case refinement decision is always false.
 * \param [in] user_data    Pointer passed to refine callback context as is.
 * \return                  NULL on success, error object otherwise.
 */
sc3_error_t        *p4est3_set_refine (p4est3_t * p3,
                                       p4est3_refine_callback_t crefine);

/** Provide a function and data to be used as coarsening contition.
 * Must be assigned to the source forest.
 * \param [in,out] p3       New forest object under construction.
 * \param [in] ccoarse      Callback function prototype to decide
 *                          for coarsening.  NULL value is possible,
 *                          in this case coarsening decision is always false.
 * \param [in] user_data    Pointer passed to coarsen callback context as is.
 * \return                  NULL on success, error object otherwise.
 */
sc3_error_t        *p4est3_set_coarsen (p4est3_t * p3,
                                        p4est3_coarsen_callback_t ccoarse);

/** Indicator to store quadrant's family within the same process.
 * \param [in,out] p3       The forest must not have been setup.
 * \param [in] is_family    True value indicates not splitting quadrant family
 *                          between ranks. Default value is false.
 * \return                  NULL on success, error object otherwise.
*/
sc3_error_t        *p4est3_set_family (p4est3_t * p3, int is_family);

/** TODO document */
sc3_error_t        *p4est3_set_user_data (p4est3_t * p3,
                                          void *user_data);

/** Finalize construction of a forest.
 * Afterwards, no more \c p4est3_set_* functions may be called.
 * \param [in,out] p3      Forest under construction will be finalized.
 * \return                 NULL on success, error object otherwise.
 */
sc3_error_t        *p4est3_setup (p4est3_t * p3);

/** Increase reference counter of a forest after setup.
 * \param [in,out] p3       Must be setup.  Increase its reference counter.
 * \return                  NULL on success, error object otherwise.
 */
sc3_error_t        *p4est3_ref (p4est3_t * p3);

/** Decrease reference counter of a forest after setup.
 * The lowest legal value for the reference counter is one, as after setup.
 * This function never destroys the object: it is not legal to unref below one.
 * The only way to deallocate a forest is \ref p4est3_destroy.
 * \param [in,out] p3       Must be setup and have reference counter greater one.
 *                          Decrease its reference counter.
 *                          When the count reaches one, nothing happens.
 * \return                  NULL on success, error object otherwise.
 */
sc3_error_t        *p4est3_unref (p4est3_t * p3);

/** Destroy a forest that must be valid, but may or may not be setup.
 * It must have a reference count of exactly one.  Otherwise we return an error.
 * Thus, all additional references must have been dropped before calling.
 * \param [in,out] pp3      Valid forest with one reference.  NULL on output.
 * \return                  NULL on success, error object otherwise.
 */
sc3_error_t        *p4est3_destroy (p4est3_t ** pp3);

/** Retrieve the connectivity registered with the forest.
 * It must be returned to the forest before the forest is destructed.
 * To this end, use \ref p4est3_restore_connectivity.
 * We allow an arbitrary number of simultaneous or staggered accesses.
 * \param [in,out] p3       Must be setup.  We increment its access count.
 * \param [out] pconn       Non-NULL reference argument.
 * \return                  NULL on success, error object otherwise.
 */
sc3_error_t        *p4est3_access_connectivity (p4est3_t * p3,
                                                p4est3_connectivity_t **
                                                pconn);

/** Release a connectivity previously obtained
 * with \ref p4est3_access_connectivity.
 * Every connectivity access must be restored before \ref p4est3_destroy.
 * \param [in,out] p3       Must be setup.  We decrement the access count.
 * \param [out] conn        Same pointer as passed to corresponding access.
 * \return                  NULL on success, error object otherwise.
 */
sc3_error_t        *p4est3_restore_connectivity (p4est3_t * p3,
                                                 p4est3_connectivity_t *
                                                 conn);

/** Query the range of processor-local trees of this forest.
 * \param [in] p3                   Initialized, valid forest.
 * \param [out] first_local_tree    The first local tree, or -1 if empty.
 *                          Pointer to this output variable must not be NULL.
 * \param [out] last_local_tree     The last local tree (inclusive), or -2.
 *                          Pointer to this output variable must not be NULL.
 * \return                  NULL on success, error object otherwise.
 */
sc3_error_t        *p4est3_get_local_num_trees (const p4est3_t * p3,
                                                p4est3_topidx *
                                                first_local_tree,
                                                p4est3_topidx *
                                                last_local_tree);

/*----------------------- accessing quadrants ------------------------*/

/* TODO: think about this interface */
sc3_error_t        *p4est3_get_quadrants (const p4est3_t * p3, char **q);

sc3_error_t        *p4est3_get_global_num_quads (const p4est3_t * p3,
                                                 p4est3_gloidx * n);
sc3_error_t        *p4est3_get_local_num_quads (const p4est3_t * p3,
                                                p4est3_locidx * n);

sc3_error_t        *p4est3_get_user_data (p4est3_t * p3,
                                          void **user_data);

/*-------------------- working with quadrants ------------------------*/

/* TODO: move quadrant functions into a new file p4est3_quadrant.{c,h} */
#if 0
int                 p4est3_quadrant_is2_valid (p4est3_t * p3,
                                               const void *q, char *reason);
#endif
#if 0

typedef struct p4est3_access_attr p4est3_access_attr_t;
typedef struct p4est3_access p4est3_access_t;

p4est3_access_t    *p4est3_access_new (p4est3_t * p3,
                                       p4est3_access_attr_t * ia);
void                p4est3_access_ref (p4est3_access_t * a3);
void                p4est3_access_unref (p4est3_access_t ** a3);
void                p4est3_access_destroy (p4est3_access_t ** a3);

p4est3_locidx_t     p4est3_access_get_length (p4est3_access_t * a3);
p4est3_quadrant_t  *p4est3_access_get_begin (p4est3_access_t * a3);
p4est3_quadrant_t  *p4est3_access_index (p4est3_access_t * a3,
                                         p4est3_locidx_t li);

#endif /* 0 */

#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif /* !P4EST3_H */
