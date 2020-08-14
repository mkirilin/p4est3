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

#ifndef P4EST3_H
#define P4EST3_H

#include <sc3_mpi.h>
#include <p4est3_base.h>
#include <p4est3_connectivity.h>
#include <p4est3_quadrant_vtable.h>

#ifdef __cplusplus
extern              "C"
{
#if 0
}
#endif
#endif

typedef struct p4est3 p4est3_t;

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

/* p4est construction parameters: connectivity, uniform level, etc. */
/* While we're not ready defining the connectivity, use abstract trees. */

int                 p4est3_is_valid (const p4est3_t * p3, char *reason);
int                 p4est3_is_new (const p4est3_t * p3, char *reason);
int                 p4est3_is_setup (const p4est3_t * p3, char *reason);

sc3_error_t        *p4est3_new (sc3_allocator_t * alloc, p4est3_t ** pp3);

/** Provide a communicator to use.
 * \param [in,out] p3       The forest must not have been setup.
 * \param [in] comm         This communicator replaces any previous one.
 *                          If it is dupd, we also set it to return errors.
 * \param [in] dup          If true, the input communicator is dupd
 *                          and set to return errors.
 */
sc3_error_t        *p4est3_set_comm (p4est3_t * p3,
                                     sc3_MPI_Comm_t comm, int dup);
sc3_error_t        *p4est3_set_connectivity (p4est3_t * p3,
                                             p4est3_connectivity_t * conn);
sc3_error_t        *p4est3_set_vtable (p4est3_t * p3,
                                       p4est3_quadrant_vtable_t * qvt);
sc3_error_t        *p4est3_set_level (p4est3_t * p3, int level);

/* TODO: document default value for all _set_ */
/** Set a way that creates quadrants in a tree in a setup p4est3 phase.
 * \param [in,out] p3       The forest must not have been setup.
 * \param [in] mode         See \ref p4est3_setup_mode_t type for
 *                          available options. Default value is
 *                          P4EST3_NEW_MORTON.
 */
sc3_error_t        *p4est3_set_setup_mode (p4est3_t * p3,
                                           p4est3_setup_mode_t mode);
sc3_error_t        *p4est3_setup (p4est3_t * p3);

sc3_error_t        *p4est3_ref (p4est3_t * p3);
sc3_error_t        *p4est3_unref (p4est3_t ** pp3);
sc3_error_t        *p4est3_destroy (p4est3_t ** pp3);

/*----------------------- accessing quadrants ------------------------*/

/* TODO: think about this interface */
sc3_error_t        *p4est3_get_quadrants (const p4est3_t * p3, char **q);

sc3_error_t        *p4est3_get_global_num_quads (const p4est3_t * p3,
                                                 p4est3_gloidx * n);
sc3_error_t        *p4est3_get_local_num_quads (const p4est3_t * p3,
                                                p4est3_locidx * n);

#if 0

p4est3_quadrant_t  *p4est3_quadrant_range (p4est3_t * p3,
                                           p4est3_topidx_t tbegin,
                                           p4est3_topidx_t tend,
                                           p4est3_locidx_t qbegin,
                                           p4est3_locidx_t qend);
void                p4est3_quadrant_restore (p4est3_t * p3,
                                             p4est3_quadrant_t * q3,
                                             p4est3_topidx_t tbegin,
                                             p4est3_topidx_t tend,
                                             p4est3_locidx_t qbegin,
                                             p4est3_locidx_t qend);

typedef struct p4est3_access_attr p4est3_access_attr_t;
typedef struct p4est3_access p4est3_access_t;

p4est3_access_t    *p4est3_iter_new (p4est3_t * p3,
                                     p4est3_access_attr_t * ia);
p4est3_quadrant_t  *p4est3_iter_end (p4est3_access_t * pi);
void                p4est3_iter_inc (p4est3_access_t * pi);

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
