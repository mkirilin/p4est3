/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2010 The University of Texas System
  Additional copyright (C) 2011 individual authors
  Written by Carsten Burstedde, Lucas C. Wilcox, and Tobin Isaac

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

/** \file p8est_dune.h
 *
 * Provide various functions to interface dune to a p4est mesh.
 *
 * For example, populate element corner and face number tables, where the
 * range for each codimension begins with 0 and is contiguous.
 *
 * We also provide a face-only variant of the iterator that visits all small
 * quadrants on one side of a hanging face from separate callbacks in order.
 *
 * \ingroup p8est
 */

#ifndef P8EST_DUNE_H
#define P8EST_DUNE_H

#include <p8est_iterate.h>

SC_EXTERN_C_BEGIN;

/** Parameters to the constructor of the dune number tables. */
typedef struct p8est_dune_numbers_params
{
  /** Determine which codimensions are numbered. */
  p8est_connect_type_t ctype;
}
p8est_dune_numbers_params_t;

/** Element corner, edge and face numbers to interface to the DUNE library.
 *
 * The element corner, edge and face arrays hold local indices that are
 * unique within this process among all faces, and likewise unique within
 * this process among all edges, and likewise all corners.  The same number
 * may occur in two or even three of the arrays, but it refers to a
 * different mesh object each instance.
 *
 * The indices are numbered starting from zero.
 *
 * The array entries are of type \ref p4est_locidx_t.
 */
typedef struct p8est_dune_numbers
{
  /** Parameters the numbers are built with are copied by value. */
  p8est_dune_numbers_params_t params;

  /** Highest number (exclusive) among all element corner entries. */
  p4est_locidx_t      num_corner_numbers;

  /** For each local element corner a process-local index. */
  sc_array_t         *element_corners;

  /** Highest number (exclusive) among all element edge entries. */
  p4est_locidx_t      num_edge_numbers;

  /** For each local element edge a process-local index. */
  sc_array_t         *element_edges;

  /** Highest number (exclusive) among all element face entries. */
  p4est_locidx_t      num_face_numbers;

  /** For each local element face a process-local index. */
  sc_array_t         *element_faces;
}
p8est_dune_numbers_t;

/** Set default parameters to pass to \ref p8est_dune_numbers_new.
 * \param [out] params  Pointer must not be NULL.
 *                      The structure is filled with default values.
 */
void                p8est_dune_numbers_params_init
  (p8est_dune_numbers_params_t * params);

/** Create lookup tables for unique corners, edges and faces.
 * Hanging corners, edges and faces are included as independent entities.
 * All numbers are process-local.
 *
 * \param [in] p8est    Required input parameter is not modified.  It must
 *                      be balanced at least to \ref P8EST_CONNECT_ALMOST.
 * \param [in] ghost    The ghost layer must have been generated with
 *                      \ref p8est_ghost_new using the same \c p8est and
 *                      the parameter \ref P8EST_CONNECT_FULL.
 * \param [in] params   Further parameters to control the mode of operation.
 *                      When passing NULL, the behavior is identical to using
 *                      defaults by \ref p8est_dune_numbers_params_init.
 * \return              A fully initialized dune node numbering.
 *                      Deallocate with \ref p8est_dune_numbers_destroy.
 */
p8est_dune_numbers_t *p8est_dune_numbers_new (p8est_t * p8est,
                                              p8est_ghost_t * ghost,
                                              const
                                              p8est_dune_numbers_params_t *
                                              params);

/** Destroy the dune element number tables previously generated.
 * \param [in] dn       Valid dune number tables from \ref
 *                      p8est_dune_numbers_new are deallocated.
 */
void                p8est_dune_numbers_destroy (p8est_dune_numbers_t * dn);

/** Execute user supplied callbacks at every local volume and face.
 *
 * Execute the user-supplied callback functions at every volume and face in
 * the local forest.  When multiple small quadrants are on the other side of
 * the face of a large quadrant, the callback is executed for every small
 * face neighbor in turn, with the same large quadrant on the other side.
 * The \a user_data pointer is not touched by the iterator, only passed to
 * each of the callbacks.  Any of the callbacks may be NULL.
 * The callback functions are interspersed with each other, i.e. some face
 * callbacks will occur between volume callbacks:
 *
 * 1. Volume callbacks occur in the sorted Morton-index order.
 * 2. A face callback is not executed until after the volume callbacks have
 *    been executed for the quadrants that share it.
 * 3. The face callbacks for the small elements at the same hanging face
 *    are executed right after each other in their Morton-index order.
 * 4. Callbacks are not executed at faces that only involve ghost
 *    quadrants, i.e. that are not adjacent in the local section of the
 *    forest.
 *
 * The convention for the context information passed to the callbacks
 * is the same as for \ref p8est_iterate with the following exception:
 * For a hanging face, we call the callback anew for every small neighbor
 * quadrant with the same large neighbor quadrant.  The hanging side data
 * (see \ref p8est_iter_face_side_hanging_t) is populated at index [0] each
 * time, and we store the number of the small quadrant within the order of
 * its face siblings in the variable is.hanging.quadid[1], values in [0, 4).
 *
 * \param[in] p4est          The forest to iterate over.
 * \param[in] ghost_layer    Required valid ghost structure.
 * \param[in,out] user_data  optional context to supply to each callback.
 * \param[in] iter_volume    callback function for every quadrant interior.
 * \param[in] iter_face      callback function for every face between.
 */
void                p8est_dune_iterate (p8est_t * p4est,
                                        p8est_ghost_t * ghost_layer,
                                        void *user_data,
                                        p8est_iter_volume_t iter_volume,
                                        p8est_iter_face_t iter_face);

SC_EXTERN_C_END;

#endif /* P8EST_DUNE_H */
