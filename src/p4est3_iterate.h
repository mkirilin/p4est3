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

/** \file p4est3_iterate.h
 * Iterate through the leaves of a forest as well as the neighbor connections.
 *
 * \ingroup p4est3
 */

#ifndef P4EST3_ITERATE_H
#define P4EST3_ITERATE_H

#include <p4est3.h>

#ifdef __cplusplus
extern              "C"
{
#if 0
}
#endif
#endif

/** Pass context information about a local element to the quadrant iteration. */
typedef struct p4est3_iterate_volume_info
{
  p4est3_t           *p3;               /**< Pointer to the forest */
  void               *user_data;        /**< Passed into \a p4est_iterate_* */
  void               *quadrant;         /**< Pointer to iterated quadrant */
  p4est3_topidx       ntree;            /**< Number of tree of quadrant */
  p4est3_locidx       nquad;            /**< Local index of quadrant */
}
p4est3_iterate_volume_info_t;

/** Callback function invoked for every local quadrant during iteration.
    Guaranteed to be called in ascending tree and then quadrant order. */
typedef sc3_error_t *(*p4est3_iterate_volume_t) (p4est3_iterate_volume_info_t
                                                 * vi);

/** Context representing each one of the quadrants at a face connection. */
typedef struct p4est3_iterate_face_side
{
  p4est3_topidx       ntree;            /**< Number of tree of quadrant */
  p4est3_locidx       nquad;            /**< Local index of quadrant */
  int                 nface;            /**< Number of face at connection */
  int                 is_ghost;         /**< Is this side a ghost quadrant */
  void               *quadrant;         /**< Pointer to quadrant itself */
}
p4est3_iterate_face_side_t;

/** Pass context information about face connection to the quadrant iteration. */
typedef struct p4est3_iterate_face_info
{
  p4est3_t           *p3;               /**< Pointer to the forest */
  void               *user_data;        /**< Passed into \a p4est_iterate_* */
  int                 orientation;      /**< Orientation of the sides relative
                                            to each other, as in the definition
                                            of p4est_connectivity_t */
  int                 tree_boundary;    /**< Boolean */
  sc3_array_t        *sides;            /**< Array of 1 or 2 elements of type
                                    p4est3_iterate_face_side_t.  At a domain
                                    boundary, it's 1, otherwise it's 2.  In the
                                    latter case ordered ascending by quadrant. */
}
p4est3_iterate_face_info_t;

/** Callback function invoked for every face connection during iteration.
 * Guaranteed to be called in ascending tree, quadrant and face order
 * as seen from the lower-ordered of the two connecting quadrants.
 * For hanging face connections, it is called separately for each small face
 * repeating the large face.
 */
typedef sc3_error_t *(*p4est3_iterate_face_t) (p4est3_iterate_face_info_t *
                                               fi);

/** Context representing each one of edge-/corner-connecting quadrants. */
typedef struct p4est3_iterate_codim_side
{
  p4est3_topidx       ntree;            /**< Number of tree of quadrant */
  p4est3_locidx       nquad;            /**< Local index of quadrant */
  int                 nbound;           /**< Number of quadrant's connecting
                                             boundary entity */
  int                 orientation;      /**< Relative orientation around the
                                            connecting entity for non-corners
                                            (3D edges), as in the definition
                                            of p4est_connectivity_t.
                                            It's 0 for corners */
  int                 is_ghost;         /**< Is this side a ghost quadrant */
  void               *quadrant;         /**< Pointer to quadrant itself */
}
p4est3_iterate_codim_side_t;

/** Pass context information about higher codimension connection to iteration. */
typedef struct p4est3_iterate_codim_info
{
  p4est3_t           *p3;               /**< Pointer to the forest */
  void               *user_data;        /**< Passed into \a p4est_iterate_* */
  int                 codimension;      /**< Codimension of connection */
  int                 tree_boundary;    /**< Boolean */
  sc3_array_t        *sides;            /**< One entry per neighbor across
                                             codimension boundary entity. */
}
p4est3_iterate_codim_info_t;

/** Callback function invoked for every edge and corner connection.
 * Guaranteed to be called in ascending tree, quadrant and edge/corner order
 * as seen from the lower-ordered of the connecting quadrants.
 * For hanging face connections, it is called separately for each small
 * boundary entity (only relevant for edges), repeating the quadrants
 * that are larger neighbors around a parent of the boundary entity.
 */
typedef sc3_error_t *(*p4est3_iterate_codim_t) (p4est3_iterate_codim_info_t *
                                                fi);

/** Iterate through the forest for volumes and face connections.
 * \param [in] p3       Forest passed for reference.
 * \param [in] cvolume  Volume callback called for every local quadrant.
 *                      Ignored if NULL.
 * \param [in] cface    Callback for every face connection involving
 *                      local quadrants.  Ignored if NULL.
 * \param [in,out] user_data        Passed through to the callbacks.
 * \return              NULL on success, error object otherwise.
 */
sc3_error_t        *p4est3_iterate_face (p4est3_t * p3,
                                         p4est3_iterate_volume_t * cvolume,
                                         p4est3_iterate_face_t * cface,
                                         void *user_data);

/** Iterate through the forest for connections of any codimension.
 * \param [in] p3       Forest passed for reference.
 * \param [in] codims   Binary OR of bit 0 (volume), 1 (face),
 *                      2 (2D: corner; 3D: edge), 3 (3D: corner).
 * \param [in] cvolume  Volume callback called for every local quadrant
 *                      when volumes enabled in \a codims.  Ignored if NULL.
 * \param [in] cface    Callback for every face connection involving
 *                      local quadrants when faces enabled in \a codims.
 *                      Ignored if NULL.
 * \param [in] ccodim   Callback for every non-face connection involving
 *                      local quadrants when enabled by \a codims.
 *                      Ignored if NULL.
 * \param [in,out] user_data        Passed through to the callbacks.
 * \return              NULL on success, error object otherwise.
 */
sc3_error_t        *p4est3_iterate_codim (p4est3_t * p3, int codims,
                                          p4est3_iterate_volume_t * cvolume,
                                          p4est3_iterate_face_t * cface,
                                          p4est3_iterate_codim_t * ccodim,
                                          void *user_data);

#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif /* !P4EST3_ITERATE_H */
