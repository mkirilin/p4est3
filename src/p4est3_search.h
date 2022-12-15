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

#ifndef P4EST_SEARCH_H
#define P4EST_SEARCH_H

#include <p4est3.h>

/** Binary search in partition array.
 * Given two targets \a my_begin and \a my_end, find offsets such that
 * `search_in[begin] >= my_begin`, `my_end <= search_in[end]`.
 * If more than one index satisfies the conditions, then the minimal index is the
 * result. If there is no index that satisfies the conditions, then \a begin
 * and \a end are tried to set equal such that `search_in[begin] >= my_end`.
 * If \a my_begin is less or equal than the smallest value of \a search_in
 * \a begin is set to 0 and if \a my_end is bigger or equal than the largest
 * value of \a search_in \a end is set to \a num_entities - 1.
 * If none of the above conditions is satisfied, the output is not well defined.
 * We require `my_begin <= my_begin'.
 * \param [in] alloc        Valid allocator to setup temporary arrays.
 *                          Must be setup.
 * \param [in] num_entities Number of entities to get the length of
 *                          \a search_in.
 * \param [in] search_in    The sorted array (ascending) in that the function
 *                          will search.
 *                          If `k` indexes search_in, then
 *                          `0 <= k < num_procs`.
 * \param [in] my_begin     The first target that defines the start of the
 *                          search window.
 * \param [in] my_end       The second target that defines the end of the
 *                          search window.
 * \param [in,out] begin    The first offset such that
 *                          `search_in[begin] >= my_begin`.
 * \param [in,out] end      The second offset such that
 *                          `my_end <= search_in[end]`.
 */
sc3_error_t        *p4est3_find_partition (sc3_allocator_t * alloc,
                                           const int num_entities,
                                           p4est3_gloidx * search_in,
                                           p4est3_gloidx my_begin,
                                           p4est3_gloidx my_end,
                                           p4est3_gloidx * begin,
                                           p4est3_gloidx * end);

#endif /* !P4EST_SEARCH_H */