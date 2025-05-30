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

/* translate 2D definitions to 3D for reuse of a single code file */
#include <p4est_to_p8est.h>

/* these definitions are private and thus not in the above header */
#define p4est_dune_numbers_t            p8est_dune_numbers_t
#define p4est_dune_numbers_params_t     p8est_dune_numbers_params_t
#define p4est_dune_numbers_params_init  p8est_dune_numbers_params_init
#define p4est_dune_numbers_new          p8est_dune_numbers_new
#define p4est_dune_numbers_destroy      p8est_dune_numbers_destroy

/* one code file to generate independent objects for 2D and 3D */
#include "p4est_dune.c"
