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

#ifndef P4EST3_BASE_H
#define P4EST3_BASE_H

#include <sc3_alloc.h>
#include <sc3_error.h>
#include <p4est3_config.h>

#define P3A_IS(f,o) SC3A_IS (f,o)
#ifndef P4EST_ENABLE_DEBUG
#define P3A_IS2(f,o,p) SC3_NOOP
#else
#define P3A_IS2(f,o,p) do {                                             \
  char _r[SC3_BUFSIZE];                                                 \
  if (!(f ((o), (p), _r))) {                                            \
    char _errmsg[SC3_BUFSIZE];                                          \
    sc3_snprintf (_errmsg, SC3_BUFSIZE,                                 \
                  "%s(%s,%s): %s", #f, #o, #p, _r);                     \
    return sc3_error_new_bug (__FILE__, __LINE__, _errmsg);             \
  }} while (0)
#endif

/* TODO: P3E_TAIL will be removed */
#define P3E_TAIL(f) do { SC3E (f); return NULL; } while (0)

/* TODO: replace with versions with different message prefix */
#define P3A_CHECK(x) SC3A_CHECK (x)
#define P3A_STACK(x) SC3A_STACK (x)
#define P3E(f) SC3E (f)
#define P3E_DEMAND(f) SC3E_DEMAND (f)

typedef int         p4est3_topidx;
#define P4EST3_TOPIDX_MAX INT_MAX

typedef int         p4est3_locidx;
#define P4EST3_LOCIDX_MAX INT_MAX
#define p4est3_loccut sc3_intcut

typedef long        p4est3_gloidx;
#define P4EST3_GLOIDX_MAX LONG_MAX
#define p4est3_glopow sc3_longpow
#define p4est3_glocut sc3_longcut

#ifdef __cplusplus
extern              "C"
{
#if 0
}
#endif
#endif

/* no function prototypes yet */

#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif /* !P4EST3_BASE_H */
