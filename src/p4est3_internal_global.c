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

#include <sc3_refcount.h>
#include <sc3_mpienv.h>
#include <p4est3_internal.h>

static int
p4est3_glopart_is_valid (const p4est3_glopart_t * m, char * reason)
{
  SC3E_TEST (m != NULL, reason);
  SC3E_IS (sc3_refcount_is_valid, &m->rc, reason);

  if (!m->setup) {
    SC3E_TEST (m->gfpos == NULL, reason);
    SC3E_TEST (m->gftree == NULL, reason);
  }
  else {
    SC3E_IS (sc3_mpienv_is_setup, &m->mpienv, reason);
    SC3E_TEST (m->gfpos != NULL, reason);
    SC3E_TEST (m->gftree != NULL, reason);
    SC3E_TEST (m->gftree > 0, reason);
  }
  SC3E_YES (reason);
}

static int
p4est3_glooffs_is_valid (const p4est3_glooffs_t * m, char * reason)
{
  SC3E_TEST (m != NULL, reason);
  SC3E_IS (sc3_refcount_is_valid, &m->rc, reason);

  if (!m->setup) {
    SC3E_TEST (m->goffset == NULL, reason);
  }
  else {
    SC3E_IS (sc3_mpienv_is_setup, &m->mpienv, reason);
    SC3E_TEST (m->goffset != NULL, reason);
  }
  SC3E_YES (reason);
}

static int
p4est3_glopart_is_new (const p4est3_glopart_t * m, char *reason)
{
  SC3E_IS (p4est3_glopart_is_valid, m, reason);
  SC3E_TEST (!m->setup, reason);
  SC3E_YES (reason);
}

static int
p4est3_glooffs_is_new (const p4est3_glooffs_t * m, char *reason)
{
  SC3E_IS (p4est3_glooffs_is_valid, m, reason);
  SC3E_TEST (!m->setup, reason);
  SC3E_YES (reason);
}

static int
p4est3_glopart_is_setup (const p4est3_glopart_t * m, char *reason)
{
  SC3E_IS (p4est3_glopart_is_valid, m, reason);
  SC3E_TEST (m->setup, reason);
  SC3E_YES (reason);
}

static int
p4est3_glooffs_is_setup (const p4est3_glooffs_t * m, char *reason)
{
  SC3E_IS (p4est3_glooffs_is_valid, m, reason);
  SC3E_TEST (m->setup, reason);
  SC3E_YES (reason);
}

sc3_error_t        *
p4est3_glopart_new (sc3_allocator_t * mator, p4est3_glopart_t ** mp)
{
  p4est3_glopart_t *m;

  SC3E_RETVAL (mp, NULL);
  SC3A_IS (sc3_allocator_is_setup, mator);

  SC3E (sc3_allocator_ref (mator));
  SC3E (sc3_allocator_calloc_one (mator, sizeof (p4est3_glopart_t), &m));
  SC3E (sc3_refcount_init (&m->rc));
  m->gfpos = NULL;
  m->gftree = NULL;
  m->qsize = 0;

  SC3A_IS (p4est3_glopart_is_new, m);
  *mp = m;
  return NULL;
}

sc3_error_t        *
p4est3_glooffs_new (sc3_allocator_t * mator, p4est3_glooffs_t ** mp)
{
  p4est3_glooffs_t *m;

  SC3E_RETVAL (mp, NULL);
  SC3A_IS (sc3_allocator_is_setup, mator);

  SC3E (sc3_allocator_ref (mator));
  SC3E (sc3_allocator_calloc_one (mator, sizeof (p4est3_glooffs_t), &m));
  SC3E (sc3_refcount_init (&m->rc));
  m->goffset = NULL;

  SC3A_IS (p4est3_glooffs_is_new, m);
  *mp = m;
  return NULL;
}

sc3_error_t *
p4est3_glopart_set_mpienv (sc3_mpienv_t * mpienv, p4est3_glopart_t * m)
{
  SC3A_IS (p4est3_glopart_is_new, m);
  SC3A_IS (sc3_mpienv_is_setup, mpienv);
  m->mpienv = mpienv;
  return NULL;
}

sc3_error_t *
p4est3_glooffs_set_mpienv (sc3_mpienv_t * mpienv, p4est3_glooffs_t * m)
{
  SC3A_IS (p4est3_glooffs_is_new, m);
  SC3A_IS (sc3_mpienv_is_setup, mpienv);
  m->mpienv = mpienv;
  return NULL;
}

sc3_error_t *
p4est3_glopart_set_qsize (p4est3_glopart_t * m, int qsize)
{
  SC3A_IS (p4est3_glopart_is_new, m);
  SC3A_CHECK (qsize > 0);
  m->qsize = qsize;
  return NULL;
}

sc3_error_t *
p4est3_glopart_setup (p4est3_glopart_t * m)
{
  int                 noderank, nodesize, mpisize;
  int                 dispunit;
  sc3_MPI_Aint_t      gftreebytes, gfposbytes, tempbytes;
  sc3_MPI_Comm_t      nodecomm;
  sc3_MPI_Info_t      info_noncontig;

  SC3A_IS (p4est3_glopart_is_new, m);

  SC3E (sc3_mpienv_get_noderank (m->mpienv, &noderank));
  SC3E (sc3_mpienv_get_nodesize (m->mpienv, &nodesize));
  SC3E (sc3_mpienv_get_mpisize (m->mpienv, mpisize));
  SC3E (sc3_mpienv_get_info_noncont (m->mpienv, &info_noncontig));
  SC3E (sc3_mpienv_get_nodecomm (m->mpienv, &nodecomm));

  SC3A_CHECK (0 <= noderank && noderank < nodesize);
  SC3A_CHECK (m->qsize > 0);

  /* create shared partition arrays */
  gftreebytes = (mpisize + 1) * sizeof (p4est3_topidx);
  SC3E (sc3_MPI_Win_allocate_shared
        (noderank == 0 ? gftreebytes : 0, sizeof (p4est3_topidx),
         info_noncontig, nodecomm, &m->gftree, &m->gftreewin));
  gfposbytes = (mpisize + 1) * m->qsize;
  SC3E (sc3_MPI_Win_allocate_shared
        (noderank == 0 ? gfposbytes : 0, m->qsize,
         info_noncontig, nodecomm, &m->gfpos, &m->gfposwin));

  if (noderank > 0) {
    SC3E (sc3_MPI_Win_shared_query
          (m->gftreewin, 0, &tempbytes, dispunit, &m->gftree));
    SC3A_CHECK (tempbytes >= gftreebytes);
    SC3A_CHECK (dispunit == sizeof (p4est3_topidx));
    SC3A_CHECK (m->gfpos != NULL);
    SC3E (sc3_MPI_Win_shared_query
          (m->gfposwin, 0, &tempbytes, &dispunit, &m->gfpos));
    SC3A_CHECK (tempbytes >= gfposbytes);
    SC3A_CHECK (dispunit == m->qsize);
    SC3A_CHECK (m->gfpos != NULL);
  }
  m->setup = 1;
  return NULL;
}

sc3_error_t *
p4est3_glooffs_setup (p4est3_glooffs_t * m)
{
  int                 noderank, nodesize, mpisize; 
  int                 dispunit;
  sc3_MPI_Aint_t      goffsetbytes, tempbytes;
  sc3_MPI_Comm_t      nodecomm;
  sc3_MPI_Info_t      info_noncontig;

  SC3A_IS (p4est3_glopart_is_new, m);

  SC3E (sc3_mpienv_get_noderank (m->mpienv, &noderank));
  SC3E (sc3_mpienv_get_nodesize (m->mpienv, &nodesize));
  SC3E (sc3_mpienv_get_mpisize (m->mpienv, mpisize));
  SC3E (sc3_mpienv_get_info_noncont (m->mpienv, &info_noncontig));
  SC3E (sc3_mpienv_get_nodecomm (m->mpienv, &nodecomm));

  SC3A_CHECK (0 <= noderank && noderank < nodesize);

  /* create shared partition arrays */
  goffsetbytes = (mpisize + 1) * sizeof (p4est3_gloidx);
  SC3E (sc3_MPI_Win_allocate_shared
        (noderank == 0 ? goffsetbytes : 0, sizeof (p4est3_gloidx),
         info_noncontig, nodecomm, &m->goffset, &m->goffsetwin));

  if (noderank > 0) {
    SC3E (sc3_MPI_Win_shared_query (m->goffsetwin, 0,
                                    &tempbytes, &dispunit, &m->goffset));
    SC3A_CHECK (tempbytes >= goffsetbytes);
    SC3A_CHECK (dispunit == (int) sizeof (p4est3_gloidx));
    SC3A_CHECK (m->goffset != NULL);
  }
  m->setup = 1;
  return NULL;
}

sc3_error_t *
p4est3_glopart_ref (p4est3_glopart_t * m)
{
  SC3E (sc3_refcount_ref (&m->rc));
  return NULL;
}

sc3_error_t *
p4est3_glooffs_ref (p4est3_glooffs_t * m)
{
  SC3E (sc3_refcount_ref (&m->rc));
  return NULL;
}

sc3_error_t *
p4est3_glopart_unref (p4est3_glopart_t ** mp)
{
  int waslast;
  sc3_allocator_t    *mator;
  p4est3_glopart_t *m;
  sc3_error_t        *leak = NULL;

  SC3E_INOUTP (mp, m);
  SC3A_IS (p4est3_glopart_is_valid, m);
  SC3E (sc3_refcount_unref (&m->rc, &waslast));
  if (waslast) {
    *mp = NULL;
    mator = m->mator;
    if (m->setup) {
      /* deallocate data created on setup here */
      SC3E (sc3_MPI_Win_free (&m->gfposwin));
      SC3E (sc3_MPI_Win_free (&m->gftreewin));
    }
    /* deallocate data knonw on setup here */
    SC3L (&leak, sc3_mpienv_unref (&m->mpienv));
    SC3E (sc3_allocator_free (mator, m));
    SC3L (&leak, sc3_allocator_unref (&mator));
  }
  return leak;
}

sc3_error_t *
p4est3_glooffs_unref (p4est3_glooffs_t ** mp)
{
  int waslast;
  p4est3_glopart_t *m;
  sc3_allocator_t    *mator;
  sc3_error_t        *leak = NULL;

  SC3E_INOUTP (mp, m);
  SC3A_IS (p4est3_glooffs_is_valid, m);
  SC3E (sc3_refcount_unref (&m->rc, &waslast));
  if (waslast) {
    *mp = NULL;
    mator = m->mator;
    if (m->setup) {
      /* deallocate data created on setup here */
      SC3E (sc3_MPI_Win_free (&m->gfposwin));
    }
    /* deallocate data knonw on setup here */
    SC3L (&leak, sc3_mpienv_unref (&m->mpienv));
    SC3E (sc3_allocator_free (mator, m));
    SC3L (&leak, sc3_allocator_unref (&mator));
  }
  return leak;
}

sc3_error_t        *
p4est3_glopart_destroy (p4est3_glopart_t ** mp)
{
  sc3_error_t        *leak = NULL;
  p4est3_glopart_t       *m;

  SC3E_INULLP (mp, m);
  SC3L_DEMAND (&leak, sc3_refcount_is_last (&m->rc, NULL));
  SC3L (&leak, sc3_mpienv_unref (&m));

  SC3A_CHECK (m == NULL || leak != NULL);
  return leak;
}

sc3_error_t        *
p4est3_glooffs_destroy (p4est3_glooffs_t ** mp)
{
  sc3_error_t        *leak = NULL;
  p4est3_glooffs_t       *m;

  SC3E_INULLP (mp, m);
  SC3L_DEMAND (&leak, sc3_refcount_is_last (&m->rc, NULL));
  SC3L (&leak, sc3_mpienv_unref (&m));

  SC3A_CHECK (m == NULL || leak != NULL);
  return leak;
}

sc3_error_t        *
p4est3_glopart_get_gftree (const p4est3_glopart_t * m, int **gftree)
{
  SC3A_IS (p4est3_glopart_is_setup, m);
  SC3A_CHECK (gftree != NULL);

  *gftree = m->gftree;
  return NULL;
}

sc3_error_t        *
p4est3_glopart_get_gfpos (const p4est3_glopart_t * m, int **gfpos)
{
  SC3A_IS (p4est3_glopart_is_setup, m);
  SC3A_CHECK (gfpos != NULL);

  *gfpos = m->gfpos;
  return NULL;
}

sc3_error_t        *
p4est3_gloffs_get_goffset (const p4est3_glooffs_t * m, int **goffset)
{
  SC3A_IS (p4est3_glooffs_is_setup, m);
  SC3A_CHECK (goffset != NULL);

  *goffset = m->goffset;
  return NULL;
}