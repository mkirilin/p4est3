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
p4est3_glotree_is_valid (const p4est3_glotree_t * m, char *reason)
{
  SC3E_TEST (m != NULL, reason);
  SC3E_IS (sc3_refcount_is_valid, &m->rc, reason);
  SC3E_IS (sc3_allocator_is_setup, m->mator, reason);

  if (!m->setup) {
    SC3E_TEST (m->gftree == NULL, reason);
  }
  else {
    SC3E_IS (sc3_mpienv_is_setup, &m->mpienv, reason);
    SC3E_TEST (m->gftree != NULL, reason);
  }
  SC3E_YES (reason);
}

static int
p4est3_glopos_is_valid (const p4est3_glopos_t * m, char *reason)
{
  SC3E_TEST (m != NULL, reason);
  SC3E_IS (sc3_refcount_is_valid, &m->rc, reason);
  SC3E_IS (sc3_allocator_is_setup, m->mator, reason);

  if (!m->setup) {
    SC3E_TEST (m->gfpos == NULL, reason);
  }
  else {
    SC3E_IS (sc3_mpienv_is_setup, &m->mpienv, reason);
    SC3E_TEST (m->gfpos != NULL, reason);
    SC3E_TEST (m->qsize > 0, reason);
  }
  SC3E_YES (reason);
}

static int
p4est3_glooffs_is_valid (const p4est3_glooffs_t * m, char *reason)
{
  SC3E_TEST (m != NULL, reason);
  SC3E_IS (sc3_refcount_is_valid, &m->rc, reason);
  SC3E_IS (sc3_allocator_is_setup, m->mator, reason);

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
p4est3_glotree_is_new (const p4est3_glotree_t * m, char *reason)
{
  SC3E_IS (p4est3_glotree_is_valid, m, reason);
  SC3E_TEST (!m->setup, reason);
  SC3E_YES (reason);
}

static int
p4est3_glopos_is_new (const p4est3_glopos_t * m, char *reason)
{
  SC3E_IS (p4est3_glopos_is_valid, m, reason);
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

sc3_error_t        *
p4est3_glotree_new (sc3_allocator_t * mator, p4est3_glotree_t ** mp)
{
  p4est3_glotree_t   *m;

  SC3E_RETVAL (mp, NULL);
  SC3A_IS (sc3_allocator_is_setup, mator);

  SC3E (sc3_allocator_ref (mator));
  SC3E (sc3_allocator_calloc_one (mator, sizeof (p4est3_glotree_t), &m));
  SC3E (sc3_refcount_init (&m->rc));
  m->gftree = NULL;

  SC3A_IS (p4est3_glotree_is_new, m);
  *mp = m;
  return NULL;
}

sc3_error_t        *
p4est3_glopos_new (sc3_allocator_t * mator, p4est3_glopos_t ** mp)
{
  p4est3_glopos_t    *m;

  SC3E_RETVAL (mp, NULL);
  SC3A_IS (sc3_allocator_is_setup, mator);

  SC3E (sc3_allocator_ref (mator));
  SC3E (sc3_allocator_calloc_one (mator, sizeof (p4est3_glopos_t), &m));
  SC3E (sc3_refcount_init (&m->rc));
  m->gfpos = NULL;
  m->qsize = 0;

  SC3A_IS (p4est3_glopos_is_new, m);
  *mp = m;
  return NULL;
}

sc3_error_t        *
p4est3_glooffs_new (sc3_allocator_t * mator, p4est3_glooffs_t ** mp)
{
  p4est3_glooffs_t   *m;

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

sc3_error_t        *
p4est3_glopartition_set_mpienv (p4est3_glotree_t * mt, p4est3_glopos_t * mp,
                                p4est3_glooffs_t * mo, sc3_mpienv_t * mpienv)
{
  SC3A_IS (sc3_mpienv_is_setup, mpienv);
  SC3E_DEMAND (mt != NULL || mp != NULL || mo != NULL,
               "At least one global partition object should not be NULL");
  if (mt != NULL) {
    SC3A_IS (p4est3_glotree_is_new, mt);
    SC3E (sc3_mpienv_ref (mpienv));
    mt->mpienv = mpienv;
  }
  if (mp != NULL) {
    SC3A_IS (p4est3_glopos_is_new, mp);
    SC3E (sc3_mpienv_ref (mpienv));
    mp->mpienv = mpienv;
  }
  if (mo != NULL) {
    SC3A_IS (p4est3_glooffs_is_new, mo);
    SC3E (sc3_mpienv_ref (mpienv));
    mo->mpienv = mpienv;
  }
  return NULL;
}

sc3_error_t        *
p4est3_glopos_set_qsize (p4est3_glopos_t * m, int qsize)
{
  SC3A_IS (p4est3_glotree_is_new, m);
  SC3A_CHECK (qsize > 0);
  m->qsize = qsize;
  return NULL;
}

sc3_error_t        *
p4est3_glopartition_setup (p4est3_glotree_t * mt, p4est3_glopos_t * mp,
                           p4est3_glooffs_t * mo)
{
  int                 noderank, nodesize, mpisize;
  int                 dispunit;
  sc3_MPI_Aint_t      gftreebytes, gfposbytes, goffsetbytes, tempbytes;
  sc3_MPI_Comm_t      nodecomm;
  sc3_MPI_Info_t      info_noncontig;
  sc3_mpienv_t       *mpienv;

  SC3E_DEMAND (mt != NULL || mp != NULL || mo != NULL,
               "At least one global partition object should not be NULL");
  if (mt != NULL) {
    mpienv = mt->mpienv;
  }
  else if (mp != NULL) {
    mpienv = mp->mpienv;
  }
  else {
    mpienv = mo->mpienv;
  }
  SC3A_IS (p4est3_glopos_is_new, mp);
  SC3A_IS (p4est3_glooffs_is_new, mo);

  SC3E (sc3_mpienv_get_noderank (mpienv, &noderank));
  SC3E (sc3_mpienv_get_nodesize (mpienv, &nodesize));
  SC3E (sc3_mpienv_get_mpisize (mpienv, mpisize));
  SC3E (sc3_mpienv_get_info_noncont (mpienv, &info_noncontig));
  SC3E (sc3_mpienv_get_nodecomm (mpienv, &nodecomm));

  SC3A_CHECK (0 <= noderank && noderank < nodesize);

  /* create shared partition arrays */
  if (mt != NULL) {
    SC3A_IS (p4est3_glotree_is_new, mt);
    gftreebytes = (mpisize + 1) * sizeof (p4est3_topidx);
    SC3E (sc3_MPI_Win_allocate_shared
          (noderank == 0 ? gftreebytes : 0, sizeof (p4est3_topidx),
           info_noncontig, nodecomm, &mt->gftree, &mt->gftreewin));
  }
  if (mp != NULL) {
    SC3A_CHECK (mp->qsize > 0);
    gfposbytes = (mpisize + 1) * mp->qsize;
    SC3E (sc3_MPI_Win_allocate_shared
          (noderank == 0 ? gfposbytes : 0, mp->qsize,
           info_noncontig, nodecomm, &mp->gfpos, &mp->gfposwin));
  }
  if (mo != NULL) {
    goffsetbytes = (mpisize + 1) * sizeof (p4est3_gloidx);
    SC3E (sc3_MPI_Win_allocate_shared
          (noderank == 0 ? goffsetbytes : 0, sizeof (p4est3_gloidx),
           info_noncontig, nodecomm, &mo->goffset, &mo->goffsetwin));
  }
  if (noderank > 0) {
    if (mt != NULL) {
      SC3E (sc3_MPI_Win_shared_query
            (mt->gftreewin, 0, &tempbytes, dispunit, &mt->gftree));
      SC3A_CHECK (tempbytes >= gftreebytes);
      SC3A_CHECK (dispunit == sizeof (p4est3_topidx));
      SC3A_CHECK (mt->gftree != NULL);
    }
    if (mp != NULL) {
      SC3E (sc3_MPI_Win_shared_query
            (mp->gfposwin, 0, &tempbytes, &dispunit, &mp->gfpos));
      SC3A_CHECK (tempbytes >= gfposbytes);
      SC3A_CHECK (dispunit == mp->qsize);
      SC3A_CHECK (mp->gfpos != NULL);
    }
    if (mo != NULL) {
      SC3E (sc3_MPI_Win_shared_query
            (mo->goffsetwin, 0, &tempbytes, &dispunit, &mo->goffset));
      SC3A_CHECK (tempbytes >= goffsetbytes);
      SC3A_CHECK (dispunit == (int) sizeof (p4est3_gloidx));
      SC3A_CHECK (mo->goffset != NULL);
    }
  }
  if (mt != NULL) {
    mt->setup = 1;
  }
  if (mp != NULL) {
    mp->setup = 1;
  }
  if (mo != NULL) {
    mt->setup = 1;
  }
  return NULL;
}

sc3_error_t        *
p4est3_glotree_ref (p4est3_glotree_t * m)
{
  SC3E (sc3_refcount_ref (&m->rc));
  return NULL;
}

sc3_error_t        *
p4est3_glopos_ref (p4est3_glopos_t * m)
{
  SC3E (sc3_refcount_ref (&m->rc));
  return NULL;
}

sc3_error_t        *
p4est3_glooffs_ref (p4est3_glooffs_t * m)
{
  SC3E (sc3_refcount_ref (&m->rc));
  return NULL;
}

sc3_error_t        *
p4est3_glotree_unref (p4est3_glotree_t ** mp)
{
  int                 waslast;
  sc3_allocator_t    *mator;
  p4est3_glotree_t   *m;
  sc3_error_t        *leak = NULL;

  SC3E_INOUTP (mp, m);
  SC3A_IS (p4est3_glotree_is_valid, m);
  SC3E (sc3_refcount_unref (&m->rc, &waslast));
  if (waslast) {
    *mp = NULL;
    mator = m->mator;
    if (m->setup) {
      /* deallocate data created on setup here */
      SC3E (sc3_MPI_Win_free (&m->gftreewin));
    }
    /* deallocate data knonw on setup here */
    SC3L (&leak, sc3_mpienv_unref (&m->mpienv));
    SC3E (sc3_allocator_free (mator, m));
    SC3L (&leak, sc3_allocator_unref (&mator));
  }
  return leak;
}

sc3_error_t        *
p4est3_glopos_unref (p4est3_glopos_t ** mp)
{
  int                 waslast;
  sc3_allocator_t    *mator;
  p4est3_glopos_t    *m;
  sc3_error_t        *leak = NULL;

  SC3E_INOUTP (mp, m);
  SC3A_IS (p4est3_glopos_is_valid, m);
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
p4est3_glooffs_unref (p4est3_glooffs_t ** mp)
{
  int                 waslast;
  p4est3_glooffs_t   *m;
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
      SC3E (sc3_MPI_Win_free (&m->goffsetwin));
    }
    /* deallocate data knonw on setup here */
    SC3L (&leak, sc3_mpienv_unref (&m->mpienv));
    SC3E (sc3_allocator_free (mator, m));
    SC3L (&leak, sc3_allocator_unref (&mator));
  }
  return leak;
}

sc3_error_t        *
p4est3_glotree_destroy (p4est3_glotree_t ** mp)
{
  sc3_error_t        *leak = NULL;
  p4est3_glotree_t   *m;

  SC3E_INULLP (mp, m);
  SC3L_DEMAND (&leak, sc3_refcount_is_last (&m->rc, NULL));
  SC3L (&leak, sc3_mpienv_unref (&m));

  SC3A_CHECK (m == NULL || leak != NULL);
  return leak;
}

sc3_error_t        *
p4est3_glopos_destroy (p4est3_glopos_t ** mp)
{
  sc3_error_t        *leak = NULL;
  p4est3_glopos_t    *m;

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
  p4est3_glooffs_t   *m;

  SC3E_INULLP (mp, m);
  SC3L_DEMAND (&leak, sc3_refcount_is_last (&m->rc, NULL));
  SC3L (&leak, sc3_mpienv_unref (&m));

  SC3A_CHECK (m == NULL || leak != NULL);
  return leak;
}
