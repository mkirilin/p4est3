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

#ifndef P4_TO_P8
#include <p4est_algorithms.h>
#include <p4est_dune.h>
#include <p4est_lnodes.h>
#else
#include <p8est_algorithms.h>
#include <p8est_dune.h>
#include <p8est_lnodes.h>
#endif

/* map element corner to element node index */
static const int    corner_enode[P4EST_CHILDREN] = {
  0, 4, 20, 24
#ifdef P4_TO_P8
    , 100, 104, 120, 124
#endif
};

#ifdef P4_TO_P8

/* map element edge to element node index */
static const int    edge_enode[P8EST_EDGES] = {
  2, 22, 102, 122, 10, 14, 110, 114, 50, 54, 70, 74,
};

#endif

/* map element face to element node index */
static const int    face_enode[P4EST_FACES] = {
#ifndef P4_TO_P8
  10, 14, 2, 22
#else
  60, 64, 52, 72, 12, 112
#endif
};

/* element node next to a corner c one step along a face of given normal i */
static int
corner_face_enode (int c, int i)
{
  int                 j;
  int                 m;
  int                 eno;

  P4EST_ASSERT (0 <= c && c < P4EST_CHILDREN);
  P4EST_ASSERT (0 <= i && i < P4EST_DIM);

  /* move by one step in each tangential direction */
  eno = corner_enode[c];
  m = 1;
  for (j = 0; j < P4EST_DIM; ++j) {
    if (j != i) {
      eno += (1 - 2 * (c & 1)) * m;
    }
    m *= 5;
    c >>= 1;
  }
  P4EST_ASSERT (c == 0);

  /* return final element node number */
  P4EST_ASSERT (0 <= eno && eno < m);
  return eno;
}

#ifdef P4_TO_P8

/* element node next to a corner c along j and two steps along k */
static int
corner_twostep_enode (int c, int j, int k)
{
  int                 i;
  int                 m;
  int                 eno;

  P4EST_ASSERT (0 <= c && c < P4EST_CHILDREN);
  P4EST_ASSERT (0 <= j && j < P4EST_DIM);
  P4EST_ASSERT ((0 <= k && k < P4EST_DIM) || k == -1);
  P4EST_ASSERT (j != k);

  /* move by one or two steps depending on direction */
  eno = corner_enode[c];
  m = 1;
  for (i = 0; i < P4EST_DIM; ++i) {
    if (i == j) {
      eno += (1 - 2 * (c & 1)) * m;
    }
    else if (i == k) {
      eno += (1 - 2 * (c & 1)) * 2 * m;
    }
    m *= 5;
    c >>= 1;
  }
  P4EST_ASSERT (c == 0);

  /* return final element node number */
  P4EST_ASSERT (0 <= eno && eno < m);
  return eno;
}

/* element node next to a corner c one step along a given edge j */
static int
corner_edge_enode (int c, int j)
{
  return corner_twostep_enode (c, j, -1);
}

#endif /* P4_TO_P8 */

/* set default parameters */
void
p4est_dune_numbers_params_init (p4est_dune_numbers_params_t * params)
{
  P4EST_ASSERT (params != NULL);

  /* set default parameters to construct a dune numbers structure */
  params->ctype = P4EST_CONNECT_FULL;
}

#if 0                           /* presently unused */

static              p4est_gloidx_t
lni_to_gni (p4est_lnodes_t * ln, p4est_locidx_t lni)
{
  p4est_locidx_t      oco;

  P4EST_ASSERT (0 <= lni);

  if (lni < (oco = ln->owned_count)) {
    return ln->global_offset + lni;
  }
  else {
    P4EST_ASSERT (lni < ln->num_local_nodes);
    return ln->nonlocal_nodes[lni - oco];
  }
}

#endif

/* Populate number arrays for real.  This function does all the work. */
static void
generate_numbers (p4est_dune_numbers_t * dn, p4est_lnodes_t * ln)
{
  p4est_lnodes_code_t fc;
  p4est_locidx_t      lne, el;
  p4est_locidx_t     *enos;
  int                 i;
  int                 cid, work;
  int                 c, f;
  int                 do_c, do_f;
  p4est_locidx_t     *geco, *gefa;
#ifdef P4_TO_P8
  int                 j, k;
  int                 e;
  int                 do_e;
  p4est_locidx_t     *geed;
#endif

  /* basic verification of input parameters */
  P4EST_ASSERT (dn != NULL);
  P4EST_ASSERT (ln != NULL);
  P4EST_ASSERT (ln->degree == 4);

  /* examine which codimensions are handled */
  do_c = dn->element_corners != NULL;
  geco = NULL;
#ifdef P4_TO_P8
  do_e = dn->element_edges != NULL;
  geed = NULL;
#endif
  do_f = dn->element_faces != NULL;
  gefa = NULL;

  /* maybe there is nothing to do */
  if (!do_c &&
#ifdef P4_TO_P8
      !do_e &&
#endif
      !do_f) {
    /* there is absolutely nothing to do */
    return;
  }

  /* prepare loop over elements to assign globally unique node indices */
  lne = ln->num_local_elements;
  enos = ln->element_nodes;
  for (el = 0; el < lne; ++el) {

    /* assign standard values for an element that is not hanging anywhere */
    if (do_c) {
      geco = (p4est_locidx_t *)
        sc_array_index (dn->element_corners, el * P4EST_CHILDREN);
      for (c = 0; c < P4EST_CHILDREN; ++c) {
        /* retrieve global number of corner node */
        geco[c] = enos[corner_enode[c]];
      }
    }
#ifdef P4_TO_P8
    if (do_e) {
      geed = (p4est_locidx_t *)
        sc_array_index (dn->element_edges, el * P8EST_EDGES);
      for (e = 0; e < P8EST_EDGES; ++e) {
        /* retrieve global number of edge node */
        geed[e] = enos[edge_enode[e]];
      }
    }
#endif
    if (do_f) {
      gefa = (p4est_locidx_t *)
        sc_array_index (dn->element_faces, el * P4EST_FACES);
      for (f = 0; f < P4EST_FACES; ++f) {
        /* retrieve global number of face node */
        gefa[f] = enos[face_enode[f]];
      }
    }

    /* modify node numbers for hanging faces and edges */
    if (0 != (fc = ln->face_code[el])) {
      /* element has child number cid and some hanging faces or edges */
      cid = (fc & (P4EST_CHILDREN - 1));
      work = (fc >> P4EST_DIM);

      /* iterate over the faces adjacent to corner cid */
      for (i = 0; i < P4EST_DIM; ++i) {
        if (work & 1) {
          /* face normal to direction i is hanging */
          f = p4est_corner_faces[cid][i];
          c = cid ^ (P4EST_CHILDREN - 1) ^ (1 << i);
          if (do_c) {
            /* assign to mid-face corner */
            geco[c] = enos[face_enode[f]];
          }
#ifdef P4_TO_P8
          if (do_e) {
            /* assign to the two mid-face edges */
            j = (i + 1) % 3;
            k = (j + 1) % 3;
            geed[p8est_corner_edges[c][j]] =
              enos[corner_twostep_enode (cid, j, k)];
            geed[p8est_corner_edges[c][k]] =
              enos[corner_twostep_enode (cid, k, j)];
          }
#endif
          if (do_f) {
            /* assign to small mid-face */
            gefa[f] = enos[corner_face_enode (cid, i)];
          }
        }
        work >>= 1;
      }
#ifdef P4_TO_P8
      /* iterate over the edges adjacent to corner cid */
      for (j = 0; j < P4EST_DIM; ++j) {
        if (work & 1) {
          /* edge in direction j is hanging */
          e = p8est_corner_edges[cid][j];
          if (do_c) {
            /* assign to mid-edge corner */
            c = cid ^ (1 << j);
            geco[c] = enos[edge_enode[e]];
          }
          if (do_e) {
            /* assign to small mid-edge */
            geed[e] = enos[corner_edge_enode (cid, j)];
          }
        }
        work >>= 1;
      }
#endif
      P4EST_ASSERT (work == 0);
    }

    /* proceed to the next element */
    enos += ln->vnodes;
  }
}

static void
consecutive_numbers (p4est_dune_numbers_t * dn,
                     p4est_locidx_t num_local_nodes)
{
  /* an array of bits, one per local node number, for codimensions */
  uint8_t            *local_node_bytes;
  uint8_t             mask_byte, node_byte;
  int                 j;
  sc_array_t         *dim_numbers[P4EST_DIM], *numbers;
  p4est_locidx_t     *dim_results[P4EST_DIM], *newnums;
  p4est_locidx_t      li, lcount, lnum, ncount;

  /* input checks */
  P4EST_ASSERT (dn != NULL);

  /* prepare temporary data */
  local_node_bytes = P4EST_ALLOC_ZERO (uint8_t, num_local_nodes);
  dim_numbers[0] = dn->element_corners;
  dim_results[0] = &dn->num_corner_numbers;
#ifdef P4_TO_P8
  dim_numbers[1] = dn->element_edges;
  dim_results[1] = &dn->num_edge_numbers;
#endif
  dim_numbers[P4EST_DIM - 1] = dn->element_faces;
  dim_results[P4EST_DIM - 1] = &dn->num_face_numbers;

  /* go through node arrays to count unique entries */
  for (j = 0; j < P4EST_DIM; ++j) {
    P4EST_ASSERT (*dim_results[j] == 0);
    if ((numbers = dim_numbers[j]) == NULL) {
      continue;
    }
    P4EST_ASSERT (numbers->elem_size == sizeof (p4est_locidx_t));
    mask_byte = 1 << j;

    /* mark each local node that occurs in any of the arrays */
    lcount = (p4est_locidx_t) numbers->elem_count;
    for (li = 0; li < lcount; ++li) {
      lnum = *(p4est_locidx_t *) sc_array_index (numbers, li);
      P4EST_ASSERT (0 <= lnum && lnum < num_local_nodes);
      node_byte = local_node_bytes[lnum];
      if (!(node_byte & mask_byte)) {
        local_node_bytes[lnum] = node_byte | mask_byte;
      }
    }
  }

  /* renumber nodes to contiguous subrange */
  newnums = P4EST_ALLOC (p4est_locidx_t, num_local_nodes);
  for (j = 0; j < P4EST_DIM; ++j) {
    P4EST_ASSERT (*dim_results[j] == 0);
    if ((numbers = dim_numbers[j]) == NULL) {
      continue;
    }
    mask_byte = 1 << j;
    memset (newnums, -1, num_local_nodes * sizeof (*newnums));
    ncount = 0;

    /* figure out which local nodes are relevant for dimension j */
    lcount = num_local_nodes;
    for (li = 0; li < lcount; ++li) {
      if (local_node_bytes[li] & mask_byte) {
        newnums[li] = ncount++;
      }
    }

    /* reassign entries of numbers array */
    lcount = (p4est_locidx_t) numbers->elem_count;
    for (li = 0; li < lcount; ++li) {
      lnum = *(p4est_locidx_t *) sc_array_index (numbers, li);
      P4EST_ASSERT (0 <= lnum && lnum < num_local_nodes);
      P4EST_ASSERT (local_node_bytes[lnum] & mask_byte);
      lnum = newnums[lnum];
      P4EST_ASSERT (0 <= lnum && lnum < ncount);
      *(p4est_locidx_t *) sc_array_index (numbers, li) = lnum;
    }
    *dim_results[j] = ncount;
  }
  P4EST_FREE (newnums);
  P4EST_FREE (local_node_bytes);
}

/* create lookup tables for unique corners and faces */
p4est_dune_numbers_t *
p4est_dune_numbers_new (p4est_t * p4est, p4est_ghost_t * ghost,
                        const p4est_dune_numbers_params_t * params)
{
  p4est_dune_numbers_t *dn;
  p4est_dune_numbers_params_t *pa;
  p4est_lnodes_t     *ln;
  p4est_locidx_t      lne, lnums;

  /* formally verify arguments */
  P4EST_ASSERT (p4est != NULL);
  P4EST_ASSERT (p4est_is_valid (p4est));
  P4EST_ASSERT (p4est_is_balanced (p4est, P4EST_CONNECT_ALMOST));
  P4EST_ASSERT (ghost != NULL);

  /* allocate numbers structure to populate */
  dn = P4EST_ALLOC_ZERO (p4est_dune_numbers_t, 1);
  pa = &dn->params;

  /* establish input parameters and remember them */
  if (params == NULL) {
    p4est_dune_numbers_params_init (pa);
  }
  else {
    *pa = *params;
    params = NULL;
  }

  /* generate temporary local ghost (if needed) and node numbering */
  ln = p4est_lnodes_new (p4est, ghost, 4);

  /* construct corner and face numbers and return them */
  lne = p4est->local_num_quadrants;
  if (pa->ctype >= P4EST_CONNECT_CORNER) {
    dn->element_corners =
      sc_array_new_count (sizeof (p4est_locidx_t), lne * P4EST_CHILDREN);
  }
#ifdef P4_TO_P8
  if (pa->ctype >= P8EST_CONNECT_EDGE) {
    dn->element_edges =
      sc_array_new_count (sizeof (p4est_locidx_t), lne * P8EST_EDGES);
  }
#endif
  if (pa->ctype >= P4EST_CONNECT_FACE) {
    dn->element_faces =
      sc_array_new_count (sizeof (p4est_locidx_t), lne * P4EST_FACES);
  }
  generate_numbers (dn, ln);

  /* free temporary data */
  lnums = ln->num_local_nodes;
  p4est_lnodes_destroy (ln);

  /* transform numbering into contiguous ranges per codimension */
  consecutive_numbers (dn, lnums);

  /* that's it! */
  return dn;
}

void
p4est_dune_numbers_destroy (p4est_dune_numbers_t * dn)
{
  P4EST_ASSERT (dn != NULL);

  if (dn->element_corners != NULL) {
    sc_array_destroy (dn->element_corners);
  }
#ifdef P4_TO_P8
  if (dn->element_edges != NULL) {
    sc_array_destroy (dn->element_edges);
  }
#endif
  if (dn->element_faces != NULL) {
    sc_array_destroy (dn->element_faces);
  }
  P4EST_FREE (dn);
}

typedef struct p4est_dune_iter
{
  /* wrapped iteration context */
  p4est_t            *p4est;
  p4est_iter_volume_t iter_volume;
  p4est_iter_face_t   iter_face;
  void               *user_data;

  /* this structure will be reused many times */
  p4est_iter_face_info_t sfinfo;
  p4est_iter_face_info_t *finfo;
}
p4est_dune_iter_t;

static void
dune_iter_volume (p4est_iter_volume_info_t * info, void *user_data)
{
  p4est_dune_iter_t  *dune_iter = (p4est_dune_iter_t *) user_data;

  /* consistency checks */
  P4EST_ASSERT (info != NULL);
  P4EST_ASSERT (dune_iter != NULL);
  P4EST_ASSERT (dune_iter->p4est == info->p4est);

  /* unwrap volume iteration callback */
  P4EST_ASSERT (dune_iter->iter_volume != NULL);
  dune_iter->iter_volume (info, dune_iter->user_data);
}

static void
dune_iter_face (p4est_iter_face_info_t * info, void *user_data)
{
  p4est_dune_iter_t  *dune_iter = (p4est_dune_iter_t *) user_data;
  p4est_iter_face_info_t *finfo;
  p4est_iter_face_side_t *s[2], *sh;
  p4est_iter_face_side_t *f[2], *fh;
  int                 i;
  int                 fside, hside;

  /* consistency checks */
  P4EST_ASSERT (info != NULL);
  P4EST_ASSERT (dune_iter != NULL);
  P4EST_ASSERT (dune_iter->p4est == info->p4est);

  /* unwrap face callback */
  P4EST_ASSERT (dune_iter->iter_face != NULL);

  /* examine face boundary situation */
  P4EST_ASSERT (info->sides.elem_count > 0);
  if (info->sides.elem_count == 1) {
    /* this face is on a physical boundary */
    dune_iter->iter_face (info, dune_iter->user_data);
    return;
  }

  /* examine hanging face situation */
  P4EST_ASSERT (info->sides.elem_count == 2);
  for (i = 0; i < 2; ++i) {
    s[i] = (p4est_iter_face_side_t *) sc_array_index_int (&info->sides, i);
    P4EST_ASSERT (0 <= s[i]->face && s[i]->face < P4EST_FACES);
  }
  if (!s[0]->is_hanging && !s[1]->is_hanging) {
    dune_iter->iter_face (info, dune_iter->user_data);
    return;
  }

  /* we have precisely one hanging face.  Reconstruct its info */
  P4EST_ASSERT (s[0]->is_hanging != s[1]->is_hanging);
  finfo = dune_iter->finfo;
  P4EST_ASSERT (finfo != NULL);
  P4EST_ASSERT (finfo->sides.elem_size == sizeof (p4est_iter_face_side_t));
  P4EST_ASSERT (finfo->sides.elem_count == 2);
  if (s[0]->is_hanging) {
    hside = 0;
    fside = 1;
  }
  else {
    hside = 1;
    fside = 0;
  }
  P4EST_ASSERT (s[hside]->is_hanging);
  P4EST_ASSERT (!s[fside]->is_hanging);
  for (i = 0; i < 2; ++i) {
    f[i] = (p4est_iter_face_side_t *) sc_array_index_int (&finfo->sides, i);
  }

  /* copy full side information as is and prepare setting the other side */
  *f[fside] = *s[fside];
  sh = s[hside];
  fh = f[hside];
  fh->treeid = sh->treeid;
  fh->face = sh->face;
  fh->is_hanging = sh->is_hanging;
  P4EST_ASSERT (fh->is_hanging);

  /* iterate over the hanging face neighbors */
  for (i = 0; i < P4EST_HALF; ++i) {
    if (s[fside]->is.full.is_ghost && sh->is.hanging.is_ghost[i]) {
      /* if both sides are remote, we do nothing for this pair */
      continue;
    }

    /* we copy the hanging face into the first position */
    fh->is.hanging.is_ghost[0] = sh->is.hanging.is_ghost[i];
    fh->is.hanging.quad[0] = sh->is.hanging.quad[i];
    fh->is.hanging.quadid[0] = sh->is.hanging.quadid[i];

    /* abuse the second position to store this face's sequence number */
    fh->is.hanging.quadid[1] = i;

    /* execute callback for this face pairing */
    dune_iter->iter_face (finfo, dune_iter->user_data);
  }
}

void
p4est_dune_iterate (p4est_t * p4est, p4est_ghost_t * ghost_layer,
                    void *user_data, p4est_iter_volume_t iter_volume,
                    p4est_iter_face_t iter_face)
{
  P4EST_ASSERT (p4est != NULL);
  P4EST_ASSERT (ghost_layer != NULL);

  if (iter_face == NULL) {
    /* no need to do anything special: pass through */
    p4est_iterate (p4est, ghost_layer, user_data, iter_volume, NULL,
#ifdef P4_TO_P8
                   NULL,
#endif
                   NULL);
  }
  else {
    /* we need to translate the face iteration calls */
    p4est_iter_face_info_t *finfo;
    p4est_dune_iter_t   sdune_iter, *dune_iter = &sdune_iter;
    memset (dune_iter, 0, sizeof (*dune_iter));

    /* remember wrapped information */
    dune_iter->p4est = p4est;
    dune_iter->iter_volume = iter_volume;
    dune_iter->iter_face = iter_face;
    dune_iter->user_data = user_data;

    /* prepare internal context */
    finfo = dune_iter->finfo = &dune_iter->sfinfo;
    finfo->p4est = p4est;
    finfo->ghost_layer = ghost_layer;
    sc_array_init_count (&finfo->sides, sizeof (p4est_iter_face_side_t), 2);
    sc_array_memset (&finfo->sides, 0);

    /* wrapped iterator call */
    p4est_iterate (p4est, ghost_layer, dune_iter,
                   iter_volume != NULL ? dune_iter_volume : NULL,
                   /* we know this is not NULL here, regardless */
                   iter_face != NULL ? dune_iter_face : NULL,
#ifdef P4_TO_P8
                   NULL,
#endif
                   NULL);

    /* free work memory */
    sc_array_reset (&finfo->sides);
  }
}
