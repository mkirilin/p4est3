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

#include <p4est3_internal.h>
#include <p8est_extended.h>
#include <p8est_search.h>
#include <p4est3_convert_p8est.h>
#include <p4est3_p8est.h>
#include <p4est3_quadrant_zyx.h>
#include <p4est3_quadrant_mort3d.h>
#include <p8est_vtk.h>

#include <sc_statistics.h>
#include <sc_flops.h>
#include <sc_options.h>
#include "model.h"
#include "tribox.h"
#include <string.h>

static int          max_ref_level = 1;

typedef struct triangle
{
  double              v0[3];
  double              v1[3];
  double              v2[3];
}
triangle_t;

/*static int
triangulation_is_vertex_inside_aabb (const double * aabb, const double * v)
{
  return ((v[0] >= aabb[0] && v[0] <= aabb[3]) &&
          (v[1] >= aabb[1] && v[1] <= aabb[4]) &&
          (v[2] >= aabb[2] && v[2] <= aabb[5]));
}*/

static int
triangulation_intersect_model (p4est_topidx_t which_tree,
                               const double aabb[6], void *model, void *point)
{
  float               boxcenter[3], boxhalfsize[3], triverts[3][3];

  triangle_t         *t =
    (triangle_t *) ((p4est_model_t *) model)->primitives;
  size_t              p = *(size_t *) point;

  /*int is_v0_in = triangulation_is_vertex_inside_aabb (aabb, t[p].v0);
     int is_v1_in = triangulation_is_vertex_inside_aabb (aabb, t[p].v1);
     int is_v2_in = triangulation_is_vertex_inside_aabb (aabb, t[p].v2);

     if (is_v0_in || is_v1_in || is_v2_in)
     {
     return 1;
     }
     else {
     return 0;
     } */
  boxcenter[0] = (aabb[3] + aabb[0]) * 0.5;
  boxcenter[1] = (aabb[4] + aabb[1]) * 0.5;
  boxcenter[2] = (aabb[5] + aabb[2]) * 0.5;

  boxhalfsize[0] = (aabb[3] - aabb[0]) * 0.5;
  boxhalfsize[1] = (aabb[4] - aabb[1]) * 0.5;
  boxhalfsize[2] = (aabb[5] - aabb[2]) * 0.5;

  triverts[0][0] = t[p].v0[0];
  triverts[0][1] = t[p].v0[1];
  triverts[0][2] = t[p].v0[2];

  triverts[1][0] = t[p].v1[0];
  triverts[1][1] = t[p].v1[1];
  triverts[1][2] = t[p].v1[2];

  triverts[2][0] = t[p].v2[0];
  triverts[2][1] = t[p].v2[1];
  triverts[2][2] = t[p].v2[2];

  return triBoxOverlap (boxcenter, boxhalfsize, triverts);
}

/*
 * Read a line from a file. Obtained from:
 * http://stackoverflow.com/questions/314401/
 * how-to-read-a-line-from-the-console-in-c/314422#314422
 *
 * Using this avoids a dependence on IEEE Std 1003.1-2008 (``POSIX.1'') for the
 * getline function.
 * 
 * Copied from p4est_connectivity.c:
 * p4est_connectivity_getline_upper (FILE * stream).
 */
static char        *
triangulation_getline_upper (FILE *stream)
{
  char               *line = P4EST_ALLOC (char, 1024), *linep = line;
  size_t              lenmax = 1024, len = lenmax;
  int                 c;

  if (line == NULL)
    return NULL;

  for (;;) {
    c = fgetc (stream);
    if (c == EOF && linep == line) {
      P4EST_FREE (linep);
      return NULL;
    }
    c = toupper (c);

    if (--len == 0) {
      char               *linen;

      len = lenmax;
      lenmax *= 2;

      linen = P4EST_REALLOC (linep, char, lenmax);
      if (linen == NULL) {
        P4EST_FREE (linep);
        return NULL;
      }

      line = linen + (line - linep);
      linep = linen;
    }
    if ((*line++ = c) == '\n')
      break;
  }
  *line = '\0';
  return linep;
}

/** This function goes through the file fin that has .off (Object File Format).
 * Object File Format (.off) files are used to represent the geometry of a
 * model by specifying the polygons of the model's surface. The polygons can
 * have any number of vertices.
 * The .off files in the conform to the following standard. OFF files are all
 * ASCII files beginning with the keyword OFF. The next line states the number
 * of vertices, the number of faces, and the number of edges. The number of
 * edges can be safely ignored.
 *
 * The vertices are listed with x, y, z coordinates, written one per line.
 * After the list of vertices, the faces are listed, with one face per line.
 *
 * ****************************************************************************
 * !!!Our model considers triangulated meshes only, so its all faces should be
 * triangles!!!
 * ****************************************************************************
 *
 * For each face, the number of vertices is specified, followed by indices into
 * the list of vertices. See the examples below.
 * Note that earlier versions of the model files had faces with -1 indices into
 * the vertex list. That was due to an error in the conversion program and
 * should be corrected now.
 * OFF
 * numVertices numFaces numEdges
 * x y z
 * x y z
 *
 * ... numVertices like above
 * NVertices v1 v2 v3 ... vN
 * MVertices v1 v2 v3 ... vM
 * ... numFaces like above
 *
 * Note that vertices are numbered starting at 0 (not starting at 1).
 * A simple example for a cube (that is not supported by this model):
 *
 * OFF
 * 8 6 0
 * -0.500000 -0.500000 0.500000
 * 0.500000 -0.500000 0.500000
 * -0.500000 0.500000 0.500000
 * 0.500000 0.500000 0.500000
 * -0.500000 0.500000 -0.500000
 * 0.500000 0.500000 -0.500000
 * -0.500000 -0.500000 -0.500000
 * 0.500000 -0.500000 -0.500000
 * 4 0 1 3 2
 * 4 2 3 5 4
 * 4 4 5 7 6
 * 4 6 7 1 0
 * 4 1 7 5 3
 * 4 6 0 2 4
*/
static int
triangulation_read_off_file_stream (p4est_model_t *m, FILE *fin)
{
  char               *line;
  int                 lines_read = 0;
  size_t              v = 0, f = 0;
  int                 retval;
  size_t              num_vertices, num_edges, v2f, vid[3];
  triangle_t         *faces = NULL;
  double             *vertices = NULL;
  double              x = 0, y = 0, z = 0;
  double              min_x = 0., min_y = 0., min_z = 0.;
  double              max_x = 0., max_y = 0., max_z = 0.;
  double              axis_scale = 1.;

  for (;;) {
    line = triangulation_getline_upper (fin);

    if (line == NULL) {
      break;
    }

    ++lines_read;

    /* check for file format */
    if (lines_read == 1) {
      if (!strstr (line, "OFF")) {
        P4EST_LERROR ("Wrong file format to read");
        P4EST_FREE (line);
        return 0;
      }
      P4EST_FREE (line);
      continue;
    }

    /* check for number of vertices and faces in the object to read */
    if (lines_read == 2) {
      retval
        =
        sscanf (line, "%lu %lu %lu", &num_vertices, &m->num_prim, &num_edges);

      if (retval != 3) {
        P4EST_LERROR ("Wrong file format to read");
        P4EST_FREE (line);
        return 0;
      }
      if (num_vertices < 3) {
        P4EST_LERROR ("Not enough vertices to build a triangle");
        P4EST_FREE (line);
        return 0;
      }
      if (m->num_prim == 0) {
        P4EST_LERROR ("No primitives in the file");
        P4EST_FREE (line);
        return 0;
      }

      /* allocate memory for vertices and face (triangles (primitives)) */
      vertices = P4EST_ALLOC (double, 3 * num_vertices);
      faces = P4EST_ALLOC (triangle_t, m->num_prim);
      m->primitives = (void *) faces;

      P4EST_FREE (line);
      continue;
    }

    if (v < num_vertices) {
      /* read vertices */
      retval = sscanf (line, "%lf %lf %lf", &x, &y, &z);
      if (retval != 3) {
        P4EST_LERROR ("Premature end of file");
        P4EST_FREE (line);
        return 0;
      }
      if (v == 0) {
        min_x = max_x = x;
        min_y = max_y = y;
        min_z = max_z = z;
      }
      min_x = SC_MIN (min_x, x);
      min_y = SC_MIN (min_y, y);
      min_z = SC_MIN (min_z, z);

      max_x = SC_MAX (max_x, x);
      max_y = SC_MAX (max_y, y);
      max_z = SC_MAX (max_z, z);

      vertices[3 * v + 0] = x;
      vertices[3 * v + 1] = y;
      vertices[3 * v + 2] = z;

      ++v;

      if (v == num_vertices) {
        axis_scale
          =
          1. / SC_MAX (SC_MAX (max_x - min_x, max_y - min_y), max_z - min_z);
      }
      P4EST_FREE (line);
      continue;
    }
    else {
      /* all vertices are read, now read faces */
      retval
        = sscanf (line, "%lu %lu %lu %lu", &v2f, &vid[0], &vid[1], &vid[2]);
      if (retval != 4) {
        P4EST_LERROR ("Premature end of file");
        P4EST_FREE (line);
        return 0;
      }
      if (v2f != 3) {
        P4EST_LERROR ("Unsupported face format in file");
        P4EST_FREE (line);
        return 0;
      }
      faces[f].v0[0] = (vertices[3 * vid[0] + 0] - min_x) * axis_scale;
      faces[f].v0[1] = (vertices[3 * vid[0] + 1] - min_y) * axis_scale;
      faces[f].v0[2] = (vertices[3 * vid[0] + 2] - min_z) * axis_scale;

      faces[f].v1[0] = (vertices[3 * vid[1] + 0] - min_x) * axis_scale;
      faces[f].v1[1] = (vertices[3 * vid[1] + 1] - min_y) * axis_scale;
      faces[f].v1[2] = (vertices[3 * vid[1] + 2] - min_z) * axis_scale;

      faces[f].v2[0] = (vertices[3 * vid[2] + 0] - min_x) * axis_scale;
      faces[f].v2[1] = (vertices[3 * vid[2] + 1] - min_y) * axis_scale;
      faces[f].v2[2] = (vertices[3 * vid[2] + 2] - min_z) * axis_scale;

      ++f;
    }
    P4EST_FREE (line);
  }
  if (lines_read >= 2) {
    P4EST_FREE (vertices);
  }
  return 1;
}

static int
triangulation_read_off_file (p4est_model_t *m, const char *filename)
{
  int                 retval = 0;
  FILE               *fin = NULL;

  P4EST_GLOBAL_PRODUCTIONF ("Reading connectivity from %s\n", filename);

  fin = fopen (filename, "rb");
  if (fin == NULL) {
    P4EST_LERRORF ("Failed to open %s\n", filename);
    return 0;
  }

  if (!triangulation_read_off_file_stream (m, fin)) {
    P4EST_LERRORF ("Failed to read %s: pass 1\n", filename);
    return 0;
  }

  retval = fclose (fin);
  if (retval != 0) {
    P4EST_LERRORF ("Failed to close %s\n", filename);
    return 0;
  }

  return 1;
}

static void
triangulation_desroy_primitives (void *primitives)
{
  P4EST_FREE (primitives);
}

static int
triangulation_setup_model (p4est_model_t **m, const char *filename,
                           const char *model_name)
{
  p4est_model_t      *model = P4EST_ALLOC_ZERO (p4est_model_t, 1);
  model->conn = p8est_connectivity_new_unitcube ();
  if (model->conn == NULL) {
    P4EST_LERROR ("Failed to create a model's connectivity");
    return 0;
  }
  model->geom = NULL;
  model->intersect = triangulation_intersect_model;
  if (!triangulation_read_off_file (model, filename)) {
    P4EST_LERRORF ("Failed to read a valid model from %s\n", filename);
    return 0;
  }
  model->output_prefix = model_name;
  model->destroy_primitives = triangulation_desroy_primitives;
  *m = model;
  return 1;
}

static sc3_error_t *
make_allocator (sc3_allocator_t *oa, sc3_allocator_t **alloc)
{
  SC3A_IS (sc3_allocator_is_setup, oa);
  SC3E (sc3_allocator_new (oa, alloc));
  SC3E (sc3_allocator_setup (*alloc));
  return NULL;
}

#ifndef P4EST_ENABLE_DEBUG
static sc3_error_t *
array_new (sc3_allocator_t *alloc, size_t esize, int ealloc,
           int ecount, sc3_array_t **arr)
{
  SC3E_RETVAL (arr, NULL);
  SC3A_IS (sc3_allocator_is_setup, alloc);
  SC3A_CHECK (ealloc >= 0);

  SC3E (sc3_array_new (alloc, arr));
  SC3E (sc3_array_set_elem_size (*arr, esize));
  SC3E (sc3_array_set_elem_alloc (*arr, ealloc));
  SC3E (sc3_array_set_elem_count (*arr, ecount));
  SC3E (sc3_array_set_initzero (*arr, 1));
  SC3E (sc3_array_setup (*arr));

  return NULL;
}

static sc3_error_t *
compare_results (sc3_allocator_t *alloc, p4est3_t *p3, p4est_t *p,
                 const p4est3_quadrant_vtable_t *qvt)
{

  char               *q3;
  int                *level, *p3level;
  p4est_quadrant_t   *q;
  p4est3_topidx       i, tt;
  size_t              nq;
  p4est3_topidx       fltree, lltree;
  p4est3_gloidx       num_glo_quads;
  p4est3_locidx       num_loc_quads;
  size_t              processed_quads, processed_quads_p3;
  p4est_tree_t       *tree;
  sc3_array_t        *p3levels, *levels;

  SC3E (p4est3_get_local_num_trees (p3, &fltree, &lltree));
  SC3E (p4est3_get_global_num_quads (p3, &num_glo_quads));
  SC3E (p4est3_get_local_num_quads (p3, &num_loc_quads));
  SC3E_DEMAND (fltree == p->first_local_tree && lltree == p->last_local_tree,
               "Different trees at processor");
  SC3E_DEMAND (num_glo_quads == p->global_num_quadrants,
               "different #global quadrants");
  SC3E_DEMAND (num_loc_quads == p->local_num_quadrants,
               "different #local quadrants");
  SC3E (array_new (alloc, sizeof (int), num_loc_quads, 0, &p3levels));
  SC3E (array_new (alloc, sizeof (int), p->local_num_quadrants, 0, &levels));
  for (tt = p->first_local_tree; tt <= p->last_local_tree; ++tt) {
    tree = p4est_tree_array_index (p->trees, tt);
    for (nq = 0; nq < tree->quadrants.elem_count; ++nq) {
      q = (p4est_quadrant_t *) sc_array_index (&tree->quadrants, nq);
      SC3E (sc3_array_push (levels, &level));
      *level = q->level;
    }
  }

  SC3E (p4est3_get_quadrants (p3, &q3));
  if (num_loc_quads > 0) {
    SC3E (sc3_array_push (p3levels, &level));
    SC3E (p4est3_quadrant_level (qvt, q3, level));
  }
  for (i = 1; i < num_loc_quads; ++i) {
    q3 += qvt->quadrant_size;
    SC3E (sc3_array_push (p3levels, &level));
    SC3E (p4est3_quadrant_level (qvt, q3, level));
  }
  SC3E (sc3_array_get_elem_count (p3levels, &processed_quads_p3));
  SC3E (sc3_array_get_elem_count (levels, &processed_quads));
  SC3E_DEMAND (processed_quads_p3 == processed_quads, "wrong #p3levels");
  SC3E_DEMAND (processed_quads == (size_t) num_loc_quads, "wrong #levels");
  for (i = 0; i < num_loc_quads; ++i) {
    SC3E (sc3_array_index (p3levels, i, &p3level));
    SC3E (sc3_array_index (levels, i, &level));
    SC3E_DEMAND (*level == *p3level, "levels mismatch");
  }

  SC3E (sc3_array_destroy (&p3levels));
  SC3E (sc3_array_destroy (&levels));
  return NULL;
}
#endif /* P4EST_ENABLE_DEBUG */

static sc3_error_t *
p4est3_new_shortcut (p4est3_t **p3, sc3_allocator_t *alloc,
                     sc3_MPI_Comm_t mpicomm, p4est3_connectivity_t *conn,
                     const p4est3_quadrant_vtable_t *qvt,
                     p4est3_t *src, int is_partition, int is_family,
                     p4est3_weight_callback_t cweight, void *user_data)
{
  SC3E (p4est3_new (alloc, p3));
  SC3E (p4est3_set_comm (*p3, mpicomm, 1));
  SC3E (p4est3_set_connectivity (*p3, conn));
  SC3E (p4est3_set_quadrant_vtable (*p3, qvt));
  SC3E (p4est3_set_setup_mode (*p3, P4EST3_NEW_RECURSIVE_CHILD));
  SC3E (p4est3_set_source (*p3, src));
  SC3E (p4est3_set_shared (*p3, 1));
  SC3E (p4est3_set_contiguous (*p3, 1));
  SC3E (p4est3_set_family (*p3, is_family));
  SC3E (p4est3_set_partition (*p3, is_partition, cweight));
  /*SC3E (p4est3_set_user_data (*p3, user_data)); */
  if ((*p3)->old != NULL) {
    (*p3)->old->user_data = user_data;
  }

  return NULL;
}

static sc3_error_t *
run_program (sc_MPI_Comm *mpicomm, p4est_model_t *model)
{
  size_t              zz;
  char                filename[BUFSIZ];
  sc_array_t         *primitives;
  p4est_t            *p4est;
  p4est3_t           *p4est3, *p3part;
  sc3_allocator_t    *alloc, *mainalloc;
  const size_t        quad_data_size = 0;
  const int           start_level = 3;
  int                 level;
  p4est3_connectivity_t *conn3;

  sc_flopinfo_t       fi, sshot_total, sshot_new, sshot_convert,
    sshot_part_p4est, sshot_part_p4est3;
  sc_statinfo_t       stats;

  mainalloc = sc3_allocator_nothread ();
  SC3E (make_allocator (mainalloc, &alloc));
  SC3E (p4est3_new (alloc, &p4est3));
  SC3E (p4est3_set_shared (p4est3, 1));
  SC3E (p4est3_set_contiguous (p4est3, 1));

  /* run mesh refinement based on data */
  P4EST_GLOBAL_PRODUCTIONF ("Setting up %lld search objects\n",
                            (long long) model->num_prim);
  primitives = sc_array_new_count (sizeof (triangle_t), model->num_prim);
  for (zz = 0; zz < model->num_prim; ++zz) {
    *(size_t *) sc_array_index (primitives, zz) = zz;
  }

  /* create mesh */
  P4EST_GLOBAL_PRODUCTION ("Create initial mesh\n");
  sc_flops_snap (&fi, &sshot_total);
  sc_flops_snap (&fi, &sshot_new);
  p4est = p4est_new_ext (*mpicomm, model->conn, 0, start_level, 1,
                         quad_data_size, p4est_model_quad_init, model);
  sc_flops_shot (&fi, &sshot_new);
  for (level = start_level; level < max_ref_level; ++level) {

    p4est_search_local (p4est, 0, NULL, p4est_model_intersect, primitives);
    p4est_refine (p4est, 0, p4est_model_refine, p4est_model_quad_init);

    if (level == max_ref_level - 1) {

      snprintf (filename, BUFSIZ, "./%s/%d/before",
                model->output_prefix, level + 1);
      p4est_vtk_write_file (p4est, model->geom, filename);

      sc_flops_snap (&fi, &sshot_convert);
      SC3E (p4est3_convert_p8est (p4est, p4est3, &conn3));
      sc_flops_shot (&fi, &sshot_convert);

      sc_flops_snap (&fi, &sshot_part_p4est);
      p4est_partition (p4est, 0, NULL);
      sc_flops_shot (&fi, &sshot_part_p4est);

      sc_flops_snap (&fi, &sshot_part_p4est3);
      SC3E (p4est3_new_shortcut (&p3part, alloc, *mpicomm, p4est3->conn,
                                 p4est3->qvt, p4est3, 1, 0, NULL, NULL));
      SC3E (p4est3_setup (p3part));
      sc_flops_shot (&fi, &sshot_part_p4est3);

      snprintf (filename, BUFSIZ, "./%s/%d/after",
                model->output_prefix, level + 1);
      p4est_vtk_write_file (p4est, model->geom, filename);
    }
  }
  sc_flops_shot (&fi, &sshot_total);
  sc_stats_set1 (&stats, sshot_total.iwtime, "Total");
  sc_stats_compute (*mpicomm, 1, &stats);
  sc_stats_print (p4est_package_id, SC_LP_ESSENTIAL, 1, &stats, 1, 1);

  sc_stats_set1 (&stats, sshot_new.iwtime, "p4est_new");
  sc_stats_compute (*mpicomm, 1, &stats);
  sc_stats_print (p4est_package_id, SC_LP_ESSENTIAL, 1, &stats, 1, 1);

  sc_stats_set1 (&stats, sshot_convert.iwtime, "Conversion p4est -> p4est3");
  sc_stats_compute (*mpicomm, 1, &stats);
  sc_stats_print (p4est_package_id, SC_LP_ESSENTIAL, 1, &stats, 1, 1);

  sc_stats_set1 (&stats, sshot_part_p4est.iwtime, "p4est_partition");
  sc_stats_compute (*mpicomm, 1, &stats);
  sc_stats_print (p4est_package_id, SC_LP_ESSENTIAL, 1, &stats, 1, 1);

  sc_stats_set1 (&stats, sshot_part_p4est3.iwtime, "p4est3_partition");
  sc_stats_compute (*mpicomm, 1, &stats);
  sc_stats_print (p4est_package_id, SC_LP_ESSENTIAL, 1, &stats, 1, 1);

  /* cleanup */
  sc_array_destroy (primitives);
  p4est_destroy (p4est);
  SC3E (p4est3_destroy (&p4est3));
  SC3E (p4est3_destroy (&p3part));
  SC3E (p4est3_connectivity_destroy (&conn3));
  SC3E (sc3_allocator_destroy (&alloc));
  return NULL;
}

static int
usagerrf (sc_options_t *opt, const char *fmt, ...)
{
  va_list             ap;
  char                msg[BUFSIZ];

  va_start (ap, fmt);
  vsnprintf (msg, BUFSIZ, fmt, ap);
  va_end (ap);

  P4EST_GLOBAL_LERROR ("ERROR/\n");
  P4EST_GLOBAL_LERRORF ("ERROR: %s\n", msg);
  P4EST_GLOBAL_LERROR ("ERROR\\\n");
  return 1;
}

static int
usagerr (sc_options_t *opt, const char *msg)
{
  return usagerrf (opt, "%s", msg);
}

int
main (int argc, char **argv)
{
  sc_MPI_Comm         mpicomm;
  int                 mpiret;
  int                 ue, fa;
  sc_options_t       *opt;
  p4est_model_t      *model = NULL;
  const char         *fn_par;
  char               *model_name;
  char                filename[BUFSIZ], filename_temp[BUFSIZ];

  /* initialize MPI */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  /* initialize global context */
  mpicomm = sc_MPI_COMM_WORLD;

  /* set global logging options for p4est */
  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  /* initialize global application state */
  opt = sc_options_new (argv[0]);
  sc_options_add_int (opt, 'L', "maxlevel", &max_ref_level, P4EST_QMAXLEVEL,
                      "Maximum refinement level");
  sc_options_add_string (opt, 'F', "filename", &fn_par, "model.off",
                         "Input file in .off format");

  /* proceed in run-once loop for cleaner error checking */
  ue = 0;
  do {
    /* parse command line and assign configuration variables */
    fa = sc_options_parse (p4est_package_id, SC_LP_DEFAULT, opt, argc, argv);
    if (fa < 0 || fa != argc) {
      ue = usagerr (opt, "invalid option format or non-option argument");
      break;
    }
    P4EST_GLOBAL_PRODUCTIONF ("Manifold dimension is %d\n", P4EST_DIM);
    sc_options_print_summary (p4est_package_id, SC_LP_PRODUCTION, opt);

    /* check consistency of parameters */
    if (max_ref_level < 1 || max_ref_level > P4EST_QMAXLEVEL) {
      ue = usagerrf (opt, "maxlevel not between 1 and %d", P4EST_QMAXLEVEL);
    }
  }
  while (0);
  if (ue) {
    sc_options_print_usage (p4est_package_id, SC_LP_ERROR, opt, NULL);
  }
  strcpy (filename, fn_par);
  strcpy (filename_temp, fn_par);
  model_name = strtok (filename_temp, ".");

  /* setup appplication model */
  if (!ue && !triangulation_setup_model (&model, filename, model_name)) {
    P4EST_ASSERT (model == NULL);
    ue = usagerr (opt, "model-specific initialization error");
  }

  /* execute application model */
  if (!ue) {
    P4EST_ASSERT (model != NULL);
    SC3X (run_program (&mpicomm, model));
  }

  /* cleanup application model */
  if (model != NULL) {
    p4est_connectivity_destroy (model->conn);
    p4est_model_destroy (model);
    P4EST_FREE (model);
  }

  /* deinit main program */
  sc_options_destroy (opt);
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return ue ? EXIT_FAILURE : EXIT_SUCCESS;
}
