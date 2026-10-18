/*
 * meshing_showcase.c - the meshing module: meshes, references, remeshing,
 * quads, tetrahedra and hexahedra
 *
 * Writes PNG figures into ./showcase/:
 *   mesh1_core.png     the corner store: a torus and its dual (the stored
 *                      adjacency), a cube's creases kept through random edits,
 *                      a linear attribute kept exact through edits
 *   mesh2_reference.png  the surface between the vertices: a coarse blob flat and
 *                      as its PN reference, a coarse can rounded without features
 *                      and sharp with its creases, curvature from the reference
 *   mesh3_remesh.png   remeshing on the reference: uniform and adaptive by
 *                      curvature, creases kept, attributes and their seams transferred
 *   mesh4_quads.png    quads by Instant Meshes on the reference: a cross field and
 *                      its singularities, quad-dominant and pure, along creases, a torus
 *
 * Each figure lives in meshing_showcase/NN_*.c, compiled as one translation
 * unit with the rasterizer in example_raster.h.
 */

#include "qaws.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
#include <direct.h>
#define MAKE_DIR(p) _mkdir(p)
#else
#include <sys/stat.h>
#define MAKE_DIR(p) mkdir(p, 0755)
#endif

#include "example_raster.h"

#define MS_PI 3.14159265358979323846

/* ---- shared helpers ---------------------------------------------------- */

static unsigned int ms_rand(unsigned int* s)
{
	*s = *s * 1664525u + 1013904223u;
	return *s >> 8;
}

static double ms_uniform(unsigned int* s)
{
	return (ms_rand(s) & 0xFFFFFF) / (double)0x1000000;
}

/* blue to red through white */
static void ms_colormap(double t, double* rgb)
{
	static double const k[5][3] = { { 59, 76, 192 }, { 141, 176, 254 }, { 221, 221, 221 }, { 244, 154, 123 }, { 180, 4, 38 } };
	double s;
	int i, j;
	t = t < 0 ? 0 : t > 1 ? 1 : t;
	s = t * 4;
	i = (int)s;
	if (i > 3)
		i = 3;
	s -= i;
	for (j = 0; j < 3; j++)
		rgb[j] = k[i][j] + s * (k[i + 1][j] - k[i][j]);
}

static qaws_mesh_3d* ms_make(double const* x, unsigned int nv, unsigned int const* f, unsigned int const* fs, unsigned int nf)
{
	qaws_mesh_3d_desc d;
	qaws_mesh_3d* m = NULL;
	qaws_scalar* p = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * (nv ? nv : 1));
	unsigned int i;
	for (i = 0; i < 3 * nv; i++)
		p[i] = (qaws_scalar)x[i];
	memset(&d, 0, sizeof d);
	d.positions = p;
	d.vertex_count = nv;
	d.face_vertices = f;
	d.face_sizes = fs;
	d.face_count = nf;
	qaws_mesh_3d_create(&d, &m);
	free(p);
	return m;
}

/* the mesh's positions as doubles (caller frees) */
static double* ms_positions(qaws_mesh_3d const* m)
{
	unsigned int n = qaws_mesh_3d_get_vertex_count(m), i;
	qaws_scalar* p = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * (n ? n : 1));
	double* x = (double*)malloc(sizeof(double) * 3 * (n ? n : 1));
	qaws_mesh_3d_get_positions(m, p);
	for (i = 0; i < 3 * n; i++)
		x[i] = (double)p[i];
	free(p);
	return x;
}

/* A mesh drawn: smooth (vertex normals) or flat, a colour per vertex
   (vrgb, 3 each) or one colour, a wireframe of its edges (wire may be NULL). */
static void ms_draw(er_canvas* c, qaws_mesh_3d const* m, double const* rgb, double const* vrgb, int smooth, double const* wire, double wire_width)
{
	unsigned int nv = qaws_mesh_3d_get_vertex_count(m), nt = qaws_mesh_3d_get_triangle_count(m), ne = qaws_mesh_3d_get_edge_count(m), i;
	double* x = ms_positions(m);
	double* n = NULL;
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 3 * (nt ? nt : 1));
	unsigned int* e = (unsigned int*)malloc(sizeof(unsigned int) * 2 * (ne ? ne : 1));
	qaws_mesh_3d_get_triangles(m, t);
	if (smooth)
	{
		qaws_scalar* ns = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * (nv ? nv : 1));
		n = (double*)malloc(sizeof(double) * 3 * (nv ? nv : 1));
		qaws_mesh_3d_get_vertex_normals(m, ns);
		for (i = 0; i < 3 * nv; i++)
			n[i] = (double)ns[i];
		free(ns);
	}
	for (i = 0; (rgb || vrgb) && i < nt; i++)   /* (no colour: the edges only) */
	{
		unsigned int a = t[3 * i], b = t[3 * i + 1], d = t[3 * i + 2];
		er_triangle(c, x + 3 * a, x + 3 * b, x + 3 * d, n ? n + 3 * a : NULL, n ? n + 3 * b : NULL, n ? n + 3 * d : NULL,
			vrgb ? vrgb + 3 * a : rgb, vrgb ? vrgb + 3 * b : rgb, vrgb ? vrgb + 3 * d : rgb);
	}
	if (wire)
	{
		qaws_mesh_3d_get_edges(m, e);
		for (i = 0; i < ne; i++)
			er_line(c, x + 3 * e[2 * i], x + 3 * e[2 * i + 1], wire, wire_width);
	}
	free(x); free(n); free(t); free(e);
}

/* a panel fitted to a mesh, with its title and a line of numbers */
static void ms_panel(er_canvas* c, int x, int y, int w, int h, qaws_mesh_3d const* m, double yaw, double pitch, char const* title,
	char const* sub)
{
	double* p = ms_positions(m);
	er_panel(c, x, y + 72, w, h - 72, p, qaws_mesh_3d_get_vertex_count(m), yaw, pitch, 0.06);
	free(p);
	er_text(c, x + 14, y + 12, title, 3, 30, 30, 30);
	if (sub)
		er_text(c, x + 14, y + 46, sub, 2, 95, 95, 95);
}

#include "meshing_showcase/01_core.c"
#include "meshing_showcase/02_reference.c"
#include "meshing_showcase/03_remesh.c"
#include "meshing_showcase/04_quads.c"

int main(void)
{
	setvbuf(stdout, NULL, _IONBF, 0);
	MAKE_DIR("showcase");
	demo_core();
	demo_reference();
	demo_remesh();
	demo_quads();
	return 0;
}
