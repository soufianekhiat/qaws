/*
 * Test 171: offsets of meshes (qaws_mesh_3d_offset), small meshes
 *
 *   - a unit cube contoured outward by 0.1: closed, its volume the rounded
 *     offset's 1 + 0.6 + 0.03 pi + 0.004 pi / 3, its vertices 0.1 away;
 *     inward by 0.1: the cube of side 0.8, corners sharp;
 *   - the cube moved (mitred): the cubes of side 1.2 and 0.8 exactly, the
 *     corners 0.1 sqrt 3 away;
 *   - a sphere's shell 0.2 thick (contoured): the wall's volume;
 *   - a unit square sheet thickened by 0.2: moved, the box exactly;
 *     contoured, the slab with rounded rims;
 *   - a U whose slot 0.4 wide closes when moved out by 0.25: its faces
 *     crossing, resolved, closed, the slot's middle inside;
 *   - a slab 0.1 thick offset inward by 0.1: nothing left;
 *   - refusals: distance 0, an unknown kind;
 *   - a distance per vertex (MOVE): the cube's top out by 0.2, its bottom by
 *     0.1, the faces between them tilted; all alike, the uniform result;
 *   - ROUND: the cube out by 0.1 the rounded cube (its volume, every vertex
 *     within the chords' depth of 0.1), in by 0.1 the cube of side 0.8
 *     exactly; the U's slot closed; a sheet thickened with rounded rims;
 *     refused for a step of 0.001.
 */

#include "test_common.h"
#include "qaws_mesh_offset.h"
#include "qaws_mesh_boolean.h"
#include "qaws_mesh_bvh.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define OT_PI 3.14159265358979323846

/* a box's 8 corners and 12 triangles, outward */
static void ot_box_arrays(qaws_scalar* x, unsigned int* t, double const* lo, double const* hi)
{
	static unsigned int const f[36] = { 0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 2, 3, 7, 2, 7, 6, 1, 2, 6, 1, 6, 5, 0, 4, 7, 0, 7, 3 };
	unsigned int i;
	for (i = 0; i < 8; i++)
	{
		x[3 * i] = (qaws_scalar)((i == 1 || i == 2 || i == 5 || i == 6) ? hi[0] : lo[0]);
		x[3 * i + 1] = (qaws_scalar)((i == 2 || i == 3 || i == 6 || i == 7) ? hi[1] : lo[1]);
		x[3 * i + 2] = (qaws_scalar)(i >= 4 ? hi[2] : lo[2]);
	}
	memcpy(t, f, sizeof f);
}

static qaws_mesh_3d* ot_box(double const* lo, double const* hi)
{
	qaws_scalar x[24];
	unsigned int t[36];
	qaws_mesh_3d_desc d;
	qaws_mesh_3d* m = NULL;
	ot_box_arrays(x, t, lo, hi);
	memset(&d, 0, sizeof d);
	d.positions = x;
	d.vertex_count = 8;
	d.face_vertices = t;
	d.face_count = 12;
	qaws_mesh_3d_create(&m, &d);
	return m;
}

/* a unit sphere: rings between single poles */
static qaws_mesh_3d* ot_sphere(unsigned int nu, unsigned int nv)
{
	unsigned int n = nu * (nv - 1) + 2, i, j, k = 0;
	qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * n);
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 6 * nu * (nv - 1));
	qaws_mesh_3d_desc d;
	qaws_mesh_3d* m = NULL;
	for (j = 1; j < nv; j++)
		for (i = 0; i < nu; i++)
		{
			double th = OT_PI * j / nv, ph = 2 * OT_PI * i / nu;
			qaws_scalar* p = x + 3 * ((j - 1) * nu + i);
			p[0] = (qaws_scalar)(sin(th) * cos(ph)); p[1] = (qaws_scalar)(sin(th) * sin(ph)); p[2] = (qaws_scalar)cos(th);
		}
	x[3 * (n - 2)] = 0; x[3 * (n - 2) + 1] = 0; x[3 * (n - 2) + 2] = 1;
	x[3 * (n - 1)] = 0; x[3 * (n - 1) + 1] = 0; x[3 * (n - 1) + 2] = -1;
	for (i = 0; i < nu; i++)
	{
		unsigned int a = i, b = (i + 1) % nu, o = (nv - 2) * nu;
		t[k++] = n - 2; t[k++] = a; t[k++] = b;
		t[k++] = n - 1; t[k++] = o + b; t[k++] = o + a;
	}
	for (j = 0; j + 2 < nv; j++)
		for (i = 0; i < nu; i++)
		{
			unsigned int a = j * nu + i, b = j * nu + (i + 1) % nu;
			t[k++] = a; t[k++] = b + nu; t[k++] = b;
			t[k++] = a; t[k++] = a + nu; t[k++] = b + nu;
		}
	memset(&d, 0, sizeof d);
	d.positions = x;
	d.vertex_count = n;
	d.face_vertices = t;
	d.face_count = k / 3;
	qaws_mesh_3d_create(&m, &d);
	free(x);
	free(t);
	return m;
}

/* a unit square sheet at z = 0, n x n cells */
static qaws_mesh_3d* ot_sheet(unsigned int n)
{
	unsigned int i, j, k = 0, nv = (n + 1) * (n + 1);
	qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 6 * n * n);
	qaws_mesh_3d_desc d;
	qaws_mesh_3d* m = NULL;
	for (j = 0; j <= n; j++)
		for (i = 0; i <= n; i++)
		{
			x[3 * (j * (n + 1) + i)] = (qaws_scalar)((double)i / n);
			x[3 * (j * (n + 1) + i) + 1] = (qaws_scalar)((double)j / n);
			x[3 * (j * (n + 1) + i) + 2] = 0;
		}
	for (j = 0; j < n; j++)
		for (i = 0; i < n; i++)
		{
			unsigned int v = j * (n + 1) + i;
			t[k++] = v; t[k++] = v + 1; t[k++] = v + n + 2;
			t[k++] = v; t[k++] = v + n + 2; t[k++] = v + n + 1;
		}
	memset(&d, 0, sizeof d);
	d.positions = x;
	d.vertex_count = nv;
	d.face_vertices = t;
	d.face_count = 2 * n * n;
	qaws_mesh_3d_create(&m, &d);
	free(x);
	free(t);
	return m;
}

/* the enclosed volume (divergence theorem over the triangles) */
static double ot_volume(qaws_mesh_3d const* m)
{
	unsigned int nv, nt, i;
	qaws_scalar* x;
	unsigned int* t;
	double v = 0;
	if (!m)
		return 0;
	nv = qaws_mesh_3d_get_vertex_count(m);
	nt = qaws_mesh_3d_get_triangle_count(m);
	x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	t = (unsigned int*)malloc(sizeof(unsigned int) * 3 * nt);
	qaws_mesh_3d_get_positions(x, m);
	qaws_mesh_3d_get_triangles(t, m);
	for (i = 0; i < nt; i++)
	{
		qaws_scalar const *a = x + 3 * t[3 * i], *b = x + 3 * t[3 * i + 1], *c = x + 3 * t[3 * i + 2];
		v += ((double)a[0] * ((double)b[1] * c[2] - (double)b[2] * c[1]) - (double)a[1] * ((double)b[0] * c[2] - (double)b[2] * c[0]) +
			(double)a[2] * ((double)b[0] * c[1] - (double)b[1] * c[0])) / 6;
	}
	free(x);
	free(t);
	return v;
}

static qaws_status ot_offset(qaws_mesh_3d** out, qaws_mesh_offset_report* r, qaws_mesh_3d const* m, qaws_mesh_offset_kind kind, qaws_mesh_offset_method method,
	double dist, unsigned int res)
{
	qaws_mesh_offset_desc d;
	memset(&d, 0, sizeof d);
	d.kind = kind;
	d.method = method;
	d.distance = (qaws_scalar)dist;
	d.resolution = res;
	return qaws_mesh_3d_offset(out, r, m, &d);
}

static void test_cube(void)
{
	static double const lo[3] = { 0, 0, 0 }, hi[3] = { 1, 1, 1 };
	qaws_mesh_3d *c = ot_box(lo, hi), *o = NULL;
	qaws_mesh_offset_report r;
	char msg[300];
	double v, ve = 1 + 0.6 + 0.3 * OT_PI * 0.1 + 4.0 / 3 * OT_PI * 0.001;
	qaws_status st;
	st = ot_offset(&o, &r, c, QAWS_MESH_OFFSET_SURFACE, QAWS_MESH_OFFSET_CONTOUR, 0.1, 40);
	v = ot_volume(o);
	sprintf(msg, "cube contoured out by 0.1: volume %.4f (the offset's %.4f), distances %.4f to %.4f (cell %.3f), closed %d; the distance evaluated at %u of %u grid points",
		v, ve, (double)r.min_distance, (double)r.max_distance, (double)r.cell, r.closed, r.evaluated, r.grid[0] * r.grid[1] * r.grid[2]);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && o && r.closed && fabs(v - ve) < 0.01 * ve && fabs((double)r.min_distance - 0.1) < 0.1 * (double)r.cell &&
		fabs((double)r.max_distance - 0.1) < 0.1 * (double)r.cell && 10 * r.evaluated < 6 * r.grid[0] * r.grid[1] * r.grid[2], msg);
	qaws_mesh_3d_destroy(o);
	o = NULL;
	st = ot_offset(&o, &r, c, QAWS_MESH_OFFSET_SURFACE, QAWS_MESH_OFFSET_CONTOUR, -0.1, 40);
	v = ot_volume(o);
	sprintf(msg, "cube contoured in by 0.1: volume %.5f (0.512), %u sharp vertices, closed %d", v, r.sharp, r.closed);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && o && r.closed && fabs(v - 0.512) < 0.002 && r.sharp >= 8, msg);
	qaws_mesh_3d_destroy(o);
	o = NULL;
	st = ot_offset(&o, &r, c, QAWS_MESH_OFFSET_SURFACE, QAWS_MESH_OFFSET_MOVE, 0.1, 0);
	v = ot_volume(o);
	sprintf(msg, "cube moved out by 0.1: volume %.12f (1.728), its corners %.6f to %.6f away (0.1 sqrt 3), closed %d", v,
		(double)r.min_distance, (double)r.max_distance, r.closed);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && o && r.closed && fabs(v - 1.728) < 1e-9 && fabs((double)r.max_distance - 0.1 * sqrt(3.0)) < 1e-9, msg);
	qaws_mesh_3d_destroy(o);
	o = NULL;
	st = ot_offset(&o, &r, c, QAWS_MESH_OFFSET_SURFACE, QAWS_MESH_OFFSET_MOVE, -0.1, 0);
	v = ot_volume(o);
	sprintf(msg, "cube moved in by 0.1: volume %.12f (0.512), closed %d", v, r.closed);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && o && r.closed && fabs(v - 0.512) < 1e-9, msg);
	qaws_mesh_3d_destroy(o);
	qaws_mesh_3d_destroy(c);
}

static void test_shell(void)
{
	qaws_mesh_3d *s = ot_sphere(32, 16), *o = NULL;
	qaws_mesh_offset_report r;
	char msg[300];
	double vs = ot_volume(s), v, ve, ri;
	qaws_status st = ot_offset(&o, &r, s, QAWS_MESH_OFFSET_SHELL, QAWS_MESH_OFFSET_CONTOUR, 0.2, 40);
	v = ot_volume(o);
	/* the inner surface 0.2 inside the faceted sphere: a sphere of its inradius-ish, about 0.8 */
	ri = pow(vs / (4.0 / 3 * OT_PI), 1.0 / 3) - 0.2;
	ve = vs - 4.0 / 3 * OT_PI * ri * ri * ri;
	sprintf(msg, "a sphere's shell 0.2 thick: volume %.4f (about %.4f), inner distances %.4f to %.4f, closed %d", v, ve,
		(double)r.min_distance, (double)r.max_distance, r.closed);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && o && r.closed && fabs(v - ve) < 0.03 * ve && fabs((double)r.min_distance - 0.2) < 0.2 * (double)r.cell, msg);
	qaws_mesh_3d_destroy(o);
	qaws_mesh_3d_destroy(s);
}

static void test_thicken(void)
{
	qaws_mesh_3d *s = ot_sheet(6), *o = NULL;
	qaws_mesh_offset_report r;
	char msg[300];
	double v, ve = 0.2 + 4 * OT_PI * 0.01 / 2 + 4.0 / 3 * OT_PI * 0.001;
	qaws_status st = ot_offset(&o, &r, s, QAWS_MESH_OFFSET_THICKEN, QAWS_MESH_OFFSET_MOVE, 0.2, 0);
	v = ot_volume(o);
	sprintf(msg, "a sheet thickened by 0.2, moved: volume %.12f (0.2), closed %d", v, r.closed);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && o && r.closed && fabs(v - 0.2) < 1e-9, msg);
	qaws_mesh_3d_destroy(o);
	o = NULL;
	st = ot_offset(&o, &r, s, QAWS_MESH_OFFSET_THICKEN, QAWS_MESH_OFFSET_CONTOUR, 0.2, 40);
	v = ot_volume(o);
	sprintf(msg, "contoured: volume %.4f (the slab with round rims %.4f), closed %d", v, ve, r.closed);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && o && r.closed && fabs(v - ve) < 0.02 * ve, msg);
	qaws_mesh_3d_destroy(o);
	qaws_mesh_3d_destroy(s);
}

/* a U: a base and two prongs, the slot between x 0.3 and 0.7 */
static qaws_mesh_3d* ot_u(void)
{
	static double const b0[3] = { 0, 0, 0 }, b1[3] = { 1, 0.3, 0.3 }, p0[3] = { 0, 0.2, 0 }, p1[3] = { 0.3, 1, 0.3 }, q0[3] = { 0.7, 0.2, 0 }, q1[3] = { 1, 1, 0.3 };
	qaws_scalar x[48];
	unsigned int t[72], i;
	qaws_mesh_3d_desc d, e;
	qaws_mesh_3d *ab = NULL, *u = NULL;
	qaws_scalar* y;
	unsigned int* s;
	unsigned int nv, nt;
	ot_box_arrays(x, t, b0, b1);
	ot_box_arrays(x + 24, t + 36, p0, p1);
	for (i = 36; i < 72; i++) t[i] += 8;
	memset(&d, 0, sizeof d);
	d.positions = x;
	d.vertex_count = 16;
	d.face_vertices = t;
	d.face_count = 24;
	qaws_mesh_boolean(&ab, NULL, &d, NULL, QAWS_BOOLEAN_UNION, NULL);
	nv = qaws_mesh_3d_get_vertex_count(ab);
	nt = qaws_mesh_3d_get_triangle_count(ab);
	y = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * (nv + 8));
	s = (unsigned int*)malloc(sizeof(unsigned int) * 3 * (nt + 12));
	qaws_mesh_3d_get_positions(y, ab);
	qaws_mesh_3d_get_triangles(s, ab);
	ot_box_arrays(y + 3 * nv, s + 3 * nt, q0, q1);
	for (i = 3 * nt; i < 3 * nt + 36; i++) s[i] += nv;
	memset(&e, 0, sizeof e);
	e.positions = y;
	e.vertex_count = nv + 8;
	e.face_vertices = s;
	e.face_count = nt + 12;
	qaws_mesh_boolean(&u, NULL, &e, NULL, QAWS_BOOLEAN_UNION, NULL);
	qaws_mesh_3d_destroy(ab);
	free(y);
	free(s);
	return u;
}

static void test_u(void)
{
	qaws_mesh_3d *u = ot_u(), *o = NULL;
	qaws_mesh_offset_report r;
	qaws_mesh_bvh* b = NULL;
	qaws_scalar w = 0, mid[3] = { 0.5, 0.7, 0.15 };
	char msg[300];
	qaws_status st = ot_offset(&o, &r, u, QAWS_MESH_OFFSET_SURFACE, QAWS_MESH_OFFSET_MOVE, 0.25, 0);
	if (o)
	{
		qaws_mesh_bvh_create_mesh_3d(&b, o, NULL);
		qaws_mesh_bvh_winding_number(&w, b, mid, 1, NULL);
	}
	sprintf(msg, "a U moved out by 0.25 across its 0.4 slot: %u pairs of faces crossing resolved, closed %d, the slot's middle's winding %.3f",
		r.intersect.pairs, r.closed, (double)w);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && u && o && r.intersect.pairs > 0 && r.closed && fabs((double)w - 1) < 1e-3, msg);
	qaws_mesh_bvh_destroy(b);
	qaws_mesh_3d_destroy(o);
	qaws_mesh_3d_destroy(u);
}

static void test_empty(void)
{
	static double const lo[3] = { 0, 0, 0 }, hi[3] = { 1, 1, 0.1 };
	qaws_mesh_3d *s = ot_box(lo, hi), *a = NULL, *b = NULL;
	qaws_mesh_offset_report r;
	qaws_mesh_offset_desc d;
	qaws_status sa = ot_offset(&a, &r, s, QAWS_MESH_OFFSET_SURFACE, QAWS_MESH_OFFSET_CONTOUR, -0.1, 32), sb = ot_offset(&b, &r, s, QAWS_MESH_OFFSET_SURFACE, QAWS_MESH_OFFSET_MOVE, -0.1, 0);
	qaws_status s0, s1;
	qaws_mesh_3d* z = NULL;
	char msg[300];
	memset(&d, 0, sizeof d);
	s0 = qaws_mesh_3d_offset(&z, NULL, s, &d);
	d.distance = 1;
	d.kind = (qaws_mesh_offset_kind)7;
	s1 = qaws_mesh_3d_offset(&z, NULL, s, &d);
	sprintf(msg, "a slab 0.1 thick offset in by 0.1: nothing left (contoured %d, moved %d); refused: distance 0 (%d), an unknown kind (%d)",
		sa == QAWS_STATUS_OK && !a, sb == QAWS_STATUS_OK && !b, (int)s0, (int)s1);
	printf("  %s\n", msg);
	TEST_ASSERT(sa == QAWS_STATUS_OK && !a && sb == QAWS_STATUS_OK && !b && s0 == QAWS_STATUS_INVALID_ARGUMENT && s1 == QAWS_STATUS_INVALID_ARGUMENT && !z, msg);
	qaws_mesh_3d_destroy(a);
	qaws_mesh_3d_destroy(b);
	qaws_mesh_3d_destroy(s);
}

/* a distance per vertex (MOVE): the cube's top corners out by 0.2, its bottom ones by 0.1 */
static void test_varying(void)
{
	static double const lo[3] = { -0.5, -0.5, -0.5 }, hi[3] = { 0.5, 0.5, 0.5 };
	qaws_mesh_3d *cube = ot_box(lo, hi), *o = NULL, *u = NULL, *z = NULL;
	qaws_mesh_offset_desc d;
	qaws_mesh_offset_report r;
	qaws_scalar dist[8], same[8];
	double zmin = 1e30, zmax = -1e30, v = 0;
	unsigned int i, nv;
	qaws_status st, s0;
	char msg[300];
	for (i = 0; i < 8; i++)
	{
		qaws_vec3 q = qaws_mesh_3d_get_vertex(cube, i);
		dist[i] = (qaws_scalar)(q.z > 0 ? 0.2 : 0.1);
		same[i] = (qaws_scalar)0.1;
	}
	memset(&d, 0, sizeof d);
	d.kind = QAWS_MESH_OFFSET_SURFACE;
	d.method = QAWS_MESH_OFFSET_MOVE;
	d.distance = (qaws_scalar)1;
	d.distances = dist;
	st = qaws_mesh_3d_offset(&o, &r, cube, &d);
	if (o)
	{
		nv = qaws_mesh_3d_get_vertex_count(o);
		for (i = 0; i < nv; i++)
		{
			qaws_vec3 q = qaws_mesh_3d_get_vertex(o, i);
			zmin = fmin(zmin, (double)q.z);
			zmax = fmax(zmax, (double)q.z);
		}
	}
	d.distances = same;
	qaws_mesh_3d_offset(&u, NULL, cube, &d);
	if (u)
		v = ot_volume(u);
	d.method = QAWS_MESH_OFFSET_CONTOUR;
	s0 = qaws_mesh_3d_offset(&z, NULL, cube, &d);
	sprintf(msg, "a distance per vertex: the cube's top out by 0.2, its bottom by 0.1: z from %.12f to %.12f; all 0.1: volume %.12f (1.728); with CONTOUR refused (%d)",
		zmin, zmax, v, (int)s0);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && o && fabs(zmin + 0.6) < 1e-12 && fabs(zmax - 0.7) < 1e-12 && fabs(v - 1.728) < 1e-9 && s0 == QAWS_STATUS_INVALID_ARGUMENT && !z, msg);
	qaws_mesh_3d_destroy(o);
	qaws_mesh_3d_destroy(u);
	qaws_mesh_3d_destroy(cube);
}

/* ROUND: the points within the distance, as a union of slabs, wedges and cones (steps of pi / 8: the chords inside by
   1 - cos(pi / 16)^2 of the distance at most) */
static void test_round(void)
{
	static double const lo[3] = { 0, 0, 0 }, hi[3] = { 1, 1, 1 };
	qaws_mesh_3d *c = ot_box(lo, hi), *u = ot_u(), *sh = ot_sheet(4), *o = NULL;
	qaws_mesh_offset_report r;
	qaws_mesh_offset_desc d;
	char msg[300];
	double v, ve = 1 + 0.6 + 0.3 * OT_PI * 0.1 + 4.0 / 3 * OT_PI * 0.001, vt = 0.2 + 4 * 0.5 * OT_PI * 0.01 + 4.0 / 3 * OT_PI * 0.001, chord = 1 - cos(OT_PI / 16) * cos(OT_PI / 16);
	qaws_status st, s2;
	st = ot_offset(&o, &r, c, QAWS_MESH_OFFSET_SURFACE, QAWS_MESH_OFFSET_ROUND, 0.1, 0);
	v = ot_volume(o);
	sprintf(msg, "cube rounded out by 0.1: volume %.5f (the offset's %.5f), its vertices %.6f to %.12f away, closed %d", v, ve,
		(double)r.min_distance, (double)r.max_distance, r.closed);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && o && r.closed && v < ve && v > ve - 0.006 && (double)r.min_distance > 0.1 * (1 - chord) - 1e-9 &&
		(double)r.max_distance < 0.1 + 1e-9, msg);
	qaws_mesh_3d_destroy(o);
	o = NULL;
	st = ot_offset(&o, &r, c, QAWS_MESH_OFFSET_SURFACE, QAWS_MESH_OFFSET_ROUND, -0.1, 0);
	v = ot_volume(o);
	sprintf(msg, "cube rounded in by 0.1: volume %.12f (0.512), closed %d", v, r.closed);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && o && r.closed && fabs(v - 0.512) < 1e-9, msg);
	qaws_mesh_3d_destroy(o);
	o = NULL;
	st = ot_offset(&o, &r, u, QAWS_MESH_OFFSET_SURFACE, QAWS_MESH_OFFSET_ROUND, 0.25, 0);
	sprintf(msg, "a U rounded out by 0.25, its slot closed: distances %.4f to %.12f (0.25, the chords and the saddles' balls inside by 5%% at most), closed %d",
		(double)r.min_distance, (double)r.max_distance, r.closed);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && o && r.closed && (double)r.min_distance > 0.95 * 0.25 && (double)r.max_distance < 0.25 + 1e-9, msg);
	qaws_mesh_3d_destroy(o);
	o = NULL;
	st = ot_offset(&o, &r, sh, QAWS_MESH_OFFSET_THICKEN, QAWS_MESH_OFFSET_ROUND, 0.2, 0);
	v = ot_volume(o);
	sprintf(msg, "a sheet thickened by 0.2, rounded: volume %.5f (%.5f), its vertices %.6f to %.12f away, closed %d", v, vt,
		(double)r.min_distance, (double)r.max_distance, r.closed);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && o && r.closed && v < vt && v > vt - 0.004 && (double)r.min_distance > 0.1 * (1 - chord) - 1e-9 &&
		(double)r.max_distance < 0.1 + 1e-9, msg);
	qaws_mesh_3d_destroy(o);
	o = NULL;
	memset(&d, 0, sizeof d);
	d.kind = QAWS_MESH_OFFSET_SURFACE;
	d.method = QAWS_MESH_OFFSET_ROUND;
	d.distance = (qaws_scalar)0.1;
	d.round_angle = (qaws_scalar)0.001;
	s2 = qaws_mesh_3d_offset(&o, NULL, c, &d);
	sprintf(msg, "ROUND refused for a step of 0.001 (%d)", (int)s2);
	printf("  %s\n", msg);
	TEST_ASSERT(s2 == QAWS_STATUS_INVALID_ARGUMENT && !o, msg);
	qaws_mesh_3d_destroy(c);
	qaws_mesh_3d_destroy(u);
	qaws_mesh_3d_destroy(sh);
}

int test_171_mesh_offset_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 171: offsets of meshes\n");
	test_cube();
	test_shell();
	test_thicken();
	test_u();
	test_empty();
	test_varying();
	test_round();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
