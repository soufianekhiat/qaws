/*
 * Test 175: tetrahedra keeping the border (qaws_tetrahedralize_desc
 * keep_border), small meshes
 *
 *   - Schonhardt's polyhedron (a triangular prism twisted by 30 degrees, its
 *     sides cut along the reflex diagonals): no tetrahedra of its own
 *     vertices fill it, so the recovery adds border points; kept, they move
 *     inside: the border's 8 triangles and 6 vertices as given, the volume
 *     exact, the cells valid;
 *   - a sphere's rings: the border kept triangle for triangle, as without;
 *   - the plain recovery, for comparison: the border split.
 */

#include "test_common.h"
#include "qaws_tetrahedralize.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define TB_PI 3.14159265358979323846

/* the twisted prism: the bottom a0 a1 a2 at z = 0, the top b0 b1 b2 at z = 1 turned by twist; the sides cut along a_i b_{i+1}
   (reflex) or a_{i+1} b_i */
static qaws_mesh_3d* tb_prism(double twist, int reflex)
{
	qaws_scalar x[18];
	unsigned int t[24], i, k = 0;
	qaws_mesh_3d_desc d;
	qaws_mesh_3d* m = NULL;
	for (i = 0; i < 3; i++)
	{
		double th = TB_PI / 2 + 2 * TB_PI * i / 3;
		x[3 * i] = (qaws_scalar)cos(th);
		x[3 * i + 1] = (qaws_scalar)sin(th);
		x[3 * i + 2] = 0;
		x[9 + 3 * i] = (qaws_scalar)cos(th + twist);
		x[9 + 3 * i + 1] = (qaws_scalar)sin(th + twist);
		x[9 + 3 * i + 2] = 1;
	}
	t[k++] = 0; t[k++] = 2; t[k++] = 1;
	t[k++] = 3; t[k++] = 4; t[k++] = 5;
	for (i = 0; i < 3; i++)
	{
		unsigned int a = i, an = (i + 1) % 3, b = 3 + i, bn = 3 + (i + 1) % 3;
		if (reflex)
		{
			t[k++] = a; t[k++] = an; t[k++] = bn;
			t[k++] = a; t[k++] = bn; t[k++] = b;
		}
		else
		{
			t[k++] = a; t[k++] = an; t[k++] = b;
			t[k++] = an; t[k++] = bn; t[k++] = b;
		}
	}
	memset(&d, 0, sizeof d);
	d.positions = x;
	d.vertex_count = 6;
	d.face_vertices = t;
	d.face_count = 8;
	qaws_mesh_3d_create(&m, &d);
	return m;
}

static double tb_surface_volume(qaws_mesh_3d const* m)
{
	unsigned int nt = qaws_mesh_3d_get_triangle_count(m), i;
	double v = 0;
	for (i = 0; i < nt; i++)
	{
		unsigned int t[3];
		qaws_vec3 a, b, c;
		qaws_mesh_3d_get_triangle(t, m, i);
		a = qaws_mesh_3d_get_vertex(m, t[0]);
		b = qaws_mesh_3d_get_vertex(m, t[1]);
		c = qaws_mesh_3d_get_vertex(m, t[2]);
		v += ((double)a.x * ((double)b.y * c.z - (double)b.z * c.y) - (double)a.y * ((double)b.x * c.z - (double)b.z * c.x) +
			(double)a.z * ((double)b.x * c.y - (double)b.y * c.x)) / 6;
	}
	return v;
}

typedef struct tb_result
{
	qaws_status st;
	qaws_mesh_report r;
	qaws_volume_quality q;
	unsigned int border_faces, border_vertices, valid;
} tb_result;

static tb_result tb_run(qaws_mesh_3d const* s, int border_only, int keep)
{
	qaws_tetrahedralize_desc d;
	qaws_volume_mesh* v = NULL;
	tb_result o;
	memset(&o, 0, sizeof o);
	memset(&d, 0, sizeof d);
	d.border_only = border_only;
	d.keep_border = keep;
	o.st = qaws_volume_mesh_tetrahedralize(&v, &o.r, s, &d);
	if (v)
	{
		qaws_mesh_3d* b = NULL;
		qaws_volume_mesh_measure_quality(&o.q, v);
		qaws_volume_mesh_get_boundary(&b, v);
		if (b)
		{
			o.border_faces = qaws_mesh_3d_get_face_count(b);
			o.border_vertices = qaws_mesh_3d_get_vertex_count(b);
			qaws_mesh_3d_destroy(b);
		}
		o.valid = qaws_volume_mesh_validate(v) == QAWS_STATUS_OK;
		qaws_volume_mesh_destroy(v);
	}
	return o;
}

static void test_schonhardt(void)
{
	qaws_mesh_3d *p = tb_prism(TB_PI / 6, 1), *q = tb_prism(TB_PI / 6, 0);
	tb_result a = tb_run(p, 1, 0), b = tb_run(p, 1, 1), c = tb_run(q, 1, 0), d = tb_run(p, 0, 1);
	double vol = tb_surface_volume(p);
	char msg[500];
	sprintf(msg, "Schonhardt's polyhedron: recovered %u border points (border %u triangles, %u vertices); kept: %u moved inside, border %u triangles %u vertices, volume %.15f (%.15f), dihedral %.1f..%.1f, valid %d; the convex twist needs none (%u); filled inside too: border %u, moved %u, valid %d",
		a.r.splits, a.border_faces, a.border_vertices, b.r.collapses, b.border_faces, b.border_vertices, (double)b.q.volume, vol, (double)b.q.min_dihedral,
		(double)b.q.max_dihedral, b.valid, c.r.splits, d.border_faces, d.r.collapses, d.valid);
	printf("  %s\n", msg);
	TEST_ASSERT(a.st == QAWS_STATUS_OK && a.r.splits > 0 && a.border_faces > 8 && b.st == QAWS_STATUS_OK && b.r.collapses == b.r.splits && b.border_faces == 8 &&
		b.border_vertices == 6 && fabs((double)b.q.volume - vol) < 1e-12 && b.valid && b.q.min_dihedral > 0 && c.r.splits == 0 && d.st == QAWS_STATUS_OK &&
		d.border_faces == 8 && d.valid, msg);
	qaws_mesh_3d_destroy(p);
	qaws_mesh_3d_destroy(q);
}

/* a sphere of radius 1: rings between single poles */
static qaws_mesh_3d* tb_sphere(unsigned int nu, unsigned int nv)
{
	unsigned int n = nu * (nv - 1) + 2, i, j, k = 0;
	qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * n);
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 6 * nu * (nv - 1));
	qaws_mesh_3d_desc d;
	qaws_mesh_3d* m = NULL;
	for (j = 1; j < nv; j++)
		for (i = 0; i < nu; i++)
		{
			double th = TB_PI * j / nv, ph = 2 * TB_PI * i / nu;
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

static void test_sphere(void)
{
	qaws_mesh_3d* s = tb_sphere(16, 8);
	tb_result a = tb_run(s, 0, 0), b = tb_run(s, 0, 1);
	char msg[300];
	sprintf(msg, "a sphere's rings (%u triangles): border %u -> %u triangles (%u points added, %u moved inside), valid %d, dihedral %.1f",
		qaws_mesh_3d_get_face_count(s), a.border_faces, b.border_faces, b.r.splits, b.r.collapses, b.valid, (double)b.q.min_dihedral);
	printf("  %s\n", msg);
	TEST_ASSERT(b.st == QAWS_STATUS_OK && b.border_faces == qaws_mesh_3d_get_face_count(s) && b.valid && b.q.min_dihedral > 0, msg);
	qaws_mesh_3d_destroy(s);
}

int test_175_tetrahedra_border_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 175: tetrahedra keeping the border\n");
	test_schonhardt();
	test_sphere();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
