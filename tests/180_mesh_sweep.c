/*
 * Test 180: swept solids (qaws_mesh_sweep.h), small meshes
 *
 *   - a unit cube translated by (1, 1, 0): its Minkowski sum with the
 *     segment, volume 3 exactly, one slab;
 *   - a thin box turned a quarter about the z axis: the annular sector's
 *     volume (pi / 4 (r1^2 - r0^2) times its height, and the boxes at both
 *     ends half outside it: 2 h), within the chords;
 *   - a small ball along a Hermite curve: a tube, its volume the ball's
 *     cross-section along the length plus the ball (within the polyhedra's
 *     own error);
 *   - a vertex motion of another vertex count refused.
 */

#include "test_common.h"
#include "qaws_mesh_sweep.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define SW_PI 3.14159265358979323846

static double sw_volume(qaws_mesh_3d const* m)
{
	unsigned int nt = m ? qaws_mesh_3d_get_triangle_count(m) : 0, i;
	double v = 0;
	for (i = 0; i < nt; i++)
	{
		unsigned int c[3];
		qaws_vec3 a, b, e;
		qaws_mesh_3d_get_triangle(c, m, i);
		a = qaws_mesh_3d_get_vertex(m, c[0]);
		b = qaws_mesh_3d_get_vertex(m, c[1]);
		e = qaws_mesh_3d_get_vertex(m, c[2]);
		v += ((double)a.x * ((double)b.y * e.z - (double)b.z * e.y) + (double)a.y * ((double)b.z * e.x - (double)b.x * e.z) +
			(double)a.z * ((double)b.x * e.y - (double)b.y * e.x)) / 6;
	}
	return v;
}

/* the box [x0, x1] x [y0, y1] x [z0, z1], outward */
static qaws_mesh_3d* sw_box(double x0, double x1, double y0, double y1, double z0, double z1)
{
	static unsigned int const f[36] = { 0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 2, 3, 7, 2, 7, 6, 1, 2, 6, 1, 6, 5, 0, 4, 7, 0, 7, 3 };
	qaws_scalar x[24];
	qaws_mesh_3d_desc d;
	qaws_mesh_3d* m = NULL;
	unsigned int i;
	for (i = 0; i < 8; i++)
	{
		x[3 * i] = (qaws_scalar)((i == 1 || i == 2 || i == 5 || i == 6) ? x1 : x0);
		x[3 * i + 1] = (qaws_scalar)((i == 2 || i == 3 || i == 6 || i == 7) ? y1 : y0);
		x[3 * i + 2] = (qaws_scalar)(i >= 4 ? z1 : z0);
	}
	memset(&d, 0, sizeof d);
	d.positions = x; d.vertex_count = 8; d.face_vertices = f; d.face_count = 12;
	qaws_mesh_3d_create(&m, &d);
	return m;
}

/* an icosahedron split once, on the sphere of radius r: 80 faces */
static qaws_mesh_3d* sw_ball(double r)
{
	static int const ic[60] = { 0, 11, 5, 0, 5, 1, 0, 1, 7, 0, 7, 10, 0, 10, 11, 1, 5, 9, 5, 11, 4, 11, 10, 2, 10, 7, 6, 7, 1, 8, 3, 9, 4, 3, 4, 2,
		3, 2, 6, 3, 6, 8, 3, 8, 9, 4, 9, 5, 2, 4, 11, 6, 2, 10, 8, 6, 7, 9, 8, 1 };
	double p = (1 + sqrt(5.0)) / 2, x[3 * 42];
	double const v0[36] = { -1, p, 0, 1, p, 0, -1, -p, 0, 1, -p, 0, 0, -1, p, 0, 1, p, 0, -1, -p, 0, 1, -p, p, 0, -1, p, 0, 1, -p, 0, -1, -p, 0, 1 };
	unsigned int t[3 * 80], nv = 12, nt = 0, i, j, mid[30][3], nm = 0;
	qaws_scalar xs[3 * 42];
	qaws_mesh_3d_desc d;
	qaws_mesh_3d* m = NULL;
	for (i = 0; i < 36; i++) x[i] = v0[i];
	for (i = 0; i < 20; i++)
	{
		unsigned int c[3], e;
		for (e = 0; e < 3; e++)
		{
			unsigned int a = (unsigned int)ic[3 * i + e], b = (unsigned int)ic[3 * i + (e + 1) % 3], lo = a < b ? a : b, hi = a < b ? b : a, s;
			for (s = 0; s < nm && !(mid[s][0] == lo && mid[s][1] == hi); s++)
				;
			if (s == nm)
			{
				mid[nm][0] = lo; mid[nm][1] = hi; mid[nm][2] = nv;
				for (j = 0; j < 3; j++) x[3 * nv + j] = 0.5 * (x[3 * a + j] + x[3 * b + j]);
				nv++;
				nm++;
			}
			c[e] = mid[s][2];
		}
		t[3 * nt] = (unsigned int)ic[3 * i]; t[3 * nt + 1] = c[0]; t[3 * nt + 2] = c[2]; nt++;
		t[3 * nt] = (unsigned int)ic[3 * i + 1]; t[3 * nt + 1] = c[1]; t[3 * nt + 2] = c[0]; nt++;
		t[3 * nt] = (unsigned int)ic[3 * i + 2]; t[3 * nt + 1] = c[2]; t[3 * nt + 2] = c[1]; nt++;
		t[3 * nt] = c[0]; t[3 * nt + 1] = c[1]; t[3 * nt + 2] = c[2]; nt++;
	}
	for (i = 0; i < nv; i++)
	{
		double l = sqrt(x[3 * i] * x[3 * i] + x[3 * i + 1] * x[3 * i + 1] + x[3 * i + 2] * x[3 * i + 2]);
		for (j = 0; j < 3; j++) xs[3 * i + j] = (qaws_scalar)(r * x[3 * i + j] / l);
	}
	memset(&d, 0, sizeof d);
	d.positions = xs; d.vertex_count = nv; d.face_vertices = t; d.face_count = nt;
	qaws_mesh_3d_create(&m, &d);
	return m;
}

static qaws_motion* sw_line(double dx, double dy, double dz)
{
	qaws_motion_pose s, e;
	qaws_motion* m = NULL;
	memset(&s, 0, sizeof s);
	memset(&e, 0, sizeof e);
	s.rotation[0] = e.rotation[0] = 1;
	e.position.x = (qaws_scalar)dx; e.position.y = (qaws_scalar)dy; e.position.z = (qaws_scalar)dz;
	s.velocity = e.velocity = e.position;
	qaws_motion_create_poses(&m, &s, &e, 0);
	return m;
}

static void test_translation(void)
{
	qaws_mesh_3d *c = sw_box(0, 1, 0, 1, 0, 1), *o = NULL;
	qaws_motion* m = sw_line(1, 1, 0);
	qaws_sweep_report r;
	qaws_status st = qaws_mesh_3d_sweep(&o, &r, c, m, NULL);
	double v = sw_volume(o);
	char msg[300];
	sprintf(msg, "a unit cube moved by (1, 1, 0): volume %.15f (3), %u slab, %u front faces, closed %d", v, r.slabs, r.front_faces, r.closed);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && o && fabs(v - 3) < 1e-12 && r.slabs == 1 && r.closed, msg);
	qaws_mesh_3d_destroy(c);
	qaws_mesh_3d_destroy(o);
	qaws_motion_destroy(m);
}

static void test_turn(void)
{
	double h = 0.01, r0 = 1, r1 = sqrt(4 + h * h), want = 0.25 * SW_PI * (r1 * r1 - r0 * r0) + 2 * h, v;   /* (the sector, and the boxes at both ends half outside it) */
	qaws_mesh_3d *b = sw_box(1, 2, -h, h, 0, 1), *o = NULL;
	qaws_motion_pose s, e;
	qaws_motion* m = NULL;
	qaws_sweep_report r;
	qaws_status st;
	char msg[300];
	memset(&s, 0, sizeof s);
	memset(&e, 0, sizeof e);
	s.rotation[0] = 1;
	e.rotation[0] = (qaws_scalar)cos(SW_PI / 4);   /* (a quarter turn: the half angle) */
	e.rotation[3] = (qaws_scalar)sin(SW_PI / 4);
	s.angular_velocity.z = e.angular_velocity.z = (qaws_scalar)(SW_PI / 2);
	qaws_motion_create_poses(&m, &s, &e, 0);
	st = qaws_mesh_3d_sweep(&o, &r, b, m, NULL);
	v = sw_volume(o);
	sprintf(msg, "a thin box turned a quarter about z: volume %.5f (the annular sector %.5f), %u slabs, chords within %.1e, closed %d", v, want, r.slabs,
		(double)r.deviation, r.closed);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && o && fabs(v - want) < 0.01 * want && r.slabs > 4 && r.closed, msg);
	qaws_mesh_3d_destroy(b);
	qaws_mesh_3d_destroy(o);
	qaws_motion_destroy(m);
}

static void test_tube(void)
{
	double rad = 0.2, L = 0, vb, v, want;
	qaws_mesh_3d *b = sw_ball(rad), *o = NULL;
	qaws_motion_pose s, e;
	qaws_motion* m = NULL;
	qaws_sweep_report r;
	qaws_sweep_desc d;
	qaws_status st;
	unsigned int k;
	char msg[300];
	memset(&s, 0, sizeof s);
	memset(&e, 0, sizeof e);
	s.rotation[0] = e.rotation[0] = 1;
	e.position.x = 2;
	s.velocity.x = 1; s.velocity.y = 2;
	e.velocity.x = 1; e.velocity.y = -2;
	qaws_motion_create_poses(&m, &s, &e, 0);
	/* the path's length, by fine chords */
	{
		qaws_motion_state q0, q1;
		qaws_motion_evaluate(&q0, m, 0);
		for (k = 1; k <= 2000; k++)
		{
			qaws_motion_evaluate(&q1, m, (qaws_scalar)(k / 2000.0));
			L += sqrt(((double)q1.position.x - q0.position.x) * ((double)q1.position.x - q0.position.x) + ((double)q1.position.y - q0.position.y) * ((double)q1.position.y - q0.position.y));
			q0 = q1;
		}
	}
	memset(&d, 0, sizeof d);
	d.tolerance = (qaws_scalar)0.004;
	st = qaws_mesh_3d_sweep(&o, &r, b, m, &d);
	v = sw_volume(o);
	vb = sw_volume(b);
	/* (the polyhedral ball's cross-section scaled as its volume: the tube of a ball of its volume) */
	want = SW_PI * rad * rad * L * pow(vb / (4.0 / 3 * SW_PI * rad * rad * rad), 2.0 / 3) + vb;
	sprintf(msg, "a ball (80 faces, r 0.2) along a Hermite curve of length %.4f: volume %.4f (the tube %.4f), %u slabs, chords within %.1e, closed %d",
		L, v, want, r.slabs, (double)r.deviation, r.closed);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && o && fabs(v - want) < 0.04 * want && r.closed, msg);
	qaws_mesh_3d_destroy(b);
	qaws_mesh_3d_destroy(o);
	qaws_motion_destroy(m);
}

static void test_refusals(void)
{
	qaws_mesh_3d *c = sw_box(0, 1, 0, 1, 0, 1), *o = NULL;
	qaws_scalar p[9] = { 0 };
	qaws_motion* m = NULL;
	qaws_status st;
	char msg[200];
	qaws_motion_create_vertices(&m, p, p, p, p, NULL, NULL, 3, 0);
	st = qaws_mesh_3d_sweep(&o, NULL, c, m, NULL);
	sprintf(msg, "refused: a vertex motion of 3 vertices for an 8-vertex mesh (%d)", (int)st);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_INVALID_ARGUMENT && !o, msg);
	qaws_mesh_3d_destroy(c);
	qaws_motion_destroy(m);
}

int test_180_mesh_sweep_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 180: swept solids\n");
	test_translation();
	test_turn();
	test_tube();
	test_refusals();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
