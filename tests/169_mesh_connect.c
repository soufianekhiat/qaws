/*
 * Test 169: the mesh connector (qaws_mesh_3d_connect), small meshes
 *
 *   - two coaxial tubes facing each other across a gap, joined in quads:
 *     one valid mesh in one piece with the two far loops left, the bridge's
 *     vertices on the tubes' cylinder (the Hermite blend leaves both along
 *     the axis), no twist;
 *   - tubes of 24 and 16 sides: joined in triangles, valid, one piece, its
 *     smallest angle reasonable; smoothed: still valid;
 *   - a sphere with its two caps cut off, its two holes joined through it: a
 *     closed torus (genus 1);
 *   - refusals: a vertex off the border, one loop named twice.
 */

#include "test_common.h"
#include "qaws_mesh_connect.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define CT_PI 3.14159265358979323846

/* an open tube of radius r about z, from z0 to z1, nu sides and nv rows, faces outward */
static qaws_mesh_3d* ct_tube(unsigned int nu, unsigned int nv, double r, double z0, double z1)
{
	unsigned int n = nu * (nv + 1), i, j, k = 0;
	qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * n);
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 6 * nu * nv);
	qaws_mesh_3d_desc d;
	qaws_mesh_3d* m = NULL;
	for (j = 0; j <= nv; j++)
		for (i = 0; i < nu; i++)
		{
			double a = 2 * CT_PI * i / nu;
			x[3 * (j * nu + i)] = (qaws_scalar)(r * cos(a));
			x[3 * (j * nu + i) + 1] = (qaws_scalar)(r * sin(a));
			x[3 * (j * nu + i) + 2] = (qaws_scalar)(z0 + (z1 - z0) * j / nv);
		}
	for (j = 0; j < nv; j++)
		for (i = 0; i < nu; i++)
		{
			unsigned int a = j * nu + i, b = j * nu + (i + 1) % nu;
			t[k++] = a; t[k++] = b; t[k++] = b + nu;
			t[k++] = a; t[k++] = b + nu; t[k++] = a + nu;
		}
	memset(&d, 0, sizeof d);
	d.positions = x;
	d.vertex_count = n;
	d.face_vertices = t;
	d.face_count = 2 * nu * nv;
	qaws_mesh_3d_create(&m, &d);
	free(x);
	free(t);
	return m;
}

static void test_tubes(void)
{
	qaws_mesh_3d *a = ct_tube(24, 4, 0.5, 0, 1), *b = ct_tube(24, 4, 0.5, 1.6, 2.6), *m = NULL;
	qaws_connect_desc d;
	qaws_connect_report r;
	qaws_mesh_topology top = { 0, 0, 0, 0, 0 };
	qaws_status st;
	double dev = 0;
	unsigned int i, nv0 = 2 * 24 * 5;
	char msg[240];
	memset(&d, 0, sizeof d);
	d.quads = 1;
	/* (a's top loop, b's bottom one) */
	st = qaws_mesh_3d_connect(&m, &r, a, 24 * 4, b, 0, &d);
	if (m)
	{
		qaws_mesh_3d_get_topology(&top, m);
		for (i = nv0; i < qaws_mesh_3d_get_vertex_count(m); i++)
		{
			qaws_vec3 p = qaws_mesh_3d_get_vertex(m, i);
			dev = fmax(dev, fabs(sqrt((double)p.x * p.x + (double)p.y * p.y) - 0.5));
		}
	}
	sprintf(msg, "two tubes across a gap of %.2f, in quads: %u rows, %u faces, valid %d, pieces %u, loops %u, bridge off the cylinder by %.1e",
		(double)r.gap, r.rows, r.bridge_faces, m && qaws_mesh_3d_validate(m) == QAWS_STATUS_OK, m ? top.components : 0, m ? top.boundary_loops : 0, dev);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && m && qaws_mesh_3d_validate(m) == QAWS_STATUS_OK && top.components == 1 && top.boundary_loops == 2
		&& r.bridge_faces == 24 * r.rows && dev < 1e-9 && fabs(r.gap - 0.6) < 1e-6, msg);
	qaws_mesh_3d_destroy(m);
	qaws_mesh_3d_destroy(a);
	qaws_mesh_3d_destroy(b);
}

static void test_counts(void)
{
	unsigned int k;
	char msg[240];
	for (k = 0; k < 2; k++)
	{
		qaws_mesh_3d *a = ct_tube(24, 3, 0.5, 0, 1), *b = ct_tube(16, 3, 0.4, 1.5, 2.5), *m = NULL;
		qaws_connect_desc d;
		qaws_connect_report r;
		qaws_mesh_topology top = { 0, 0, 0, 0, 0 };
		qaws_status st;
		memset(&d, 0, sizeof d);
		d.quads = 1;
		d.flags = k ? QAWS_CONNECT_FAIR : 0;
		st = qaws_mesh_3d_connect(&m, &r, a, 24 * 3, b, 0, &d);
		if (m)
			qaws_mesh_3d_get_topology(&top, m);
		sprintf(msg, "tubes of 24 and 16 sides%s: %u rows, %u bridge faces, valid %d, pieces %u, loops %u, smallest angle %.1f degrees",
			k ? ", smoothed" : "", r.rows, r.bridge_faces, m && qaws_mesh_3d_validate(m) == QAWS_STATUS_OK, m ? top.components : 0,
			m ? top.boundary_loops : 0, (double)r.min_angle * 180 / CT_PI);
		printf("  %s\n", msg);
		TEST_ASSERT(st == QAWS_STATUS_OK && m && qaws_mesh_3d_validate(m) == QAWS_STATUS_OK && top.components == 1 && top.boundary_loops == 2
			&& r.min_angle > 15 * CT_PI / 180, msg);
		qaws_mesh_3d_destroy(m);
		qaws_mesh_3d_destroy(a);
		qaws_mesh_3d_destroy(b);
	}
}

static void test_handle(void)
{
	/* a sphere's band between latitudes: rings 1 to nb - 1 of a UV sphere, its two caps cut off */
	unsigned int na = 24, nb = 12, rings = nb - 1, n = na * rings, i, j, k = 0;
	qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * n);
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 6 * na * (rings - 1));
	qaws_mesh_3d_desc d;
	qaws_mesh_3d *s = NULL, *m = NULL;
	qaws_connect_report r;
	qaws_mesh_topology top = { 0, 0, 0, 0, 0 };
	qaws_status st;
	char msg[240];
	for (j = 0; j < rings; j++)
		for (i = 0; i < na; i++)
		{
			double th = CT_PI * (j + 1) / nb, ph = 2 * CT_PI * i / na;
			x[3 * (j * na + i)] = (qaws_scalar)(sin(th) * cos(ph));
			x[3 * (j * na + i) + 1] = (qaws_scalar)(sin(th) * sin(ph));
			x[3 * (j * na + i) + 2] = (qaws_scalar)cos(th);
		}
	for (j = 0; j + 1 < rings; j++)
		for (i = 0; i < na; i++)
		{
			unsigned int a = j * na + i, b = j * na + (i + 1) % na;
			/* (outward: down the rings is up the latitude) */
			t[k++] = a; t[k++] = a + na; t[k++] = b + na;
			t[k++] = a; t[k++] = b + na; t[k++] = b;
		}
	memset(&d, 0, sizeof d);
	d.positions = x;
	d.vertex_count = n;
	d.face_vertices = t;
	d.face_count = k / 3;
	qaws_mesh_3d_create(&s, &d);
	st = qaws_mesh_3d_connect(&m, &r, s, 0, NULL, (rings - 1) * na, NULL);
	if (m)
		qaws_mesh_3d_get_topology(&top, m);
	sprintf(msg, "a sphere without its caps, its two holes joined through it: %u rows, valid %d, closed %d, genus %u", r.rows,
		m && qaws_mesh_3d_validate(m) == QAWS_STATUS_OK, m ? top.closed : 0, m ? top.genus : 0);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && m && qaws_mesh_3d_validate(m) == QAWS_STATUS_OK && top.closed && top.genus == 1, msg);
	qaws_mesh_3d_destroy(m);
	m = NULL;
	TEST_ASSERT(qaws_mesh_3d_connect(&m, &r, s, na + 3, NULL, 0, NULL) == QAWS_STATUS_INVALID_ARGUMENT, "a vertex off the border refused");
	TEST_ASSERT(qaws_mesh_3d_connect(&m, &r, s, 0, NULL, 5, NULL) == QAWS_STATUS_INVALID_ARGUMENT && !m, "one loop named twice refused");
	qaws_mesh_3d_destroy(s);
	free(x);
	free(t);
}

int test_169_mesh_connect_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 169: the mesh connector\n");
	test_tubes();
	test_counts();
	test_handle();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
