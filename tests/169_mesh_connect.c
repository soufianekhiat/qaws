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
 *   - refusals: a vertex off the border, one loop named twice;
 *   - two unit squares joined by a strip between their facing sides (chains):
 *     the area 2.5 exact, square quads; the same as planar meshes;
 *   - attributes and features carried: a vertex attribute 1 on one tube and
 *     3 on the other kept, between them on the bridge; the features kept;
 *   - volumes: two blocks of 2 x 2 hexahedra a unit apart joined by
 *     hexahedra (valid, cubes: scaled Jacobian 1); patches that differ
 *     refused; two prisms joined by prisms.
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

/* a unit square of n x n cells at x offset ox, counter-clockwise faces (z = 0) */
static void ct_square(unsigned int n, double ox, qaws_scalar* x, unsigned int* t)
{
	unsigned int i, j, k = 0;
	for (j = 0; j <= n; j++)
		for (i = 0; i <= n; i++)
		{
			x[3 * (j * (n + 1) + i)] = (qaws_scalar)(ox + (double)i / n);
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
}

static double ct_area(qaws_mesh_3d const* m)
{
	unsigned int nt = qaws_mesh_3d_get_triangle_count(m), i;
	double s = 0;
	for (i = 0; i < nt; i++)
	{
		unsigned int c[3];
		qaws_vec3 a, b, e;
		qaws_mesh_3d_get_triangle(c, m, i);
		a = qaws_mesh_3d_get_vertex(m, c[0]);
		b = qaws_mesh_3d_get_vertex(m, c[1]);
		e = qaws_mesh_3d_get_vertex(m, c[2]);
		s += 0.5 * (((double)b.x - a.x) * ((double)e.y - a.y) - ((double)b.y - a.y) * ((double)e.x - a.x));
	}
	return s;
}

static void test_strip(void)
{
	unsigned int n = 4, nv = (n + 1) * (n + 1);
	qaws_scalar xa[3 * 25], xb[3 * 25];
	unsigned int ta[6 * 16], tb[6 * 16];
	qaws_mesh_3d_desc d;
	qaws_mesh_3d *a = NULL, *b = NULL, *m = NULL;
	qaws_connect_desc cd;
	qaws_connect_report r;
	qaws_mesh_topology top = { 0, 0, 0, 0, 0 };
	qaws_status st;
	char msg[240];
	ct_square(n, 0, xa, ta);
	ct_square(n, 1.5, xb, tb);
	memset(&d, 0, sizeof d);
	d.positions = xa;
	d.vertex_count = nv;
	d.face_vertices = ta;
	d.face_count = 2 * n * n;
	qaws_mesh_3d_create(&a, &d);
	d.positions = xb;
	d.face_vertices = tb;
	qaws_mesh_3d_create(&b, &d);
	memset(&cd, 0, sizeof cd);
	cd.quads = 1;
	cd.flags = QAWS_CONNECT_CHAINS;
	/* a's right side up (along its border), b's left side up (against its border) */
	cd.a_end = n * (n + 1) + n;
	cd.b_end = n * (n + 1);
	st = qaws_mesh_3d_connect(&m, &r, a, n, b, 0, &cd);
	if (m)
		qaws_mesh_3d_get_topology(&top, m);
	sprintf(msg, "two squares joined by a strip between their facing sides: %u rows, %u faces, valid %d, pieces %u, loops %u, area %.12f",
		r.rows, r.bridge_faces, m && qaws_mesh_3d_validate(m) == QAWS_STATUS_OK, top.components, top.boundary_loops, m ? ct_area(m) : 0.0);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && m && qaws_mesh_3d_validate(m) == QAWS_STATUS_OK && top.components == 1 && top.boundary_loops == 1
		&& fabs(ct_area(m) - 2.5) < 1e-12 && fabs(r.min_scaled_jacobian - 1) < 1e-9, msg);
	qaws_mesh_3d_destroy(m);
	{
		/* the same as planar meshes */
		qaws_scalar pa[2 * 25], pb[2 * 25];
		qaws_mesh_2d_desc d2;
		qaws_mesh_2d *a2 = NULL, *b2 = NULL, *m2 = NULL;
		unsigned int i;
		for (i = 0; i < nv; i++)
		{
			pa[2 * i] = xa[3 * i];
			pa[2 * i + 1] = xa[3 * i + 1];
			pb[2 * i] = xb[3 * i];
			pb[2 * i + 1] = xb[3 * i + 1];
		}
		memset(&d2, 0, sizeof d2);
		d2.positions = pa;
		d2.vertex_count = nv;
		d2.face_vertices = ta;
		d2.face_count = 2 * n * n;
		qaws_mesh_2d_create(&a2, &d2);
		d2.positions = pb;
		d2.face_vertices = tb;
		qaws_mesh_2d_create(&b2, &d2);
		st = qaws_mesh_2d_connect(&m2, &r, a2, n, b2, 0, &cd);
		sprintf(msg, "the same as planar meshes: area %.12f", m2 ? ct_area(qaws_mesh_2d_get_mesh(m2)) : 0.0);
		printf("  %s\n", msg);
		TEST_ASSERT(st == QAWS_STATUS_OK && m2 && fabs(ct_area(qaws_mesh_2d_get_mesh(m2)) - 2.5) < 1e-12, msg);
		qaws_mesh_2d_destroy(m2);
		qaws_mesh_2d_destroy(a2);
		qaws_mesh_2d_destroy(b2);
	}
	qaws_mesh_3d_destroy(a);
	qaws_mesh_3d_destroy(b);
}

static void test_carried(void)
{
	qaws_mesh_3d *a = ct_tube(16, 2, 0.5, 0, 1), *b = ct_tube(16, 2, 0.5, 1.6, 2.6), *m = NULL;
	qaws_scalar va[48], vb[48], *vo;
	unsigned int ia = 0, ib = 0, io = 0, i, e[2 * 16], ok = 1, lo = 1, nv;
	qaws_attribute_desc ad;
	qaws_connect_desc cd;
	qaws_connect_report r;
	char msg[200];
	for (i = 0; i < 48; i++)
	{
		va[i] = 1;
		vb[i] = 3;
	}
	memset(&ad, 0, sizeof ad);
	ad.name = "heat";
	ad.element = QAWS_MESH_VERTEX;
	ad.components = 1;
	qaws_mesh_3d_add_attribute(&ia, a, va, &ad);
	qaws_mesh_3d_add_attribute(&ib, b, vb, &ad);
	for (i = 0; i < 16; i++)
	{
		e[2 * i] = i;
		e[2 * i + 1] = (i + 1) % 16;
	}
	qaws_mesh_3d_set_feature_edges(a, e, 16);
	memset(&cd, 0, sizeof cd);
	cd.quads = 1;
	qaws_mesh_3d_connect(&m, &r, a, 16 * 2, b, 0, &cd);
	if (m && qaws_mesh_3d_find_attribute(&io, m, "heat") == QAWS_STATUS_OK)
	{
		nv = qaws_mesh_3d_get_vertex_count(m);
		vo = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nv);
		qaws_mesh_3d_get_attribute(vo, m, io);
		for (i = 0; i < nv; i++)
		{
			if (i < 48)
				ok &= vo[i] == 1;
			else if (i < 96)
				ok &= vo[i] == 3;
			else
				lo &= vo[i] > 1 && vo[i] < 3;
		}
		free(vo);
	}
	else
		ok = 0;
	sprintf(msg, "tubes with a vertex attribute 1 and 3 and features: kept on each %d, the bridge between them %d, %u feature edges", ok, lo,
		m ? qaws_mesh_3d_get_feature_edge_count(m) : 0);
	printf("  %s\n", msg);
	TEST_ASSERT(m && ok && lo && qaws_mesh_3d_get_feature_edge_count(m) == 16, msg);
	qaws_mesh_3d_destroy(m);
	qaws_mesh_3d_destroy(a);
	qaws_mesh_3d_destroy(b);
}

/* an n x n x 1 block of unit hexahedra at height z0 */
static qaws_volume_mesh* ct_block(unsigned int n, double z0)
{
	unsigned int s = n + 1, i, j, k, c = 0;
	qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * s * s * 2);
	unsigned int* cv = (unsigned int*)malloc(sizeof(unsigned int) * 8 * n * n);
	qaws_cell_type* ty = (qaws_cell_type*)malloc(sizeof(qaws_cell_type) * n * n);
	qaws_volume_mesh_desc d;
	qaws_volume_mesh* m = NULL;
	for (k = 0; k < 2; k++)
		for (j = 0; j < s; j++)
			for (i = 0; i < s; i++)
			{
				unsigned int v = (k * s + j) * s + i;
				x[3 * v] = (qaws_scalar)i;
				x[3 * v + 1] = (qaws_scalar)j;
				x[3 * v + 2] = (qaws_scalar)(z0 + k);
			}
	for (j = 0; j < n; j++)
		for (i = 0; i < n; i++, c++)
		{
			unsigned int b0 = j * s + i, t0 = s * s + b0;
			unsigned int const q[8] = { b0, b0 + 1, b0 + s + 1, b0 + s, t0, t0 + 1, t0 + s + 1, t0 + s };
			memcpy(cv + 8 * c, q, sizeof q);
			ty[c] = QAWS_CELL_HEXAHEDRON;
		}
	memset(&d, 0, sizeof d);
	d.positions = x;
	d.vertex_count = 2 * s * s;
	d.cell_vertices = cv;
	d.cell_types = ty;
	d.cell_count = n * n;
	qaws_volume_mesh_create(&m, &d);
	free(x);
	free(cv);
	free(ty);
	return m;
}

static void test_volumes(void)
{
	qaws_volume_mesh *a = ct_block(2, 0), *b = ct_block(2, 2), *c3 = ct_block(3, 2), *m = NULL;
	unsigned int top[4], bot[4], bot3[9], i;
	qaws_connect_report r;
	qaws_status st;
	char msg[240];
	/* (a's top facets: each hexahedron's facet 1; b's bottom ones: facet 0; six facets per hexahedron) */
	for (i = 0; i < 4; i++)
	{
		top[i] = 6 * i + 1;
		bot[i] = 6 * i;
	}
	for (i = 0; i < 9; i++)
		bot3[i] = 6 * i;
	st = qaws_volume_mesh_connect(&m, &r, a, top, b, bot, 4, 4, NULL);
	sprintf(msg, "two blocks of 2 x 2 hexahedra a unit apart, joined: %u layers, %u cells, valid %d, gap %.3f, scaled Jacobian %.6f", r.rows,
		m ? qaws_volume_mesh_get_cell_count(m) : 0, m && qaws_volume_mesh_validate(m) == QAWS_STATUS_OK, (double)r.gap, (double)r.min_scaled_jacobian);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && m && qaws_volume_mesh_validate(m) == QAWS_STATUS_OK && qaws_volume_mesh_get_cell_count(m) == 8 + 4 * r.rows
		&& fabs(r.gap - 1) < 1e-12 && fabs(r.min_scaled_jacobian - 1) < 1e-9, msg);
	qaws_volume_mesh_destroy(m);
	m = NULL;
	TEST_ASSERT(qaws_volume_mesh_connect(&m, &r, a, top, c3, bot3, 4, 9, NULL) == QAWS_STATUS_INVALID_ARGUMENT && !m,
		"patches of 2 x 2 and 3 x 3 facets refused");
	{
		/* a prism under another two apart: triangles joined by prisms */
		qaws_scalar const px[36] = { 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 2, 1, 0, 2, 0, 1, 2, 0, 0, 3, 1, 0, 3, 0, 1, 3 };
		unsigned int const pv[12] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
		qaws_cell_type const pt[2] = { QAWS_CELL_PRISM, QAWS_CELL_PRISM };
		qaws_volume_mesh_desc d;
		qaws_volume_mesh* p = NULL;
		unsigned int f0 = 1, f1 = 5;   /* (the first prism's top facet; the second's bottom one: 5 facets per prism) */
		memset(&d, 0, sizeof d);
		d.positions = px;
		d.vertex_count = 12;
		d.cell_vertices = pv;
		d.cell_types = pt;
		d.cell_count = 2;
		qaws_volume_mesh_create(&p, &d);
		st = qaws_volume_mesh_connect(&m, &r, p, &f0, NULL, &f1, 1, 1, NULL);
		sprintf(msg, "two prisms of one mesh a unit apart, joined: %u cells, valid %d", m ? qaws_volume_mesh_get_cell_count(m) : 0,
			m && qaws_volume_mesh_validate(m) == QAWS_STATUS_OK);
		printf("  %s\n", msg);
		TEST_ASSERT(st == QAWS_STATUS_OK && m && qaws_volume_mesh_validate(m) == QAWS_STATUS_OK && qaws_volume_mesh_get_cell_count(m) == 2 + r.rows, msg);
		qaws_volume_mesh_destroy(m);
		qaws_volume_mesh_destroy(p);
	}
	qaws_volume_mesh_destroy(a);
	qaws_volume_mesh_destroy(b);
	qaws_volume_mesh_destroy(c3);
}

int test_169_mesh_connect_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 169: the mesh connector\n");
	test_tubes();
	test_counts();
	test_handle();
	test_strip();
	test_carried();
	test_volumes();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
