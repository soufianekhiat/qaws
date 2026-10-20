/*
 * Test 123: size fields (qaws_mesh_field.h)
 *
 *   - uniform and callback fields, and the smaller of two;
 *   - from curvature on a sphere of radius 2 with its normals: the chord
 *     tolerance's 2 sqrt(e (2R - e)) and the angle tolerance's 2R sin(t / 2),
 *     the smaller of the two, the clamps;
 *   - gradation: on a sphere with a sharp bump, sizes no faster growing than
 *     the gradation along every edge, and the bump the smallest;
 *   - refusals.
 */

#include "test_common.h"
#include "qaws_mesh_field.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define HM_NONE_TEST 0xFFFFFFFFu

/* an icosphere: the icosahedron split `level` times, its points on the unit sphere */
static qaws_mesh_3d* sf_icosphere(unsigned int level, double radius, qaws_scalar** out_normals)
{
	static unsigned int const f0[60] = { 0, 11, 5, 0, 5, 1, 0, 1, 7, 0, 7, 10, 0, 10, 11, 1, 5, 9, 5, 11, 4, 11, 10, 2, 10, 7, 6, 7, 1, 8,
		3, 9, 4, 3, 4, 2, 3, 2, 6, 3, 6, 8, 3, 8, 9, 4, 9, 5, 2, 4, 11, 6, 2, 10, 8, 6, 7, 9, 8, 1 };
	double g = (1 + sqrt(5.0)) / 2;
	double const v0[36] = { -1, g, 0, 1, g, 0, -1, -g, 0, 1, -g, 0, 0, -1, g, 0, 1, g, 0, -1, -g, 0, 1, -g, g, 0, -1, g, 0, 1, -g, 0, -1, -g, 0, 1 };
	unsigned int nv = 12, nf = 20, l, i, cap = 12;
	double* x;
	unsigned int* f;
	qaws_mesh_3d_desc d;
	qaws_mesh_3d* m = NULL;
	qaws_scalar* p;
	for (l = 0; l < level; l++)
		cap = cap * 4;
	x = (double*)malloc(sizeof(double) * 3 * (cap + 64));
	f = (unsigned int*)malloc(sizeof(unsigned int) * 3 * 20 * cap);
	memcpy(x, v0, sizeof v0);
	memcpy(f, f0, sizeof f0);
	for (l = 0; l < level; l++)
	{
		unsigned int* keys = (unsigned int*)malloc(sizeof(unsigned int) * 9 * nf), nk = 0, onf = nf, *old = (unsigned int*)malloc(sizeof(unsigned int) * 3 * nf);
		memcpy(old, f, sizeof(unsigned int) * 3 * nf);
		nf = 0;
		for (i = 0; i < onf; i++)
		{
			unsigned int e[3], j, q;
			for (j = 0; j < 3; j++)
			{
				unsigned int a = old[3 * i + j], b = old[3 * i + (j + 1) % 3], lo = a < b ? a : b, hi = a < b ? b : a;
				e[j] = HM_NONE_TEST;
				for (q = 0; q < nk; q++)
					if (keys[3 * q] == lo && keys[3 * q + 1] == hi)
						e[j] = keys[3 * q + 2];
				if (e[j] == HM_NONE_TEST)
				{
					unsigned int k;
					for (k = 0; k < 3; k++)
						x[3 * nv + k] = 0.5 * (x[3 * a + k] + x[3 * b + k]);
					keys[3 * nk] = lo; keys[3 * nk + 1] = hi; keys[3 * nk + 2] = nv;
					e[j] = nv++;
					nk++;
				}
			}
			{
				unsigned int a = old[3 * i], b = old[3 * i + 1], c = old[3 * i + 2];
				unsigned int t[12] = { a, e[0], e[2], e[0], b, e[1], e[2], e[1], c, e[0], e[1], e[2] };
				memcpy(f + 3 * nf, t, sizeof t);
				nf += 4;
			}
		}
		free(keys);
		free(old);
	}
	p = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	if (out_normals)
		*out_normals = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	for (i = 0; i < nv; i++)
	{
		double r = sqrt(x[3 * i] * x[3 * i] + x[3 * i + 1] * x[3 * i + 1] + x[3 * i + 2] * x[3 * i + 2]);
		unsigned int k;
		for (k = 0; k < 3; k++)
		{
			p[3 * i + k] = (qaws_scalar)(radius * x[3 * i + k] / r);
			if (out_normals)
				(*out_normals)[3 * i + k] = (qaws_scalar)(x[3 * i + k] / r);
		}
	}
	memset(&d, 0, sizeof d);
	d.positions = p;
	d.vertex_count = nv;
	d.face_vertices = f;
	d.face_count = nf;
	qaws_mesh_3d_create(&m, &d);
	free(x);
	free(f);
	free(p);
	return m;
}

static qaws_scalar sf_by_x(void* user, qaws_scalar const* p)
{
	(void)user;
	return (qaws_scalar)(0.1 + 0.05 * p[0]);
}

static void test_simple(void)
{
	qaws_size_field *u = NULL, *cb = NULL, *mn = NULL, *bad = NULL;
	qaws_scalar pts[6] = { 0, 0, 0, 2, 0, 0 }, out[2] = { 0, 0 };
	qaws_status s1;
	char msg[300];
	qaws_size_field_create_uniform(&u, (qaws_scalar)0.15);
	qaws_size_field_create_callback(&cb, sf_by_x, NULL);
	qaws_size_field_create_min(&mn, u, cb);
	s1 = qaws_size_field_create_uniform(&bad, 0);
	if (mn)
		qaws_size_field_evaluate(out, mn, NULL, pts, 2);
	sprintf(msg, "uniform 0.15 and 0.1 + 0.05 x: their min %.3g at x = 0, %.3g at x = 2; a zero length %s", (double)out[0], (double)out[1],
		qaws_status_to_string(s1));
	TEST_ASSERT(u && cb && mn && fabs(out[0] - 0.1) < 1e-6 && fabs(out[1] - 0.15) < 1e-6 && s1 == QAWS_STATUS_INVALID_ARGUMENT && !bad, msg);
	qaws_size_field_destroy(mn);
	qaws_size_field_destroy(cb);
	qaws_size_field_destroy(u);
}

static unsigned int sf_rand(unsigned int* s)
{
	*s = *s * 1664525u + 1013904223u;
	return *s >> 8;
}

static void test_curvature(void)
{
	qaws_scalar* nrm = NULL;
	qaws_mesh_3d* m = sf_icosphere(3, 2, &nrm);
	qaws_mesh_reference_desc rd;
	qaws_mesh_reference* r = NULL;
	qaws_curvature_size_desc cd;
	qaws_size_field *fc = NULL, *fa = NULL, *fb = NULL, *fl = NULL;
	unsigned int i, s = 3, nt;
	double R = 2, e = 0.01, t = 0.2, wc = 2 * sqrt(e * (2 * R - e)), wa = 2 * R * sin(t / 2), ec = 0, ea = 0, eb = 0, elo = 0;
	char msg[400];
	memset(&rd, 0, sizeof rd);
	rd.normals = nrm;
	qaws_mesh_reference_create(&r, m, &rd);
	if (!r)
	{
		TEST_ASSERT(0, "curvature: the reference");
		goto done;
	}
	memset(&cd, 0, sizeof cd);
	cd.chord_tolerance = (qaws_scalar)e;
	cd.max_length = 10;
	qaws_size_field_create_curvature(&fc, r, &cd);
	cd.chord_tolerance = 0;
	cd.angle_tolerance = (qaws_scalar)t;
	qaws_size_field_create_curvature(&fa, r, &cd);
	cd.chord_tolerance = (qaws_scalar)(e / 4);
	qaws_size_field_create_curvature(&fb, r, &cd);
	cd.min_length = (qaws_scalar)0.5;
	qaws_size_field_create_curvature(&fl, r, &cd);
	nt = qaws_mesh_reference_get_triangle_count(r);
	for (i = 0; fc && fa && fb && fl && i < 400; i++)
	{
		qaws_mesh_location l;
		qaws_scalar x[3], v[4];
		double a = (sf_rand(&s) & 0xFFFF) / 65536.0, b = (sf_rand(&s) & 0xFFFF) / 65536.0;
		if (a + b > 1) { a = 1 - a; b = 1 - b; }
		l.patch = sf_rand(&s) % nt; l.a = (qaws_scalar)a; l.b = (qaws_scalar)b;
		qaws_mesh_reference_evaluate(x, NULL, NULL, r, &l, 1);
		qaws_size_field_evaluate(v, fc, &l, x, 1);
		qaws_size_field_evaluate(v + 1, fa, &l, x, 1);
		qaws_size_field_evaluate(v + 2, fb, &l, x, 1);
		qaws_size_field_evaluate(v + 3, fl, &l, x, 1);
		ec = fmax(ec, fabs(v[0] - wc) / wc);
		ea = fmax(ea, fabs(v[1] - wa) / wa);
		eb = fmax(eb, fabs(v[2] - fmin(wa, 2 * sqrt(e / 4 * (2 * R - e / 4)))) / wa);
		elo = fmax(elo, fabs(v[3] - 0.5));
	}
	sprintf(msg, "sphere R 2: chord 0.01 gives %.4g (off %.2g), angle 0.2 gives %.4g (off %.2g), both the smaller (off %.2g), clamped to 0.5 (off %.2g)",
		wc, ec, wa, ea, eb, elo);
	TEST_ASSERT(fc && fa && fb && fl && ec < 0.06 && ea < 0.06 && eb < 0.06 && elo < 1e-9, msg);
done:
	qaws_size_field_destroy(fc);
	qaws_size_field_destroy(fa);
	qaws_size_field_destroy(fb);
	qaws_size_field_destroy(fl);
	qaws_mesh_reference_destroy(r);
	qaws_mesh_3d_destroy(m);
	free(nrm);
}

static void test_gradation(void)
{
	qaws_mesh_3d* m = sf_icosphere(3, 1, NULL);
	qaws_mesh_reference* r = NULL;
	qaws_curvature_size_desc cd;
	qaws_size_field *f = NULL, *g = NULL;
	unsigned int nv, ne, i, *e;
	double beta = 0.3, worst = 0, hmin = 1e300, hbump = 0, hfar = 0;
	qaws_scalar *x, *h;
	char msg[300];
	nv = qaws_mesh_3d_get_vertex_count(m);
	x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	h = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nv);
	qaws_mesh_3d_get_positions(x, m);
	for (i = 0; i < nv; i++)
	{
		/* a sharp bump round (0, 0, 1) */
		double d = (x[3 * i]) * (x[3 * i]) + (x[3 * i + 1]) * (x[3 * i + 1]) + (x[3 * i + 2] - 1) * (x[3 * i + 2] - 1), s = 1 + 0.5 * exp(-25 * d);
		x[3 * i] = (qaws_scalar)(x[3 * i] * s);
		x[3 * i + 1] = (qaws_scalar)(x[3 * i + 1] * s);
		x[3 * i + 2] = (qaws_scalar)(x[3 * i + 2] * s);
	}
	qaws_mesh_3d_set_positions(m, x);
	qaws_mesh_reference_create(&r, m, NULL);
	memset(&cd, 0, sizeof cd);
	cd.angle_tolerance = (qaws_scalar)0.15;
	cd.min_length = (qaws_scalar)0.005;
	cd.max_length = (qaws_scalar)1;
	cd.gradation = (qaws_scalar)beta;
	if (r)
		qaws_size_field_create_curvature(&f, r, &cd);
	cd.gradation = 0;
	if (r)
		qaws_size_field_create_curvature(&g, r, &cd);
	if (!f || !g)
	{
		TEST_ASSERT(0, "gradation: the fields");
		goto done;
	}
	for (i = 0; i < nv; i++)
	{
		qaws_mesh_location l;
		qaws_mesh_reference_get_vertex_location(&l, r, i);
		qaws_size_field_evaluate(h + i, f, &l, x + 3 * i, 1);
		if (h[i] < hmin)
			hmin = h[i];
		if (x[3 * i + 2] > 1.3)
			hbump = fmax(hbump, h[i]);
		if (x[3 * i + 2] < -0.5)
			hfar = fmax(hfar, h[i]);
	}
	ne = qaws_mesh_3d_get_edge_count(m);
	e = (unsigned int*)malloc(sizeof(unsigned int) * 2 * ne);
	qaws_mesh_3d_get_edges(e, m);
	for (i = 0; i < ne; i++)
	{
		unsigned int a = e[2 * i], b = e[2 * i + 1];
		double d = sqrt((x[3 * a] - x[3 * b]) * (x[3 * a] - x[3 * b]) + (x[3 * a + 1] - x[3 * b + 1]) * (x[3 * a + 1] - x[3 * b + 1])
			+ (x[3 * a + 2] - x[3 * b + 2]) * (x[3 * a + 2] - x[3 * b + 2]));
		worst = fmax(worst, (fabs(h[a] - h[b]) - beta * d) / d);
	}
	free(e);
	sprintf(msg, "gradation 0.3 on a bumped sphere: growth past it along the edges %.2g; sizes %.3g on the bump, %.3g far from it, smallest %.3g",
		worst, hbump, hfar, hmin);
	TEST_ASSERT(worst < TEST_BOUND(1e-9, 32) && hbump < 0.5 * hfar && hmin < hbump + TEST_BOUND(1e-12, 32), msg);
done:
	qaws_size_field_destroy(f);
	qaws_size_field_destroy(g);
	qaws_mesh_reference_destroy(r);
	qaws_mesh_3d_destroy(m);
	free(x);
	free(h);
}

static void test_refusals(void)
{
	qaws_mesh_3d *a = sf_icosphere(1, 1, NULL), *b = sf_icosphere(1, 1, NULL);
	qaws_mesh_reference *ra = NULL, *rb = NULL;
	qaws_curvature_size_desc cd;
	qaws_size_field *fa = NULL, *fb = NULL, *mn = NULL, *none = NULL;
	qaws_status s1, s2, s3;
	qaws_scalar p[3] = { 0, 0, 1 }, v;
	char msg[300];
	qaws_mesh_reference_create(&ra, a, NULL);
	qaws_mesh_reference_create(&rb, b, NULL);
	memset(&cd, 0, sizeof cd);
	s1 = qaws_size_field_create_curvature(&none, ra, &cd);
	cd.chord_tolerance = (qaws_scalar)0.01;
	qaws_size_field_create_curvature(&fa, ra, &cd);
	qaws_size_field_create_curvature(&fb, rb, &cd);
	s2 = qaws_size_field_create_min(&mn, fa, fb);
	s3 = fa ? qaws_size_field_evaluate(&v, fa, NULL, p, 1) : QAWS_STATUS_OK;
	sprintf(msg, "refusals: no tolerance %s, min of fields on two references %s, a reference's field without locations %s",
		qaws_status_to_string(s1), qaws_status_to_string(s2), qaws_status_to_string(s3));
	TEST_ASSERT(s1 == QAWS_STATUS_INVALID_ARGUMENT && !none && s2 == QAWS_STATUS_INVALID_ARGUMENT && !mn && s3 == QAWS_STATUS_INVALID_ARGUMENT, msg);
	qaws_size_field_destroy(fa);
	qaws_size_field_destroy(fb);
	qaws_mesh_reference_destroy(ra);
	qaws_mesh_reference_destroy(rb);
	qaws_mesh_3d_destroy(a);
	qaws_mesh_3d_destroy(b);
}

int test_123_size_field_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 123: size fields\n");
	test_simple();
	test_curvature();
	test_gradation();
	test_refusals();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
