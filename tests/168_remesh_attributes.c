/*
 * Test 168: remeshing from attributes (qaws_mesh_3d_cut_isolines, metric
 * intersection and gradation), small meshes
 *
 *   - a flat grid with a linear attribute: its isolines at two levels cut in
 *     as feature edges, their vertices exactly on the levels (value and
 *     position), all faces triangles, the mesh valid;
 *   - an icosphere with the attribute z, its isoline z = 0.3 cut in, then
 *     remeshed (its own reference: the line a crease of it): the line stays,
 *     its vertices still at f = 0.3;
 *   - metric intersection: an isotropic and an anisotropic field give the
 *     smaller size in each direction; a rotated one: the result's ball
 *     inside both; a field with itself: unchanged;
 *   - gradation: a size of 0.01 at one corner in 0.2 elsewhere grows from it
 *     no faster than the gradation allows, and stays 0.2 far away;
 *   - a cube's attributes remeshed (its own reference): per-side UV charts
 *     (a seam on every edge) exact at every corner, a value per triangle (a
 *     seam along every side's diagonal, no crease) with no face across a
 *     diagonal (each face's corners on its own triangle), the edge values'
 *     lengths kept, the features the cube's edges alone.
 */

#include "test_common.h"
#include "qaws_remesh.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

static qaws_mesh_3d* ra_grid(unsigned int n)
{
	unsigned int nv = (n + 1) * (n + 1), i, j, k = 0;
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

static qaws_mesh_3d* ra_icosphere(void)
{
	static int const ic[60] = { 0, 11, 5, 0, 5, 1, 0, 1, 7, 0, 7, 10, 0, 10, 11, 1, 5, 9, 5, 11, 4, 11, 10, 2, 10, 7, 6, 7, 1, 8, 3, 9, 4, 3, 4, 2,
		3, 2, 6, 3, 6, 8, 3, 8, 9, 4, 9, 5, 2, 4, 11, 6, 2, 10, 8, 6, 7, 9, 8, 1 };
	double p = (1 + sqrt(5.0)) / 2;
	double const v0[36] = { -1, p, 0, 1, p, 0, -1, -p, 0, 1, -p, 0, 0, -1, p, 0, 1, p, 0, -1, -p, 0, 1, -p, p, 0, -1, p, 0, 1, -p, 0, -1, -p, 0, 1 };
	unsigned int nt = 20, nv = 12, i, j, lv, nm;
	double* x = (double*)malloc(sizeof(double) * 3 * 700);
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 3 * 1280);
	unsigned int* t2 = (unsigned int*)malloc(sizeof(unsigned int) * 3 * 1280);
	unsigned int(*mid)[3] = (unsigned int(*)[3])malloc(sizeof(unsigned int) * 3 * 2000);
	qaws_scalar* xs;
	qaws_mesh_3d_desc d;
	qaws_mesh_3d* m = NULL;
	for (i = 0; i < 12; i++)
	{
		double l = sqrt(v0[3 * i] * v0[3 * i] + v0[3 * i + 1] * v0[3 * i + 1] + v0[3 * i + 2] * v0[3 * i + 2]);
		for (j = 0; j < 3; j++)
			x[3 * i + j] = v0[3 * i + j] / l;
	}
	for (i = 0; i < 60; i++)
		t[i] = (unsigned int)ic[i];
	for (lv = 0; lv < 3; lv++)
	{
		nm = 0;
		for (i = 0; i < nt; i++)
		{
			unsigned int c[3], e;
			for (e = 0; e < 3; e++)
			{
				unsigned int a = t[3 * i + e], b = t[3 * i + (e + 1) % 3], lo = a < b ? a : b, hi = a < b ? b : a, s;
				double l;
				for (s = 0; s < nm && !(mid[s][0] == lo && mid[s][1] == hi); s++)
					;
				if (s == nm)
				{
					mid[nm][0] = lo;
					mid[nm][1] = hi;
					mid[nm++][2] = nv;
					for (j = 0; j < 3; j++)
						x[3 * nv + j] = 0.5 * (x[3 * a + j] + x[3 * b + j]);
					l = sqrt(x[3 * nv] * x[3 * nv] + x[3 * nv + 1] * x[3 * nv + 1] + x[3 * nv + 2] * x[3 * nv + 2]);
					for (j = 0; j < 3; j++)
						x[3 * nv + j] /= l;
					nv++;
				}
				c[e] = mid[s][2];
			}
			{
				unsigned int const q[12] = { t[3 * i], c[0], c[2], c[0], t[3 * i + 1], c[1], c[2], c[1], t[3 * i + 2], c[0], c[1], c[2] };
				memcpy(t2 + 12 * i, q, sizeof q);
			}
		}
		nt *= 4;
		memcpy(t, t2, sizeof(unsigned int) * 3 * nt);
	}
	xs = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	for (i = 0; i < 3 * nv; i++)
		xs[i] = (qaws_scalar)x[i];
	memset(&d, 0, sizeof d);
	d.positions = xs;
	d.vertex_count = nv;
	d.face_vertices = t;
	d.face_count = nt;
	qaws_mesh_3d_create(&m, &d);
	free(x); free(t); free(t2); free(mid); free(xs);
	return m;
}

/* a vertex attribute of one component from f(x, y, z) = a x + b y + c z */
static unsigned int ra_linear_attribute(qaws_mesh_3d* m, double a, double b, double c)
{
	unsigned int nv = qaws_mesh_3d_get_vertex_count(m), i, id = 0;
	qaws_scalar* v = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nv);
	qaws_attribute_desc d;
	for (i = 0; i < nv; i++)
	{
		qaws_vec3 p = qaws_mesh_3d_get_vertex(m, i);
		v[i] = (qaws_scalar)(a * p.x + b * p.y + c * p.z);
	}
	memset(&d, 0, sizeof d);
	d.name = "f";
	d.element = QAWS_MESH_VERTEX;
	d.rule = QAWS_ATTRIBUTE_LINEAR;
	d.components = 1;
	qaws_mesh_3d_add_attribute(&id, m, v, &d);
	free(v);
	return id;
}

static void test_isolines_flat(void)
{
	qaws_mesh_3d* m = ra_grid(10);
	unsigned int a = ra_linear_attribute(m, 1, 0.3, 0), ne = 0, nv, nfe, i, ok = 1;
	qaws_scalar const levels[2] = { (qaws_scalar)0.25, (qaws_scalar)0.55 };
	qaws_scalar* val;
	unsigned int* fe;
	qaws_status st = qaws_mesh_3d_cut_isolines(&ne, m, a, levels, 2);
	char msg[200];
	nv = qaws_mesh_3d_get_vertex_count(m);
	nfe = qaws_mesh_3d_get_feature_edge_count(m);
	val = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nv);
	fe = (unsigned int*)malloc(sizeof(unsigned int) * 2 * (nfe + 1));
	qaws_mesh_3d_get_attribute(val, m, a);
	qaws_mesh_3d_get_feature_edges(fe, m);
	for (i = 0; i < 2 * nfe; i++)
	{
		qaws_vec3 p = qaws_mesh_3d_get_vertex(m, fe[i]);
		double f = val[fe[i]];
		ok &= (f == levels[0] || f == levels[1]) && fabs(p.x + 0.3 * p.y - f) < 1e-12;
	}
	sprintf(msg, "grid with f = x + 0.3 y, isolines 0.25 and 0.55: %u feature edges, on their levels %d, %u faces all triangles %d, valid %d", ne, ok,
		qaws_mesh_3d_get_face_count(m), qaws_mesh_3d_get_face_count(m) == qaws_mesh_3d_get_triangle_count(m),
		qaws_mesh_3d_validate(m) == QAWS_STATUS_OK);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && ne > 10 && nfe == ne && ok && qaws_mesh_3d_get_face_count(m) == qaws_mesh_3d_get_triangle_count(m)
		&& qaws_mesh_3d_validate(m) == QAWS_STATUS_OK, msg);
	TEST_ASSERT(qaws_mesh_3d_cut_isolines(NULL, m, 99, levels, 1) == QAWS_STATUS_INVALID_ARGUMENT, "an attribute that is not there refused");
	qaws_mesh_3d_destroy(m);
	free(val);
	free(fe);
}

static void test_isoline_remesh(void)
{
	qaws_mesh_3d* sphere = ra_icosphere();
	qaws_mesh_3d *cut = NULL, *out = NULL;
	qaws_mesh_reference* ref = NULL;
	qaws_remesh_desc d;
	qaws_scalar const level = (qaws_scalar)0.3;
	unsigned int a = ra_linear_attribute(sphere, 0, 0, 1), ne = 0, nfe = 0, i, ok = 1, oa = 0;
	qaws_status st;
	char msg[200];
	qaws_mesh_reference_create(&ref, sphere, NULL);
	qaws_mesh_3d_clone(&cut, sphere);
	st = qaws_mesh_3d_cut_isolines(&ne, cut, a, &level, 1);
	memset(&d, 0, sizeof d);
	d.target_length = (qaws_scalar)0.15;
	if (st == QAWS_STATUS_OK)
		st = qaws_mesh_remesh(&out, NULL, cut, &d);
	if (out)
	{
		unsigned int nv = qaws_mesh_3d_get_vertex_count(out);
		qaws_scalar* val = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nv);
		unsigned int* fe;
		nfe = qaws_mesh_3d_get_feature_edge_count(out);
		fe = (unsigned int*)malloc(sizeof(unsigned int) * 2 * (nfe + 1));
		qaws_mesh_3d_find_attribute(&oa, out, "f");
		qaws_mesh_3d_get_attribute(val, out, oa);
		qaws_mesh_3d_get_feature_edges(fe, out);
		for (i = 0; i < 2 * nfe; i++)
			ok &= fabs(val[fe[i]] - level) < 1e-9;
		free(val);
		free(fe);
	}
	sprintf(msg, "icosphere cut at z = 0.3 (%u edges), remeshed: %u feature edges, their vertices at f = 0.3 %d", ne, nfe, ok);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && out && nfe >= 10 && ok, msg);
	{
		/* against the uncut sphere's reference: its line not one of the reference's creases, kept across its smooth surface */
		qaws_mesh_3d* out2 = NULL;
		unsigned int nfe2 = 0, ok2 = 1, oa2 = 0;
		d.reference = ref;
		st = qaws_mesh_remesh(&out2, NULL, cut, &d);
		if (out2)
		{
			unsigned int nv = qaws_mesh_3d_get_vertex_count(out2);
			qaws_scalar* val = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nv);
			unsigned int* fe;
			nfe2 = qaws_mesh_3d_get_feature_edge_count(out2);
			fe = (unsigned int*)malloc(sizeof(unsigned int) * 2 * (nfe2 + 1));
			qaws_mesh_3d_find_attribute(&oa2, out2, "f");
			qaws_mesh_3d_get_attribute(val, out2, oa2);
			qaws_mesh_3d_get_feature_edges(fe, out2);
			/* (the cut polyline's chords projected onto the curved uncut surface: some 1e-4 off the level) */
			for (i = 0; i < 2 * nfe2; i++)
			{
				ok2 &= fabs(val[fe[i]] - level) < 5e-4;
			}
			free(val);
			free(fe);
		}
		sprintf(msg, "the same against the uncut sphere: %u faces (%u with its own), %u feature edges at f = 0.3 (5e-4) %d", out2 ? qaws_mesh_3d_get_face_count(out2) : 0,
			out ? qaws_mesh_3d_get_face_count(out) : 0, nfe2, ok2);
		printf("  %s\n", msg);
		TEST_ASSERT(st == QAWS_STATUS_OK && out2 && out && qaws_mesh_3d_get_face_count(out2) < 1.3 * qaws_mesh_3d_get_face_count(out) && nfe2 >= 10 && ok2, msg);
		qaws_mesh_3d_destroy(out2);
	}
	qaws_mesh_3d_destroy(out);
	qaws_mesh_3d_destroy(cut);
	qaws_mesh_reference_destroy(ref);
	qaws_mesh_3d_destroy(sphere);
}

/* a constant metric field on the reference: R diag(1/h^2) R^T with R a turn by th about z */
static qaws_metric_field* ra_metric(qaws_mesh_reference const* ref, double hx, double hy, double hz, double th)
{
	unsigned int nv = qaws_mesh_3d_get_vertex_count(qaws_mesh_reference_get_mesh(ref)), i;
	qaws_scalar* M = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 6 * nv);
	qaws_metric_field* f = NULL;
	double c = cos(th), s = sin(th), a = 1 / (hx * hx), b = 1 / (hy * hy);
	for (i = 0; i < nv; i++)
	{
		M[6 * i] = (qaws_scalar)(a * c * c + b * s * s);
		M[6 * i + 1] = (qaws_scalar)((a - b) * c * s);
		M[6 * i + 2] = 0;
		M[6 * i + 3] = (qaws_scalar)(a * s * s + b * c * c);
		M[6 * i + 4] = 0;
		M[6 * i + 5] = (qaws_scalar)(1 / (hz * hz));
	}
	qaws_metric_field_create_vertices(&f, ref, M);
	free(M);
	return f;
}

static double ra_q(double const* M, double const* u)
{
	return M[0] * u[0] * u[0] + 2 * M[1] * u[0] * u[1] + 2 * M[2] * u[0] * u[2] + M[3] * u[1] * u[1] + 2 * M[4] * u[1] * u[2] + M[5] * u[2] * u[2];
}

static void test_intersection(void)
{
	qaws_mesh_3d* g = ra_grid(4);
	qaws_mesh_reference* ref = NULL;
	qaws_metric_field *iso, *ani, *rot, *I = NULL, *J = NULL, *K = NULL;
	qaws_mesh_location l;
	qaws_scalar mi[6], ma[6], mr[6], mj[6], mk[6];
	double e = 0, inside = 1e300, tight = 0;
	unsigned int s, k;
	char msg[200];
	qaws_mesh_reference_create(&ref, g, NULL);
	iso = ra_metric(ref, 0.2, 0.2, 0.2, 0);
	ani = ra_metric(ref, 0.05, 0.4, 0.4, 0);
	rot = ra_metric(ref, 0.05, 0.4, 0.4, 0.5);
	{
		qaws_metric_field const* p[2] = { iso, ani };
		qaws_metric_field_create_intersection(&I, p, 2);
		p[1] = rot;
		qaws_metric_field_create_intersection(&J, p, 2);
		p[1] = iso;
		qaws_metric_field_create_intersection(&K, p, 2);
	}
	qaws_mesh_reference_get_vertex_location(&l, ref, 7);
	qaws_metric_field_evaluate(mi, I, &l, 1);
	qaws_metric_field_evaluate(mj, J, &l, 1);
	qaws_metric_field_evaluate(mk, K, &l, 1);
	qaws_metric_field_evaluate(ma, iso, &l, 1);
	qaws_metric_field_evaluate(mr, rot, &l, 1);
	{
		/* expected: diag(1 / 0.05^2, 1 / 0.2^2, 1 / 0.2^2) */
		double const want[6] = { 400, 0, 0, 25, 0, 25 };
		for (k = 0; k < 6; k++)
			e = fmax(e, fabs(mi[k] - want[k]) / 400 + fabs(mk[k] - ma[k]) / 400);
	}
	for (s = 0; s < 64; s++)
	{
		/* the turned case: in every direction at least as fine as both, and as fine as one of them in some */
		double th = 3.14159265358979323846 * s / 64, u[3] = { cos(th), sin(th), 0 }, dj[6], da[6], dr[6], q, qa, qr;
		for (k = 0; k < 6; k++)
		{
			dj[k] = mj[k];
			da[k] = ma[k];
			dr[k] = mr[k];
		}
		q = ra_q(dj, u);
		qa = ra_q(da, u);
		qr = ra_q(dr, u);
		inside = fmin(inside, q / fmax(qa, qr));
		tight = fmax(tight, fmax(qa, qr) / q);
	}
	sprintf(msg, "metric intersection: isotropic with anisotropic off by %.1e, with itself unchanged; turned: q / max(qa, qr) >= %.6f, tight to %.4f", e,
		inside, tight);
	printf("  %s\n", msg);
	TEST_ASSERT(e < 1e-9 && inside > 1 - 1e-9 && tight > 0.999, msg);
	qaws_metric_field_destroy(iso); qaws_metric_field_destroy(ani); qaws_metric_field_destroy(rot);
	qaws_metric_field_destroy(I); qaws_metric_field_destroy(J); qaws_metric_field_destroy(K);
	qaws_mesh_reference_destroy(ref);
	qaws_mesh_3d_destroy(g);
}

static void test_gradation(void)
{
	unsigned int n = 20, nv = (n + 1) * (n + 1), i;
	qaws_mesh_3d* g = ra_grid(n);
	qaws_mesh_reference_desc rd;
	qaws_mesh_reference* ref = NULL;
	qaws_scalar* M = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 6 * nv);
	qaws_metric_field *f = NULL, *gr = NULL;
	double worst = 0, far = 0, lb = log(1.3);
	char msg[200];
	qaws_status st;
	memset(&rd, 0, sizeof rd);
	rd.interpolant = QAWS_INTERPOLANT_LINEAR;
	qaws_mesh_reference_create(&ref, g, &rd);
	for (i = 0; i < nv; i++)
	{
		double h = i == 0 ? 0.01 : 0.2;
		M[6 * i] = M[6 * i + 3] = M[6 * i + 5] = (qaws_scalar)(1 / (h * h));
		M[6 * i + 1] = M[6 * i + 2] = M[6 * i + 4] = 0;
	}
	qaws_metric_field_create_vertices(&f, ref, M);
	st = qaws_metric_field_create_graded(&gr, f, (qaws_scalar)1.3);
	for (i = 0; i < nv && st == QAWS_STATUS_OK; i++)
	{
		qaws_mesh_location l;
		qaws_scalar m[6];
		qaws_vec3 p = qaws_mesh_3d_get_vertex(g, i);
		double h, d = sqrt((double)p.x * p.x + (double)p.y * p.y);
		qaws_mesh_reference_get_vertex_location(&l, ref, i);
		qaws_metric_field_evaluate(m, gr, &l, 1);
		h = 1 / sqrt((double)m[0]);
		/* (growth along the way out: h' = h ln(1.3) per unit of metric length, so dh/dd = ln(1.3), linear) */
		worst = fmax(worst, h / (0.01 + d * lb) - 1);
		if (d > 1)
			far = fmax(far, fabs(h - 0.2));
	}
	sprintf(msg, "gradation 1.3 from a size 0.01 at a corner in 0.2: at most %.3f over 0.01 + d ln 1.3, %.1e off 0.2 past distance 1", worst, far);
	printf("  %s\n", msg);
	/* (the growth runs along the mesh's edges, whose paths are up to some percent longer than the straight distance) */
	TEST_ASSERT(st == QAWS_STATUS_OK && worst < 0.1 && far < 1e-9, msg);
	{
		qaws_metric_field* bad = NULL;
		TEST_ASSERT(qaws_metric_field_create_graded(&bad, f, (qaws_scalar)1) == QAWS_STATUS_INVALID_ARGUMENT && !bad, "a gradation of 1 refused");
	}
	qaws_metric_field_destroy(f);
	qaws_metric_field_destroy(gr);
	qaws_mesh_reference_destroy(ref);
	qaws_mesh_3d_destroy(g);
	free(M);
}

static void test_seams(void)
{
	static unsigned int const F[36] = { 0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 2, 3, 7, 2, 7, 6, 1, 2, 6, 1, 6, 5, 0, 4, 7, 0, 7, 3 };
	qaws_scalar x[24], uv[72], tri[12], *e;
	unsigned int i, k, id, ne, nt = 0, nfe = 0, onfe = 0;
	qaws_mesh_3d_desc d;
	qaws_mesh_3d *m = NULL, *o = NULL;
	qaws_attribute_desc ad;
	qaws_remesh_desc rd;
	double wu = -1, lin[8] = { 0 }, lout[8] = { 0 }, wl = 0;
	unsigned int across = 1000;
	char msg[300];
	for (i = 0; i < 8; i++)
	{
		x[3 * i] = (qaws_scalar)((i == 1 || i == 2 || i == 5 || i == 6) ? 1 : 0);
		x[3 * i + 1] = (qaws_scalar)((i == 2 || i == 3 || i == 6 || i == 7) ? 1 : 0);
		x[3 * i + 2] = (qaws_scalar)(i >= 4 ? 1 : 0);
	}
	memset(&d, 0, sizeof d);
	d.positions = x; d.vertex_count = 8; d.face_vertices = F; d.face_count = 12;
	qaws_mesh_3d_create(&m, &d);
	qaws_mesh_3d_detect_features(m, (qaws_scalar)0.5, 0, 0);
	nfe = qaws_mesh_3d_get_feature_edge_count(m);
	for (i = 0; i < 12; i++)
	{
		unsigned int ax = i < 4 ? 2 : i < 8 ? 1 : 0;   /* (the sides: z-, z+, y-, y+, x+, x-) */
		tri[i] = (qaws_scalar)i;
		for (k = 0; k < 3; k++)
		{
			unsigned int v = F[3 * i + k];
			uv[6 * i + 2 * k] = (qaws_scalar)((ax == 0 ? x[3 * v + 1] : x[3 * v]) + 2 * (i / 2));
			uv[6 * i + 2 * k + 1] = (qaws_scalar)(ax == 2 ? x[3 * v + 1] : x[3 * v + 2]);
		}
	}
	memset(&ad, 0, sizeof ad);
	ad.name = "uv"; ad.element = QAWS_MESH_CORNER; ad.components = 2; ad.rule = QAWS_ATTRIBUTE_LINEAR;
	qaws_mesh_3d_add_attribute(&id, m, uv, &ad);
	ad.name = "tri"; ad.element = QAWS_MESH_FACE; ad.components = 1; ad.rule = QAWS_ATTRIBUTE_NEAREST;
	qaws_mesh_3d_add_attribute(&id, m, tri, &ad);
	ne = qaws_mesh_3d_get_edge_count(m);
	e = (qaws_scalar*)malloc(sizeof(qaws_scalar) * ne);
	{
		unsigned int* ev = (unsigned int*)malloc(sizeof(unsigned int) * 2 * ne);
		qaws_mesh_3d_get_edges(ev, m);
		for (i = 0; i < ne; i++)
		{
			qaws_vec3 a = qaws_mesh_3d_get_vertex(m, ev[2 * i]), b = qaws_mesh_3d_get_vertex(m, ev[2 * i + 1]);
			int diff = (a.x != b.x) + (a.y != b.y) + (a.z != b.z);
			e[i] = (qaws_scalar)(diff == 1 ? 1 + (ev[2 * i] + ev[2 * i + 1]) % 7 : 0);
			if (e[i] > 0)
				lin[(int)e[i]] += 1;
		}
		free(ev);
	}
	ad.name = "e"; ad.element = QAWS_MESH_EDGE; ad.components = 1; ad.rule = QAWS_ATTRIBUTE_LINEAR;
	qaws_mesh_3d_add_attribute(&id, m, e, &ad);
	memset(&rd, 0, sizeof rd);
	rd.target_length = (qaws_scalar)0.2;
	if (qaws_mesh_remesh(&o, NULL, m, &rd) == QAWS_STATUS_OK && o)
	{
		unsigned int c[3], iu, it, ie, ne2 = qaws_mesh_3d_get_edge_count(o), *ev2 = (unsigned int*)malloc(sizeof(unsigned int) * 2 * (ne2 + 1)), *fe;
		qaws_scalar *ouv, *ot, *oe = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (ne2 + 1));
		nt = qaws_mesh_3d_get_triangle_count(o);
		ouv = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 6 * nt);
		ot = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nt);
		qaws_mesh_3d_find_attribute(&iu, o, "uv"); qaws_mesh_3d_get_attribute(ouv, o, iu);
		qaws_mesh_3d_find_attribute(&it, o, "tri"); qaws_mesh_3d_get_attribute(ot, o, it);
		qaws_mesh_3d_find_attribute(&ie, o, "e"); qaws_mesh_3d_get_attribute(oe, o, ie);
		wu = 0;
		across = 0;
		for (i = 0; i < nt; i++)
		{
			unsigned int j = (unsigned int)lround((double)ot[i]), ax = j < 4 ? 2 : j < 8 ? 1 : 0, out = 0;
			qaws_mesh_3d_get_triangle(c, o, i);
			for (k = 0; k < 3; k++)
			{
				qaws_vec3 pp = qaws_mesh_3d_get_vertex(o, c[k]);
				double P[3] = { pp.x, pp.y, pp.z }, A[3], B[3], C[3], u[3], w[3], q[3], d00 = 0, d01 = 0, d11 = 0, d20 = 0, d21 = 0, den, bu, bv;
				unsigned int z;
				wu = fmax(wu, fabs((double)ouv[6 * i + 2 * k] - ((ax == 0 ? P[1] : P[0]) + 2 * (j / 2))) + fabs((double)ouv[6 * i + 2 * k + 1] - (ax == 2 ? P[1] : P[2])));
				/* (on triangle j: its barycentric weights not under -1e-9) */
				for (z = 0; z < 3; z++) { A[z] = x[3 * F[3 * j] + z]; B[z] = x[3 * F[3 * j + 1] + z]; C[z] = x[3 * F[3 * j + 2] + z]; }
				for (z = 0; z < 3; z++) { u[z] = B[z] - A[z]; w[z] = C[z] - A[z]; q[z] = P[z] - A[z]; d00 += u[z] * u[z]; d01 += u[z] * w[z]; d11 += w[z] * w[z]; d20 += q[z] * u[z]; d21 += q[z] * w[z]; }
				den = d00 * d11 - d01 * d01; bu = (d11 * d20 - d01 * d21) / den; bv = (d00 * d21 - d01 * d20) / den;
				out |= bu < -1e-9 || bv < -1e-9 || 1 - bu - bv < -1e-9;
			}
			across += out;
		}
		qaws_mesh_3d_get_edges(ev2, o);
		for (i = 0; i < ne2; i++)
		{
			qaws_vec3 a = qaws_mesh_3d_get_vertex(o, ev2[2 * i]), b = qaws_mesh_3d_get_vertex(o, ev2[2 * i + 1]);
			int v = (int)lround((double)oe[i]);
			if (v > 0 && v < 8)
				lout[v] += sqrt(((double)a.x - b.x) * ((double)a.x - b.x) + ((double)a.y - b.y) * ((double)a.y - b.y) + ((double)a.z - b.z) * ((double)a.z - b.z));
		}
		for (i = 1; i < 8; i++)
			wl = fmax(wl, fabs(lin[i] - lout[i]));
		/* the features on the cube's edges alone (two coordinates of their middles 0 or 1) */
		onfe = qaws_mesh_3d_get_feature_edge_count(o);
		fe = (unsigned int*)malloc(sizeof(unsigned int) * 2 * (onfe + 1));
		qaws_mesh_3d_get_feature_edges(fe, o);
		for (i = 0; i < onfe; i++)
		{
			qaws_vec3 a = qaws_mesh_3d_get_vertex(o, fe[2 * i]), b = qaws_mesh_3d_get_vertex(o, fe[2 * i + 1]);
			double mid[3] = { 0.5 * ((double)a.x + b.x), 0.5 * ((double)a.y + b.y), 0.5 * ((double)a.z + b.z) };
			unsigned int on = 0;
			for (k = 0; k < 3; k++)
				on += fabs(mid[k]) < 1e-9 || fabs(mid[k] - 1) < 1e-9;
			if (on < 2)
				onfe = 0;   /* (a feature off the cube's edges) */
		}
		free(fe); free(ev2); free(oe); free(ouv); free(ot);
	}
	sprintf(msg, "a cube's attributes remeshed: %u faces, UVs off their side's chart by %.1e, %u faces across a seam, edge values' lengths off by %.1e, %u features on the cube's edges (input %u)",
		nt, wu, across, wl, onfe, nfe);
	printf("  %s\n", msg);
	TEST_ASSERT(o && wu >= 0 && wu < 1e-9 && across == 0 && wl < 1e-9 && onfe >= nfe, msg);
	free(e);
	qaws_mesh_3d_destroy(m);
	qaws_mesh_3d_destroy(o);
}

int test_168_remesh_attributes_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 168: remeshing from attributes\n");
	test_isolines_flat();
	test_isoline_remesh();
	test_intersection();
	test_gradation();
	test_seams();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
