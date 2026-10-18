/*
 * Test 132: CVT remeshing (qaws_mesh_remesh, QAWS_REMESH_CVT), against
 * LOCAL on the same inputs
 *
 *   - an icosphere with its normals at 0.12: edges on target, mostly of
 *     valence 6 and none outside 5..7, angles in 30..100, edge lengths less
 *     spread than LOCAL's, the vertices on the reference;
 *   - the showcase blob at 0.11: few irregular vertices, angles above 30;
 *   - a cube with its creases: corners kept, vertices on the cube, the
 *     crease edges along its edges (length 24), a linear attribute exact;
 *   - adaptive: test 124's curvature field on a bumped sphere, edges near
 *     the size (closer than LOCAL's), short on the bump;
 *   - convergence: the energy no higher after L-BFGS than after Lloyd,
 *     edges no more spread after 30 steps than after 1;
 *   - a closed slab thinner than the target: a valid closed genus 0 mesh;
 *     an open bumpy square: its border kept, its corners, angles above 25;
 *   - determinism: two runs, the same arrays;
 *   - about 2.4k seeds: the time printed.
 */

#include "test_common.h"
#include "qaws_remesh.h"
#include "mesh/qaws_remesh_cvt.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

static qaws_mesh_3d* ct_make(double const* x, unsigned int nv, unsigned int const* f, unsigned int const* fs, unsigned int nf)
{
	qaws_mesh_3d_desc d;
	qaws_mesh_3d* m = NULL;
	qaws_scalar* p = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
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

/* the icosahedron split once on the unit sphere (normals: radial, may be NULL); blob: pushed in and out as the showcase's */
static qaws_mesh_3d* ct_icosphere(qaws_scalar* normals, int blob)
{
	static unsigned int const f0[60] = { 0, 11, 5, 0, 5, 1, 0, 1, 7, 0, 7, 10, 0, 10, 11, 1, 5, 9, 5, 11, 4, 11, 10, 2, 10, 7, 6, 7, 1, 8,
		3, 9, 4, 3, 4, 2, 3, 2, 6, 3, 6, 8, 3, 8, 9, 4, 9, 5, 2, 4, 11, 6, 2, 10, 8, 6, 7, 9, 8, 1 };
	double g = (1 + sqrt(5.0)) / 2, x[3 * 42];
	double const v0[36] = { -1, g, 0, 1, g, 0, -1, -g, 0, 1, -g, 0, 0, -1, g, 0, 1, g, 0, -1, -g, 0, 1, -g, g, 0, -1, g, 0, 1, -g, 0, -1, -g, 0, 1 };
	unsigned int f[240], keys[3 * 30], nk = 0, nv = 12, nf = 0, i, j, q, k;
	memcpy(x, v0, sizeof v0);
	for (i = 0; i < 20; i++)
	{
		unsigned int e[3];
		for (j = 0; j < 3; j++)
		{
			unsigned int a = f0[3 * i + j], b = f0[3 * i + (j + 1) % 3], lo = a < b ? a : b, hi = a < b ? b : a;
			e[j] = 0xFFFFFFFFu;
			for (q = 0; q < nk; q++)
				if (keys[3 * q] == lo && keys[3 * q + 1] == hi)
					e[j] = keys[3 * q + 2];
			if (e[j] == 0xFFFFFFFFu)
			{
				for (k = 0; k < 3; k++)
					x[3 * nv + k] = 0.5 * (x[3 * a + k] + x[3 * b + k]);
				keys[3 * nk] = lo; keys[3 * nk + 1] = hi; keys[3 * nk + 2] = nv;
				e[j] = nv++;
				nk++;
			}
		}
		{
			unsigned int a = f0[3 * i], b = f0[3 * i + 1], c = f0[3 * i + 2];
			unsigned int t[12] = { a, e[0], e[2], e[0], b, e[1], e[2], e[1], c, e[0], e[1], e[2] };
			memcpy(f + 3 * nf, t, sizeof t);
			nf += 4;
		}
	}
	for (i = 0; i < nv; i++)
	{
		double r = sqrt(x[3 * i] * x[3 * i] + x[3 * i + 1] * x[3 * i + 1] + x[3 * i + 2] * x[3 * i + 2]), s;
		for (k = 0; k < 3; k++)
		{
			x[3 * i + k] /= r;
			if (normals)
				normals[3 * i + k] = (qaws_scalar)x[3 * i + k];
		}
		if (!blob)
			continue;
		s = 1 + 0.32 * sin(2.2 * x[3 * i] + 0.4) * cos(1.7 * x[3 * i + 1]) + 0.2 * x[3 * i + 2] * x[3 * i + 2];
		x[3 * i] *= 1.25 * s;
		x[3 * i + 1] *= s;
		x[3 * i + 2] *= 0.9 * s;
	}
	return ct_make(x, nv, f, NULL, nf);
}

static qaws_mesh_3d* ct_remesh(qaws_mesh_3d const* m, qaws_mesh_reference const* r, qaws_remesh_method method, double L, qaws_size_field const* sf,
	unsigned int iterations, qaws_mesh_report* rep, qaws_status* st, double* seconds)
{
	qaws_remesh_desc d;
	qaws_mesh_3d* out = NULL;
	clock_t t0 = clock();
	memset(&d, 0, sizeof d);
	d.method = method;
	d.reference = r;
	d.target_length = (qaws_scalar)L;
	d.size = sf;
	d.iterations = iterations;
	*st = qaws_mesh_remesh(m, &d, &out, rep);
	if (seconds)
		*seconds = (double)(clock() - t0) / CLOCKS_PER_SEC;
	return out;
}

/* what a result's shape is: interior vertices by valence, edge lengths' mean and coefficient of variation */
typedef struct ct_shape
{
	unsigned int nv, interior, v6, v567;
	double mean, cv;
} ct_shape;

static void ct_measure(qaws_mesh_3d const* m, ct_shape* s)
{
	unsigned int ne = qaws_mesh_3d_get_edge_count(m), i, *e = (unsigned int*)malloc(sizeof(unsigned int) * 2 * (ne + 1));
	qaws_scalar* p;
	double sum = 0, sq = 0;
	memset(s, 0, sizeof *s);
	s->nv = qaws_mesh_3d_get_vertex_count(m);
	p = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * (s->nv + 1));
	qaws_mesh_3d_get_positions(m, p);
	qaws_mesh_3d_get_edges(m, e);
	for (i = 0; i < ne; i++)
	{
		qaws_scalar const *a = p + 3 * e[2 * i], *b = p + 3 * e[2 * i + 1];
		double l = sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]));
		sum += l;
		sq += l * l;
	}
	s->mean = ne ? sum / ne : 0;
	s->cv = ne && s->mean > 0 ? sqrt(fmax(0, sq / ne - s->mean * s->mean)) / s->mean : 0;
	for (i = 0; i < s->nv; i++)
	{
		unsigned int ring[64], n = 0;
		if (qaws_mesh_3d_is_boundary_vertex(m, i))
			continue;
		qaws_mesh_3d_get_vertex_ring(m, i, ring, 64, &n);
		s->interior++;
		s->v6 += n == 6;
		s->v567 += n >= 5 && n <= 7;
	}
	free(p);
	free(e);
}

/* the largest gap between the vertices and the reference at their footpoints, and off the unit sphere */
static void ct_feet(qaws_mesh_3d const* m, qaws_mesh_reference const* r, double* foot, double* off)
{
	unsigned int nv = qaws_mesh_3d_get_vertex_count(m), i;
	qaws_mesh_location* l = (qaws_mesh_location*)malloc(sizeof(qaws_mesh_location) * (nv + 1));
	qaws_scalar *x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * (nv + 1)), *y = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * (nv + 1));
	*foot = *off = 0;
	qaws_mesh_3d_get_positions(m, x);
	qaws_mesh_3d_get_footpoints(m, l);
	qaws_mesh_reference_evaluate(r, l, nv, y, NULL, NULL);
	for (i = 0; i < nv; i++)
	{
		*foot = fmax(*foot, fabs(x[3 * i] - y[3 * i]) + fabs(x[3 * i + 1] - y[3 * i + 1]) + fabs(x[3 * i + 2] - y[3 * i + 2]));
		*off = fmax(*off, fabs(sqrt(x[3 * i] * x[3 * i] + x[3 * i + 1] * x[3 * i + 1] + x[3 * i + 2] * x[3 * i + 2]) - 1));
	}
	free(l); free(x); free(y);
}

static void test_sphere(void)
{
	qaws_scalar nrm[3 * 42];
	qaws_mesh_3d *m = ct_icosphere(nrm, 0), *out = NULL, *loc = NULL;
	qaws_mesh_reference_desc rd;
	qaws_mesh_reference* r = NULL;
	qaws_mesh_report rep, lrep;
	qaws_mesh_topology t;
	qaws_status st = QAWS_STATUS_INTERNAL_ERROR, ls = QAWS_STATUS_INTERNAL_ERROR;
	ct_shape s, ls_;
	double foot = 0, off = 0, sec = 0;
	char msg[600];
	memset(&rd, 0, sizeof rd);
	rd.normals = nrm;
	qaws_mesh_reference_create(m, &rd, &r);
	if (r)
	{
		out = ct_remesh(m, r, QAWS_REMESH_CVT, 0.12, NULL, 0, &rep, &st, &sec);
		loc = ct_remesh(m, r, QAWS_REMESH_LOCAL, 0.12, NULL, 0, &lrep, &ls, NULL);
	}
	if (!out || !loc)
	{
		TEST_ASSERT(0, "sphere: remeshed");
		goto done;
	}
	ct_measure(out, &s);
	ct_measure(loc, &ls_);
	ct_feet(out, r, &foot, &off);
	qaws_mesh_3d_get_topology(out, &t);
	sprintf(msg, "sphere by CVT at 0.12: %u vertices in %.2f s, valid %d, closed %d, genus %u; mean edge %.4f (LOCAL %.4f), valence 6 at %.1f%% (LOCAL %.1f%%), 5..7 at %.1f%%, angles %.1f..%.1f (LOCAL %.1f..%.1f), edge CV %.3f (LOCAL %.3f), deviation %.2g, footpoints %.2g, off the sphere %.2g",
		s.nv, sec, qaws_mesh_3d_validate(out) == QAWS_STATUS_OK, t.closed, t.genus, s.mean, ls_.mean, 100.0 * s.v6 / s.interior, 100.0 * ls_.v6 / ls_.interior,
		100.0 * s.v567 / s.interior, (double)rep.quality.min_angle, (double)rep.quality.max_angle, (double)lrep.quality.min_angle, (double)lrep.quality.max_angle,
		s.cv, ls_.cv, (double)rep.max_deviation, foot, off);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && ls == QAWS_STATUS_OK && qaws_mesh_3d_validate(out) == QAWS_STATUS_OK && t.closed && t.genus == 0 && t.components == 1
		&& fabs(s.mean / 0.12 - 1) <= 0.03 && s.v6 >= 0.92 * s.interior && s.v567 == s.interior && rep.quality.min_angle >= 30
		&& rep.quality.max_angle <= 100 && s.cv <= 0.85 * ls_.cv && rep.max_deviation < 0.01 && foot < 1e-9 && off < 0.01, msg);
done:
	qaws_mesh_reference_destroy(r);
	qaws_mesh_3d_destroy(m);
	qaws_mesh_3d_destroy(out);
	qaws_mesh_3d_destroy(loc);
}

static void test_blob(void)
{
	qaws_mesh_3d *m = ct_icosphere(NULL, 1), *out = NULL, *loc = NULL;
	qaws_mesh_reference* r = NULL;
	qaws_mesh_report rep, lrep;
	qaws_status st = QAWS_STATUS_INTERNAL_ERROR, ls = QAWS_STATUS_INTERNAL_ERROR;
	ct_shape s, l;
	double sec = 0;
	char msg[500];
	qaws_mesh_reference_create(m, NULL, &r);
	if (r)
	{
		out = ct_remesh(m, r, QAWS_REMESH_CVT, 0.11, NULL, 0, &rep, &st, &sec);
		loc = ct_remesh(m, r, QAWS_REMESH_LOCAL, 0.11, NULL, 0, &lrep, &ls, NULL);
	}
	if (!out || !loc)
	{
		TEST_ASSERT(0, "blob: remeshed");
		goto done;
	}
	ct_measure(out, &s);
	ct_measure(loc, &l);
	sprintf(msg, "blob by CVT at 0.11: %u vertices in %.2f s, irregular %.1f%% (LOCAL %.1f%%), mean edge %.4f (LOCAL %.4f), angles %.1f..%.1f (LOCAL %.1f..%.1f), edge CV %.3f (LOCAL %.3f), valid %d",
		s.nv, sec, 100.0 * (s.interior - s.v6) / s.nv, 100.0 * (l.interior - l.v6) / l.nv, s.mean, l.mean, (double)rep.quality.min_angle,
		(double)rep.quality.max_angle, (double)lrep.quality.min_angle, (double)lrep.quality.max_angle, s.cv, l.cv, qaws_mesh_3d_validate(out) == QAWS_STATUS_OK);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && qaws_mesh_3d_validate(out) == QAWS_STATUS_OK && s.interior - s.v6 <= 0.10 * s.nv && rep.quality.min_angle >= 30, msg);
done:
	qaws_mesh_reference_destroy(r);
	qaws_mesh_3d_destroy(m);
	qaws_mesh_3d_destroy(out);
	qaws_mesh_3d_destroy(loc);
}

static void test_cube(void)
{
	static double const x[24] = { -1, -1, -1, 1, -1, -1, 1, 1, -1, -1, 1, -1, -1, -1, 1, 1, -1, 1, 1, 1, 1, -1, 1, 1 };
	static unsigned int const f[24] = { 0, 3, 2, 1, 4, 5, 6, 7, 0, 1, 5, 4, 1, 2, 6, 5, 2, 3, 7, 6, 3, 0, 4, 7 };
	static unsigned int const fs[6] = { 4, 4, 4, 4, 4, 4 };
	qaws_mesh_3d *m = ct_make(x, 8, f, fs, 6), *out = NULL;
	qaws_mesh_report rep;
	qaws_status st = QAWS_STATUS_INTERNAL_ERROR;
	qaws_attribute_desc ad;
	qaws_scalar vals[8];
	unsigned int i, va = 0, nv, nfe, corners = 0, found[8] = { 0 }, a2 = 0;
	double onc = 0, flen = 0, fon = 0, ev = 0;
	char msg[500];
	for (i = 0; i < 8; i++)
		vals[i] = (qaws_scalar)(x[3 * i] + 2 * x[3 * i + 1] + 3 * x[3 * i + 2]);
	memset(&ad, 0, sizeof ad);
	ad.name = "f"; ad.element = QAWS_MESH_VERTEX; ad.components = 1;
	qaws_mesh_3d_add_attribute(m, &ad, vals, &va);
	qaws_mesh_3d_detect_features(m, (qaws_scalar)(M_PI / 6), 0, 0);
	out = ct_remesh(m, NULL, QAWS_REMESH_CVT, 0.25, NULL, 0, &rep, &st, NULL);
	if (!out)
	{
		TEST_ASSERT(0, "cube: remeshed");
		goto done;
	}
	nv = qaws_mesh_3d_get_vertex_count(out);
	nfe = qaws_mesh_3d_get_feature_edge_count(out);
	{
		qaws_scalar *p = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv), *fv = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nv);
		unsigned int* fe = (unsigned int*)malloc(sizeof(unsigned int) * 2 * (nfe + 1)), k;
		qaws_mesh_3d_get_positions(out, p);
		qaws_mesh_3d_find_attribute(out, "f", &a2);
		qaws_mesh_3d_get_attribute(out, a2, fv);
		for (i = 0; i < nv; i++)
		{
			double mx = fmax(fabs(p[3 * i]), fmax(fabs(p[3 * i + 1]), fabs(p[3 * i + 2])));
			onc = fmax(onc, fabs(mx - 1));
			ev = fmax(ev, fabs(fv[i] - (p[3 * i] + 2 * p[3 * i + 1] + 3 * p[3 * i + 2])));
			for (k = 0; k < 8; k++)
				if (fabs(p[3 * i] - x[3 * k]) + fabs(p[3 * i + 1] - x[3 * k + 1]) + fabs(p[3 * i + 2] - x[3 * k + 2]) < 1e-9)
					found[k] = 1;
			corners += qaws_mesh_3d_is_feature_corner(out, i);
		}
		qaws_mesh_3d_get_feature_edges(out, fe);
		for (i = 0; i < nfe; i++)
		{
			qaws_scalar const *a = p + 3 * fe[2 * i], *b = p + 3 * fe[2 * i + 1];
			unsigned int on = 0;
			flen += sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]));
			for (k = 0; k < 3; k++)
				on += fabs(fabs(a[k]) - 1) < 1e-9 && fabs(fabs(b[k]) - 1) < 1e-9;
			fon += on < 2;
		}
		for (k = 0; k < 8; k++)
			corners = corners >= 8 && found[k] ? corners : 0;
		free(p); free(fv); free(fe);
	}
	sprintf(msg, "cube by CVT at 0.25: %u vertices, valid %d, on the cube %.2g, corners %u, %u crease edges of length %.12g (%.0f off the cube's edges), linear attribute %.2g, angles %.1f..%.1f",
		nv, qaws_mesh_3d_validate(out) == QAWS_STATUS_OK, onc, corners, nfe, flen, fon, ev, (double)rep.quality.min_angle, (double)rep.quality.max_angle);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && qaws_mesh_3d_validate(out) == QAWS_STATUS_OK && onc < 1e-9 && corners >= 8 && fabs(flen - 24) < 1e-9 && fon == 0
		&& ev < 1e-9 && rep.quality.min_angle >= 30, msg);
done:
	qaws_mesh_3d_destroy(m);
	qaws_mesh_3d_destroy(out);
}

/* the mean of |edge / size - 1| (the size the mean of the ends'), and the mean edge on the bump and far from it */
static void ct_sizes(qaws_mesh_3d const* m, qaws_size_field const* sf, double* dev, double* bump, double* far_)
{
	unsigned int nv = qaws_mesh_3d_get_vertex_count(m), ne = qaws_mesh_3d_get_edge_count(m), i, nb = 0, nf = 0;
	unsigned int* e = (unsigned int*)malloc(sizeof(unsigned int) * 2 * (ne + 1));
	qaws_scalar *p = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * (nv + 1)), *h = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (nv + 1));
	qaws_mesh_location* l = (qaws_mesh_location*)malloc(sizeof(qaws_mesh_location) * (nv + 1));
	double lb = 0, lf = 0, sd = 0;
	qaws_mesh_3d_get_positions(m, p);
	qaws_mesh_3d_get_footpoints(m, l);
	qaws_mesh_3d_get_edges(m, e);
	qaws_size_field_evaluate(sf, l, p, nv, h);
	for (i = 0; i < ne; i++)
	{
		qaws_scalar const *a = p + 3 * e[2 * i], *b = p + 3 * e[2 * i + 1];
		double len = sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]));
		sd += fabs(len / (0.5 * (h[e[2 * i]] + h[e[2 * i + 1]])) - 1);
		if (a[2] > 1.25 && b[2] > 1.25)
		{
			lb += len;
			nb++;
		}
		else if (a[2] < -0.3 && b[2] < -0.3)
		{
			lf += len;
			nf++;
		}
	}
	*dev = ne ? sd / ne : 0;
	*bump = nb ? lb / nb : 0;
	*far_ = nf ? lf / nf : 1;
	free(e); free(p); free(h); free(l);
}

static void test_adaptive(void)
{
	qaws_mesh_3d *m = ct_icosphere(NULL, 0), *out = NULL, *loc = NULL;
	qaws_mesh_reference* r = NULL;
	qaws_curvature_size_desc cd;
	qaws_size_field* sf = NULL;
	qaws_mesh_report rep, lrep;
	qaws_status st = QAWS_STATUS_INTERNAL_ERROR, ls = QAWS_STATUS_INTERNAL_ERROR;
	unsigned int nv, i;
	double dev = 0, ldev = 0, bump = 0, far_ = 1, lbump = 0, lfar = 1, sec = 0;
	qaws_scalar* p;
	char msg[500];
	nv = qaws_mesh_3d_get_vertex_count(m);
	p = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	qaws_mesh_3d_get_positions(m, p);
	for (i = 0; i < nv; i++)
	{
		double d = p[3 * i] * p[3 * i] + p[3 * i + 1] * p[3 * i + 1] + (p[3 * i + 2] - 1) * (p[3 * i + 2] - 1), s = 1 + 0.6 * exp(-6 * d);
		p[3 * i] = (qaws_scalar)(p[3 * i] * s);
		p[3 * i + 1] = (qaws_scalar)(p[3 * i + 1] * s);
		p[3 * i + 2] = (qaws_scalar)(p[3 * i + 2] * s);
	}
	qaws_mesh_3d_set_positions(m, p);
	free(p);
	qaws_mesh_reference_create(m, NULL, &r);
	memset(&cd, 0, sizeof cd);
	cd.angle_tolerance = (qaws_scalar)0.2;
	cd.min_length = (qaws_scalar)0.02;
	cd.max_length = (qaws_scalar)0.4;
	cd.gradation = (qaws_scalar)0.4;
	if (r)
		qaws_size_field_create_curvature(r, &cd, &sf);
	if (sf)
	{
		out = ct_remesh(m, r, QAWS_REMESH_CVT, 0, sf, 0, &rep, &st, &sec);
		loc = ct_remesh(m, r, QAWS_REMESH_LOCAL, 0, sf, 0, &lrep, &ls, NULL);
	}
	if (!out || !loc)
	{
		TEST_ASSERT(0, "adaptive: remeshed");
		goto done;
	}
	ct_sizes(out, sf, &dev, &bump, &far_);
	ct_sizes(loc, sf, &ldev, &lbump, &lfar);
	sprintf(msg, "adaptive by CVT: %u vertices in %.2f s, valid %d, mean |edge / size - 1| %.3f (LOCAL %.3f), mean edge on the bump over far from it %.3f (LOCAL %.3f), angles %.1f..%.1f (LOCAL %.1f..%.1f)",
		qaws_mesh_3d_get_vertex_count(out), sec, qaws_mesh_3d_validate(out) == QAWS_STATUS_OK, dev, ldev, bump / far_, lbump / lfar,
		(double)rep.quality.min_angle, (double)rep.quality.max_angle, (double)lrep.quality.min_angle, (double)lrep.quality.max_angle);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && qaws_mesh_3d_validate(out) == QAWS_STATUS_OK && dev <= 0.15 && dev < ldev && bump > 0 && bump / far_ < 0.6, msg);
done:
	qaws_size_field_destroy(sf);
	qaws_mesh_reference_destroy(r);
	qaws_mesh_3d_destroy(m);
	qaws_mesh_3d_destroy(out);
	qaws_mesh_3d_destroy(loc);
}

static void test_convergence(void)
{
	qaws_mesh_3d *m = ct_icosphere(NULL, 1), *a = NULL, *b = NULL;
	qaws_mesh_reference* r = NULL;
	qaws_mesh_report rep;
	qaws_status s1 = QAWS_STATUS_INTERNAL_ERROR, s2 = QAWS_STATUS_INTERNAL_ERROR;
	rv_probe one, all;
	ct_shape sa, sb;
	char msg[400];
	memset(&one, 0, sizeof one);
	memset(&all, 0, sizeof all);
	qaws_mesh_reference_create(m, NULL, &r);
	if (r)
	{
		rv_trace = &one;
		a = ct_remesh(m, r, QAWS_REMESH_CVT, 0.11, NULL, 1, &rep, &s1, NULL);
		rv_trace = &all;
		b = ct_remesh(m, r, QAWS_REMESH_CVT, 0.11, NULL, 30, &rep, &s2, NULL);
		rv_trace = NULL;
	}
	if (!a || !b)
	{
		TEST_ASSERT(0, "convergence: remeshed");
		goto done;
	}
	ct_measure(a, &sa);
	ct_measure(b, &sb);
	sprintf(msg, "blob by CVT: energy after Lloyd %.6g, after 1 L-BFGS step %.6g, after 30 %.6g (%u evaluations); edge CV after 1 step %.4f, after 30 %.4f",
		all.lloyd_energy, one.energy, all.energy, all.evaluations, sa.cv, sb.cv);
	printf("  %s\n", msg);
	TEST_ASSERT(s1 == QAWS_STATUS_OK && s2 == QAWS_STATUS_OK && all.energy > 0 && all.energy <= all.lloyd_energy && one.energy <= one.lloyd_energy
		&& sb.cv <= sa.cv, msg);
done:
	qaws_mesh_reference_destroy(r);
	qaws_mesh_3d_destroy(m);
	qaws_mesh_3d_destroy(a);
	qaws_mesh_3d_destroy(b);
}

static void test_slab(void)
{
	static double const x[24] = { -0.5, -0.5, -0.01, 0.5, -0.5, -0.01, 0.5, 0.5, -0.01, -0.5, 0.5, -0.01, -0.5, -0.5, 0.01, 0.5, -0.5, 0.01, 0.5, 0.5, 0.01,
		-0.5, 0.5, 0.01 };
	static unsigned int const f[24] = { 0, 3, 2, 1, 4, 5, 6, 7, 0, 1, 5, 4, 1, 2, 6, 5, 2, 3, 7, 6, 3, 0, 4, 7 };
	static unsigned int const fs[6] = { 4, 4, 4, 4, 4, 4 };
	qaws_mesh_3d *m = ct_make(x, 8, f, fs, 6), *out = NULL;
	qaws_mesh_report rep;
	qaws_mesh_topology t;
	qaws_status st = QAWS_STATUS_INTERNAL_ERROR;
	char msg[300];
	qaws_mesh_3d_detect_features(m, (qaws_scalar)(M_PI / 6), 0, 0);
	out = ct_remesh(m, NULL, QAWS_REMESH_CVT, 0.05, NULL, 0, &rep, &st, NULL);
	if (!out)
	{
		TEST_ASSERT(0, "slab: remeshed");
		goto done;
	}
	qaws_mesh_3d_get_topology(out, &t);
	sprintf(msg, "a slab 0.02 thick by CVT at 0.05: %u vertices, valid %d, closed %d, genus %u, angles %.1f..%.1f", qaws_mesh_3d_get_vertex_count(out),
		qaws_mesh_3d_validate(out) == QAWS_STATUS_OK, t.closed, t.genus, (double)rep.quality.min_angle, (double)rep.quality.max_angle);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && qaws_mesh_3d_validate(out) == QAWS_STATUS_OK && t.closed && t.genus == 0 && t.components == 1, msg);
done:
	qaws_mesh_3d_destroy(m);
	qaws_mesh_3d_destroy(out);
}

/* an open bumpy square (test 124's): its border a locked line like a crease */
static void test_open(void)
{
	unsigned int n = 6, i, j, k = 0, nb = 0, okb = 1, found = 0, f[6 * 36];
	double x[3 * 49];
	qaws_mesh_3d *m, *out = NULL;
	qaws_mesh_report rep;
	qaws_mesh_topology t;
	qaws_status st = QAWS_STATUS_INTERNAL_ERROR;
	char msg[300];
	for (j = 0; j <= n; j++)
		for (i = 0; i <= n; i++)
		{
			double u = (double)i / n, v = (double)j / n;
			x[3 * (j * (n + 1) + i)] = u;
			x[3 * (j * (n + 1) + i) + 1] = v;
			x[3 * (j * (n + 1) + i) + 2] = 0.15 * sin(3 * u) * cos(2 * v);
		}
	for (j = 0; j < n; j++)
		for (i = 0; i < n; i++)
		{
			unsigned int a = j * (n + 1) + i;
			f[k++] = a; f[k++] = a + 1; f[k++] = a + n + 2;
			f[k++] = a; f[k++] = a + n + 2; f[k++] = a + n + 1;
		}
	m = ct_make(x, (n + 1) * (n + 1), f, NULL, 2 * n * n);
	out = ct_remesh(m, NULL, QAWS_REMESH_CVT, 0.06, NULL, 0, &rep, &st, NULL);
	if (!out)
	{
		TEST_ASSERT(0, "open square: remeshed");
		goto done;
	}
	qaws_mesh_3d_get_topology(out, &t);
	{
		unsigned int nv = qaws_mesh_3d_get_vertex_count(out);
		qaws_scalar* p = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
		qaws_mesh_3d_get_positions(out, p);
		for (i = 0; i < nv; i++)
		{
			double u = p[3 * i], v = p[3 * i + 1];
			found += (fabs(u) < 1e-9 || fabs(u - 1) < 1e-9) && (fabs(v) < 1e-9 || fabs(v - 1) < 1e-9);
			if (!qaws_mesh_3d_is_boundary_vertex(out, i))
				continue;
			nb++;
			okb &= fabs(u) < 1e-9 || fabs(u - 1) < 1e-9 || fabs(v) < 1e-9 || fabs(v - 1) < 1e-9;
		}
		free(p);
	}
	sprintf(msg, "an open bumpy square by CVT at 0.06: valid %d, %u border loop, %u border vertices on the square's border %u, corners kept %u, angles %.1f..%.1f",
		qaws_mesh_3d_validate(out) == QAWS_STATUS_OK, t.boundary_loops, nb, okb, found, (double)rep.quality.min_angle, (double)rep.quality.max_angle);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && qaws_mesh_3d_validate(out) == QAWS_STATUS_OK && t.boundary_loops == 1 && okb && found == 4 && nb > 40
		&& rep.quality.min_angle >= 25, msg);   /* (measured 28.8: the triangles at the locked corners) */
done:
	qaws_mesh_3d_destroy(m);
	qaws_mesh_3d_destroy(out);
}

static void test_determinism(void)
{
	qaws_mesh_3d *m = ct_icosphere(NULL, 1), *a = NULL, *b = NULL;
	qaws_mesh_report rep;
	qaws_status s1 = QAWS_STATUS_INTERNAL_ERROR, s2 = QAWS_STATUS_INTERNAL_ERROR;
	int same = 0;
	char msg[200];
	a = ct_remesh(m, NULL, QAWS_REMESH_CVT, 0.15, NULL, 10, &rep, &s1, NULL);
	b = ct_remesh(m, NULL, QAWS_REMESH_CVT, 0.15, NULL, 10, &rep, &s2, NULL);
	if (a && b && qaws_mesh_3d_get_vertex_count(a) == qaws_mesh_3d_get_vertex_count(b) && qaws_mesh_3d_get_triangle_count(a) == qaws_mesh_3d_get_triangle_count(b))
	{
		unsigned int nv = qaws_mesh_3d_get_vertex_count(a), nt = qaws_mesh_3d_get_triangle_count(a);
		qaws_scalar *pa = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv), *pb = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
		unsigned int *ta = (unsigned int*)malloc(sizeof(unsigned int) * 3 * nt), *tb = (unsigned int*)malloc(sizeof(unsigned int) * 3 * nt);
		qaws_mesh_3d_get_positions(a, pa);
		qaws_mesh_3d_get_positions(b, pb);
		qaws_mesh_3d_get_triangles(a, ta);
		qaws_mesh_3d_get_triangles(b, tb);
		same = !memcmp(pa, pb, sizeof(qaws_scalar) * 3 * nv) && !memcmp(ta, tb, sizeof(unsigned int) * 3 * nt);
		free(pa); free(pb); free(ta); free(tb);
	}
	sprintf(msg, "CVT twice on the blob: the same vertices and faces %d", same);
	printf("  %s\n", msg);
	TEST_ASSERT(s1 == QAWS_STATUS_OK && s2 == QAWS_STATUS_OK && same, msg);
	qaws_mesh_3d_destroy(m);
	qaws_mesh_3d_destroy(a);
	qaws_mesh_3d_destroy(b);
}

/* about 2k seeds: the time printed */
static void test_timing(void)
{
	qaws_mesh_3d *m = ct_icosphere(NULL, 1), *out = NULL;
	qaws_mesh_reference* r = NULL;
	qaws_mesh_report rep;
	qaws_status st = QAWS_STATUS_INTERNAL_ERROR;
	rv_probe pr;
	double sec = 0;
	memset(&pr, 0, sizeof pr);
	qaws_mesh_reference_create(m, NULL, &r);
	rv_trace = &pr;
	if (r)
		out = ct_remesh(m, r, QAWS_REMESH_CVT, 0.09, NULL, 0, &rep, &st, &sec);
	rv_trace = NULL;
	printf("  blob by CVT at 0.09: %u vertices, %u energy evaluations, %.2f s\n", out ? qaws_mesh_3d_get_vertex_count(out) : 0, pr.evaluations, sec);
	TEST_ASSERT(st == QAWS_STATUS_OK && out, "timing run");
	qaws_mesh_reference_destroy(r);
	qaws_mesh_3d_destroy(m);
	qaws_mesh_3d_destroy(out);
}

int test_132_remesh_cvt_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 132: CVT remeshing\n");
	test_sphere();
	test_blob();
	test_cube();
	test_adaptive();
	test_convergence();
	test_slab();
	test_open();
	test_determinism();
	test_timing();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
