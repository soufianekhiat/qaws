/*
 * Test 131: remeshing quality (LOCAL's defaults, the angle-bounded pass, the
 * distance guard)
 *
 *   - default: the can of test 125 (8 sides, creases at 50 degrees) within
 *     [30, 110] degrees; the bumpy blob of showcase figure 2 with its mean
 *     edge within 4% of the target, 80% of its vertices of valence 6, 99%
 *     of valence 5 to 7, 99% of its edges Delaunay; the adaptive bumped
 *     sphere of test 124 twice as fine on the bump as far from it;
 *   - angles bounded to [35, 86] degrees (Wang et al. 2019's defaults) on an
 *     icosphere, the blob and the adaptive sphere, 98% of valence 5 to 7,
 *     the mean edge within 10% of the target (the uniform ones); [30, 90]
 *     and [35, 86] on the can, corners out of bounds only in exempt
 *     sectors; a cube at 0.08, 0.15, 0.25 and 0.6 (its 90-degree corner
 *     sectors start with one face at the coarse sizes) and an open square
 *     (one face at each border corner), every corner not exempt inside;
 *   - tolerances of 2e-3, 5e-4 and 2e-4 on the blob, and 2e-3 with
 *     [35, 86] on the blob and the adaptive sphere: every sampled point
 *     (face centres, edge midpoints and quarter points) within it, no
 *     corner past 170 degrees;
 *   - in every mode: input corners kept where they were, output corners
 *     at input vertices, crease vertices on creases (the border counted),
 *     footpoints evaluated, valid, closed, genus kept;
 *   - open meshes from qaws_mesh_generate_surface (the showcase's wavy
 *     patch, a bilinear square, a cylinder piece; their footpoints on the
 *     surfaces, not on their own reference) at 0.1 and 0.06: area and
 *     border within 2%, faces near area / (0.43 L^2) (before: 4 to 23
 *     faces, the input's (face, u, v) read as reference locations);
 *   - refusals; timings printed.
 *
 * By default the blob's tolerance is 5e-3, and the cube at 0.08 and 0.15,
 * the can in [35, 86] and the smaller tolerances and those with the angle
 * bounds run with QAWS_TEST_HEAVY=1 only (about 80 s of the 100).
 *
 * Measured when written (all targets met): the can 33.9..108.0 degrees
 * (27.2..117.4 before); the blob's mean edge +1.5% (-7%), valence 6 at
 * 82.5% (78%), Delaunay 100%; every bounded case 0 corners out (the cube
 * at 0.15 kept one at 86.45 with 8 relocation directions, at 0.6 13 with
 * 2k attempts a round); the tolerance held at every sample, 2e-3 with 66%
 * more vertices, 2e-4 with angles 9.8..144.9 (0 to 180 before); the angle
 * pass 1.1 times the default remesh's time on the blob.
 */

#include "test_common.h"
#include "qaws_remesh.h"
#include "qaws_mesh_generate.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

static qaws_mesh_3d* rq_make(double const* x, unsigned int nv, unsigned int const* f, unsigned int const* fs, unsigned int nf)
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
	qaws_mesh_3d_create(&m, &d);
	free(p);
	return m;
}

/* test 124's icosahedron split once, on the unit sphere, with its radial normals */
static qaws_mesh_3d* rq_icosphere(qaws_scalar* normals)
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
		double r = sqrt(x[3 * i] * x[3 * i] + x[3 * i + 1] * x[3 * i + 1] + x[3 * i + 2] * x[3 * i + 2]);
		for (k = 0; k < 3; k++)
		{
			x[3 * i + k] /= r;
			if (normals)
				normals[3 * i + k] = (qaws_scalar)x[3 * i + k];
		}
	}
	return rq_make(x, nv, f, NULL, nf);
}

/* test 124's icosphere bumped at its pole */
static qaws_mesh_3d* rq_bumped(void)
{
	qaws_mesh_3d* m = rq_icosphere(NULL);
	unsigned int nv = qaws_mesh_3d_get_vertex_count(m), i;
	qaws_scalar* p = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	qaws_mesh_3d_get_positions(p, m);
	for (i = 0; i < nv; i++)
	{
		double d = p[3 * i] * p[3 * i] + p[3 * i + 1] * p[3 * i + 1] + (p[3 * i + 2] - 1) * (p[3 * i + 2] - 1), s = 1 + 0.6 * exp(-6 * d);
		p[3 * i] = (qaws_scalar)(p[3 * i] * s);
		p[3 * i + 1] = (qaws_scalar)(p[3 * i + 1] * s);
		p[3 * i + 2] = (qaws_scalar)(p[3 * i + 2] * s);
	}
	qaws_mesh_3d_set_positions(m, p);
	free(p);
	return m;
}

/* showcase figure 2's blob: the icosahedron split once by edge splits and flips, pushed in and out smoothly */
static qaws_mesh_3d* rq_blob(void)
{
	static unsigned int const f[60] = { 0, 11, 5, 0, 5, 1, 0, 1, 7, 0, 7, 10, 0, 10, 11, 1, 5, 9, 5, 11, 4, 11, 10, 2, 10, 7, 6, 7, 1, 8,
		3, 9, 4, 3, 4, 2, 3, 2, 6, 3, 6, 8, 3, 8, 9, 4, 9, 5, 2, 4, 11, 6, 2, 10, 8, 6, 7, 9, 8, 1 };
	double g = (1 + sqrt(5.0)) / 2, x[36];
	double const v[36] = { -1, g, 0, 1, g, 0, -1, -g, 0, 1, -g, 0, 0, -1, g, 0, 1, g, 0, -1, -g, 0, 1, -g, g, 0, -1, g, 0, 1, -g, 0, -1, -g, 0, 1 };
	qaws_mesh_3d* m;
	unsigned int i, n, ne, nv0, vn;
	unsigned int *e, *src;
	qaws_scalar* p;
	for (i = 0; i < 36; i++)
		x[i] = v[i] / sqrt(1 + g * g);
	m = rq_make(x, 12, f, NULL, 20);
	ne = qaws_mesh_3d_get_edge_count(m);
	nv0 = qaws_mesh_3d_get_vertex_count(m);
	e = (unsigned int*)malloc(sizeof(unsigned int) * 2 * ne);
	src = (unsigned int*)malloc(sizeof(unsigned int) * 2 * ne);
	qaws_mesh_3d_get_edges(e, m);
	for (i = 0; i < ne; i++)
		if (qaws_mesh_3d_split_edge(&vn, m, NULL, NULL, e[2 * i], e[2 * i + 1]) == QAWS_STATUS_OK)
		{
			src[2 * (vn - nv0)] = e[2 * i];
			src[2 * (vn - nv0) + 1] = e[2 * i + 1];
		}
	free(e);
	ne = qaws_mesh_3d_get_edge_count(m);
	e = (unsigned int*)malloc(sizeof(unsigned int) * 2 * ne);
	qaws_mesh_3d_get_edges(e, m);
	for (i = 0; i < ne; i++)
	{
		unsigned int a = e[2 * i], b = e[2 * i + 1], nw, old;
		if ((a >= nv0) == (b >= nv0))
			continue;
		nw = a >= nv0 ? a : b;
		old = a >= nv0 ? b : a;
		if (old != src[2 * (nw - nv0)] && old != src[2 * (nw - nv0) + 1] && qaws_mesh_3d_can_flip_edge(m, a, b))
			qaws_mesh_3d_flip_edge(m, a, b);
	}
	free(e);
	free(src);
	n = qaws_mesh_3d_get_vertex_count(m);
	p = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * n);
	qaws_mesh_3d_get_positions(p, m);
	for (i = 0; i < n; i++)
	{
		double q[3], l, s;
		q[0] = p[3 * i]; q[1] = p[3 * i + 1]; q[2] = p[3 * i + 2];
		l = sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2]);
		q[0] /= l; q[1] /= l; q[2] /= l;
		s = 1 + 0.32 * sin(2.2 * q[0] + 0.4) * cos(1.7 * q[1]) + 0.2 * q[2] * q[2];
		p[3 * i] = (qaws_scalar)(q[0] * 1.25 * s);
		p[3 * i + 1] = (qaws_scalar)(q[1] * s);
		p[3 * i + 2] = (qaws_scalar)(q[2] * 0.9 * s);
	}
	qaws_mesh_3d_set_positions(m, p);
	free(p);
	return m;
}

/* test 125's can: 8 sides, its rims creases (the sides meet at 45 degrees, under the 50 detected) */
static qaws_mesh_3d* rq_can(void)
{
	unsigned int n = 8, i, k = 0;
	double x[3 * 18];
	unsigned int fv[3 * 32];
	qaws_mesh_3d* can;
	for (i = 0; i < n; i++)
	{
		double a = 2 * M_PI * i / n;
		x[3 * i] = cos(a); x[3 * i + 1] = sin(a); x[3 * i + 2] = -0.8;
		x[3 * (n + i)] = cos(a); x[3 * (n + i) + 1] = sin(a); x[3 * (n + i) + 2] = 0.8;
	}
	x[6 * n] = 0; x[6 * n + 1] = 0; x[6 * n + 2] = -0.8;
	x[6 * n + 3] = 0; x[6 * n + 4] = 0; x[6 * n + 5] = 0.8;
	for (i = 0; i < n; i++)
	{
		unsigned int a = i, b = (i + 1) % n;
		fv[k++] = a; fv[k++] = b; fv[k++] = n + b;
		fv[k++] = a; fv[k++] = n + b; fv[k++] = n + a;
		fv[k++] = 2 * n; fv[k++] = b; fv[k++] = a;
		fv[k++] = 2 * n + 1; fv[k++] = n + a; fv[k++] = n + b;
	}
	can = rq_make(x, 2 * n + 2, fv, NULL, 4 * n);
	qaws_mesh_3d_detect_features(can, (qaws_scalar)(50 * M_PI / 180), 0, 0);
	return can;
}

/* test 124's cube of quads with its creases */
static qaws_mesh_3d* rq_cube(void)
{
	static double const x[24] = { -1, -1, -1, 1, -1, -1, 1, 1, -1, -1, 1, -1, -1, -1, 1, 1, -1, 1, 1, 1, 1, -1, 1, 1 };
	static unsigned int const f[24] = { 0, 3, 2, 1, 4, 5, 6, 7, 0, 1, 5, 4, 1, 2, 6, 5, 2, 3, 7, 6, 3, 0, 4, 7 };
	static unsigned int const fs[6] = { 4, 4, 4, 4, 4, 4 };
	qaws_mesh_3d* m = rq_make(x, 8, f, fs, 6);
	qaws_mesh_3d_detect_features(m, (qaws_scalar)(M_PI / 6), 0, 0);
	return m;
}

/* an open flat square [-1, 1]^2 of 4 x 4 cells, its border a feature with a 90-degree corner (one face) at each end */
static qaws_mesh_3d* rq_square(void)
{
	double x[3 * 25];
	unsigned int f[3 * 32], i, j, k = 0;
	for (j = 0; j < 5; j++)
		for (i = 0; i < 5; i++)
		{
			x[3 * (5 * j + i)] = -1 + 0.5 * i;
			x[3 * (5 * j + i) + 1] = -1 + 0.5 * j;
			x[3 * (5 * j + i) + 2] = 0;
		}
	for (j = 0; j < 4; j++)
		for (i = 0; i < 4; i++)
		{
			unsigned int a = 5 * j + i, b = a + 1, c = a + 6, d = a + 5;
			f[k++] = a; f[k++] = b; f[k++] = c;
			f[k++] = a; f[k++] = c; f[k++] = d;
		}
	return rq_make(x, 25, f, NULL, 32);
}

static qaws_mesh_3d* rq_remesh(qaws_mesh_3d const* m, qaws_mesh_reference const* r, double L, qaws_size_field const* sf, double lo, double hi,
	double tol, qaws_mesh_report* rep, double* seconds)
{
	qaws_remesh_desc d;
	qaws_mesh_3d* out = NULL;
	clock_t t0 = clock();
	memset(&d, 0, sizeof d);
	d.reference = r;
	d.target_length = (qaws_scalar)L;
	d.size = sf;
	d.min_angle = (qaws_scalar)(lo * M_PI / 180);
	d.max_angle = (qaws_scalar)(hi * M_PI / 180);
	d.tolerance = (qaws_scalar)tol;
	if (qaws_mesh_remesh(&out, rep, m, &d) != QAWS_STATUS_OK)
		out = NULL;
	if (seconds)
		*seconds = (double)(clock() - t0) / CLOCKS_PER_SEC;
	return out;
}

/* ---- measures ---------------------------------------------------------- */

typedef struct rq_meas
{
	double amin, amax;              /* every corner, degrees */
	double emin, emax;              /* the corners not exempt */
	unsigned int corners, inside;   /* all corners; inside the bounds */
	unsigned int out_free;          /* out of the bounds and not exempt */
	double v6, v567;                /* shares of the vertices off features and the border */
	double mean_edge;
	double delaunay;                /* share of the inner edges off the features with their opposite angles at most pi */
	unsigned int nv, nf;
} rq_meas;

static int rq_key_cmp(void const* a, void const* b)
{
	unsigned long long x = *(unsigned long long const*)a, y = *(unsigned long long const*)b;
	return x < y ? -1 : x > y;
}

static unsigned long long rq_key(unsigned int a, unsigned int b)
{
	return a < b ? ((unsigned long long)a << 32) | b : ((unsigned long long)b << 32) | a;
}

static double rq_angle(qaws_scalar const* p, unsigned int v, unsigned int a, unsigned int b)
{
	double u[3], t[3], uu = 0, tt = 0, ut = 0, c;
	unsigned int k;
	for (k = 0; k < 3; k++)
	{
		u[k] = p[3 * a + k] - p[3 * v + k];
		t[k] = p[3 * b + k] - p[3 * v + k];
		uu += u[k] * u[k]; tt += t[k] * t[k]; ut += u[k] * t[k];
	}
	c = ut / sqrt(uu * tt);
	return acos(c < -1 ? -1 : c > 1 ? 1 : c) * 180 / M_PI;
}

/* the angles (exempt as the header says: a sector of a feature or border vertex no whole number of corners within [lo, hi]
   divides), valences, edges, the Delaunay share */
static void rq_measure(qaws_mesh_3d const* m, double lo, double hi, rq_meas* q)
{
	unsigned int nv = qaws_mesh_3d_get_vertex_count(m), nf = qaws_mesh_3d_get_triangle_count(m), ne = qaws_mesh_3d_get_edge_count(m);
	unsigned int nfe = qaws_mesh_3d_get_feature_edge_count(m), i, j, k, *t, *fe, *e, *first, *list, *grp, *cnt, *onf, ninner = 0, n6 = 0, n567 = 0;
	unsigned int ndel = 0, nin = 0;
	unsigned long long* keys;
	double *ang, *sec, len = 0;
	qaws_scalar* p = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	t = (unsigned int*)malloc(sizeof(unsigned int) * 3 * nf);
	fe = (unsigned int*)malloc(sizeof(unsigned int) * 2 * (nfe + 1));
	e = (unsigned int*)malloc(sizeof(unsigned int) * 2 * ne);
	keys = (unsigned long long*)malloc(sizeof(unsigned long long) * (nfe + 1));
	first = (unsigned int*)calloc(nv + 1, sizeof(unsigned int));
	list = (unsigned int*)malloc(sizeof(unsigned int) * 3 * nf);
	grp = (unsigned int*)malloc(sizeof(unsigned int) * 3 * nf);
	cnt = (unsigned int*)calloc(nv + 1, sizeof(unsigned int));
	onf = (unsigned int*)calloc(nv + 1, sizeof(unsigned int));
	ang = (double*)malloc(sizeof(double) * 3 * nf);
	sec = (double*)malloc(sizeof(double) * 3 * nf);
	qaws_mesh_3d_get_positions(p, m);
	qaws_mesh_3d_get_triangles(t, m);
	qaws_mesh_3d_get_feature_edges(fe, m);
	qaws_mesh_3d_get_edges(e, m);
	for (i = 0; i < nfe; i++)
	{
		keys[i] = rq_key(fe[2 * i], fe[2 * i + 1]);
		onf[fe[2 * i]] = onf[fe[2 * i + 1]] = 1;
	}
	qsort(keys, nfe, sizeof *keys, rq_key_cmp);
	for (i = 0; i < nv; i++)
		onf[i] |= qaws_mesh_3d_is_boundary_vertex(m, i);
	/* corners grouped per vertex, then into sectors across the edges that are not features */
	for (i = 0; i < 3 * nf; i++)
		first[t[i] + 1]++;
	for (i = 0; i < nv; i++)
		first[i + 1] += first[i];
	for (i = 0; i < 3 * nf; i++)
	{
		list[first[t[i]] + cnt[t[i]]++] = i;
		ang[i] = rq_angle(p, t[i], t[3 * (i / 3) + (i + 1) % 3], t[3 * (i / 3) + (i + 2) % 3]);
		grp[i] = i;
	}
	for (i = 0; i < nv; i++)
		for (j = first[i]; j < first[i + 1]; j++)
			for (k = j + 1; k < first[i + 1]; k++)
			{
				unsigned int a = list[j], b = list[k], s, x = 0xFFFFFFFFu;
				/* a shared neighbour: an edge (i, x) both corners hold */
				for (s = 1; s < 3 && x == 0xFFFFFFFFu; s++)
				{
					unsigned int na = t[3 * (a / 3) + (a + s) % 3], r;
					for (r = 1; r < 3; r++)
						if (t[3 * (b / 3) + (b + r) % 3] == na)
							x = na;
				}
				if (x == 0xFFFFFFFFu)
				{
					continue;
				}
				{
					unsigned long long key = rq_key(i, x);
					if (nfe && bsearch(&key, keys, nfe, sizeof *keys, rq_key_cmp))
						continue;
				}
				/* union: relabel b's group as a's */
				{
					unsigned int ga = grp[a], gb = grp[b], r;
					if (ga == gb)
						continue;
					for (r = first[i]; r < first[i + 1]; r++)
						if (grp[list[r]] == gb)
							grp[list[r]] = ga;
				}
			}
	for (i = 0; i < 3 * nf; i++)
		sec[i] = 0;
	for (i = 0; i < 3 * nf; i++)
		sec[grp[i]] += ang[i];
	q->amin = q->emin = 1e9;
	q->amax = q->emax = 0;
	q->corners = 3 * nf;
	q->inside = 0;
	q->out_free = 0;
	for (i = 0; i < 3 * nf; i++)
	{
		int ex = 0, in = (lo <= 0 || ang[i] >= lo) && (hi <= 0 || ang[i] <= hi);
		if (onf[t[i]] && lo > 0)
		{
			double th = sec[grp[i]], kk = floor(th / lo);
			ex = kk < 1 || (hi > 0 && th / kk > hi);
		}
		q->amin = fmin(q->amin, ang[i]);
		q->amax = fmax(q->amax, ang[i]);
		q->inside += in;
		if (!ex)
		{
			q->emin = fmin(q->emin, ang[i]);
			q->emax = fmax(q->emax, ang[i]);
			q->out_free += !in;
		}
	}
	for (i = 0; i < nv; i++)
	{
		unsigned int val = first[i + 1] - first[i];
		if (onf[i])
			continue;
		ninner++;
		n6 += val == 6;
		n567 += val >= 5 && val <= 7;
	}
	q->v6 = ninner ? (double)n6 / ninner : 1;
	q->v567 = ninner ? (double)n567 / ninner : 1;
	/* the inner edges off the features: the angles across them */
	for (i = 0; i < 3 * nf; i++)
	{
		unsigned int a = t[i], b = t[3 * (i / 3) + (i + 1) % 3], j2;
		unsigned long long key = rq_key(a, b);
		if (a > b || (nfe && bsearch(&key, keys, nfe, sizeof *keys, rq_key_cmp)))
			continue;
		/* the other face holding (b, a) */
		for (j2 = first[a]; j2 < first[a + 1]; j2++)
		{
			unsigned int c = list[j2];
			if (t[3 * (c / 3) + (c + 2) % 3] == b)
			{
				double s = ang[3 * (i / 3) + (i + 2) % 3] + ang[3 * (c / 3) + (c + 1) % 3];
				nin++;
				ndel += s <= 180 + 1e-6 * 180 / M_PI;
			}
		}
	}
	q->delaunay = nin ? (double)ndel / nin : 1;
	for (i = 0; i < ne; i++)
	{
		qaws_scalar const *a = p + 3 * e[2 * i], *b = p + 3 * e[2 * i + 1];
		len += sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]));
	}
	q->mean_edge = ne ? len / ne : 0;
	q->nv = nv;
	q->nf = nf;
	free(p); free(t); free(fe); free(e); free(keys); free(first); free(list); free(grp); free(cnt); free(onf); free(ang); free(sec);
}

/* the largest distance to the reference of the faces' centres, edge midpoints and quarter points */
static double rq_deviation(qaws_mesh_3d const* m, qaws_mesh_reference const* r)
{
	unsigned int nv = qaws_mesh_3d_get_vertex_count(m), nf = qaws_mesh_3d_get_triangle_count(m), i, j, k, n = 0;
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 3 * nf);
	qaws_scalar *p = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv), *s = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * 10 * nf);
	qaws_scalar* d = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 10 * nf);
	qaws_mesh_location* l = (qaws_mesh_location*)malloc(sizeof(qaws_mesh_location) * 10 * nf);
	double dmax = 0;
	qaws_mesh_3d_get_positions(p, m);
	qaws_mesh_3d_get_triangles(t, m);
	for (i = 0; i < nf; i++)
	{
		for (k = 0; k < 3; k++)
			s[3 * n + k] = (p[3 * t[3 * i] + k] + p[3 * t[3 * i + 1] + k] + p[3 * t[3 * i + 2] + k]) / 3;
		n++;
		for (j = 0; j < 3; j++)
		{
			qaws_scalar const *a = p + 3 * t[3 * i + j], *b = p + 3 * t[3 * i + (j + 1) % 3];
			double w[3] = { 0.25, 0.5, 0.75 };
			unsigned int h;
			for (h = 0; h < 3; h++, n++)
				for (k = 0; k < 3; k++)
					s[3 * n + k] = (qaws_scalar)((1 - w[h]) * a[k] + w[h] * b[k]);
		}
	}
	qaws_mesh_reference_project(l, d, r, s, NULL, n);
	for (i = 0; i < n; i++)
		dmax = fmax(dmax, d[i]);
	free(t); free(p); free(s); free(d); free(l);
	return dmax;
}

/* the invariants of every mode: valid, closedness and genus kept, footpoints evaluated, feature corners at input vertices,
   vertices on features with their footpoints on the reference's creases */
static int rq_invariants(qaws_mesh_3d const* in, qaws_mesh_3d const* out, qaws_mesh_reference const* r, char* why)
{
	unsigned int nv = qaws_mesh_3d_get_vertex_count(out), nvi = qaws_mesh_3d_get_vertex_count(in), nfe = qaws_mesh_3d_get_feature_edge_count(out);
	unsigned int nfi = qaws_mesh_3d_get_feature_edge_count(in), i, j, k, *fe, *fi, corners = 1, crease = 1;
	unsigned int nt = qaws_mesh_reference_get_triangle_count(r), nk, *tv = (unsigned int*)malloc(sizeof(unsigned int) * 3 * (nt + 1));
	unsigned long long* keys;
	unsigned char* onf;
	qaws_scalar *x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv), *y = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	qaws_scalar* xi = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nvi);
	qaws_mesh_location* l = (qaws_mesh_location*)malloc(sizeof(qaws_mesh_location) * nv);
	qaws_mesh_topology t0, t1;
	double foot = 0;
	int ok;
	fe = (unsigned int*)malloc(sizeof(unsigned int) * 2 * (nfe + 1));
	fi = (unsigned int*)malloc(sizeof(unsigned int) * 2 * (nfi + 1));
	keys = (unsigned long long*)malloc(sizeof(unsigned long long) * (nfi + 6 * nt + 1));
	onf = (unsigned char*)calloc(nvi + 1, 1);
	qaws_mesh_3d_get_positions(x, out);
	qaws_mesh_3d_get_positions(xi, in);
	qaws_mesh_3d_get_footpoints(l, out);
	qaws_mesh_reference_evaluate(y, NULL, NULL, r, l, nv);
	qaws_mesh_3d_get_feature_edges(fe, out);
	qaws_mesh_3d_get_feature_edges(fi, in);
	for (i = 0; i < nfi; i++)
	{
		keys[i] = rq_key(fi[2 * i], fi[2 * i + 1]);
		onf[fi[2 * i]] = onf[fi[2 * i + 1]] = 1;
	}
	/* the reference's border is a crease too: its edges in one triangle only */
	for (i = 0; i < nt; i++)
	{
		qaws_mesh_reference_get_triangle(tv + 3 * i, r, i);
		for (k = 0; k < 3; k++)
			keys[nfi + 3 * nt + 3 * i + k] = rq_key(tv[3 * i + k], tv[3 * i + (k + 1) % 3]);
	}
	qsort(keys + nfi + 3 * nt, 3 * nt, sizeof *keys, rq_key_cmp);
	for (i = 0, nk = nfi; i < 3 * nt; i++)
	{
		unsigned long long e = keys[nfi + 3 * nt + i];
		if ((i > 0 && keys[nfi + 3 * nt + i - 1] == e) || (i + 1 < 3 * nt && keys[nfi + 3 * nt + i + 1] == e))
			continue;
		keys[nk++] = e;
		onf[(unsigned int)(e >> 32)] = onf[(unsigned int)(e & 0xFFFFFFFFu)] = 1;
	}
	qsort(keys, nk, sizeof *keys, rq_key_cmp);
	for (i = 0; i < nv; i++)
	{
		foot = fmax(foot, fabs(x[3 * i] - y[3 * i]) + fabs(x[3 * i + 1] - y[3 * i + 1]) + fabs(x[3 * i + 2] - y[3 * i + 2]));
		if (qaws_mesh_3d_is_feature_corner(out, i))
		{
			int found = 0;
			for (j = 0; j < nvi && !found; j++)
				found = fabs(x[3 * i] - xi[3 * j]) + fabs(x[3 * i + 1] - xi[3 * j + 1]) + fabs(x[3 * i + 2] - xi[3 * j + 2]) < TEST_BOUND(1e-12, 32);
			corners &= found;
		}
	}
	/* and every input corner still one, where it was */
	for (j = 0; j < nvi; j++)
		if (qaws_mesh_3d_is_feature_corner(in, j))
		{
			int found = 0;
			for (i = 0; i < nv && !found; i++)
				found = qaws_mesh_3d_is_feature_corner(out, i)
					&& fabs(x[3 * i] - xi[3 * j]) + fabs(x[3 * i + 1] - xi[3 * j + 1]) + fabs(x[3 * i + 2] - xi[3 * j + 2]) < TEST_BOUND(1e-12, 32);
			corners &= found;
		}
	/* each end of a feature edge: its footpoint on a crease edge of its reference triangle (or at a vertex on one) */
	for (i = 0; i < 2 * nfe; i++)
	{
		qaws_mesh_location const* f = &l[fe[i]];
		unsigned int tri[3], on = 0;
		double b[3];
		qaws_mesh_reference_get_triangle(tri, r, f->patch);
		b[0] = 1 - f->a - f->b;
		b[1] = f->a;
		b[2] = f->b;
		for (k = 0; k < 3; k++)
		{
			unsigned long long key = rq_key(tri[k], tri[(k + 1) % 3]);
			if (b[(k + 2) % 3] <= TEST_BOUND(1e-9, 32) && nk && bsearch(&key, keys, nk, sizeof *keys, rq_key_cmp))
				on = 1;
			if (b[k] >= 1 - TEST_BOUND(1e-9, 32) && onf[tri[k]])
				on = 1;
		}
		crease &= on;
	}
	qaws_mesh_3d_get_topology(&t0, in);
	qaws_mesh_3d_get_topology(&t1, out);
	ok = qaws_mesh_3d_validate(out) == QAWS_STATUS_OK && t0.closed == t1.closed && t0.genus == t1.genus && foot <= TEST_BOUND(1e-9, 32) && corners && crease;
	sprintf(why, "valid %d, closed %d genus %u (input %d, %u), footpoints %.1e, corners kept %d, crease vertices on creases %d",
		qaws_mesh_3d_validate(out) == QAWS_STATUS_OK, t1.closed, t1.genus, t0.closed, t0.genus, foot, corners, crease);
	free(x); free(y); free(xi); free(l); free(fe); free(fi); free(keys); free(onf); free(tv);
	return ok;
}

/* the adaptive mesh's mean edges on the bump and far from it */
static double rq_bump_ratio(qaws_mesh_3d const* out)
{
	unsigned int nv = qaws_mesh_3d_get_vertex_count(out), ne = qaws_mesh_3d_get_edge_count(out), i, nb = 0, nfar = 0;
	qaws_scalar* p = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	unsigned int* e = (unsigned int*)malloc(sizeof(unsigned int) * 2 * ne);
	double lb = 0, lf = 0;
	qaws_mesh_3d_get_positions(p, out);
	qaws_mesh_3d_get_edges(e, out);
	for (i = 0; i < ne; i++)
	{
		qaws_scalar const *a = p + 3 * e[2 * i], *b = p + 3 * e[2 * i + 1];
		double l = sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]));
		if (a[2] > 1.25 && b[2] > 1.25)
		{
			lb += l;
			nb++;
		}
		else if (a[2] < -0.3 && b[2] < -0.3)
		{
			lf += l;
			nfar++;
		}
	}
	free(p);
	free(e);
	return nb > 10 && nfar > 10 ? (lb / nb) / (lf / nfar) : 1;
}

/* ---- the cases --------------------------------------------------------- */

typedef struct rq_shapes
{
	qaws_mesh_3d *ico, *blob, *bump, *can, *cube, *square;
	qaws_mesh_reference *rico, *rblob, *rbump, *rcan, *rcube, *rsquare;
	qaws_size_field* sf;
} rq_shapes;

static void rq_shapes_make(rq_shapes* s)
{
	qaws_scalar nrm[3 * 42];
	qaws_mesh_reference_desc rd;
	qaws_curvature_size_desc cd;
	memset(s, 0, sizeof *s);
	s->ico = rq_icosphere(nrm);
	memset(&rd, 0, sizeof rd);
	rd.normals = nrm;
	qaws_mesh_reference_create(&s->rico, s->ico, &rd);
	s->blob = rq_blob();
	qaws_mesh_reference_create(&s->rblob, s->blob, NULL);
	s->bump = rq_bumped();
	qaws_mesh_reference_create(&s->rbump, s->bump, NULL);
	memset(&cd, 0, sizeof cd);
	cd.angle_tolerance = (qaws_scalar)0.2;
	cd.min_length = (qaws_scalar)0.02;
	cd.max_length = (qaws_scalar)0.4;
	cd.gradation = (qaws_scalar)0.4;
	if (s->rbump)
		qaws_size_field_create_curvature(&s->sf, s->rbump, &cd);
	s->can = rq_can();
	qaws_mesh_reference_create(&s->rcan, s->can, NULL);
	s->cube = rq_cube();
	qaws_mesh_reference_create(&s->rcube, s->cube, NULL);
	s->square = rq_square();
	qaws_mesh_reference_create(&s->rsquare, s->square, NULL);
}

static void rq_shapes_free(rq_shapes* s)
{
	qaws_size_field_destroy(s->sf);
	qaws_mesh_reference_destroy(s->rico);
	qaws_mesh_reference_destroy(s->rblob);
	qaws_mesh_reference_destroy(s->rbump);
	qaws_mesh_reference_destroy(s->rcan);
	qaws_mesh_reference_destroy(s->rcube);
	qaws_mesh_reference_destroy(s->rsquare);
	qaws_mesh_3d_destroy(s->ico);
	qaws_mesh_3d_destroy(s->blob);
	qaws_mesh_3d_destroy(s->bump);
	qaws_mesh_3d_destroy(s->can);
	qaws_mesh_3d_destroy(s->cube);
	qaws_mesh_3d_destroy(s->square);
}

static double g_blob_seconds;   /* the default remesh of the blob, for the angle pass's timing */

static void test_default(rq_shapes const* s)
{
	qaws_mesh_3d *can, *blob, *bump;
	qaws_mesh_report rep;
	rq_meas q;
	char msg[600], inv[300];
	double sec, ratio;
	/* the can at 0.12: before this round its angles ran 27.2..117.4 */
	can = rq_remesh(s->can, s->rcan, 0.12, NULL, 0, 0, 0, &rep, &sec);
	if (!can)
	{
		TEST_ASSERT(0, "default: the can remeshed");
		return;
	}
	rq_measure(can, 0, 0, &q);
	sprintf(msg, "default, the can at 0.12: %u triangles, angles %.1f..%.1f (want [30, 110]); %.3f s", q.nf, q.amin, q.amax, sec);
	printf("  %s\n", msg);
	TEST_ASSERT(q.amin >= 30 && q.amax <= 110, msg);
	TEST_ASSERT(rq_invariants(s->can, can, s->rcan, inv), inv);
	/* the blob at 0.11: before, its mean edge 7% short, 78% of valence 6 */
	blob = rq_remesh(s->blob, s->rblob, 0.11, NULL, 0, 0, 0, &rep, &g_blob_seconds);
	if (blob)
	{
		rq_measure(blob, 0, 0, &q);
		sprintf(msg, "default, the blob at 0.11: %u triangles, mean edge %.4f (%.1f%% off), valence 6 at %.1f%%, 5..7 at %.1f%%, Delaunay %.2f%%, angles %.1f..%.1f, deviation %.2e; %.3f s",
			q.nf, q.mean_edge, 100 * (q.mean_edge / 0.11 - 1), 100 * q.v6, 100 * q.v567, 100 * q.delaunay, q.amin, q.amax, (double)rep.max_deviation, g_blob_seconds);
		printf("  %s\n", msg);
		TEST_ASSERT(fabs(q.mean_edge / 0.11 - 1) <= 0.04 && q.v6 >= 0.8 && q.v567 >= 0.99 && q.delaunay >= 0.99, msg);
		TEST_ASSERT(rq_invariants(s->blob, blob, s->rblob, inv), inv);
	}
	else
		TEST_ASSERT(0, "default: the blob remeshed");
	/* adaptive: test 124 asserts the bump's edges under 0.6 of those far from it */
	bump = s->sf ? rq_remesh(s->bump, s->rbump, 0, s->sf, 0, 0, 0, &rep, &sec) : NULL;
	if (bump)
	{
		ratio = rq_bump_ratio(bump);
		rq_measure(bump, 0, 0, &q);
		sprintf(msg, "default, adaptive bumped sphere: %u triangles, mean edge on the bump over far from it %.3f (want < 0.5), angles %.1f..%.1f; %.3f s",
			q.nf, ratio, q.amin, q.amax, sec);
		printf("  %s\n", msg);
		TEST_ASSERT(ratio < 0.5, msg);
		TEST_ASSERT(rq_invariants(s->bump, bump, s->rbump, inv), inv);
	}
	else
		TEST_ASSERT(0, "default: the bumped sphere remeshed");
	qaws_mesh_3d_destroy(can);
	qaws_mesh_3d_destroy(blob);
	qaws_mesh_3d_destroy(bump);
}

static void test_angles(rq_shapes const* s)
{
	struct
	{
		char const* name;
		qaws_mesh_3d const* m;
		qaws_mesh_reference const* r;
		double L;
		qaws_size_field const* sf;
	} const c[3] = { { "icosphere at 0.12", s->ico, s->rico, 0.12, NULL }, { "blob at 0.11", s->blob, s->rblob, 0.11, NULL },
		{ "adaptive bumped sphere", s->bump, s->rbump, 0, s->sf } };
	unsigned int i;
	char msg[600], inv[300];
	for (i = 0; i < 3; i++)
	{
		qaws_mesh_report rep;
		rq_meas q;
		double sec = 0;
		qaws_mesh_3d* out = c[i].m && c[i].r ? rq_remesh(c[i].m, c[i].r, c[i].L, c[i].sf, 35, 86, 0, &rep, &sec) : NULL;
		if (!out)
		{
			TEST_ASSERT(0, "angles: remeshed");
			continue;
		}
		rq_measure(out, 35, 86, &q);
		sprintf(msg, "angles in [35, 86], the %s: %u triangles, angles %.2f..%.2f, %.2f%% of the corners inside (%u out), valence 5..7 at %.1f%%, mean edge %.4f; %.3f s",
			c[i].name, q.nf, q.emin, q.emax, 100.0 * q.inside / q.corners, q.out_free, 100 * q.v567, q.mean_edge, sec);
		printf("  %s\n", msg);
		TEST_ASSERT(q.out_free == 0 && q.v567 >= 0.98 && (c[i].L == 0 || fabs(q.mean_edge / c[i].L - 1) <= 0.1), msg);
		TEST_ASSERT(rq_invariants(c[i].m, out, c[i].r, inv), inv);
		if (i == 1)
			printf("  the angle pass on the blob: %.3f s against %.3f s by default (%.2fx)\n", sec, g_blob_seconds, g_blob_seconds > 0 ? sec / g_blob_seconds : 0);
		qaws_mesh_3d_destroy(out);
	}
	/* the can: [30, 90] every corner not exempt, [35, 86] (heavy only, 1.5 s) 99.5% of them */
	{
		qaws_mesh_report rep;
		rq_meas q, q2;
		double sec, sec2 = 0;
		qaws_mesh_3d *a = rq_remesh(s->can, s->rcan, 0.12, NULL, 30, 90, 0, &rep, &sec);
		qaws_mesh_3d *b = test_heavy() ? rq_remesh(s->can, s->rcan, 0.12, NULL, 35, 86, 0, &rep, &sec2) : NULL;
		if (a && (b || !test_heavy()))
		{
			rq_measure(a, 30, 90, &q);
			sprintf(msg, "angles on the can at 0.12 in [30, 90]: %.2f..%.2f, %u out not exempt; %.3f s", q.emin, q.emax, q.out_free, sec);
			printf("  %s\n", msg);
			TEST_ASSERT(q.out_free == 0, msg);
			TEST_ASSERT(rq_invariants(s->can, a, s->rcan, inv), inv);
			if (b)
			{
				rq_measure(b, 35, 86, &q2);
				sprintf(msg, "angles on the can at 0.12 in [35, 86]: %.2f..%.2f, %.2f%% inside, %u out not exempt; %.3f s", q2.emin, q2.emax,
					100.0 * q2.inside / q2.corners, q2.out_free, sec2);
				printf("  %s\n", msg);
				TEST_ASSERT(q2.inside >= 0.995 * q2.corners && q2.out_free == 0, msg);
				TEST_ASSERT(rq_invariants(s->can, b, s->rcan, inv), inv);
			}
		}
		else
			TEST_ASSERT(0, "angles: the can remeshed");
		qaws_mesh_3d_destroy(a);
		qaws_mesh_3d_destroy(b);
	}
	/* a cube at several sizes (its 90-degree corner sectors need 2 faces each: at the coarse sizes the pass must add them) and
	   the open square (its border corners the same), corners and creases kept through the pass */
	{
		static double const L[5] = { 0.08, 0.15, 0.25, 0.6, 0.15 };
		for (i = 0; i < 5; i++)
		{
			if (i < 2 && !test_heavy())
				continue;   /* (the cube at 0.08 and 0.15, 5 s and 2.5 s: heavy only; 0.25 and 0.6 are the coarse sizes) */
			qaws_mesh_report rep;
			rq_meas q;
			double sec, sec0;
			qaws_mesh_3d const* in = i < 4 ? s->cube : s->square;
			qaws_mesh_reference const* r = i < 4 ? s->rcube : s->rsquare;
			qaws_mesh_3d *a = rq_remesh(in, r, L[i], NULL, 35, 86, 0, &rep, &sec), *a0 = rq_remesh(in, r, L[i], NULL, 0, 0, 0, NULL, &sec0);
			if (a && a0)
			{
				rq_measure(a, 35, 86, &q);
				sprintf(msg, "angles on the %s at %.2f in [35, 86]: %.2f..%.2f, %u out not exempt; %.3f s (%.1fx the default remesh)", i < 4 ? "cube" : "open square",
					L[i], q.emin, q.emax, q.out_free, sec, sec0 > 0 ? sec / sec0 : 0);
				printf("  %s\n", msg);
				TEST_ASSERT(q.out_free == 0, msg);
				TEST_ASSERT(rq_invariants(in, a, r, inv), inv);
			}
			else
				TEST_ASSERT(0, "angles: the cube and the square remeshed");
			qaws_mesh_3d_destroy(a);
			qaws_mesh_3d_destroy(a0);
		}
	}
}

static void test_tolerance(rq_shapes const* s)
{
	qaws_mesh_report rep, rep0;
	double sec, d0, d1;
	/* (by default 5e-3, the guard's cost growing as the tolerance shrinks, and the default remesh to compare with heavy only: its
	   sampled deviation 1.5 s; measured 4.99e-3 against 6.63e-3, 1587 vertices against 1581, in 2 s against 5) */
	double const tol = test_heavy() ? 2e-3 : 5e-3;
	qaws_mesh_3d *a = test_heavy() ? rq_remesh(s->blob, s->rblob, 0.11, NULL, 0, 0, 0, &rep0, NULL) : NULL;
	qaws_mesh_3d *b = rq_remesh(s->blob, s->rblob, 0.11, NULL, 0, 0, tol, &rep, &sec);
	char msg[400], inv[300];
	if ((!a && test_heavy()) || !b)
	{
		TEST_ASSERT(0, "tolerance: remeshed");
		qaws_mesh_3d_destroy(a);
		qaws_mesh_3d_destroy(b);
		return;
	}
	d0 = a ? rq_deviation(a, s->rblob) : 0;
	d1 = rq_deviation(b, s->rblob);
	sprintf(msg, "tolerance %.0e on the blob at 0.11: sampled deviation %.3e (by default %.3e), %u vertices against %u (+%.0f%%), min angle %.1f; %.3f s",
		tol, d1, d0, qaws_mesh_3d_get_vertex_count(b), a ? qaws_mesh_3d_get_vertex_count(a) : 0,
		a ? 100.0 * qaws_mesh_3d_get_vertex_count(b) / qaws_mesh_3d_get_vertex_count(a) - 100 : 0, (double)rep.quality.min_angle, sec);
	printf("  %s\n", msg);
	TEST_ASSERT(d1 <= tol + 1e-9, msg);
	TEST_ASSERT(rq_invariants(s->blob, b, s->rblob, inv), inv);
	qaws_mesh_3d_destroy(a);
	qaws_mesh_3d_destroy(b);
	/* smaller tolerances, and the tolerance with the angle bounds: every sample within it, no sliver (none past 170 degrees); heavy
	   only (9 s, 46 s, 6 s and 7 s) */
	if (test_heavy())
	{
		struct
		{
			char const* name;
			qaws_mesh_3d const* m;
			qaws_mesh_reference const* r;
			double L, tol, lo, hi;
			qaws_size_field const* sf;
		} const c[4] = { { "blob", s->blob, s->rblob, 0.11, 5e-4, 0, 0, NULL }, { "blob", s->blob, s->rblob, 0.11, 2e-4, 0, 0, NULL },
			{ "blob", s->blob, s->rblob, 0.11, 2e-3, 35, 86, NULL }, { "adaptive bumped sphere", s->bump, s->rbump, 0, 2e-3, 35, 86, s->sf } };
		unsigned int i;
		for (i = 0; i < 4; i++)
		{
			rq_meas q;
			qaws_mesh_3d* o = rq_remesh(c[i].m, c[i].r, c[i].L, c[i].sf, c[i].lo, c[i].hi, c[i].tol, &rep, &sec);
			if (!o)
			{
				TEST_ASSERT(0, "tolerance: remeshed");
				continue;
			}
			d1 = rq_deviation(o, c[i].r);
			rq_measure(o, c[i].lo, c[i].hi, &q);
			sprintf(msg, "tolerance %.0e on the %s%s: sampled deviation %.3e, %u vertices, angles %.2f..%.2f (%u out of the bounds not exempt); %.3f s",
				c[i].tol, c[i].name, c[i].lo > 0 ? " with [35, 86]" : "", d1, q.nv, q.amin, q.amax, q.out_free, sec);
			printf("  %s\n", msg);
			TEST_ASSERT(d1 <= c[i].tol + 1e-9 && q.amax <= 170, msg);
			TEST_ASSERT(rq_invariants(c[i].m, o, c[i].r, inv), inv);
			qaws_mesh_3d_destroy(o);
		}
	}
}


/* ---- open generated meshes --------------------------------------------- */

/* the showcase's wavy patch (figure 5, 0), a bilinear square (1) and a quarter-cylinder piece (2), meshed by qaws_mesh_generate_surface */
static qaws_mesh_3d* rq_generated(unsigned int which, double L)
{
	qaws_vec3 cp[49];
	qaws_surface_bspline_desc d;
	qaws_surface_bilinear_desc b;
	qaws_generate_desc g;
	qaws_surface* srf = NULL;
	qaws_mesh_3d* out = NULL;
	unsigned int i, j;
	memset(&d, 0, sizeof d);
	if (which == 0)
	{
		for (i = 0; i < 7; i++)
			for (j = 0; j < 7; j++)
			{
				cp[i * 7 + j].x = (qaws_scalar)(i / 3.0 - 1);
				cp[i * 7 + j].y = (qaws_scalar)(j / 3.0 - 1);
				cp[i * 7 + j].z = (qaws_scalar)(0.55 * sin(1.7 * i) * cos(1.3 * j) + (i == 3 && j == 3 ? 0.6 : 0));
			}
		d.u_degree = d.v_degree = 3;
		d.u_point_count = d.v_point_count = 7;
	}
	else if (which == 2)
	{
		static double const cx[3] = { 1, 1, 0 }, cy[3] = { 0, 1, 1 };
		for (i = 0; i < 3; i++)
			for (j = 0; j < 2; j++)
			{
				cp[i * 2 + j].x = (qaws_scalar)cx[i];
				cp[i * 2 + j].y = (qaws_scalar)cy[i];
				cp[i * 2 + j].z = (qaws_scalar)(1.5 * j);
			}
		d.u_degree = 2;
		d.v_degree = 1;
		d.u_point_count = 3;
		d.v_point_count = 2;
	}
	d.control_points = cp;
	memset(&b, 0, sizeof b);
	b.p10.x = b.p11.x = b.p01.y = b.p11.y = 2;
	if ((which == 1 ? qaws_surface_create_bilinear(&srf, &b) : qaws_surface_create_bspline(&srf, &d)) != QAWS_STATUS_OK)
		return NULL;
	memset(&g, 0, sizeof g);
	g.elements = QAWS_ELEMENTS_TRIANGLES;
	g.target_length = (qaws_scalar)L;
	if (qaws_mesh_generate_surface(&out, NULL, srf, &g) != QAWS_STATUS_OK)
		out = NULL;
	qaws_surface_destroy(srf);
	return out;
}

/* the area and the border's length (its edges in one triangle) */
static void rq_area_border(qaws_mesh_3d const* m, double* area, double* border)
{
	unsigned int nv = qaws_mesh_3d_get_vertex_count(m), nf = qaws_mesh_3d_get_triangle_count(m), i, k, n;
	qaws_scalar* p = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 3 * nf);
	unsigned long long* keys = (unsigned long long*)malloc(sizeof(unsigned long long) * 3 * nf);
	qaws_mesh_3d_get_positions(p, m);
	qaws_mesh_3d_get_triangles(t, m);
	*area = *border = 0;
	for (i = 0; i < nf; i++)
	{
		qaws_scalar const *a = p + 3 * t[3 * i], *b = p + 3 * t[3 * i + 1], *c = p + 3 * t[3 * i + 2];
		double u[3], w[3], x, y, z;
		for (k = 0; k < 3; k++)
		{
			u[k] = b[k] - a[k];
			w[k] = c[k] - a[k];
			keys[3 * i + k] = rq_key(t[3 * i + k], t[3 * i + (k + 1) % 3]);
		}
		x = u[1] * w[2] - u[2] * w[1];
		y = u[2] * w[0] - u[0] * w[2];
		z = u[0] * w[1] - u[1] * w[0];
		*area += 0.5 * sqrt(x * x + y * y + z * z);
	}
	n = 3 * nf;
	qsort(keys, n, sizeof *keys, rq_key_cmp);
	for (i = 0; i < n; i++)
	{
		unsigned int a, b;
		if ((i > 0 && keys[i - 1] == keys[i]) || (i + 1 < n && keys[i + 1] == keys[i]))
			continue;
		a = (unsigned int)(keys[i] >> 32);
		b = (unsigned int)(keys[i] & 0xFFFFFFFFu);
		*border += sqrt((p[3 * a] - p[3 * b]) * (p[3 * a] - p[3 * b]) + (p[3 * a + 1] - p[3 * b + 1]) * (p[3 * a + 1] - p[3 * b + 1])
			+ (p[3 * a + 2] - p[3 * b + 2]) * (p[3 * a + 2] - p[3 * b + 2]));
	}
	free(p); free(t); free(keys);
}

/* LOCAL on the meshes qaws_mesh_generate_surface makes of open patches: their area and border kept, faces near area / (0.43 L^2) */
static void test_open_generated(void)
{
	static char const* const name[3] = { "wavy patch", "bilinear square", "cylinder piece" };
	static double const sizes[2] = { 0.1, 0.06 };
	unsigned int w, k;
	for (w = 0; w < 3; w++)
	{
		qaws_mesh_3d* g = rq_generated(w, 0.1);
		qaws_mesh_reference* r = NULL;
		double a0, b0;
		if (!g || qaws_mesh_reference_create(&r, g, NULL) != QAWS_STATUS_OK)
		{
			TEST_ASSERT(0, "open generated: the patch meshed");
			qaws_mesh_3d_destroy(g);
			continue;
		}
		rq_area_border(g, &a0, &b0);
		for (k = 0; k < 2; k++)
		{
			qaws_mesh_report rep;
			qaws_mesh_3d* out = rq_remesh(g, r, sizes[k], NULL, 0, 0, 0, &rep, NULL);
			char msg[400], inv[300];
			double a, b, want = a0 / (0.43 * sizes[k] * sizes[k]);
			unsigned int nf;
			if (!out)
			{
				TEST_ASSERT(0, "open generated: remeshed");
				continue;
			}
			rq_area_border(out, &a, &b);
			nf = qaws_mesh_3d_get_triangle_count(out);
			sprintf(msg, "open generated, %s (%u triangles) at %.2f: %u triangles (about %.0f), area %.4f (%.2f%% off), border %.4f (%.2f%% off)",
				name[w], qaws_mesh_3d_get_triangle_count(g), sizes[k], nf, want, a, 100 * (a / a0 - 1), b, 100 * (b / b0 - 1));
			printf("  %s\n", msg);
			TEST_ASSERT(fabs(a / a0 - 1) <= 0.02 && fabs(b / b0 - 1) <= 0.02 && nf > 0.75 * want && nf < 1.33 * want, msg);
			TEST_ASSERT(rq_invariants(g, out, r, inv), inv);
			qaws_mesh_3d_destroy(out);
		}
		qaws_mesh_reference_destroy(r);
		qaws_mesh_3d_destroy(g);
	}
}

static void test_refusals(rq_shapes const* s)
{
	static double const bad[5][3] = { { 0.7, 0.6, 0 }, { -0.1, 0, 0 }, { M_PI / 3, 0, 0 }, { 0, 1.0, 0 }, { 0, 0, -1 } };
	qaws_status st[5];
	qaws_mesh_3d* out = NULL;
	unsigned int i, ok = 1;
	char msg[400];
	for (i = 0; i < 5; i++)
	{
		qaws_remesh_desc d;
		memset(&d, 0, sizeof d);
		d.min_angle = (qaws_scalar)bad[i][0];
		d.max_angle = (qaws_scalar)bad[i][1];
		d.tolerance = (qaws_scalar)bad[i][2];
		st[i] = qaws_mesh_remesh(&out, NULL, s->ico, &d);
		ok &= st[i] == QAWS_STATUS_INVALID_ARGUMENT && !out;
	}
	sprintf(msg, "refusals: min 0.7 over max 0.6 %s, min -0.1 %s, min pi/3 %s, max 1.0 %s, tolerance -1 %s", qaws_status_to_string(st[0]),
		qaws_status_to_string(st[1]), qaws_status_to_string(st[2]), qaws_status_to_string(st[3]), qaws_status_to_string(st[4]));
	printf("  %s\n", msg);
	TEST_ASSERT(ok, msg);
}

int test_131_remesh_quality_main(void)
{
	rq_shapes s;
	g_pass = 0;
	g_fail = 0;
	printf("Test 131: remeshing quality\n");
	rq_shapes_make(&s);
	if (!s.rico || !s.rblob || !s.rbump || !s.rcan || !s.rcube || !s.rsquare || !s.sf)
		TEST_ASSERT(0, "shapes and references made");
	else
	{
		test_default(&s);
		test_angles(&s);
		test_tolerance(&s);
		test_open_generated();
		test_refusals(&s);
	}
	rq_shapes_free(&s);
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
