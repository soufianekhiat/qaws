/*
 * Test 130: quality measures (qaws_mesh_measure.h)
 *
 *   - analytic shapes: an equilateral triangle (Q 1, angles 60), a 30-60-90
 *     one (Q (3 - sqrt 3) / 2), a unit square (Jacobian 1, AD and ED 0), a
 *     60-degree rhombus (Jacobian sin 60, AD 30);
 *   - ideal valences: a creased quad cube is all regular (against 8 from
 *     qaws_mesh_quality), its corners want 3 and its edges 4; an icosphere
 *     has its 12 valence-5 vertices; an open quad grid is all regular, its
 *     corners want 2 and its border 3;
 *   - deviation, exact: a flat grid against itself lifted by 0.01, both
 *     ways 0.01;
 *   - deviation against brute force: a flat square against a PN grid with a
 *     raised centre, each direction's bounds bracketing dense sampling;
 *   - two-sidedness: a bump of 0.2 a flat square does not see from above;
 *   - a PN saddle against its own triangles (the closest place jumping
 *     across concave edges, the patch projection stalling on sides), its
 *     witnesses checked apart from the library; a pentagon fanned into a
 *     zero-area triangle;
 *   - consistency with the remesher's one-sided report;
 *   - refusals, and the timing on a 10k-face remesh (QAWS_TEST_HEAVY=1 only).
 */

#include "test_common.h"
#include "qaws_remesh.h"
#include <time.h>

static qaws_mesh_3d* mm_make(double const* x, unsigned int nv, unsigned int const* f, unsigned int const* fs, unsigned int nf)
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
	qaws_mesh_3d_create(&m, &d);
	free(p);
	return m;
}

/* an n x n grid on [0, 1]^2 at height z(x, y): triangles (quads 0) or quads */
static qaws_mesh_3d* mm_grid(unsigned int n, int quads, double (*z)(double, double))
{
	unsigned int nv = (n + 1) * (n + 1), i, j, k = 0, *f = (unsigned int*)malloc(sizeof(unsigned int) * 6 * n * n), *fs = NULL;
	double* x = (double*)malloc(sizeof(double) * 3 * nv);
	qaws_mesh_3d* m;
	for (j = 0; j <= n; j++)
		for (i = 0; i <= n; i++)
		{
			double u = (double)i / n, v = (double)j / n;
			x[3 * (j * (n + 1) + i)] = u;
			x[3 * (j * (n + 1) + i) + 1] = v;
			x[3 * (j * (n + 1) + i) + 2] = z ? z(u, v) : 0;
		}
	if (quads)
		fs = (unsigned int*)malloc(sizeof(unsigned int) * n * n);
	for (j = 0; j < n; j++)
		for (i = 0; i < n; i++)
		{
			unsigned int a = j * (n + 1) + i, b = a + 1, c = a + n + 2, d = a + n + 1;
			if (quads)
			{
				fs[j * n + i] = 4;
				f[k++] = a; f[k++] = b; f[k++] = c; f[k++] = d;
			}
			else
			{
				f[k++] = a; f[k++] = b; f[k++] = c;
				f[k++] = a; f[k++] = c; f[k++] = d;
			}
		}
	m = mm_make(x, nv, f, fs, quads ? n * n : 2 * n * n);
	free(x);
	free(f);
	free(fs);
	return m;
}

/* the icosahedron split once, on the unit sphere (test 124's) */
static qaws_mesh_3d* mm_icosphere(void)
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
			x[3 * i + k] /= r;
	}
	return mm_make(x, nv, f, NULL, nf);
}

/* test 124's bumped icosphere: a bump of 0.6 at the north pole */
static qaws_mesh_3d* mm_bumped(void)
{
	qaws_mesh_3d* m = mm_icosphere();
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

/* the closest point of triangle (A, B, C) to p (Ericson), its squared distance */
static double mm_closest(double const* A, double const* B, double const* C, double const* p)
{
	double ab[3], ac[3], ap[3], bp[3], cp[3], d1, d2, d3, d4, d5, d6, va, vb, vc, v = 0, w = 0, x[3], s = 0;
	int k;
	for (k = 0; k < 3; k++)
	{
		ab[k] = B[k] - A[k]; ac[k] = C[k] - A[k];
		ap[k] = p[k] - A[k]; bp[k] = p[k] - B[k]; cp[k] = p[k] - C[k];
	}
#define MM_DOT(a, b) ((a)[0] * (b)[0] + (a)[1] * (b)[1] + (a)[2] * (b)[2])
	d1 = MM_DOT(ab, ap); d2 = MM_DOT(ac, ap); d3 = MM_DOT(ab, bp); d4 = MM_DOT(ac, bp); d5 = MM_DOT(ab, cp); d6 = MM_DOT(ac, cp);
#undef MM_DOT
	vc = d1 * d4 - d3 * d2; vb = d5 * d2 - d1 * d6; va = d3 * d6 - d5 * d4;
	if (d1 <= 0 && d2 <= 0) { v = 0; w = 0; }
	else if (d3 >= 0 && d4 <= d3) { v = 1; w = 0; }
	else if (vc <= 0 && d1 >= 0 && d3 <= 0) { v = d1 / (d1 - d3); w = 0; }
	else if (d6 >= 0 && d5 <= d6) { v = 0; w = 1; }
	else if (vb <= 0 && d2 >= 0 && d6 <= 0) { v = 0; w = d2 / (d2 - d6); }
	else if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0) { w = (d4 - d3) / ((d4 - d3) + (d5 - d6)); v = 1 - w; }
	else { double den = 1 / (va + vb + vc); v = vb * den; w = vc * den; }
	for (k = 0; k < 3; k++)
	{
		x[k] = A[k] + v * ab[k] + w * ac[k] - p[k];
		s += x[k] * x[k];
	}
	return s;
}

static void test_shapes(void)
{
	double const s3 = sqrt(3.0);
	double const eq[9] = { 0, 0, 0, 1, 0, 0, 0.5, s3 / 2, 0 }, half[9] = { 0, 0, 0, s3, 0, 0, 0, 1, 0 };
	double const sq[12] = { 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0 }, rh[12] = { 0, 0, 0, 1, 0, 0, 1.5, s3 / 2, 0, 0.5, s3 / 2, 0 };
	unsigned int const t[3] = { 0, 1, 2 }, q[4] = { 0, 1, 2, 3 }, four = 4;
	qaws_mesh_3d *a = mm_make(eq, 3, t, NULL, 1), *b = mm_make(half, 3, t, NULL, 1), *c = mm_make(sq, 4, q, &four, 1), *d = mm_make(rh, 4, q, &four, 1);
	qaws_shape_quality s;
	char msg[300];
	qaws_mesh_3d_measure_shape(&s, NULL, a);
	sprintf(msg, "equilateral: Q %.15g / %.15g, under 30 %g, mean smallest angle %.15g", (double)s.min_quality, (double)s.mean_quality,
		(double)s.below_30, (double)s.mean_min_angle);
	TEST_ASSERT(fabs(s.min_quality - 1) < 1e-12 && fabs(s.mean_quality - 1) < 1e-12 && s.below_30 == 0 && fabs(s.mean_min_angle - 60) < 1e-12, msg);
	qaws_mesh_3d_measure_shape(&s, NULL, b);
	sprintf(msg, "30-60-90: Q %.15g (want %.15g), angles %.12g..%.12g, a share over 90 %g", (double)s.min_quality, (3 - s3) / 2,
		(double)s.min_angle, (double)s.max_angle, (double)s.above_90);
	TEST_ASSERT(fabs(s.min_quality - (3 - s3) / 2) < TEST_BOUND(1e-12, 64) && fabs(s.min_angle - 30) < TEST_BOUND(1e-12, 64 * 90) && fabs(s.max_angle - 90) < TEST_BOUND(1e-12, 64 * 90), msg);
	qaws_mesh_3d_measure_shape(&s, NULL, c);
	sprintf(msg, "unit square: Jacobian %.15g / %.15g, AD %g, ED %g", (double)s.min_jacobian, (double)s.mean_jacobian, (double)s.angle_deviation,
		(double)s.edge_deviation);
	TEST_ASSERT(fabs(s.min_jacobian - 1) < 1e-12 && fabs(s.mean_jacobian - 1) < 1e-12 && fabs(s.angle_deviation) < 1e-12 && fabs(s.edge_deviation) < 1e-12, msg);
	qaws_mesh_3d_measure_shape(&s, NULL, d);
	sprintf(msg, "60-degree rhombus: Jacobian %.15g (want %.15g), AD %.15g", (double)s.min_jacobian, s3 / 2, (double)s.angle_deviation);
	TEST_ASSERT(fabs(s.min_jacobian - s3 / 2) < TEST_BOUND(1e-12, 64) && fabs(s.angle_deviation - 30) < TEST_BOUND(1e-12, 64 * 90), msg);
	qaws_mesh_3d_destroy(a);
	qaws_mesh_3d_destroy(b);
	qaws_mesh_3d_destroy(c);
	qaws_mesh_3d_destroy(d);
}

/* a closed cube of n x n quads per face on [-1, 1]^3; out_corner, out_edge: a corner vertex, a vertex inside an edge */
static qaws_mesh_3d* mm_cube(unsigned int n, unsigned int* out_corner, unsigned int* out_edge)
{
	/* per face: the lattice origin and the two axes (u x v outward) */
	static int const face[6][3][3] = { { { 0, 0, 0 }, { 0, 0, 1 }, { 0, 1, 0 } }, { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } },
		{ { 0, 0, 0 }, { 1, 0, 0 }, { 0, 0, 1 } }, { { 0, 1, 0 }, { 0, 0, 1 }, { 1, 0, 0 } }, { { 0, 0, 0 }, { 0, 1, 0 }, { 1, 0, 0 } },
		{ { 0, 0, 1 }, { 1, 0, 0 }, { 0, 1, 0 } } };
	unsigned int s = n + 1, *id = (unsigned int*)malloc(sizeof(unsigned int) * s * s * s), *f = (unsigned int*)malloc(sizeof(unsigned int) * 24 * n * n),
		*fs = (unsigned int*)malloc(sizeof(unsigned int) * 6 * n * n), nv = 0, i, j, k, a, b, nf = 0;
	double* x = (double*)malloc(sizeof(double) * 3 * s * s * s);
	qaws_mesh_3d* m;
	for (k = 0; k < s; k++)
		for (j = 0; j < s; j++)
			for (i = 0; i < s; i++)
			{
				id[(k * s + j) * s + i] = 0xFFFFFFFFu;
				if (i == 0 || j == 0 || k == 0 || i == n || j == n || k == n)
				{
					x[3 * nv] = 2.0 * i / n - 1;
					x[3 * nv + 1] = 2.0 * j / n - 1;
					x[3 * nv + 2] = 2.0 * k / n - 1;
					id[(k * s + j) * s + i] = nv++;
				}
			}
	for (i = 0; i < 6; i++)
		for (b = 0; b < n; b++)
			for (a = 0; a < n; a++)
			{
				static int const da[4] = { 0, 1, 1, 0 }, db[4] = { 0, 0, 1, 1 };
				unsigned int c;
				for (c = 0; c < 4; c++)
				{
					int p[3], z;
					for (z = 0; z < 3; z++)
						p[z] = face[i][0][z] * (int)n + (int)(a + da[c]) * face[i][1][z] + (int)(b + db[c]) * face[i][2][z];
					f[4 * nf + c] = id[((unsigned int)p[2] * s + (unsigned int)p[1]) * s + (unsigned int)p[0]];
				}
				fs[nf++] = 4;
			}
	*out_corner = id[0];
	*out_edge = id[1];
	m = mm_make(x, nv, f, fs, nf);
	free(id); free(f); free(fs); free(x);
	return m;
}

static void test_valences(void)
{
	qaws_shape_quality s;
	qaws_mesh_quality mq;
	unsigned int corner, edge, *ideal;
	char msg[300];
	/* the creased cube: its corners and edges regular */
	{
		qaws_mesh_3d* m = mm_cube(3, &corner, &edge);
		ideal = (unsigned int*)malloc(sizeof(unsigned int) * qaws_mesh_3d_get_vertex_count(m));
		qaws_mesh_3d_detect_features(m, (qaws_scalar)(30 * M_PI / 180), (qaws_scalar)(60 * M_PI / 180), 0);
		qaws_mesh_3d_measure_shape(&s, ideal, m);
		qaws_mesh_3d_measure_quality(&mq, m);
		sprintf(msg, "creased quad cube (%u vertices, %u quads): irregular %u (qaws_mesh_quality: %u), ideal %u at a corner, %u on an edge, excess %u, share %g",
			qaws_mesh_3d_get_vertex_count(m), qaws_mesh_3d_get_face_count(m), s.irregular, mq.irregular, ideal[corner], ideal[edge], s.valence_excess,
			(double)s.regular_share);
		TEST_ASSERT(s.irregular == 0 && mq.irregular == 8 && ideal[corner] == 3 && ideal[edge] == 4 && s.valence_excess == 0 && s.regular_share == 1, msg);
		free(ideal);
		qaws_mesh_3d_destroy(m);
	}
	/* the icosphere: its 12 valence-5 vertices */
	{
		qaws_mesh_3d* m = mm_icosphere();
		qaws_mesh_3d_measure_shape(&s, NULL, m);
		sprintf(msg, "icosphere (%u vertices): irregular %u, V567 %g", qaws_mesh_3d_get_vertex_count(m), s.irregular, (double)s.regular_share);
		TEST_ASSERT(qaws_mesh_3d_get_vertex_count(m) == 42 && s.irregular == 12 && s.regular_share == 1, msg);
		qaws_mesh_3d_destroy(m);
	}
	/* the open quad grid: its border regular */
	{
		qaws_mesh_3d* m = mm_grid(10, 1, NULL);
		ideal = (unsigned int*)malloc(sizeof(unsigned int) * qaws_mesh_3d_get_vertex_count(m));
		qaws_mesh_3d_measure_shape(&s, ideal, m);
		sprintf(msg, "open 10 x 10 quad grid: irregular %u, ideal %u at a corner, %u on the border, %u inside", s.irregular, ideal[0], ideal[5], ideal[5 * 11 + 5]);
		TEST_ASSERT(s.irregular == 0 && ideal[0] == 2 && ideal[5] == 3 && ideal[5 * 11 + 5] == 4, msg);
		free(ideal);
		qaws_mesh_3d_destroy(m);
	}
}

static double mm_lift(double x, double y)
{
	(void)x;
	(void)y;
	return 0.01;
}

/* the 8 x 8 grid's centre raised 0.3 (its vertex at (0.5, 0.5)) */
static double mm_peak(double x, double y)
{
	return fabs(x - 0.5) < 1e-9 && fabs(y - 0.5) < 1e-9 ? 0.3 : 0;
}

static double mm_bump(double x, double y)
{
	return 0.2 * exp(-((x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5)) / 0.02);
}

static void mm_print(char const* name, qaws_mesh_deviation const* d)
{
	printf("  %s: to the reference [%.9g, %.9g], from it [%.9g, %.9g], rms %.3g, diagonal %.4g, converged %d\n", name, (double)d->to_reference[0],
		(double)d->to_reference[1], (double)d->from_reference[0], (double)d->from_reference[1], (double)d->rms, (double)d->diagonal, d->converged);
}

static void test_exact(void)
{
	qaws_mesh_3d *flat = mm_grid(10, 0, NULL), *up = mm_grid(10, 0, mm_lift);
	qaws_mesh_reference_desc rd;
	qaws_mesh_reference* r = NULL;
	qaws_mesh_deviation d;
	qaws_status st;
	double tol;
	char msg[300];
	memset(&rd, 0, sizeof rd);
	rd.interpolant = QAWS_INTERPOLANT_LINEAR;
	qaws_mesh_reference_create(&r, flat, &rd);
	st = qaws_mesh_3d_measure_deviation(&d, r, up, NULL);
	mm_print("flat grid lifted by 0.01", &d);
	tol = 1e-4 * d.diagonal;
	sprintf(msg, "a flat grid lifted by 0.01 against its LINEAR reference: to [%.12g, %.12g], from [%.12g, %.12g], converged %d", (double)d.to_reference[0],
		(double)d.to_reference[1], (double)d.from_reference[0], (double)d.from_reference[1], d.converged);
	TEST_ASSERT(st == QAWS_STATUS_OK && d.converged == 1 && d.to_reference[0] >= 0.01 - TEST_BOUND(1e-9, 32) && d.to_reference[1] <= 0.01 + tol
		&& d.from_reference[0] >= 0.01 - TEST_BOUND(1e-9, 32) && d.from_reference[1] <= 0.01 + tol && fabs(d.rms - 0.01) < TEST_BOUND(1e-9, 32), msg);
	qaws_mesh_reference_destroy(r);
	qaws_mesh_3d_destroy(flat);
	qaws_mesh_3d_destroy(up);
}

static void test_brute(void)
{
	static double const sq[12] = { 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0 };
	static unsigned int const tri[6] = { 0, 1, 2, 0, 2, 3 };
	qaws_mesh_3d *g = mm_grid(8, 0, mm_peak), *m = mm_make(sq, 4, tri, NULL, 2);
	qaws_mesh_reference* r = NULL;
	qaws_mesh_deviation d;
	qaws_status st;
	unsigned int n = 200, f, i, j, k, np, nt;
	double to = 0, from = 0, tol;
	qaws_scalar *pts, *dist;
	qaws_mesh_location* loc;
	char msg[400];
	qaws_mesh_reference_create(&r, g, NULL);
	st = qaws_mesh_3d_measure_deviation(&d, r, m, NULL);
	mm_print("flat square against a raised PN grid", &d);
	tol = 1e-4 * d.diagonal;
	/* the square sampled 200 x 200 per face, projected unhinted */
	np = (n + 1) * (n + 2) / 2;
	pts = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * np);
	dist = (qaws_scalar*)malloc(sizeof(qaws_scalar) * np);
	loc = (qaws_mesh_location*)malloc(sizeof(qaws_mesh_location) * (np > 51 * 52 / 2 ? np : 51 * 52 / 2));
	for (f = 0; f < 2; f++)
	{
		double const *A = sq + 3 * tri[3 * f], *B = sq + 3 * tri[3 * f + 1], *C = sq + 3 * tri[3 * f + 2];
		unsigned int q = 0;
		for (j = 0; j <= n; j++)
			for (i = 0; i + j <= n; i++, q++)
				for (k = 0; k < 3; k++)
					pts[3 * q + k] = (qaws_scalar)(A[k] + (B[k] - A[k]) * i / n + (C[k] - A[k]) * j / n);
		qaws_mesh_reference_project(loc, dist, r, pts, NULL, np);
		for (q = 0; q < np; q++)
			to = dist[q] > to ? dist[q] : to;
	}
	/* the reference sampled 50 x 50 per patch, against both triangles */
	nt = qaws_mesh_reference_get_triangle_count(r);
	n = 50;
	np = (n + 1) * (n + 2) / 2;
	for (f = 0; f < nt; f++)
	{
		unsigned int q = 0;
		for (j = 0; j <= n; j++)
			for (i = 0; i + j <= n; i++, q++)
			{
				loc[q].patch = f;
				loc[q].a = (qaws_scalar)((double)i / n);
				loc[q].b = (qaws_scalar)((double)j / n);
			}
		qaws_mesh_reference_evaluate(pts, NULL, NULL, r, loc, np);
		for (q = 0; q < np; q++)
		{
			double p[3], e;
			for (k = 0; k < 3; k++)
				p[k] = pts[3 * q + k];
			e = fmin(mm_closest(sq, sq + 3, sq + 6, p), mm_closest(sq, sq + 6, sq + 9, p));
			from = sqrt(e) > from ? sqrt(e) : from;
		}
	}
	sprintf(msg, "mesh to a raised PN grid: brute force %.12g in [%.12g, %.12g], gap %.3g (tol %.3g)", to, (double)d.to_reference[0], (double)d.to_reference[1],
		(double)(d.to_reference[1] - d.to_reference[0]), tol);
	TEST_ASSERT(st == QAWS_STATUS_OK && to >= d.to_reference[0] - 2 * tol && to <= d.to_reference[1] + 1e-12 && d.to_reference[1] - d.to_reference[0] <= tol, msg);
	sprintf(msg, "the raised PN grid to the mesh: brute force %.12g in [%.12g, %.12g], gap %.3g, converged %d", from, (double)d.from_reference[0],
		(double)d.from_reference[1], (double)(d.from_reference[1] - d.from_reference[0]), d.converged);
	TEST_ASSERT(st == QAWS_STATUS_OK && from >= d.from_reference[0] - 1e-12 && from <= d.from_reference[1] + 1e-12
		&& d.from_reference[1] - d.from_reference[0] <= tol && d.converged == 1, msg);
	free(pts);
	free(dist);
	free(loc);
	qaws_mesh_reference_destroy(r);
	qaws_mesh_3d_destroy(g);
	qaws_mesh_3d_destroy(m);
}

static void test_two_sided(void)
{
	static double const sq[12] = { 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0 };
	static unsigned int const tri[6] = { 0, 1, 2, 0, 2, 3 };
	qaws_mesh_3d *g = mm_grid(16, 0, mm_bump), *m = mm_make(sq, 4, tri, NULL, 2);
	qaws_mesh_reference* r = NULL;
	qaws_mesh_deviation d;
	qaws_status st;
	char msg[300];
	qaws_mesh_reference_create(&r, g, NULL);
	st = qaws_mesh_3d_measure_deviation(&d, r, m, NULL);
	mm_print("flat square under a bump of 0.2", &d);
	sprintf(msg, "a bump of 0.2 over a flat square: from the reference [%.9g, %.9g] (to it only %.9g), Hausdorff [%.9g, %.9g]", (double)d.from_reference[0],
		(double)d.from_reference[1], (double)d.to_reference[1], (double)d.hausdorff[0], (double)d.hausdorff[1]);
	TEST_ASSERT(st == QAWS_STATUS_OK && d.from_reference[0] >= 0.2 - 1e-4 * d.diagonal && d.to_reference[1] < 0.8 * d.from_reference[0] && d.hausdorff[0] == d.from_reference[0], msg);
	qaws_mesh_reference_destroy(r);
	qaws_mesh_3d_destroy(g);
	qaws_mesh_3d_destroy(m);
}

static double mm_saddle(double x, double y)
{
	return 0.25 * sin(7 * x) * sin(7 * y);
}

/* the distance from p to a reference by sampling every patch n x n: at least the true one, apart from the library's projection */
static double mm_dense(qaws_mesh_reference const* r, double const* p, unsigned int n)
{
	unsigned int nt = qaws_mesh_reference_get_triangle_count(r), np = (n + 1) * (n + 2) / 2, f, i, j, q, k;
	qaws_mesh_location* loc = (qaws_mesh_location*)malloc(sizeof(qaws_mesh_location) * np);
	qaws_scalar* pts = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * np);
	double best = HUGE_VAL;
	for (f = 0; f < nt; f++)
	{
		for (j = 0, q = 0; j <= n; j++)
			for (i = 0; i + j <= n; i++, q++)
			{
				loc[q].patch = f;
				loc[q].a = (qaws_scalar)((double)i / n);
				loc[q].b = (qaws_scalar)((double)j / n);
			}
		qaws_mesh_reference_evaluate(pts, NULL, NULL, r, loc, np);
		for (q = 0; q < np; q++)
		{
			double e = 0;
			for (k = 0; k < 3; k++)
				e += (pts[3 * q + k] - p[k]) * (pts[3 * q + k] - p[k]);
			best = e < best ? e : best;
		}
	}
	free(loc);
	free(pts);
	return sqrt(best);
}

/* A PN saddle against its own triangles: the reference passes over its mesh's valleys, where the closest place on the mesh jumps
   across a concave edge (the matched cubic maps do not shrink there, the Lipschitz bound through the net does), and its patch
   projection stalls on the patches' sides (the lower bound polished past it). Both witnesses checked apart from the library:
   the one on the mesh by sampling the reference densely, the one on the reference exactly against the mesh's triangles; the
   reference sampled against the mesh under the upper bound. */
static void test_saddle(void)
{
	qaws_mesh_3d* g = mm_grid(12, 0, mm_saddle);
	qaws_mesh_reference *r = NULL, *lin = NULL;
	qaws_mesh_reference_desc rd;
	qaws_mesh_deviation d;
	qaws_mesh_location l, *loc;
	qaws_scalar at[3], dist, *pts, *dd;
	qaws_status st;
	unsigned int nt, n = 24, np = (n + 1) * (n + 2) / 2, f, i, j, q, k;
	double tol, wt[3], at_to, at_from, sup = 0;
	char msg[500];
	memset(&rd, 0, sizeof rd);
	rd.interpolant = QAWS_INTERPOLANT_LINEAR;
	qaws_mesh_reference_create(&r, g, NULL);
	qaws_mesh_reference_create(&lin, g, &rd);
	st = qaws_mesh_3d_measure_deviation(&d, r, g, NULL);
	mm_print("a PN saddle against its own triangles", &d);
	tol = 1e-4 * d.diagonal;
	for (k = 0; k < 3; k++)
	{
		wt[k] = (double)d.to_reference_at[k];
		at[k] = d.from_reference_at[k];
	}
	at_to = mm_dense(r, wt, 100);
	qaws_mesh_reference_project(&l, &dist, lin, at, NULL, 1);
	at_from = (double)dist;
	nt = qaws_mesh_reference_get_triangle_count(r);
	loc = (qaws_mesh_location*)malloc(sizeof(qaws_mesh_location) * np);
	pts = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * np);
	dd = (qaws_scalar*)malloc(sizeof(qaws_scalar) * np);
	for (f = 0; f < nt; f++)
	{
		for (j = 0, q = 0; j <= n; j++)
			for (i = 0; i + j <= n; i++, q++)
			{
				loc[q].patch = f;
				loc[q].a = (qaws_scalar)((double)i / n);
				loc[q].b = (qaws_scalar)((double)j / n);
			}
		qaws_mesh_reference_evaluate(pts, NULL, NULL, r, loc, np);
		qaws_mesh_reference_project(loc, dd, lin, pts, NULL, np);
		for (q = 0; q < np; q++)
			sup = dd[q] > sup ? dd[q] : sup;
	}
	sprintf(msg, "PN saddle against its own triangles: to [%.9g, %.9g] (D at its witness %.9g, densely), from [%.9g, %.9g] (D at its witness %.9g, "
		"sampled sup %.9g), tol %.3g, converged %d", (double)d.to_reference[0], (double)d.to_reference[1], at_to, (double)d.from_reference[0],
		(double)d.from_reference[1], at_from, sup, tol, d.converged);
	TEST_ASSERT(st == QAWS_STATUS_OK && d.converged == 1 && d.to_reference[1] - d.to_reference[0] <= tol && d.from_reference[1] - d.from_reference[0] <= tol
		&& d.to_reference[0] <= at_to + 1e-12 && d.from_reference[0] <= at_from + 1e-12 && sup <= d.from_reference[1] + 1e-12, msg);
	free(loc);
	free(pts);
	free(dd);
	qaws_mesh_reference_destroy(r);
	qaws_mesh_reference_destroy(lin);
	qaws_mesh_3d_destroy(g);
}

/* A pentagon with a straight corner, fanned from its first vertex into a zero-area triangle first: the reference's bound to it
   does not stick on that triangle's wedge (every point mapped to its centre), in either face order */
static void test_pentagon(void)
{
	static double const x[15] = { 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0.5, 0, 0 };
	static unsigned int const f[5] = { 0, 4, 1, 2, 3 }, five = 5;
	qaws_mesh_3d *m = mm_make(x, 5, f, &five, 1), *g = mm_grid(8, 0, mm_saddle);
	qaws_mesh_reference* r = NULL;
	qaws_mesh_deviation d;
	qaws_status st;
	double tol;
	char msg[300];
	qaws_mesh_reference_create(&r, g, NULL);
	st = qaws_mesh_3d_measure_deviation(&d, r, m, NULL);
	mm_print("a pentagon with a straight corner against a PN saddle", &d);
	tol = 1e-4 * d.diagonal;
	sprintf(msg, "a pentagon with a straight corner: from the reference [%.9g, %.9g], to it [%.9g, %.9g], tol %.3g, converged %d", (double)d.from_reference[0],
		(double)d.from_reference[1], (double)d.to_reference[0], (double)d.to_reference[1], tol, d.converged);
	TEST_ASSERT(st == QAWS_STATUS_OK && d.converged == 1 && d.from_reference[1] - d.from_reference[0] <= tol && d.to_reference[1] - d.to_reference[0] <= tol, msg);
	qaws_mesh_reference_destroy(r);
	qaws_mesh_3d_destroy(g);
	qaws_mesh_3d_destroy(m);
}

static qaws_mesh_3d* mm_remesh(qaws_mesh_3d const* m, qaws_mesh_reference const* r, double L, qaws_mesh_report* rep)
{
	qaws_remesh_desc d;
	qaws_mesh_3d* out = NULL;
	memset(&d, 0, sizeof d);
	d.reference = r;
	d.target_length = (qaws_scalar)L;
	qaws_mesh_remesh(&out, rep, m, &d);
	return out;
}

static void test_consistency_and_timing(void)
{
	qaws_mesh_3d *m = mm_bumped(), *out = NULL, *big = NULL;
	qaws_mesh_reference* r = NULL;
	qaws_mesh_report rep;
	qaws_mesh_deviation d;
	qaws_shape_quality s;
	qaws_status st = QAWS_STATUS_INTERNAL_ERROR;
	char msg[400];
	qaws_mesh_reference_create(&r, m, NULL);
	out = r ? mm_remesh(m, r, 0.12, &rep) : NULL;
	if (out)
		st = qaws_mesh_3d_measure_deviation(&d, r, out, NULL);
	if (out)
		mm_print("the bumped sphere remeshed at 0.12", &d);
	sprintf(msg, "remeshed at 0.12: the certified bound %.9g over the report's face-centre deviation %.9g, converged %d",
		out ? (double)d.to_reference[1] : 0, out ? (double)rep.max_deviation : 0, out ? d.converged : 0);
	TEST_ASSERT(out && st == QAWS_STATUS_OK && d.to_reference[1] >= rep.max_deviation - 1e-12 && d.to_reference[0] <= d.to_reference[1] && d.converged == 1, msg);
	/* timing on about 10k faces (heavy only: its remesh and measure are most of this test's time) */
	big = r && test_heavy() ? mm_remesh(m, r, 0.058, &rep) : NULL;
	if (big)
	{
		clock_t c0 = clock(), c1, c2;
		unsigned int i, rounds = 20;
		for (i = 0; i < rounds; i++)
			qaws_mesh_3d_measure_shape(&s, NULL, big);
		c1 = clock();
		st = qaws_mesh_3d_measure_deviation(&d, r, big, NULL);
		c2 = clock();
		printf("  %u faces: shape %.2f ms (V567 %.3f, Q %.3f / %.3f, angles %.1f..%.1f), deviation %.3f s (to [%.3g, %.3g], from [%.3g, %.3g], rms %.3g, converged %d)\n",
			qaws_mesh_3d_get_face_count(big), 1000.0 * (double)(c1 - c0) / CLOCKS_PER_SEC / rounds, (double)s.regular_share, (double)s.min_quality,
			(double)s.mean_quality, (double)s.min_angle, (double)s.max_angle, (double)(c2 - c1) / CLOCKS_PER_SEC, (double)d.to_reference[0],
			(double)d.to_reference[1], (double)d.from_reference[0], (double)d.from_reference[1], (double)d.rms, d.converged);
		sprintf(msg, "timing on %u faces: shape %.2f ms (target 10), deviation %.3f s (target 1 in Release, 5 at most)", qaws_mesh_3d_get_face_count(big),
			1000.0 * (double)(c1 - c0) / CLOCKS_PER_SEC / rounds, (double)(c2 - c1) / CLOCKS_PER_SEC);
		TEST_ASSERT(st == QAWS_STATUS_OK && (double)(c2 - c1) / CLOCKS_PER_SEC < 5 && qaws_mesh_3d_get_face_count(big) > 8000, msg);
	}
	else if (test_heavy())
		TEST_ASSERT(0, "timing: remeshed");
	qaws_mesh_reference_destroy(r);
	qaws_mesh_3d_destroy(m);
	qaws_mesh_3d_destroy(out);
	qaws_mesh_3d_destroy(big);
}

static void test_refusals(void)
{
	static double const x[9] = { 0, 0, 0, 1, 0, 0, 0, 1, 0 };
	static unsigned int const t[3] = { 0, 1, 2 };
	qaws_mesh_3d *m = mm_make(x, 3, t, NULL, 1), *empty = mm_make(x, 3, NULL, NULL, 0);
	qaws_mesh_reference *r = NULL, *re = NULL;
	qaws_mesh_deviation d;
	qaws_deviation_desc dd;
	qaws_shape_quality s;
	qaws_status s1, s2, s3, s4, s5, s6;
	char msg[300];
	qaws_mesh_reference_create(&r, m, NULL);
	qaws_mesh_reference_create(&re, empty, NULL);
	memset(&dd, 0, sizeof dd);
	dd.precision = (qaws_scalar)-1;
	s1 = qaws_mesh_3d_measure_deviation(&d, r, NULL, NULL);
	s2 = qaws_mesh_3d_measure_deviation(NULL, r, m, NULL);
	s3 = re ? qaws_mesh_3d_measure_deviation(&d, re, m, NULL) : QAWS_STATUS_INVALID_ARGUMENT;
	s4 = qaws_mesh_3d_measure_deviation(&d, r, m, &dd);
	s5 = qaws_mesh_3d_measure_shape(&s, NULL, NULL);
	s6 = qaws_mesh_3d_measure_shape(NULL, NULL, m);
	sprintf(msg, "refusals: NULL mesh %d, NULL out %d, a reference without triangles %d, a negative precision %d, shape NULL mesh %d, NULL out %d",
		(int)s1, (int)s2, (int)s3, (int)s4, (int)s5, (int)s6);
	TEST_ASSERT(s1 == QAWS_STATUS_INVALID_ARGUMENT && s2 == QAWS_STATUS_INVALID_ARGUMENT && s3 == QAWS_STATUS_INVALID_ARGUMENT
		&& s4 == QAWS_STATUS_INVALID_ARGUMENT && s5 == QAWS_STATUS_INVALID_ARGUMENT && s6 == QAWS_STATUS_INVALID_ARGUMENT, msg);
	qaws_mesh_reference_destroy(r);
	qaws_mesh_reference_destroy(re);
	qaws_mesh_3d_destroy(m);
	qaws_mesh_3d_destroy(empty);
}

int test_130_mesh_measure_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 130: quality measures\n");
	test_shapes();
	test_valences();
	test_exact();
	test_brute();
	test_two_sided();
	test_consistency_and_timing();
	test_saddle();
	test_pentagon();
	test_refusals();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
