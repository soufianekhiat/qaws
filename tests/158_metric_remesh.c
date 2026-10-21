/*
 * Test 158: anisotropic remeshing to a metric field (qaws_metric_field_*,
 * qaws_remesh_desc.metric)
 *
 *   - a unit square, a constant metric asking for edges 0.02 across x and
 *     0.15 along y: most edges near unit length in the metric (in
 *     [1/sqrt 2, sqrt 2]), triangles well shaped in it, their count near
 *     the square's metric area over a unit triangle's, the edges' extent
 *     along y several times their extent along x;
 *   - the Hessian metric of f = tanh(12 (x - 1/2)) on a fine square: the
 *     recovered metric along x (2/9) |f''| / error to a few percent, along y
 *     max_length to a percent (f is linear there); the remeshed square
 *     interpolates f within twice the error, its elements stretched along
 *     the layer, with fewer vertices than an isotropic mesh at the
 *     smallest size would need;
 *   - an open cylinder, a constant metric across the axis 0.08 and along it
 *     0.4: edges near unit length in the metric, the vertices on the
 *     cylinder, the elements longer along the axis;
 *   - CVT with a metric against LOCAL on a unit square, 0.2 along y and
 *     across x 0.02 or graded 0.02 to 0.08: edges near unit length in the
 *     metric (their mean and spread printed), triangles well shaped in it,
 *     their count near the metric area's, stretched along y;
 *   - refusals: a metric on another reference, a metric
 *     not positive definite, a Hessian desc's error 0 or sizes crossed, a
 *     location off the reference.
 */

#include "test_common.h"
#include "qaws_remesh.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

/* an n x n grid of the unit square in z = 0, each cell two triangles */
static qaws_mesh_3d* mt_square(unsigned int n)
{
	unsigned int nv = (n + 1) * (n + 1), i, j, t = 0;
	qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	unsigned int* f = (unsigned int*)malloc(sizeof(unsigned int) * 6 * n * n);
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
			unsigned int a = j * (n + 1) + i, b = a + 1, c = a + n + 2, e = a + n + 1;
			f[t++] = a; f[t++] = b; f[t++] = c;
			f[t++] = a; f[t++] = c; f[t++] = e;
		}
	memset(&d, 0, sizeof d);
	d.positions = x;
	d.vertex_count = nv;
	d.face_vertices = f;
	d.face_count = 2 * n * n;
	qaws_mesh_3d_create(&m, &d);
	free(x);
	free(f);
	return m;
}

/* a cylinder of radius 1 and height 2 about z, na around and nz up, open; its radial normals */
static qaws_mesh_3d* mt_cylinder(unsigned int na, unsigned int nz, qaws_scalar* nrm)
{
	unsigned int nv = na * (nz + 1), i, j, t = 0;
	qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	unsigned int* f = (unsigned int*)malloc(sizeof(unsigned int) * 6 * na * nz);
	qaws_mesh_3d_desc d;
	qaws_mesh_3d* m = NULL;
	for (j = 0; j <= nz; j++)
		for (i = 0; i < na; i++)
		{
			double a = 2 * 3.14159265358979323846 * i / na;
			unsigned int v = j * na + i;
			x[3 * v] = (qaws_scalar)cos(a);
			x[3 * v + 1] = (qaws_scalar)sin(a);
			x[3 * v + 2] = (qaws_scalar)(2.0 * j / nz);
			nrm[3 * v] = (qaws_scalar)cos(a);
			nrm[3 * v + 1] = (qaws_scalar)sin(a);
			nrm[3 * v + 2] = 0;
		}
	for (j = 0; j < nz; j++)
		for (i = 0; i < na; i++)
		{
			unsigned int a = j * na + i, b = j * na + (i + 1) % na, c = b + na, e = a + na;
			f[t++] = a; f[t++] = b; f[t++] = c;
			f[t++] = a; f[t++] = c; f[t++] = e;
		}
	memset(&d, 0, sizeof d);
	d.positions = x;
	d.vertex_count = nv;
	d.face_vertices = f;
	d.face_count = 2 * na * nz;
	qaws_mesh_3d_create(&m, &d);
	free(x);
	free(f);
	return m;
}

static double mt_mlen(double const* M, double const* e)
{
	return sqrt(e[0] * (M[0] * e[0] + M[1] * e[1] + M[2] * e[2]) + e[1] * (M[1] * e[0] + M[3] * e[1] + M[4] * e[2])
		+ e[2] * (M[2] * e[0] + M[4] * e[1] + M[5] * e[2]));
}

/* a remeshed mesh's measures in the metric: the share of edges in [1/sqrt 2, sqrt 2], the edges' mean length and its
   standard deviation, the triangles' mean shape, their smallest angle and the share of angles under 30 degrees, the mean
   extents of the edges along axis `along` and across it (the other two axes' length) */
typedef struct mt_stats
{
	unsigned int nv, nt;
	double unit, shorter, longer, quality, along, across, mean, spread, min_angle, under30;
} mt_stats;

static int mt_measure(mt_stats* s, qaws_mesh_3d const* m, qaws_metric_field const* mf, unsigned int along)
{
	unsigned int nv = qaws_mesh_3d_get_vertex_count(m), ne = qaws_mesh_3d_get_edge_count(m), nt = qaws_mesh_3d_get_triangle_count(m), i, k, in = 0;
	qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	qaws_scalar* M = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 6 * nv);
	qaws_mesh_location* l = (qaws_mesh_location*)malloc(sizeof(qaws_mesh_location) * nv);
	unsigned int* e = (unsigned int*)malloc(sizeof(unsigned int) * 2 * ne);
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 3 * nt);
	int ok = x && M && l && e && t && qaws_mesh_3d_get_positions(x, m) == QAWS_STATUS_OK && qaws_mesh_3d_get_footpoints(l, m) == QAWS_STATUS_OK
		&& qaws_metric_field_evaluate(M, mf, l, nv) == QAWS_STATUS_OK && qaws_mesh_3d_get_edges(e, m) == QAWS_STATUS_OK
		&& qaws_mesh_3d_get_triangles(t, m) == QAWS_STATUS_OK;
	memset(s, 0, sizeof *s);
	s->nv = nv;
	s->nt = nt;
	s->min_angle = 180;
	for (i = 0; ok && i < ne; i++)
	{
		double d[3], A[6], len;
		unsigned int a = e[2 * i], b = e[2 * i + 1];
		for (k = 0; k < 3; k++)
			d[k] = (double)x[3 * b + k] - (double)x[3 * a + k];
		for (k = 0; k < 6; k++)
			A[k] = 0.5 * ((double)M[6 * a + k] + (double)M[6 * b + k]);
		len = mt_mlen(A, d);
		in += len >= sqrt(0.5) && len <= sqrt(2.0);
		s->shorter += len < sqrt(0.5);
		s->longer += len > sqrt(2.0);
		s->mean += len;
		s->spread += len * len;
		s->along += fabs(d[along]);
		s->across += sqrt(d[(along + 1) % 3] * d[(along + 1) % 3] + d[(along + 2) % 3] * d[(along + 2) % 3]);
	}
	for (i = 0; ok && i < nt; i++)
	{
		double u[3], v[3], w[3], A[6], uu, vv, ww, uv;
		unsigned int a = t[3 * i], b = t[3 * i + 1], c = t[3 * i + 2];
		for (k = 0; k < 3; k++)
		{
			u[k] = (double)x[3 * b + k] - (double)x[3 * a + k];
			v[k] = (double)x[3 * c + k] - (double)x[3 * a + k];
			w[k] = v[k] - u[k];
		}
		for (k = 0; k < 6; k++)
			A[k] = ((double)M[6 * a + k] + (double)M[6 * b + k] + (double)M[6 * c + k]) / 3;
		uu = mt_mlen(A, u) * mt_mlen(A, u);
		vv = mt_mlen(A, v) * mt_mlen(A, v);
		ww = mt_mlen(A, w) * mt_mlen(A, w);
		uv = 0.5 * (uu + vv - ww);
		s->quality += 2 * sqrt(3.0) * sqrt(fmax(uu * vv - uv * uv, 0)) / (uu + vv + ww);
		{
			/* the angles by the law of cosines on the metric lengths */
			double sq[3], cs;
			sq[0] = uu; sq[1] = vv; sq[2] = ww;
			for (k = 0; k < 3; k++)
			{
				double a2 = sq[k], b2 = sq[(k + 1) % 3], c2 = sq[(k + 2) % 3], ang;
				cs = (b2 + c2 - a2) / (2 * sqrt(b2 * c2));
				ang = acos(cs < -1 ? -1 : cs > 1 ? 1 : cs) * 180 / 3.14159265358979323846;
				s->min_angle = fmin(s->min_angle, ang);
				s->under30 += ang < 30;
			}
		}
	}
	if (ok && ne && nt)
	{
		s->unit = (double)in / ne;
		s->mean /= ne;
		s->spread = sqrt(fmax(s->spread / ne - s->mean * s->mean, 0));
		s->under30 /= 3.0 * nt;
		s->shorter /= ne;
		s->longer /= ne;
		s->along /= ne;
		s->across /= ne;
		s->quality /= nt;
	}
	free(x); free(M); free(l); free(e); free(t);
	return ok;
}

static qaws_mesh_3d* mt_remesh(qaws_mesh_3d const* m, qaws_mesh_reference const* r, qaws_metric_field const* mf, qaws_remesh_method method, qaws_status* st)
{
	qaws_remesh_desc d;
	qaws_mesh_report rep;
	qaws_mesh_3d* out = NULL;
	memset(&d, 0, sizeof d);
	memset(&rep, 0, sizeof rep);
	d.method = method;
	d.reference = r;
	d.metric = mf;
	*st = qaws_mesh_remesh(&out, &rep, m, &d);
	printf("    (%u iterations, %u splits, %u collapses, %u flips, %u moves)\n", rep.iterations, rep.splits, rep.collapses, rep.flips, rep.moves);
	return out;
}

static void test_constant(void)
{
	qaws_mesh_3d *m = mt_square(20), *out = NULL;
	qaws_mesh_reference* r = NULL;
	qaws_metric_field* mf = NULL;
	qaws_scalar* M = NULL;
	qaws_status st = QAWS_STATUS_INTERNAL_ERROR;
	unsigned int nv = m ? qaws_mesh_3d_get_vertex_count(m) : 0, i;
	double hx = 0.02, hy = 0.15, expect = (1 / (hx * hy)) / (sqrt(3.0) / 4);
	mt_stats s;
	char msg[400];
	qaws_mesh_reference_create(&r, m, NULL);
	M = (qaws_scalar*)calloc(6 * (size_t)nv + 1, sizeof(qaws_scalar));
	for (i = 0; M && i < nv; i++)
	{
		M[6 * i] = (qaws_scalar)(1 / (hx * hx));
		M[6 * i + 3] = (qaws_scalar)(1 / (hy * hy));
		M[6 * i + 5] = 1;
	}
	if (r && M)
		qaws_metric_field_create_vertices(&mf, r, M);
	if (mf)
		out = mt_remesh(m, r, mf, QAWS_REMESH_LOCAL, &st);
	if (!out || !mt_measure(&s, out, mf, 1))
	{
		TEST_ASSERT(0, "constant metric: remeshed");
		goto done;
	}
	sprintf(msg, "constant metric (0.02 x 0.15): %u vertices, %u triangles (metric area over a unit triangle: %.0f), %.3f of the edges of unit "
		"length, mean shape %.3f, mean extent along y %.4f, along x %.4f", s.nv, s.nt, expect, s.unit, s.quality, s.along, s.across);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && s.unit >= 0.85, msg);
	TEST_ASSERT(s.quality >= 0.8, "constant metric: triangles well shaped in the metric");
	TEST_ASSERT(s.nt >= 0.6 * expect && s.nt <= 1.4 * expect, "constant metric: triangle count near the metric area's");
	TEST_ASSERT(s.along >= 3 * s.across, "constant metric: edges stretched along y");
done:
	qaws_mesh_3d_destroy(out);
	qaws_metric_field_destroy(mf);
	qaws_mesh_reference_destroy(r);
	qaws_mesh_3d_destroy(m);
	free(M);
}

/* CVT against LOCAL on the unit square, 0.2 along y and across x 0.02 (constant) or 0.02 + 0.06 x (graded) */
static void test_cvt(void)
{
	unsigned int c, method;
	for (c = 0; c < 2; c++)
		for (method = 0; method < 2; method++)
		{
			qaws_mesh_3d *m = mt_square(20), *out = NULL;
			qaws_mesh_reference* r = NULL;
			qaws_metric_field* mf = NULL;
			qaws_scalar *M = NULL, *x = NULL;
			qaws_status st = QAWS_STATUS_INTERNAL_ERROR;
			unsigned int nv = m ? qaws_mesh_3d_get_vertex_count(m) : 0, i;
			double hy = 0.2, expect = 0, t0, secs;
			mt_stats s;
			char msg[600];
			qaws_mesh_reference_create(&r, m, NULL);
			M = (qaws_scalar*)calloc(6 * (size_t)nv + 1, sizeof(qaws_scalar));
			x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * (nv + 1));
			if (!r || !M || !x || qaws_mesh_3d_get_positions(x, m) != QAWS_STATUS_OK)
			{
				TEST_ASSERT(0, "cvt: set up");
				goto next;
			}
			for (i = 0; i < nv; i++)
			{
				double hx = c ? 0.02 + 0.06 * (double)x[3 * i] : 0.02;
				M[6 * i] = (qaws_scalar)(1 / (hx * hx));
				M[6 * i + 3] = (qaws_scalar)(1 / (hy * hy));
				M[6 * i + 5] = 1;
			}
			/* the metric area over a unit triangle's: the integral of 1 / (hx hy) */
			expect = (c ? log(0.08 / 0.02) / 0.06 : 1 / 0.02) / hy / (sqrt(3.0) / 4);
			qaws_metric_field_create_vertices(&mf, r, M);
			t0 = (double)clock();
			if (mf)
				out = mt_remesh(m, r, mf, method ? QAWS_REMESH_CVT : QAWS_REMESH_LOCAL, &st);
			secs = ((double)clock() - t0) / CLOCKS_PER_SEC;
			if (!out || !mt_measure(&s, out, mf, 1))
			{
				TEST_ASSERT(0, "cvt: remeshed");
				goto next;
			}
			sprintf(msg, "%s, %s metric (0.02%s x 0.2): %u triangles (metric area: %.0f), edges in the metric %.3f +- %.3f (%.3f in [1/sqrt 2, "
				"sqrt 2]), mean shape %.3f, smallest angle %.1f, %.3f of the angles under 30, extent along y over across %.2f, %.2f s", method ? "CVT" : "LOCAL",
				c ? "graded" : "constant", c ? "..0.08" : "", s.nt, expect, s.mean, s.spread, s.unit, s.quality, s.min_angle, s.under30, s.along / s.across, secs);
			printf("  %s\n", msg);
			TEST_ASSERT(st == QAWS_STATUS_OK && s.unit >= 0.85 && s.quality >= 0.8, msg);
			TEST_ASSERT(s.nt >= 0.7 * expect && s.nt <= 1.3 * expect, "cvt: triangle count near the metric area's");
			TEST_ASSERT(s.along >= 3 * s.across, "cvt: edges stretched along y");
			TEST_ASSERT(!method || (fabs(s.mean - 1) <= 0.05 && s.spread <= 0.15), "cvt: CVT's edges within 5% of unit length on average, spread under 0.15");
		next:
			qaws_mesh_3d_destroy(out);
			qaws_metric_field_destroy(mf);
			qaws_mesh_reference_destroy(r);
			qaws_mesh_3d_destroy(m);
			free(M);
			free(x);
		}
}

static double mt_f(double x)
{
	return tanh(12 * (x - 0.5));
}

static double mt_f2(double x)
{
	double t = tanh(12 * (x - 0.5));
	return -2 * 144 * t * (1 - t * t);
}

static void test_hessian(void)
{
	qaws_mesh_3d *m = mt_square(64), *out = NULL;
	qaws_mesh_reference* r = NULL;
	qaws_metric_field* mf = NULL;
	qaws_hessian_metric_desc hd;
	qaws_scalar *f = NULL, *x = NULL;
	qaws_status st = QAWS_STATUS_INTERNAL_ERROR;
	unsigned int nv = m ? qaws_mesh_3d_get_vertex_count(m) : 0, i, k, n = 0;
	double eps = 0.01, hmin = 0.01, hmax = 0.2, worst_xx = 0, worst_yy = 0, err = 0, layer_along = 0, layer_across = 0;
	mt_stats s;
	char msg[400];
	qaws_mesh_reference_create(&r, m, NULL);
	f = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (nv + 1));
	x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * (nv + 1));
	if (!r || !f || !x || qaws_mesh_3d_get_positions(x, m) != QAWS_STATUS_OK)
	{
		TEST_ASSERT(0, "hessian: set up");
		goto done;
	}
	for (i = 0; i < nv; i++)
		f[i] = (qaws_scalar)mt_f((double)x[3 * i]);
	memset(&hd, 0, sizeof hd);
	hd.error = (qaws_scalar)eps;
	hd.min_length = (qaws_scalar)hmin;
	hd.max_length = (qaws_scalar)hmax;
	st = qaws_metric_field_create_hessian(&mf, r, f, &hd);
	TEST_ASSERT(st == QAWS_STATUS_OK && mf, "hessian: field created");
	if (!mf)
		goto done;
	/* the recovered metric at the vertices two cells or more inside, against (2/9) |f''| / error clamped */
	for (i = 0; i < nv; i++)
	{
		qaws_mesh_location l;
		qaws_scalar M[6];
		double px = (double)x[3 * i], py = (double)x[3 * i + 1], want;
		if (px < 0.04 || px > 0.96 || py < 0.04 || py > 0.96 || qaws_mesh_reference_get_vertex_location(&l, r, i) != QAWS_STATUS_OK
			|| qaws_metric_field_evaluate(M, mf, &l, 1) != QAWS_STATUS_OK)
			continue;
		want = fmin(fmax((2.0 / 9) * fabs(mt_f2(px)) / eps, 1 / (hmax * hmax)), 1 / (hmin * hmin));
		worst_xx = fmax(worst_xx, fabs((double)M[0] - want) / want);
		worst_yy = fmax(worst_yy, fabs((double)M[3] * hmax * hmax - 1));
	}
	sprintf(msg, "hessian: recovered metric along x within %.4f of (2/9)|f''|/error (relative), along y within %.2e of 1/max_length^2", worst_xx,
		worst_yy);
	printf("  %s\n", msg);
	TEST_ASSERT(worst_xx <= 0.08, msg);
	/* (the fit's cross term tilts the axes by a few thousandths: the large x term leaks that, squared, into y) */
	TEST_ASSERT(worst_yy <= 0.01, "hessian: along y, the largest size to a percent (f linear there)");
	out = mt_remesh(m, r, mf, QAWS_REMESH_LOCAL, &st);
	if (!out || !mt_measure(&s, out, mf, 1))
	{
		TEST_ASSERT(0, "hessian: remeshed");
		goto done;
	}
	/* f's linear interpolation on the new mesh against f at each triangle's centre and edge midpoints */
	{
		unsigned int on = qaws_mesh_3d_get_vertex_count(out), nt = qaws_mesh_3d_get_triangle_count(out);
		qaws_scalar* y = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * on);
		unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 3 * nt);
		if (y && t && qaws_mesh_3d_get_positions(y, out) == QAWS_STATUS_OK && qaws_mesh_3d_get_triangles(t, out) == QAWS_STATUS_OK)
			for (i = 0; i < nt; i++)
			{
				static double const bw[4][3] = { { 1.0 / 3, 1.0 / 3, 1.0 / 3 }, { 0.5, 0.5, 0 }, { 0, 0.5, 0.5 }, { 0.5, 0, 0.5 } };
				double cx = 0, lo = 1, hi = 0;
				for (k = 0; k < 4; k++)
				{
					double px = 0, lin = 0;
					unsigned int j;
					for (j = 0; j < 3; j++)
					{
						px += bw[k][j] * (double)y[3 * t[3 * i + j]];
						lin += bw[k][j] * mt_f((double)y[3 * t[3 * i + j]]);
					}
					err = fmax(err, fabs(lin - mt_f(px)));
				}
				for (k = 0; k < 3; k++)
				{
					double vx = (double)y[3 * t[3 * i + k]];
					cx += vx / 3;
					lo = fmin(lo, vx);
					hi = fmax(hi, vx);
				}
				/* the elements across the layer's steepest part: their extent along it over across it */
				if (fabs(cx - 0.5) < 0.1)
				{
					double ylo = 2, yhi = -1;
					for (k = 0; k < 3; k++)
					{
						ylo = fmin(ylo, (double)y[3 * t[3 * i + k] + 1]);
						yhi = fmax(yhi, (double)y[3 * t[3 * i + k] + 1]);
					}
					layer_along += yhi - ylo;
					layer_across += hi - lo;
					n++;
				}
			}
		free(y);
		free(t);
	}
	sprintf(msg, "hessian: %u vertices (isotropic at the smallest size: about %.0f), %.3f of the edges of unit length (%.3f shorter, %.3f "
		"longer), mean shape %.3f, f interpolated within %.4f (error %.3f), %u triangles in the layer %.2f times longer along it than across", s.nv,
		2 / (sqrt(3.0) * hmin * hmin) * 0.2, s.unit, s.shorter, s.longer, s.quality, err, eps, n, n ? layer_along / layer_across : 0);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && err <= 2 * eps, msg);
	TEST_ASSERT(s.unit >= 0.75 && s.quality >= 0.7, "hessian: edges near unit length, triangles well shaped in the metric");
	TEST_ASSERT(n > 0 && layer_along >= 3 * layer_across, "hessian: elements stretched along the layer");
	TEST_ASSERT(s.nv < 2 / (sqrt(3.0) * hmin * hmin) * 0.2 / 4, "hessian: a quarter of the vertices of isotropic refinement over the layer");
done:
	qaws_mesh_3d_destroy(out);
	qaws_metric_field_destroy(mf);
	qaws_mesh_reference_destroy(r);
	qaws_mesh_3d_destroy(m);
	free(f);
	free(x);
}

static void test_cylinder(void)
{
	qaws_scalar* nrm = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * 64 * 17);
	qaws_mesh_3d *m = nrm ? mt_cylinder(64, 16, nrm) : NULL, *out = NULL;
	qaws_mesh_reference_desc rd;
	qaws_mesh_reference* r = NULL;
	qaws_metric_field* mf = NULL;
	qaws_scalar* M = NULL;
	qaws_status st = QAWS_STATUS_INTERNAL_ERROR;
	unsigned int nv = m ? qaws_mesh_3d_get_vertex_count(m) : 0, i;
	double off = 0;
	mt_stats s;
	char msg[400];
	memset(&rd, 0, sizeof rd);
	rd.normals = nrm;
	qaws_mesh_reference_create(&r, m, &rd);
	M = (qaws_scalar*)calloc(6 * (size_t)nv + 1, sizeof(qaws_scalar));
	for (i = 0; M && i < nv; i++)
	{
		M[6 * i] = M[6 * i + 3] = (qaws_scalar)(1 / (0.08 * 0.08));
		M[6 * i + 5] = (qaws_scalar)(1 / (0.4 * 0.4));
	}
	if (r && M)
		qaws_metric_field_create_vertices(&mf, r, M);
	if (mf)
		out = mt_remesh(m, r, mf, QAWS_REMESH_LOCAL, &st);
	if (!out || !mt_measure(&s, out, mf, 2))
	{
		TEST_ASSERT(0, "cylinder: remeshed");
		goto done;
	}
	{
		unsigned int on = qaws_mesh_3d_get_vertex_count(out);
		qaws_scalar* y = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * on);
		if (y && qaws_mesh_3d_get_positions(y, out) == QAWS_STATUS_OK)
			for (i = 0; i < on; i++)
				off = fmax(off, fabs(sqrt((double)y[3 * i] * y[3 * i] + (double)y[3 * i + 1] * y[3 * i + 1]) - 1));
		free(y);
	}
	sprintf(msg, "cylinder (0.08 around, 0.4 along): %u vertices, %.3f of the edges of unit length, mean shape %.3f, mean extent along the "
		"axis %.4f, across %.4f, %.1e off the cylinder", s.nv, s.unit, s.quality, s.along, s.across, off);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && s.unit >= 0.8 && s.quality >= 0.75, msg);
	TEST_ASSERT(s.along >= 2 * s.across, "cylinder: edges stretched along the axis");
	/* (PN between the 64 sides: the reference's own distance from the cylinder) */
	TEST_ASSERT(off <= 2e-4, "cylinder: vertices on the cylinder");
done:
	qaws_mesh_3d_destroy(out);
	qaws_metric_field_destroy(mf);
	qaws_mesh_reference_destroy(r);
	qaws_mesh_3d_destroy(m);
	free(M);
	free(nrm);
}

static void test_refusals(void)
{
	qaws_mesh_3d *m = mt_square(4), *out = NULL;
	qaws_mesh_reference *r = NULL, *r2 = NULL;
	qaws_metric_field* mf = NULL;
	qaws_metric_field* bad = NULL;
	qaws_hessian_metric_desc hd;
	qaws_remesh_desc d;
	qaws_scalar M[6 * 25], f[25];
	qaws_mesh_location l;
	unsigned int i;
	qaws_mesh_reference_create(&r, m, NULL);
	qaws_mesh_reference_create(&r2, m, NULL);
	for (i = 0; i < 25; i++)
	{
		M[6 * i] = 4; M[6 * i + 1] = 0; M[6 * i + 2] = 0; M[6 * i + 3] = 4; M[6 * i + 4] = 0; M[6 * i + 5] = 4;
		f[i] = (qaws_scalar)i;
	}
	TEST_ASSERT(qaws_metric_field_create_vertices(&mf, r, M) == QAWS_STATUS_OK, "refusals: a positive definite metric taken");
	memset(&d, 0, sizeof d);
	d.reference = r;
	d.metric = mf;
	d.reference = r2;
	TEST_ASSERT(qaws_mesh_remesh(&out, NULL, m, &d) == QAWS_STATUS_INVALID_ARGUMENT && !out, "refusals: a metric on another reference");
	M[6 * 7 + 1] = 5;   /* (xy 5 against xx = yy = 4: indefinite) */
	TEST_ASSERT(qaws_metric_field_create_vertices(&bad, r, M) == QAWS_STATUS_INVALID_ARGUMENT && !bad, "refusals: a metric not positive definite");
	memset(&hd, 0, sizeof hd);
	TEST_ASSERT(qaws_metric_field_create_hessian(&bad, r, f, &hd) == QAWS_STATUS_INVALID_ARGUMENT && !bad, "refusals: error 0");
	hd.error = 0.1f;
	hd.min_length = 0.5f;
	hd.max_length = 0.25f;
	TEST_ASSERT(qaws_metric_field_create_hessian(&bad, r, f, &hd) == QAWS_STATUS_INVALID_ARGUMENT && !bad, "refusals: min_length past max_length");
	l.patch = 1000;
	l.a = l.b = 0;
	TEST_ASSERT(qaws_metric_field_evaluate(M, mf, &l, 1) == QAWS_STATUS_INVALID_ARGUMENT, "refusals: a location off the reference");
	qaws_metric_field_destroy(mf);
	qaws_mesh_reference_destroy(r);
	qaws_mesh_reference_destroy(r2);
	qaws_mesh_3d_destroy(m);
}

int test_158_metric_remesh_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 158: metric remeshing\n");
	test_constant();
	test_hessian();
	test_cylinder();
	test_cvt();
	test_refusals();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
