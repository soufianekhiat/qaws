/*
 * Quality reports (qaws_quality_report.h): per-element measures of faces
 * and cells from their corner positions, and a summary over any array of
 * values. Everything goes through the meshes' public getters.
 */

#include "../internal/qaws_internal_memory.h"
#include "../qaws_quality_report.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define QR_PI 3.14159265358979323846

static double qr_dot(double const* a, double const* b)
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static void qr_sub(double const* a, double const* b, double* o)
{
	o[0] = a[0] - b[0];
	o[1] = a[1] - b[1];
	o[2] = a[2] - b[2];
}

static void qr_cross(double const* a, double const* b, double* o)
{
	o[0] = a[1] * b[2] - a[2] * b[1];
	o[1] = a[2] * b[0] - a[0] * b[2];
	o[2] = a[0] * b[1] - a[1] * b[0];
}

static double qr_len(double const* a)
{
	return sqrt(qr_dot(a, a));
}

static double qr_det(double const* a, double const* b, double const* c)
{
	double n[3];
	qr_cross(b, c, n);
	return qr_dot(a, n);
}

/* the Newell normal of a polygon (twice its vector area) */
static void qr_newell(double const* const* p, unsigned int n, double* N)
{
	unsigned int i;
	N[0] = N[1] = N[2] = 0;
	for (i = 0; i < n; i++)
	{
		double const *a = p[i], *b = p[(i + 1) % n];
		N[0] += (a[1] - b[1]) * (a[2] + b[2]);
		N[1] += (a[2] - b[2]) * (a[0] + b[0]);
		N[2] += (a[0] - b[0]) * (a[1] + b[1]);
	}
}

static double qr_edge_ratio(double const* const* p, unsigned int n)
{
	double lo = HUGE_VAL, hi = 0, e[3], l;
	unsigned int i;
	for (i = 0; i < n; i++)
	{
		qr_sub(p[(i + 1) % n], p[i], e);
		l = qr_len(e);
		lo = l < lo ? l : lo;
		hi = l > hi ? l : hi;
	}
	return lo > 0 ? hi / lo : HUGE_VAL;
}

/* a face of n corners; planar: measured in z = 0, signed by its turn */
static double qr_face(double const* const* p, unsigned int n, int planar, qaws_element_measure m)
{
	double N[3], un[3], area, len;
	unsigned int i;
	if (m == QAWS_MEASURE_EDGE_RATIO)
		return qr_edge_ratio(p, n);
	qr_newell(p, n, N);
	len = qr_len(N);
	if (planar)
	{
		un[0] = 0;
		un[1] = 0;
		un[2] = 1;
		area = 0.5 * N[2];
	}
	else
	{
		if (len > 0)
		{
			un[0] = N[0] / len;
			un[1] = N[1] / len;
			un[2] = N[2] / len;
		}
		else
			un[0] = un[1] = un[2] = 0;
		area = 0.5 * len;
	}
	if (m == QAWS_MEASURE_SIZE)
		return area;
	if (m == QAWS_MEASURE_MIN_ANGLE || m == QAWS_MEASURE_MAX_ANGLE || m == QAWS_MEASURE_SCALED_JACOBIAN)
	{
		double best = m == QAWS_MEASURE_MAX_ANGLE ? -HUGE_VAL : HUGE_VAL;
		if (m == QAWS_MEASURE_SCALED_JACOBIAN && n > 4)
			return NAN;
		for (i = 0; i < n; i++)
		{
			double e1[3], e2[3], c[3], v;
			qr_sub(p[(i + 1) % n], p[i], e1);
			qr_sub(p[(i + n - 1) % n], p[i], e2);
			qr_cross(e1, e2, c);
			if (m == QAWS_MEASURE_SCALED_JACOBIAN)
			{
				double l = qr_len(e1) * qr_len(e2);
				v = l > 0 ? qr_dot(c, un) / l : 0;
				if (n == 3)
					v /= sqrt(3.0) / 2;
			}
			else
			{
				/* the corner angle, past 180 where the corner turns against the normal (reflex or inverted) */
				v = atan2(qr_dot(c, un) < 0 ? -qr_len(c) : qr_len(c), qr_dot(e1, e2));
				if (v < 0)
					v += 2 * QR_PI;
				v *= 180 / QR_PI;
			}
			if (m == QAWS_MEASURE_MAX_ANGLE ? v > best : v < best)
				best = v;
		}
		return best;
	}
	if (n != 3)
		return NAN;
	{
		double a[3], b[3], c[3], la, lb, lc, s2 = 0;
		qr_sub(p[1], p[0], a);
		qr_sub(p[2], p[1], b);
		qr_sub(p[0], p[2], c);
		la = qr_len(a);
		lb = qr_len(b);
		lc = qr_len(c);
		s2 = la * la + lb * lb + lc * lc;
		if (m == QAWS_MEASURE_SHAPE)
			return s2 > 0 ? 4 * sqrt(3.0) * area / s2 : 0;
		/* ASPECT_RATIO: R / 2r = abc s / (8 A^2) */
		return area != 0 ? la * lb * lc * 0.5 * (la + lb + lc) / (8 * area * area) : HUGE_VAL;
	}
}

/* cell types: corner counts, facets (as in qaws_volume_mesh.h), corners with three edges in right-handed order and their scale */
static unsigned int const qr_cell_n[4] = { 4, 8, 6, 5 };
static unsigned int const qr_facet_n[4] = { 4, 6, 5, 5 };
static unsigned int const qr_facets[4][6][4] = {
	{ { 1, 2, 3, 9 }, { 0, 3, 2, 9 }, { 0, 1, 3, 9 }, { 0, 2, 1, 9 } },
	{ { 0, 3, 2, 1 }, { 4, 5, 6, 7 }, { 0, 1, 5, 4 }, { 1, 2, 6, 5 }, { 2, 3, 7, 6 }, { 3, 0, 4, 7 } },
	{ { 0, 2, 1, 9 }, { 3, 4, 5, 9 }, { 0, 1, 4, 3 }, { 1, 2, 5, 4 }, { 2, 0, 3, 5 } },
	{ { 0, 3, 2, 1 }, { 0, 1, 4, 9 }, { 1, 2, 4, 9 }, { 2, 3, 4, 9 }, { 3, 0, 4, 9 } }
};
static unsigned int const qr_corner_n[4] = { 4, 8, 6, 4 };
static unsigned int const qr_corners[4][8][4] = {
	{ { 0, 1, 2, 3 }, { 1, 2, 0, 3 }, { 2, 0, 1, 3 }, { 3, 0, 2, 1 } },
	{ { 0, 1, 3, 4 }, { 1, 2, 0, 5 }, { 2, 3, 1, 6 }, { 3, 0, 2, 7 }, { 4, 7, 5, 0 }, { 5, 4, 6, 1 }, { 6, 5, 7, 2 }, { 7, 6, 4, 3 } },
	{ { 0, 1, 2, 3 }, { 1, 2, 0, 4 }, { 2, 0, 1, 5 }, { 3, 5, 4, 0 }, { 4, 3, 5, 1 }, { 5, 4, 3, 2 } },
	{ { 0, 1, 3, 4 }, { 1, 2, 0, 4 }, { 2, 3, 1, 4 }, { 3, 0, 2, 4 } }
};
/* the regular cell's corner value: sqrt 2 / 2 (tetrahedron), 1 (cube), sqrt 3 / 2 (prism), sqrt 2 / 2 (pyramid of equal edges) */
static double const qr_corner_ideal[4] = { 0.70710678118654752, 1.0, 0.86602540378443865, 0.70710678118654752 };

static unsigned int qr_facet_size(unsigned int t, unsigned int f)
{
	return qr_facets[t][f][3] == 9 ? 3 : 4;
}

static double qr_volume(double const* const* p, unsigned int t)
{
	double v = 0;
	unsigned int f, j;
	for (f = 0; f < qr_facet_n[t]; f++)
	{
		unsigned int const* F = qr_facets[t][f];
		for (j = 1; j + 1 < qr_facet_size(t, f); j++)
		{
			double a[3], b[3], c[3];
			qr_sub(p[F[0]], p[0], a);
			qr_sub(p[F[j]], p[0], b);
			qr_sub(p[F[j + 1]], p[0], c);
			v += qr_det(a, b, c) / 6;
		}
	}
	return v;
}

static double qr_cell(double const* const* p, unsigned int t, qaws_element_measure m)
{
	unsigned int f, g, i, j;
	if (m == QAWS_MEASURE_SIZE)
		return qr_volume(p, t);
	if (m == QAWS_MEASURE_EDGE_RATIO)
	{
		double lo = HUGE_VAL, hi = 0;
		for (f = 0; f < qr_facet_n[t]; f++)
			for (i = 0; i < qr_facet_size(t, f); i++)
			{
				double e[3], l;
				qr_sub(p[qr_facets[t][f][(i + 1) % qr_facet_size(t, f)]], p[qr_facets[t][f][i]], e);
				l = qr_len(e);
				lo = l < lo ? l : lo;
				hi = l > hi ? l : hi;
			}
		return lo > 0 ? hi / lo : HUGE_VAL;
	}
	if (m == QAWS_MEASURE_MIN_ANGLE || m == QAWS_MEASURE_MAX_ANGLE)
	{
		/* dihedral angles: facets f and g sharing an edge (a, b in f, b, a in g) */
		double N[6][3], best = m == QAWS_MEASURE_MAX_ANGLE ? -HUGE_VAL : HUGE_VAL;
		for (f = 0; f < qr_facet_n[t]; f++)
		{
			double const* q[4];
			for (i = 0; i < qr_facet_size(t, f); i++)
				q[i] = p[qr_facets[t][f][i]];
			qr_newell(q, qr_facet_size(t, f), N[f]);
		}
		for (f = 0; f < qr_facet_n[t]; f++)
			for (i = 0; i < qr_facet_size(t, f); i++)
			{
				unsigned int a = qr_facets[t][f][i], b = qr_facets[t][f][(i + 1) % qr_facet_size(t, f)];
				for (g = f + 1; g < qr_facet_n[t]; g++)
					for (j = 0; j < qr_facet_size(t, g); j++)
						if (qr_facets[t][g][j] == b && qr_facets[t][g][(j + 1) % qr_facet_size(t, g)] == a)
						{
							double ln = qr_len(N[f]) * qr_len(N[g]), c = ln > 0 ? qr_dot(N[f], N[g]) / ln : 1, v;
							c = c > 1 ? 1 : c < -1 ? -1 : c;
							v = 180 - acos(c) * 180 / QR_PI;
							if (m == QAWS_MEASURE_MAX_ANGLE ? v > best : v < best)
								best = v;
						}
			}
		return best;
	}
	if (m == QAWS_MEASURE_SCALED_JACOBIAN)
	{
		double best = HUGE_VAL;
		for (i = 0; i < qr_corner_n[t]; i++)
		{
			unsigned int const* C = qr_corners[t][i];
			double e1[3], e2[3], e3[3], l, v;
			qr_sub(p[C[1]], p[C[0]], e1);
			qr_sub(p[C[2]], p[C[0]], e2);
			qr_sub(p[C[3]], p[C[0]], e3);
			l = qr_len(e1) * qr_len(e2) * qr_len(e3);
			v = l > 0 ? qr_det(e1, e2, e3) / l / qr_corner_ideal[t] : 0;
			best = v < best ? v : best;
		}
		return best;
	}
	if (t != 0)
		return NAN;
	{
		double a[3], b[3], c[3], s2 = 0, V = qr_volume(p, t);
		for (i = 0; i < 4; i++)
			for (j = i + 1; j < 4; j++)
			{
				double e[3];
				qr_sub(p[j], p[i], e);
				s2 += qr_dot(e, e);
			}
		if (m == QAWS_MEASURE_SHAPE)
			return s2 > 0 ? (V < 0 ? -1 : 1) * 12 * pow(3 * fabs(V), 2.0 / 3) / s2 : 0;
		/* ASPECT_RATIO: R / 3r, r = 3 V / (facet areas), R = |a^2 (b x c) + b^2 (c x a) + c^2 (a x b)| / (12 V) */
		{
			double bc[3], ca[3], ab[3], w[3], area = 0, r, R;
			qr_sub(p[1], p[0], a);
			qr_sub(p[2], p[0], b);
			qr_sub(p[3], p[0], c);
			if (V == 0)
				return HUGE_VAL;
			for (f = 0; f < 4; f++)
			{
				double const* q[3] = { p[qr_facets[0][f][0]], p[qr_facets[0][f][1]], p[qr_facets[0][f][2]] };
				double N[3];
				qr_newell(q, 3, N);
				area += 0.5 * qr_len(N);
			}
			qr_cross(b, c, bc);
			qr_cross(c, a, ca);
			qr_cross(a, b, ab);
			for (i = 0; i < 3; i++)
				w[i] = qr_dot(a, a) * bc[i] + qr_dot(b, b) * ca[i] + qr_dot(c, c) * ab[i];
			r = 3 * fabs(V) / area;
			R = qr_len(w) / (12 * fabs(V));
			return R / (3 * r);
		}
	}
}

static int qr_measure_ok(qaws_element_measure m)
{
	return (int)m >= 0 && m <= QAWS_MEASURE_SIZE;
}

qaws_status qaws_mesh_3d_measure_elements(qaws_scalar* out, qaws_mesh_3d const* mesh, qaws_element_measure m)
{
	unsigned int nv, nf, nc, f, c = 0, j, *fv, *fs;
	qaws_scalar* ps;
	double* x;
	qaws_status st;
	if (!mesh || !qr_measure_ok(m))
		return QAWS_STATUS_INVALID_ARGUMENT;
	nv = qaws_mesh_3d_get_vertex_count(mesh);
	nf = qaws_mesh_3d_get_face_count(mesh);
	nc = qaws_mesh_3d_get_corner_count(mesh);
	if (nf && !out)
		return QAWS_STATUS_INVALID_ARGUMENT;
	ps = (qaws_scalar*)qaws_internal_scratch_malloc(qaws_internal_mul3((size_t)nv + 1, 3, sizeof(qaws_scalar)));
	x = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3((size_t)nv + 1, 3, sizeof(double)));
	fv = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2((size_t)nc + 1, sizeof(unsigned int)));
	fs = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2((size_t)nf + 1, sizeof(unsigned int)));
	st = !ps || !x || !fv || !fs ? QAWS_STATUS_ALLOCATION_FAILURE : QAWS_STATUS_OK;
	if (st == QAWS_STATUS_OK)
		st = qaws_mesh_3d_get_positions(ps, mesh);
	if (st == QAWS_STATUS_OK)
		st = qaws_mesh_3d_get_faces(fv, fs, mesh);
	if (st == QAWS_STATUS_OK)
	{
		for (j = 0; j < 3 * nv; j++)
			x[j] = (double)ps[j];
		for (f = 0; f < nf; c += fs[f++])
		{
			double const* q[64];
			unsigned int n = fs[f] < 64 ? fs[f] : 64;
			for (j = 0; j < n; j++)
				q[j] = x + 3 * (size_t)fv[c + j];
			out[f] = (qaws_scalar)qr_face(q, n, 0, m);
		}
	}
	qaws_internal_free(ps);
	qaws_internal_free(x);
	qaws_internal_free(fv);
	qaws_internal_free(fs);
	return st;
}

qaws_status qaws_mesh_2d_measure_elements(qaws_scalar* out, qaws_mesh_2d const* mesh, qaws_element_measure m)
{
	unsigned int nt, nq, i, j, v[4];
	if (!mesh || !qr_measure_ok(m))
		return QAWS_STATUS_INVALID_ARGUMENT;
	nt = qaws_mesh_2d_get_triangle_count(mesh);
	nq = qaws_mesh_2d_get_quad_count(mesh);
	if (nt + nq && !out)
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (i = 0; i < nt + nq; i++)
	{
		double x[4][3];
		double const* q[4] = { x[0], x[1], x[2], x[3] };
		unsigned int n = i < nt ? 3 : 4;
		if (i < nt)
			qaws_mesh_2d_get_triangle(v, mesh, i);
		else
			qaws_mesh_2d_get_quad(v, mesh, i - nt);
		for (j = 0; j < n; j++)
		{
			qaws_vec2 p = qaws_mesh_2d_get_vertex(mesh, v[j]);
			x[j][0] = p.x;
			x[j][1] = p.y;
			x[j][2] = 0;
		}
		out[i] = (qaws_scalar)qr_face(q, n, 1, m);
	}
	return QAWS_STATUS_OK;
}

qaws_status qaws_volume_mesh_measure_elements(qaws_scalar* out, qaws_volume_mesh const* mesh, qaws_element_measure m)
{
	unsigned int nv, ncell, ncorner, i, j, c = 0, *cv;
	qaws_cell_type* ct;
	qaws_scalar* ps;
	double* x;
	qaws_status st;
	if (!mesh || !qr_measure_ok(m))
		return QAWS_STATUS_INVALID_ARGUMENT;
	nv = qaws_volume_mesh_get_vertex_count(mesh);
	ncell = qaws_volume_mesh_get_cell_count(mesh);
	ncorner = qaws_volume_mesh_get_corner_count(mesh);
	if (ncell && !out)
		return QAWS_STATUS_INVALID_ARGUMENT;
	ps = (qaws_scalar*)qaws_internal_scratch_malloc(qaws_internal_mul3((size_t)nv + 1, 3, sizeof(qaws_scalar)));
	x = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3((size_t)nv + 1, 3, sizeof(double)));
	cv = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2((size_t)ncorner + 1, sizeof(unsigned int)));
	ct = (qaws_cell_type*)qaws_internal_scratch_malloc(qaws_internal_mul2((size_t)ncell + 1, sizeof(qaws_cell_type)));
	st = !ps || !x || !cv || !ct ? QAWS_STATUS_ALLOCATION_FAILURE : QAWS_STATUS_OK;
	if (st == QAWS_STATUS_OK)
		st = qaws_volume_mesh_get_positions(ps, mesh);
	if (st == QAWS_STATUS_OK)
		st = qaws_volume_mesh_get_cells(cv, ct, mesh);
	if (st == QAWS_STATUS_OK)
	{
		for (j = 0; j < 3 * nv; j++)
			x[j] = (double)ps[j];
		for (i = 0; i < ncell; c += qr_cell_n[ct[i++]])
		{
			double const* q[8];
			for (j = 0; j < qr_cell_n[ct[i]]; j++)
				q[j] = x + 3 * (size_t)cv[c + j];
			out[i] = (qaws_scalar)qr_cell(q, (unsigned int)ct[i], m);
		}
	}
	qaws_internal_free(ps);
	qaws_internal_free(x);
	qaws_internal_free(cv);
	qaws_internal_free(ct);
	return st;
}

/* ---- the summary ---- */

typedef struct qr_pair
{
	double v;
	unsigned int i;
} qr_pair;

static int qr_pair_cmp(void const* a, void const* b)
{
	qr_pair const *x = (qr_pair const*)a, *y = (qr_pair const*)b;
	if (x->v != y->v)
		return x->v < y->v ? -1 : 1;
	return x->i < y->i ? -1 : x->i > y->i;
}

qaws_status qaws_quality_summarize(qaws_quality_summary* out, unsigned int* hist, unsigned int* worst, qaws_scalar const* values, unsigned int count,
	qaws_quality_summary_desc const* d)
{
	qr_pair* s;
	unsigned int n = 0, i, k;
	double sum = 0, lo, hi;
	int higher_worse;
	if (!out || !d || !qr_measure_ok(d->measure) || (count && !values) || (d->bins && !hist) || (d->worst_count && !worst))
		return QAWS_STATUS_INVALID_ARGUMENT;
	memset(out, 0, sizeof *out);
	higher_worse = d->measure == QAWS_MEASURE_MAX_ANGLE || d->measure == QAWS_MEASURE_ASPECT_RATIO || d->measure == QAWS_MEASURE_EDGE_RATIO;
	s = (qr_pair*)qaws_internal_scratch_malloc(qaws_internal_mul2((size_t)count + 1, sizeof(qr_pair)));
	if (!s)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; i < count; i++)
	{
		double v = (double)values[i];
		if (v != v)
		{
			out->skipped++;
			continue;
		}
		s[n].v = v;
		s[n++].i = i;
		sum += v;
		if (d->check && (higher_worse ? v > (double)d->threshold : v < (double)d->threshold))
			out->failing++;
	}
	out->count = n;
	if (d->bins)
		memset(hist, 0, sizeof(unsigned int) * d->bins);
	if (!n)
	{
		qaws_internal_free(s);
		return QAWS_STATUS_OK;
	}
	qsort(s, n, sizeof *s, qr_pair_cmp);
	out->min = (qaws_scalar)s[0].v;
	out->max = (qaws_scalar)s[n - 1].v;
	out->mean = (qaws_scalar)(sum / n);
	out->median = (qaws_scalar)(n & 1 ? s[n / 2].v : 0.5 * (s[n / 2 - 1].v + s[n / 2].v));
	lo = d->range[0] != d->range[1] ? (double)d->range[0] : s[0].v;
	hi = d->range[0] != d->range[1] ? (double)d->range[1] : s[n - 1].v;
	out->range[0] = (qaws_scalar)lo;
	out->range[1] = (qaws_scalar)hi;
	for (i = 0; i < n && d->bins; i++)
	{
		double t = hi > lo ? (s[i].v - lo) / (hi - lo) * d->bins : 0;
		unsigned int b = t <= 0 || t != t ? 0 : t >= d->bins ? d->bins - 1 : (unsigned int)t;
		hist[b]++;
	}
	k = d->worst_count < n ? d->worst_count : n;
	for (i = 0; i < k; i++)
		worst[i] = higher_worse ? s[n - 1 - i].i : s[i].i;
	out->worst_found = k;
	qaws_internal_free(s);
	return QAWS_STATUS_OK;
}
