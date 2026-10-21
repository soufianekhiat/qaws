/*
 * Parameterization (qaws_mesh_param.h): the mesh cut along its seams into
 * wedges (a vertex's corners between seams or the border), LSCM as the
 * conformal energy's sparse system, ARAP's local-global rounds on the
 * wedges' cotangent Laplacian (qaws_mesh_smooth.c's sparse helpers).
 */

#include "../internal/qaws_internal_memory.h"
#include "../qaws_mesh_param.h"
#include "qaws_mesh_internal.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define PM_NONE 0xFFFFFFFFu

static int pm_cmp_u64(void const* p, void const* q)
{
	unsigned long long a = *(unsigned long long const*)p, b = *(unsigned long long const*)q;
	return a < b ? -1 : a > b;
}

static unsigned int pm_find(unsigned int* par, unsigned int i)
{
	while (par[i] != i)
	{
		par[i] = par[par[i]];
		i = par[i];
	}
	return i;
}

static void pm_join(unsigned int* par, unsigned int a, unsigned int b)
{
	a = pm_find(par, a);
	b = pm_find(par, b);
	if (a != b)
		par[a < b ? b : a] = a < b ? a : b;
}

static unsigned long long pm_key(unsigned int a, unsigned int b)
{
	return a < b ? ((unsigned long long)a << 32 | b) : ((unsigned long long)b << 32 | a);
}

/* a sparse symmetric matrix from triplets (row, column, value; merged), rows sorted */
typedef struct pm_trip
{
	unsigned int r, c;
	double v;
} pm_trip;

static int pm_cmp_trip(void const* p, void const* q)
{
	pm_trip const *a = (pm_trip const*)p, *b = (pm_trip const*)q;
	return a->r != b->r ? (a->r < b->r ? -1 : 1) : a->c != b->c ? (a->c < b->c ? -1 : 1) : 0;
}

static int pm_csr(sm_csr* A, pm_trip* T, unsigned int nt, unsigned int n)
{
	unsigned int i, m = 0;
	memset(A, 0, sizeof *A);
	qsort(T, nt, sizeof *T, pm_cmp_trip);
	A->n = n;
	A->rp = (unsigned int*)qaws_internal_zalloc((size_t)n + 1, sizeof(unsigned int));
	A->ci = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)nt + 1));
	A->v = (double*)qaws_internal_malloc(sizeof(double) * ((size_t)nt + 1));
	if (!A->rp || !A->ci || !A->v)
	{
		sm_csr_free(A);
		return 0;
	}
	for (i = 0; i < nt; i++)
	{
		if (m && T[i].r == T[i - 1].r && T[i].c == T[i - 1].c)
		{
			A->v[m - 1] += T[i].v;
			continue;
		}
		A->ci[m] = T[i].c;
		A->v[m++] = T[i].v;
		A->rp[T[i].r + 1]++;
	}
	for (i = 0; i < n; i++)
		A->rp[i + 1] += A->rp[i];
	return 1;
}

/* triangle f laid out isometrically in the plane: X (3 corners x 2) */
static void pm_layout(double const* x, unsigned int const* t, unsigned int f, double (*X)[2])
{
	double const *a = x + 3 * (size_t)t[3 * (size_t)f], *b = x + 3 * (size_t)t[3 * (size_t)f + 1], *c = x + 3 * (size_t)t[3 * (size_t)f + 2];
	double e[3], g[3], l, px, n2;
	unsigned int k;
	for (k = 0; k < 3; k++)
	{
		e[k] = b[k] - a[k];
		g[k] = c[k] - a[k];
	}
	l = sqrt(e[0] * e[0] + e[1] * e[1] + e[2] * e[2]);
	px = l > 0 ? (e[0] * g[0] + e[1] * g[1] + e[2] * g[2]) / l : 0;
	n2 = g[0] * g[0] + g[1] * g[1] + g[2] * g[2] - px * px;
	X[0][0] = X[0][1] = 0;
	X[1][0] = l;
	X[1][1] = 0;
	X[2][0] = px;
	X[2][1] = n2 > 0 ? sqrt(n2) : 0;
}

/* the cotangent at corner e of a triangle laid out in the plane (P: 3 corners x 2) */
static double pm_cot2(double const* P, unsigned int e)
{
	double u[2] = { P[2 * ((e + 1) % 3)] - P[2 * e], P[2 * ((e + 1) % 3) + 1] - P[2 * e + 1] };
	double w[2] = { P[2 * ((e + 2) % 3)] - P[2 * e], P[2 * ((e + 2) % 3) + 1] - P[2 * e + 1] };
	double cr = fabs(u[0] * w[1] - u[1] * w[0]);
	return cr > 0 ? (u[0] * w[0] + u[1] * w[1]) / cr : 0;
}

static double pm_angle(double const* o, double const* p, double const* q, unsigned int dim)
{
	double u[3] = { 0, 0, 0 }, w[3] = { 0, 0, 0 }, c[3];
	unsigned int k;
	for (k = 0; k < dim; k++)
	{
		u[k] = p[k] - o[k];
		w[k] = q[k] - o[k];
	}
	c[0] = u[1] * w[2] - u[2] * w[1];
	c[1] = u[2] * w[0] - u[0] * w[2];
	c[2] = u[0] * w[1] - u[1] * w[0];
	return atan2(sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]), u[0] * w[0] + u[1] * w[1] + u[2] * w[2]);
}

qaws_status qaws_mesh_3d_parameterize(qaws_scalar* out_uv, qaws_param_report* report, qaws_mesh_3d const* mesh, unsigned int const* seams,
	unsigned int ns, qaws_param_desc const* desc)
{
	qaws_param_method method = desc ? desc->method : QAWS_PARAM_LSCM;
	unsigned int iters = desc && desc->iterations ? desc->iterations : 10, nv, nt, n, i, j, k, W = 0, C = 0, nb = 0, ntr = 0;
	double *x = NULL, *xw = NULL, *uv = NULL, *mass = NULL, *rhs = NULL;
	unsigned int *t = NULL, *wedge = NULL, *par = NULL, *tpar = NULL, *tw = NULL, *chart = NULL, *pin_a = NULL, *pin_b = NULL;
	unsigned long long *he = NULL, *sk = NULL;
	unsigned char *bord = NULL, *move = NULL;
	qaws_scalar* pos = NULL;
	pm_trip* T = NULL;
	sm_csr K, Q;
	sm_system S;
	qaws_param_report local;
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	memset(&K, 0, sizeof K);
	memset(&Q, 0, sizeof Q);
	memset(&S, 0, sizeof S);
	if (!report)
		report = &local;
	memset(report, 0, sizeof *report);
	if (!out_uv || !mesh || (ns && !seams) || (unsigned int)method > (unsigned int)QAWS_PARAM_BFF)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (method == QAWS_PARAM_BFF)
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	nv = qaws_mesh_3d_get_vertex_count(mesh);
	nt = qaws_mesh_3d_get_triangle_count(mesh);
	n = 3 * nt;
	pos = (qaws_scalar*)qaws_internal_malloc(qaws_internal_mul3(sizeof(qaws_scalar), 3, (size_t)nv + 1));
	x = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 3, (size_t)nv + 1));
	t = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)n + 1));
	he = (unsigned long long*)qaws_internal_malloc(qaws_internal_mul2(2 * sizeof(unsigned long long), (size_t)n + 1));
	sk = (unsigned long long*)qaws_internal_malloc(sizeof(unsigned long long) * ((size_t)ns + 1));
	par = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)n + 1));
	tpar = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)nt + 1));
	bord = (unsigned char*)qaws_internal_zalloc((size_t)n + 1, 1);
	wedge = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)n + 1));
	if (!pos || !x || !t || !he || !sk || !par || !tpar || !bord || !wedge)
		goto done;
	if ((st = qaws_mesh_3d_get_positions(pos, mesh)) != QAWS_STATUS_OK || (st = qaws_mesh_3d_get_triangles(t, mesh)) != QAWS_STATUS_OK)
		goto done;
	st = QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; i < 3 * nv; i++)
		x[i] = (double)pos[i];
	/* half-edges by their edge; the seams' keys */
	for (i = 0; i < n; i++)
	{
		he[2 * (size_t)i] = pm_key(t[i], t[3 * (i / 3) + (i % 3 + 1) % 3]);
		he[2 * (size_t)i + 1] = i;
		par[i] = i;
	}
	for (i = 0; i < nt; i++)
		tpar[i] = i;
	qsort(he, n, 2 * sizeof(unsigned long long), pm_cmp_u64);
	for (i = 0; i < ns; i++)
	{
		unsigned long long key = pm_key(seams[2 * (size_t)i], seams[2 * (size_t)i + 1]);
		unsigned int lo = 0, hi = n;
		while (lo < hi)
		{
			unsigned int mid = (lo + hi) / 2;
			if (he[2 * (size_t)mid] < key)
				lo = mid + 1;
			else
				hi = mid;
		}
		if (seams[2 * (size_t)i] >= nv || seams[2 * (size_t)i + 1] >= nv || lo == n || he[2 * (size_t)lo] != key)
		{
			st = QAWS_STATUS_INVALID_ARGUMENT;
			goto done;
		}
		sk[i] = key;
	}
	qsort(sk, ns, sizeof *sk, pm_cmp_u64);
	/* edges of exactly two triangles turned opposite ways, not seams, join the corners at their ends and the triangles; the
	   others are the cut mesh's border */
	for (i = 0; i < n; i = j)
	{
		for (j = i + 1; j < n && he[2 * (size_t)j] == he[2 * (size_t)i]; j++)
			;
		{
			unsigned long long key = he[2 * (size_t)i];
			unsigned int lo = 0, hi = ns;
			int seam;
			while (lo < hi)
			{
				unsigned int mid = (lo + hi) / 2;
				if (sk[mid] < key)
					lo = mid + 1;
				else
					hi = mid;
			}
			seam = lo < ns && sk[lo] == key;
			if (j == i + 2 && !seam)
			{
				unsigned int c0 = (unsigned int)he[2 * (size_t)i + 1], c1 = (unsigned int)he[2 * (size_t)i + 3];
				unsigned int f = c0 / 3, g = c1 / 3, a0 = c0, b0 = 3 * f + (c0 % 3 + 1) % 3, a1 = c1, b1 = 3 * g + (c1 % 3 + 1) % 3;
				if (t[a0] == t[b1] && t[b0] == t[a1])
				{
					pm_join(par, a0, b1);
					pm_join(par, b0, a1);
					pm_join(tpar, f, g);
					continue;
				}
			}
			for (k = i; k < j; k++)
				bord[(unsigned int)he[2 * (size_t)k + 1]] = 1;
		}
	}
	qaws_internal_free(he);
	he = NULL;
	/* wedges (the corners' classes) and charts (the triangles') */
	for (i = 0; i < n; i++)
		wedge[i] = pm_find(par, i) == i ? W++ : PM_NONE;
	for (i = 0; i < n; i++)
		wedge[i] = wedge[pm_find(par, i)];
	chart = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)nt + 1));
	tw = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)nt + 1));
	xw = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 3, (size_t)W + 1));
	uv = (double*)qaws_internal_zalloc(qaws_internal_mul2(2, (size_t)W + 1), sizeof(double));
	if (!chart || !tw || !xw || !uv)
		goto done;
	for (i = 0; i < nt; i++)
		tw[i] = pm_find(tpar, i) == i ? C++ : PM_NONE;
	for (i = 0; i < nt; i++)
		chart[i] = tw[pm_find(tpar, i)];
	for (i = 0; i < n; i++)
		memcpy(xw + 3 * (size_t)wedge[i], x + 3 * (size_t)t[i], 3 * sizeof(double));
	report->charts = C;
	/* pins: per chart a border wedge and the border wedge farthest from it */
	pin_a = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)C + 1));
	pin_b = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)C + 1));
	if (!pin_a || !pin_b)
		goto done;
	for (i = 0; i < C; i++)
		pin_a[i] = pin_b[i] = PM_NONE;
	for (i = 0; i < n; i++)
		if (bord[i] && pin_a[chart[i / 3]] == PM_NONE)
			pin_a[chart[i / 3]] = wedge[i];
	for (i = 0; i < C; i++)
		if (pin_a[i] == PM_NONE)
		{
			st = QAWS_STATUS_INVALID_TOPOLOGY;   /* (a closed chart: nothing to open it) */
			goto done;
		}
	{
		double* far = (double*)qaws_internal_zalloc((size_t)C + 1, sizeof(double));
		if (!far)
			goto done;
		for (i = 0; i < n; i++)
			if (bord[i])
			{
				unsigned int c = chart[i / 3], w = wedge[i];
				double const *a = xw + 3 * (size_t)pin_a[c], *b = xw + 3 * (size_t)w;
				double d = (a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]);
				if (d > far[c])
				{
					far[c] = d;
					pin_b[c] = w;
				}
			}
		for (i = 0; i < C; i++)
		{
			uv[2 * (size_t)pin_a[i]] = uv[2 * (size_t)pin_a[i] + 1] = 0;
			if (pin_b[i] != PM_NONE)
			{
				uv[2 * (size_t)pin_b[i]] = sqrt(far[i]);
				uv[2 * (size_t)pin_b[i] + 1] = 0;
			}
		}
		qaws_internal_free(far);
	}
	/* LSCM: Q = (K on u and on v) - B, 1/2 x^T B x the signed area the border encloses; u and v of a wedge adjacent rows */
	mass = (double*)qaws_internal_malloc(sizeof(double) * ((size_t)W + 1));
	if (!mass)
		goto done;
	{
		unsigned int* twc = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)n + 1));
		if (!twc)
			goto done;
		for (i = 0; i < n; i++)
			twc[i] = wedge[i];
		k = (unsigned int)sm_laplacian(&K, mass, xw, W, twc, nt, 0, 0);
		qaws_internal_free(twc);
		if (!k)
			goto done;
	}
	for (i = 0; i < n; i++)
		nb += bord[i];
	T = (pm_trip*)qaws_internal_malloc(sizeof(pm_trip) * (2 * (size_t)K.rp[W] + 4 * (size_t)nb + 1));
	move = (unsigned char*)qaws_internal_malloc(2 * (size_t)W + 1);
	rhs = (double*)qaws_internal_zalloc(qaws_internal_mul2(2, (size_t)W + 1), sizeof(double));
	if (!T || !move || !rhs)
		goto done;
	for (i = 0; i < W; i++)
	{
		unsigned int p;
		for (p = K.rp[i]; p < K.rp[i + 1]; p++)
		{
			T[ntr].r = 2 * i;
			T[ntr].c = 2 * K.ci[p];
			T[ntr++].v = K.v[p];
			T[ntr].r = 2 * i + 1;
			T[ntr].c = 2 * K.ci[p] + 1;
			T[ntr++].v = K.v[p];
		}
	}
	for (i = 0; i < n; i++)
		if (bord[i])
		{
			unsigned int a = wedge[i], b = wedge[3 * (i / 3) + (i % 3 + 1) % 3];
			pm_trip const e[4] = { { 2 * a, 2 * b + 1, -0.5 }, { 2 * b + 1, 2 * a, -0.5 }, { 2 * b, 2 * a + 1, 0.5 }, { 2 * a + 1, 2 * b, 0.5 } };
			memcpy(T + ntr, e, sizeof e);
			ntr += 4;
		}
	if (!pm_csr(&Q, T, ntr, 2 * W))
		goto done;
	{
		double* pxyz = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 6, (size_t)W + 1));
		if (!pxyz)
			goto done;
		for (i = 0; i < W; i++)
			for (k = 0; k < 3; k++)
				pxyz[6 * (size_t)i + k] = pxyz[6 * (size_t)i + 3 + k] = xw[3 * (size_t)i + k];   /* (u and v of a wedge stay together) */
		for (i = 0; i < 2 * W; i++)
			move[i] = 1;
		for (i = 0; i < C; i++)
		{
			move[2 * (size_t)pin_a[i]] = move[2 * (size_t)pin_a[i] + 1] = 0;
			if (pin_b[i] != PM_NONE)
				move[2 * (size_t)pin_b[i]] = move[2 * (size_t)pin_b[i] + 1] = 0;
		}
		st = sm_system_factor(&S, &Q, pxyz, move);
		qaws_internal_free(pxyz);
		if (st == QAWS_STATUS_OK)
			st = sm_system_solve(&S, uv, rhs, 1);
		sm_system_free(&S);
		if (st != QAWS_STATUS_OK)
			goto done;
	}
	if (method == QAWS_PARAM_ARAP)
	{
		/* each triangle's map nearest a rotation of its own layout; the wedges' Laplacian with one pin per chart, factored once */
		double* X = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 6, (size_t)nt + 1));
		unsigned int it;
		if (!X)
			goto done;
		for (i = 0; i < nt; i++)
			pm_layout(x, t, i, (double(*)[2])(X + 6 * (size_t)i));
		for (i = 0; i < W; i++)
			move[i] = 1;
		for (i = 0; i < C; i++)
			move[pin_a[i]] = 0;
		st = sm_system_factor(&S, &K, xw, move);
		for (it = 0; it < iters && st == QAWS_STATUS_OK; it++)
		{
			memset(rhs, 0, sizeof(double) * 2 * W);
			for (i = 0; i < nt; i++)
			{
				double const* P = X + 6 * (size_t)i;
				double s00 = 0, s01 = 0, s10 = 0, s11 = 0, th, c, s;
				unsigned int e;
				/* the covariance over the edges (opposite corner e: from e + 1 to e + 2), weighted by half the cotangent */
				for (e = 0; e < 3; e++)
				{
					unsigned int p = (e + 1) % 3, r = (e + 2) % 3;
					double const *up = uv + 2 * (size_t)wedge[3 * i + p], *ur = uv + 2 * (size_t)wedge[3 * i + r];
					double w = 0.5 * pm_cot2(P, e), du[2] = { up[0] - ur[0], up[1] - ur[1] }, dx[2] = { P[2 * p] - P[2 * r], P[2 * p + 1] - P[2 * r + 1] };
					s00 += w * du[0] * dx[0];
					s01 += w * du[0] * dx[1];
					s10 += w * du[1] * dx[0];
					s11 += w * du[1] * dx[1];
				}
				th = atan2(s10 - s01, s00 + s11);
				c = cos(th);
				s = sin(th);
				for (e = 0; e < 3; e++)
				{
					unsigned int p = (e + 1) % 3, r = (e + 2) % 3, a = wedge[3 * i + p], b = wedge[3 * i + r];
					double w = 0.5 * pm_cot2(P, e), dx[2] = { P[2 * p] - P[2 * r], P[2 * p + 1] - P[2 * r + 1] };
					double rx = c * dx[0] - s * dx[1], ry = s * dx[0] + c * dx[1];
					rhs[2 * (size_t)a] += w * rx;
					rhs[2 * (size_t)a + 1] += w * ry;
					rhs[2 * (size_t)b] -= w * rx;
					rhs[2 * (size_t)b + 1] -= w * ry;
				}
			}
			st = sm_system_solve(&S, uv, rhs, 2);
		}
		sm_system_free(&S);
		qaws_internal_free(X);
		if (st != QAWS_STATUS_OK)
			goto done;
	}
	{
		/* the charts side by side along u; the report */
		double *lo = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 4, (size_t)C + 1)), off = 0, a3 = 0, a2 = 0, dang = 0, gap;
		double bmin[3] = { HUGE_VAL, HUGE_VAL, HUGE_VAL }, bmax[3] = { -HUGE_VAL, -HUGE_VAL, -HUGE_VAL };
		if (!lo)
			goto done;
		for (i = 0; i < C; i++)
		{
			lo[4 * (size_t)i] = lo[4 * (size_t)i + 1] = HUGE_VAL;
			lo[4 * (size_t)i + 2] = lo[4 * (size_t)i + 3] = -HUGE_VAL;
		}
		for (i = 0; i < n; i++)
		{
			unsigned int c = chart[i / 3];
			double const* p = uv + 2 * (size_t)wedge[i];
			lo[4 * (size_t)c] = fmin(lo[4 * (size_t)c], p[0]);
			lo[4 * (size_t)c + 1] = fmin(lo[4 * (size_t)c + 1], p[1]);
			lo[4 * (size_t)c + 2] = fmax(lo[4 * (size_t)c + 2], p[0]);
			for (k = 0; k < 3; k++)
			{
				bmin[k] = fmin(bmin[k], x[3 * (size_t)t[i] + k]);
				bmax[k] = fmax(bmax[k], x[3 * (size_t)t[i] + k]);
			}
		}
		gap = 0.02 * sqrt((bmax[0] - bmin[0]) * (bmax[0] - bmin[0]) + (bmax[1] - bmin[1]) * (bmax[1] - bmin[1]) + (bmax[2] - bmin[2]) * (bmax[2] - bmin[2]));
		for (i = 0; i < C; i++)
		{
			double w = lo[4 * (size_t)i + 2] - lo[4 * (size_t)i];
			lo[4 * (size_t)i + 2] = off - lo[4 * (size_t)i];   /* (the chart's shift in u) */
			lo[4 * (size_t)i + 3] = -lo[4 * (size_t)i + 1];
			off += w + gap;
		}
		for (i = 0; i < n; i++)
		{
			unsigned int c = chart[i / 3];
			out_uv[2 * (size_t)i] = (qaws_scalar)(uv[2 * (size_t)wedge[i]] + lo[4 * (size_t)c + 2]);
			out_uv[2 * (size_t)i + 1] = (qaws_scalar)(uv[2 * (size_t)wedge[i] + 1] + lo[4 * (size_t)c + 3]);
		}
		qaws_internal_free(lo);
		/* distortion: areas' shares, corners' angles */
		for (i = 0; i < nt; i++)
		{
			double const *p0 = uv + 2 * (size_t)wedge[3 * i], *p1 = uv + 2 * (size_t)wedge[3 * i + 1], *p2 = uv + 2 * (size_t)wedge[3 * i + 2];
			double X[3][2], au = 0.5 * ((p1[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (p1[1] - p0[1]));
			pm_layout(x, t, i, X);
			a3 += 0.5 * X[1][0] * X[2][1];
			a2 += au;
			report->flipped += au <= 0;
		}
		for (i = 0; i < nt; i++)
		{
			double const *p0 = uv + 2 * (size_t)wedge[3 * i], *p1 = uv + 2 * (size_t)wedge[3 * i + 1], *p2 = uv + 2 * (size_t)wedge[3 * i + 2];
			double X[3][2], au = 0.5 * ((p1[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (p1[1] - p0[1]));
			double const* q[3] = { p0, p1, p2 };
			pm_layout(x, t, i, X);
			if (a2 != 0 && a3 > 0)
				report->area_distortion += (qaws_scalar)fabs(au / a2 - 0.5 * X[1][0] * X[2][1] / a3);
			for (k = 0; k < 3; k++)
				dang += fabs(pm_angle(q[k], q[(k + 1) % 3], q[(k + 2) % 3], 2) - pm_angle(X[k], X[(k + 1) % 3], X[(k + 2) % 3], 2));
		}
		report->angle_distortion = (qaws_scalar)(n ? dang / n : 0);
	}
	st = QAWS_STATUS_OK;
done:
	sm_system_free(&S);
	sm_csr_free(&K);
	sm_csr_free(&Q);
	qaws_internal_free(pos); qaws_internal_free(x); qaws_internal_free(t); qaws_internal_free(he); qaws_internal_free(sk);
	qaws_internal_free(par); qaws_internal_free(tpar); qaws_internal_free(bord); qaws_internal_free(wedge); qaws_internal_free(chart);
	qaws_internal_free(tw); qaws_internal_free(xw); qaws_internal_free(uv); qaws_internal_free(mass); qaws_internal_free(rhs);
	qaws_internal_free(pin_a); qaws_internal_free(pin_b); qaws_internal_free(T); qaws_internal_free(move);
	return st;
}
