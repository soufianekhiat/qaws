/*
 * Smoothing and fairing (qaws_mesh_smooth.h): the cotangent (or uniform)
 * Laplacian assembled as a sparse matrix, explicit steps on it, or the
 * implicit flows and fairing systems restricted to the moving vertices
 * and solved by the sparse LDL^T.
 */

#include "../internal/qaws_internal_memory.h"
#include "../internal/qaws_internal_ldlt.h"
#include "../qaws_mesh_smooth.h"
#include "qaws_mesh_internal.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define SM_NONE 0xFFFFFFFFu

void sm_csr_free(sm_csr* A)
{
	qaws_internal_free(A->rp);
	qaws_internal_free(A->ci);
	qaws_internal_free(A->v);
	memset(A, 0, sizeof *A);
}

static int sm_cmp_key(void const* p, void const* q)
{
	unsigned long long a = *(unsigned long long const*)p, b = *(unsigned long long const*)q;
	return a < b ? -1 : a > b;
}

/* K (symmetric, positive semi-definite: K_ii the sum of the weights, K_ij minus the weight) and the lumped masses (a third of
   the areas round each vertex; 1 for uniform weights). Cotangent weights are half the cotangents facing each edge, clamped
   at 0 with clamp; uniform ones 1 per edge */
int sm_laplacian(sm_csr* K, double* mass, double const* x, unsigned int nv, unsigned int const* t, unsigned int nt, int uniform,
	int clamp)
{
	unsigned long long* he = (unsigned long long*)qaws_internal_malloc(qaws_internal_mul3(2 * sizeof(unsigned long long), 3, (size_t)nt + 1));
	unsigned int *cnt = NULL, i, j, k, f, ne = 0, n = 3 * nt;
	memset(K, 0, sizeof *K);
	if (!he)
		return 0;
	for (i = 0; i < nv; i++)
		mass[i] = uniform ? 1 : 0;
	for (f = 0; f < nt; f++)
	{
		double const *a = x + 3 * (size_t)t[3 * (size_t)f], *b = x + 3 * (size_t)t[3 * (size_t)f + 1], *c = x + 3 * (size_t)t[3 * (size_t)f + 2];
		double u[3], w[3], q[3], a2;
		for (k = 0; k < 3; k++)
		{
			u[k] = b[k] - a[k];
			w[k] = c[k] - a[k];
		}
		q[0] = u[1] * w[2] - u[2] * w[1];
		q[1] = u[2] * w[0] - u[0] * w[2];
		q[2] = u[0] * w[1] - u[1] * w[0];
		a2 = sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2]);
		for (k = 0; k < 3; k++)
		{
			unsigned long long p = t[3 * (size_t)f + (k + 1) % 3], r = t[3 * (size_t)f + (k + 2) % 3];
			double const *o = x + 3 * (size_t)t[3 * (size_t)f + k], *e1 = x + 3 * (size_t)p, *e2 = x + 3 * (size_t)r;
			double cw;
			double d1[3] = { e1[0] - o[0], e1[1] - o[1], e1[2] - o[2] }, d2[3] = { e2[0] - o[0], e2[1] - o[1], e2[2] - o[2] };
			/* (the corner's cotangent weighs the edge facing it) */
			cw = a2 > 0 ? 0.5 * (d1[0] * d2[0] + d1[1] * d2[1] + d1[2] * d2[2]) / a2 : 0;
			he[2 * (3 * (size_t)f + k)] = p < r ? (p << 32 | r) : (r << 32 | p);
			memcpy(he + 2 * (3 * (size_t)f + k) + 1, &cw, sizeof cw);
			if (!uniform)
				mass[t[3 * (size_t)f + k]] += a2 / 6;
		}
	}
	qsort(he, n, 2 * sizeof(unsigned long long), sm_cmp_key);
	/* edges merged in place: key, weight */
	for (i = 0; i < n; i = j)
	{
		double w = 0;
		for (j = i; j < n && he[2 * (size_t)j] == he[2 * (size_t)i]; j++)
		{
			double wj;
			memcpy(&wj, he + 2 * (size_t)j + 1, sizeof wj);
			w += wj;
		}
		w = uniform ? 1 : clamp && w < 0 ? 0 : w;
		he[2 * (size_t)ne] = he[2 * (size_t)i];
		memcpy(he + 2 * (size_t)ne + 1, &w, sizeof w);
		ne++;
	}
	for (i = 0; i < nv; i++)
		if (!(mass[i] > 0))
			mass[i] = 1;   /* (a vertex on no triangle: its row is the identity's) */
	K->n = nv;
	K->rp = (unsigned int*)qaws_internal_zalloc((size_t)nv + 2, sizeof(unsigned int));
	cnt = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)nv + 1));
	if (!K->rp || !cnt)
		goto fail;
	for (i = 0; i < nv; i++)
		K->rp[i + 1] = 1;
	for (i = 0; i < ne; i++)
	{
		K->rp[(unsigned int)(he[2 * (size_t)i] >> 32) + 1]++;
		K->rp[(unsigned int)he[2 * (size_t)i] + 1]++;
	}
	for (i = 0; i < nv; i++)
		K->rp[i + 1] += K->rp[i];
	K->ci = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)K->rp[nv] + 1));
	K->v = (double*)qaws_internal_zalloc((size_t)K->rp[nv] + 1, sizeof(double));
	if (!K->ci || !K->v)
		goto fail;
	for (i = 0; i < nv; i++)
	{
		cnt[i] = K->rp[i] + 1;
		K->ci[K->rp[i]] = i;
	}
	for (i = 0; i < ne; i++)
	{
		unsigned int a = (unsigned int)(he[2 * (size_t)i] >> 32), b = (unsigned int)he[2 * (size_t)i];
		double w;
		memcpy(&w, he + 2 * (size_t)i + 1, sizeof w);
		K->ci[cnt[a]] = b;
		K->v[cnt[a]++] = -w;
		K->ci[cnt[b]] = a;
		K->v[cnt[b]++] = -w;
		K->v[K->rp[a]] += w;
		K->v[K->rp[b]] += w;
	}
	/* columns sorted in each row (the diagonal first was a placeholder) */
	for (i = 0; i < nv; i++)
	{
		unsigned int p, q;
		for (p = K->rp[i] + 1; p < K->rp[i + 1]; p++)
		{
			unsigned int c = K->ci[p];
			double w = K->v[p];
			for (q = p; q > K->rp[i] && K->ci[q - 1] > c; q--)
			{
				K->ci[q] = K->ci[q - 1];
				K->v[q] = K->v[q - 1];
			}
			K->ci[q] = c;
			K->v[q] = w;
		}
	}
	qaws_internal_free(he);
	qaws_internal_free(cnt);
	return 1;
fail:
	qaws_internal_free(he);
	qaws_internal_free(cnt);
	sm_csr_free(K);
	return 0;
}

/* C = A diag(d) B (all n x n; d NULL: the identity), columns sorted */
int sm_mul(sm_csr* C, sm_csr const* A, double const* d, sm_csr const* B)
{
	unsigned int n = A->n, i, p, q, cap = 0, nnz = 0;
	unsigned int* mark = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)n + 1));
	double* acc = (double*)qaws_internal_zalloc((size_t)n + 1, sizeof(double));
	memset(C, 0, sizeof *C);
	C->n = n;
	C->rp = (unsigned int*)qaws_internal_zalloc((size_t)n + 2, sizeof(unsigned int));
	if (!mark || !acc || !C->rp)
		goto fail;
	for (i = 0; i < n; i++)
		mark[i] = SM_NONE;
	for (i = 0; i < n; i++)
	{
		unsigned int start = nnz;
		for (p = A->rp[i]; p < A->rp[i + 1]; p++)
		{
			unsigned int j = A->ci[p];
			double a = A->v[p] * (d ? d[j] : 1);
			for (q = B->rp[j]; q < B->rp[j + 1]; q++)
			{
				unsigned int k = B->ci[q];
				if (mark[k] != i)
				{
					mark[k] = i;
					acc[k] = 0;
					if (nnz == cap)
					{
						unsigned int nc = cap ? 2 * cap : 16 * n + 16;
						unsigned int* ci = (unsigned int*)qaws_internal_realloc(C->ci, qaws_internal_mul2(sizeof(unsigned int), nc));
						double* cv;
						if (!ci)
							goto fail;
						C->ci = ci;
						cv = (double*)qaws_internal_realloc(C->v, qaws_internal_mul2(sizeof(double), nc));
						if (!cv)
							goto fail;
						C->v = cv;
						cap = nc;
					}
					C->ci[nnz++] = k;
				}
				acc[k] += a * B->v[q];
			}
		}
		/* (sorted: a row has a few rings' worth of columns) */
		for (p = start + 1; p < nnz; p++)
		{
			unsigned int c = C->ci[p];
			for (q = p; q > start && C->ci[q - 1] > c; q--)
				C->ci[q] = C->ci[q - 1];
			C->ci[q] = c;
		}
		C->rp[i + 1] = nnz;
		for (p = start; p < nnz; p++)
			C->v[p] = acc[C->ci[p]];
	}
	qaws_internal_free(mark);
	qaws_internal_free(acc);
	return 1;
fail:
	qaws_internal_free(mark);
	qaws_internal_free(acc);
	sm_csr_free(C);
	return 0;
}

qaws_status sm_system_factor(sm_system* S, sm_csr const* A, double const* xyz, unsigned char const* move)
{
	unsigned int n = A->n, i, p, k, *rp = NULL, *ci = NULL;
	double *v = NULL, *pt = NULL;
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	memset(S, 0, sizeof *S);
	S->A = A;
	S->row = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)n + 1));
	if (!S->row)
		return st;
	for (i = 0; i < n; i++)
		S->row[i] = move[i] ? S->m++ : SM_NONE;
	if (!S->m)
		return QAWS_STATUS_OK;
	rp = (unsigned int*)qaws_internal_zalloc((size_t)S->m + 1, sizeof(unsigned int));
	ci = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)A->rp[n] + 1));
	v = (double*)qaws_internal_malloc(sizeof(double) * ((size_t)A->rp[n] + 1));
	pt = xyz ? (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 3, (size_t)S->m)) : NULL;
	if (!rp || !ci || !v || (xyz && !pt))
		goto done;
	for (i = 0; i < n; i++)
	{
		unsigned int r = S->row[i];
		if (r == SM_NONE)
			continue;
		rp[r + 1] = rp[r];
		if (pt)
			for (k = 0; k < 3; k++)
				pt[3 * (size_t)r + k] = xyz[3 * (size_t)i + k];
		for (p = A->rp[i]; p < A->rp[i + 1]; p++)
		{
			unsigned int c = S->row[A->ci[p]];
			if (c != SM_NONE)
			{
				ci[rp[r + 1]] = c;
				v[rp[r + 1]++] = A->v[p];
			}
		}
	}
	st = qaws_internal_ldlt_factor(&S->L, rp, ci, v, pt, S->m, 0);
done:
	qaws_internal_free(rp); qaws_internal_free(ci); qaws_internal_free(v); qaws_internal_free(pt);
	return st;
}

qaws_status sm_system_solve(sm_system const* S, double* x, double const* rhs, unsigned int dim)
{
	sm_csr const* A = S->A;
	qaws_internal_ldlt local = S->L;
	unsigned int n = A->n, i, p, k;
	double* b;
	if (!S->m)
		return QAWS_STATUS_OK;
	b = (double*)qaws_internal_malloc(sizeof(double) * ((size_t)S->m + 1));
	local.w = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), 2 * (size_t)S->m + 1));
	if (!b || !local.w)
	{
		qaws_internal_free(b);
		qaws_internal_free(local.w);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (k = 0; k < dim; k++)
	{
		for (i = 0; i < n; i++)
		{
			unsigned int r = S->row[i];
			if (r == SM_NONE)
				continue;
			b[r] = rhs[(size_t)dim * i + k];
			for (p = A->rp[i]; p < A->rp[i + 1]; p++)
				if (S->row[A->ci[p]] == SM_NONE)
					b[r] -= A->v[p] * x[(size_t)dim * A->ci[p] + k];
		}
		qaws_internal_ldlt_solve(b, &local, b);
		for (i = 0; i < n; i++)
			if (S->row[i] != SM_NONE)
				x[(size_t)dim * i + k] = b[S->row[i]];
	}
	qaws_internal_free(b);
	qaws_internal_free(local.w);
	return QAWS_STATUS_OK;
}

void sm_system_free(sm_system* S)
{
	qaws_internal_ldlt_free(&S->L);
	qaws_internal_free(S->row);
	memset(S, 0, sizeof *S);
}

/* A x = rhs on the rows where move (3 coordinates), factored for this one solve */
static qaws_status sm_solve(sm_csr const* A, double* x, double const* rhs, unsigned char const* move)
{
	sm_system S;
	qaws_status st = sm_system_factor(&S, A, x, move);
	if (st == QAWS_STATUS_OK)
		st = sm_system_solve(&S, x, rhs, 3);
	sm_system_free(&S);
	return st;
}

static double sm_volume(double const* x, unsigned int const* t, unsigned int nt)
{
	double vol = 0;
	unsigned int f;
	for (f = 0; f < nt; f++)
	{
		double const *a = x + 3 * (size_t)t[3 * (size_t)f], *b = x + 3 * (size_t)t[3 * (size_t)f + 1], *c = x + 3 * (size_t)t[3 * (size_t)f + 2];
		vol += a[0] * (b[1] * c[2] - b[2] * c[1]) + a[1] * (b[2] * c[0] - b[0] * c[2]) + a[2] * (b[0] * c[1] - b[1] * c[0]);
	}
	return vol / 6;
}

/* the mesh read as doubles and triangles (caller frees) */
static qaws_status sm_read(qaws_mesh_3d const* mesh, double** x, unsigned int** t, unsigned int* nv, unsigned int* nt)
{
	qaws_scalar* p;
	unsigned int i;
	qaws_status st;
	*nv = qaws_mesh_3d_get_vertex_count(mesh);
	*nt = qaws_mesh_3d_get_triangle_count(mesh);
	p = (qaws_scalar*)qaws_internal_malloc(qaws_internal_mul3(sizeof(qaws_scalar), 3, (size_t)*nv + 1));
	*x = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 3, (size_t)*nv + 1));
	*t = (unsigned int*)qaws_internal_malloc(qaws_internal_mul3(sizeof(unsigned int), 3, (size_t)*nt + 1));
	if (!p || !*x || !*t)
	{
		qaws_internal_free(p);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	st = qaws_mesh_3d_get_positions(p, mesh);
	if (st == QAWS_STATUS_OK)
		st = qaws_mesh_3d_get_triangles(*t, mesh);
	for (i = 0; i < 3 * *nv; i++)
		(*x)[i] = (double)p[i];
	qaws_internal_free(p);
	return st;
}

static qaws_status sm_write(qaws_mesh_3d* mesh, double const* x, unsigned int nv)
{
	qaws_scalar* p = (qaws_scalar*)qaws_internal_malloc(qaws_internal_mul3(sizeof(qaws_scalar), 3, (size_t)nv + 1));
	unsigned int i;
	qaws_status st;
	if (!p)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; i < 3 * nv; i++)
		p[i] = (qaws_scalar)x[i];
	st = qaws_mesh_3d_set_positions(mesh, p);
	qaws_internal_free(p);
	return st;
}

/* the line Laplacian: uniform, over the line edges (two line neighbours per line vertex, ln) */
static int sm_line_laplacian(sm_csr* K, unsigned int const* ln, unsigned char const* cls, unsigned int nv)
{
	unsigned int i, k, p;
	memset(K, 0, sizeof *K);
	K->n = nv;
	K->rp = (unsigned int*)qaws_internal_zalloc((size_t)nv + 1, sizeof(unsigned int));
	K->ci = (unsigned int*)qaws_internal_malloc(qaws_internal_mul3(sizeof(unsigned int), 3, (size_t)nv + 1));
	K->v = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 3, (size_t)nv + 1));
	if (!K->rp || !K->ci || !K->v)
	{
		sm_csr_free(K);
		return 0;
	}
	for (i = 0, p = 0; i < nv; i++)
	{
		unsigned int c[3];
		double w[3];
		c[0] = i;
		w[0] = cls[i] == 1 ? 2 : 0;
		for (k = 0; k < 2; k++)
		{
			c[1 + k] = cls[i] == 1 ? ln[2 * (size_t)i + k] : i;
			w[1 + k] = cls[i] == 1 ? -1 : 0;
		}
		/* (three entries, sorted; a held vertex's row is empty: it is never solved for) */
		for (k = 0; k < 3; k++)
		{
			unsigned int m, best = k;
			for (m = k + 1; m < 3; m++)
				if (c[m] < c[best])
					best = m;
			if (best != k)
			{
				unsigned int tc = c[k];
				double tw = w[k];
				c[k] = c[best];
				w[k] = w[best];
				c[best] = tc;
				w[best] = tw;
			}
		}
		for (k = 0; k < 3; k++)
			if (w[k] != 0)
			{
				K->ci[p] = c[k];
				K->v[p++] = w[k];
			}
		K->rp[i + 1] = p;
	}
	return 1;
}

/* BILATERAL: the triangles' normals filtered (Zheng, Fu, Au, Tai 2011), each the mean of the normals of the triangles
   sharing a vertex with it weighted by their areas, by a Gaussian of the distance between centroids (sigma: their mean
   spacing) and by a Gaussian of the normals' difference (sigma_r: the range), iters times; then the vertices moved to the
   filtered triangles' planes, 20 times: each to the point nearest them in least squares (where they span; Sun et al.
   2007 moved by their mean, which rounds corners). Held
   vertices stay, a line vertex moves along its line. 0 on allocation failure. */
#define SM_PLANES 0.05   /* BILATERAL: a direction the planes fix with a twentieth of the strongest's weight */

static int sm_bilateral(double* x, unsigned int nv, unsigned int const* t, unsigned int nt, unsigned char const* cls, unsigned int const* ln,
	unsigned int iters, double sigma_r)
{
	unsigned int *v0 = (unsigned int*)qaws_internal_zalloc((size_t)nv + 2, sizeof(unsigned int)), *vf = (unsigned int*)qaws_internal_malloc(qaws_internal_mul3(sizeof(unsigned int), 3, (size_t)nt + 1));
	unsigned int* stamp = (unsigned int*)qaws_internal_zalloc((size_t)nt + 1, sizeof(unsigned int));
	double *n = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 3, (size_t)nt + 1)), *n2 = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 3, (size_t)nt + 1));
	double *c = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 3, (size_t)nt + 1)), *A = (double*)qaws_internal_malloc(sizeof(double) * ((size_t)nt + 1));
	double ss = 0;
	unsigned int i, k, it, cnt = 0, mark = 0;
	int ok = v0 && vf && stamp && n && n2 && c && A;
	if (!ok)
		goto done;
	for (i = 0; i < 3 * nt; i++) v0[t[i] + 1]++;
	for (i = 0; i < nv; i++) v0[i + 1] += v0[i];
	{
		unsigned int* fill = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)nv + 1));
		if (!fill)
		{
			ok = 0;
			goto done;
		}
		memcpy(fill, v0, sizeof(unsigned int) * nv);
		for (i = 0; i < 3 * nt; i++) vf[fill[t[i]]++] = i / 3;
		qaws_internal_free(fill);
	}
	/* the triangles' normals, centroids and areas */
	for (i = 0; i < nt; i++)
	{
		double const *a = x + 3 * (size_t)t[3 * i], *b = x + 3 * (size_t)t[3 * i + 1], *d = x + 3 * (size_t)t[3 * i + 2];
		double u[3], w[3], l;
		for (k = 0; k < 3; k++) { u[k] = b[k] - a[k]; w[k] = d[k] - a[k]; c[3 * (size_t)i + k] = (a[k] + b[k] + d[k]) / 3; }
		n[3 * (size_t)i] = u[1] * w[2] - u[2] * w[1];
		n[3 * (size_t)i + 1] = u[2] * w[0] - u[0] * w[2];
		n[3 * (size_t)i + 2] = u[0] * w[1] - u[1] * w[0];
		l = sqrt(n[3 * (size_t)i] * n[3 * (size_t)i] + n[3 * (size_t)i + 1] * n[3 * (size_t)i + 1] + n[3 * (size_t)i + 2] * n[3 * (size_t)i + 2]);
		A[i] = 0.5 * l;
		for (k = 0; k < 3; k++) n[3 * (size_t)i + k] = l > 0 ? n[3 * (size_t)i + k] / l : 0;
	}
	/* sigma: the mean distance between the centroids of triangles sharing a vertex */
	for (i = 0; i < nt; i++)
	{
		unsigned int q, j;
		mark++;
		stamp[i] = mark;
		for (q = 0; q < 3; q++)
			for (j = v0[t[3 * i + q]]; j < v0[t[3 * i + q] + 1]; j++)
			{
				unsigned int g = vf[j];
				if (stamp[g] == mark)
					continue;
				stamp[g] = mark;
				ss += sqrt((c[3 * (size_t)g] - c[3 * (size_t)i]) * (c[3 * (size_t)g] - c[3 * (size_t)i]) + (c[3 * (size_t)g + 1] - c[3 * (size_t)i + 1]) * (c[3 * (size_t)g + 1] - c[3 * (size_t)i + 1]) +
					(c[3 * (size_t)g + 2] - c[3 * (size_t)i + 2]) * (c[3 * (size_t)g + 2] - c[3 * (size_t)i + 2]));
				cnt++;
			}
	}
	ss = cnt ? ss / cnt : 1;
	for (it = 0; it < iters; it++)
	{
		for (i = 0; i < nt; i++)
		{
			double s[3] = { 0, 0, 0 }, l;
			unsigned int q, j;
			mark++;
			for (q = 0; q < 3; q++)
				for (j = v0[t[3 * i + q]]; j < v0[t[3 * i + q] + 1]; j++)
				{
					unsigned int g = vf[j];
					double dc = 0, dn = 0, w;
					if (stamp[g] == mark)
						continue;
					stamp[g] = mark;
					for (k = 0; k < 3; k++)
					{
						dc += (c[3 * (size_t)g + k] - c[3 * (size_t)i + k]) * (c[3 * (size_t)g + k] - c[3 * (size_t)i + k]);
						dn += (n[3 * (size_t)g + k] - n[3 * (size_t)i + k]) * (n[3 * (size_t)g + k] - n[3 * (size_t)i + k]);
					}
					w = A[g] * exp(-dc / (2 * ss * ss)) * exp(-dn / (2 * sigma_r * sigma_r));
					for (k = 0; k < 3; k++) s[k] += w * n[3 * (size_t)g + k];
				}
			l = sqrt(s[0] * s[0] + s[1] * s[1] + s[2] * s[2]);
			for (k = 0; k < 3; k++) n2[3 * (size_t)i + k] = l > 0 ? s[k] / l : n[3 * (size_t)i + k];
		}
		memcpy(n, n2, sizeof(double) * 3 * nt);
	}
	/* the vertices onto the filtered planes */
	for (it = 0; it < 20; it++)
	{
		for (i = 0; i < nt; i++)
			for (k = 0; k < 3; k++)
				c[3 * (size_t)i + k] = (x[3 * (size_t)t[3 * i] + k] + x[3 * (size_t)t[3 * i + 1] + k] + x[3 * (size_t)t[3 * i + 2] + k]) / 3;
		for (i = 0; i < nv; i++)
		{
			double d[3] = { 0, 0, 0 }, M[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 }, V[9], b[3] = { 0, 0, 0 }, lmax = 0;
			unsigned int j, a2;
			if (cls[i] == 2 || v0[i] == v0[i + 1])
				continue;
			/* the point nearest the filtered planes in least squares (sum n n^T) d = sum n (n . (c - x)): solved where the
			   planes span (a pseudo-inverse), so a flat vertex moves along its normal, an edge's onto the planes' line and a
			   corner's to where they meet, kept sharp */
			for (j = v0[i]; j < v0[i + 1]; j++)
			{
				unsigned int g = vf[j], a1;
				double const* nn = n + 3 * (size_t)g;
				double e = 0;
				for (k = 0; k < 3; k++) e += nn[k] * (c[3 * (size_t)g + k] - x[3 * (size_t)i + k]);
				for (k = 0; k < 3; k++)
				{
					b[k] += nn[k] * e;
					for (a1 = 0; a1 < 3; a1++) M[3 * k + a1] += nn[k] * nn[a1];
				}
			}
			wt_sym_eigen3(M, V);
			for (k = 0; k < 3; k++) lmax = fmax(lmax, fabs(M[4 * k]));
			for (a2 = 0; a2 < 3; a2++)
			{
				double lk = M[4 * a2], pr;
				if (!(fabs(lk) > SM_PLANES * lmax))
					continue;   /* (a direction the planes do not fix: not moved along) */
				pr = (V[a2] * b[0] + V[3 + a2] * b[1] + V[6 + a2] * b[2]) / lk;
				for (k = 0; k < 3; k++) d[k] += pr * V[3 * k + a2];
			}
			if (cls[i] == 1)
			{
				/* along the line: the move's part along the line's chord through the vertex */
				double u[3], l = 0, dd = 0;
				for (k = 0; k < 3; k++) { u[k] = x[3 * (size_t)ln[2 * (size_t)i] + k] - x[3 * (size_t)ln[2 * (size_t)i + 1] + k]; l += u[k] * u[k]; }
				l = sqrt(l);
				for (k = 0; k < 3; k++) { u[k] = l > 0 ? u[k] / l : 0; dd += d[k] * u[k]; }
				for (k = 0; k < 3; k++) d[k] = dd * u[k];
			}
			for (k = 0; k < 3; k++) x[3 * (size_t)i + k] += d[k];
		}
	}
done:
	qaws_internal_free(v0); qaws_internal_free(vf); qaws_internal_free(stamp); qaws_internal_free(n); qaws_internal_free(n2); qaws_internal_free(c); qaws_internal_free(A);
	return ok;
}

qaws_status qaws_mesh_3d_smooth(qaws_mesh_3d* mesh, qaws_smooth_desc const* desc)
{
	qaws_smooth_desc D;
	double *x = NULL, *y = NULL, *mass = NULL, *rhs = NULL, h = 0, vol0 = 0, step;
	unsigned int *t = NULL, *ln = NULL, *fe = NULL, nv, nt, nfe, i, k, it, iters;
	unsigned char *cls = NULL, *move = NULL;
	int nothing_held = 1;
	sm_csr K, A, B;
	qaws_status st;
	memset(&K, 0, sizeof K);
	memset(&A, 0, sizeof A);
	memset(&B, 0, sizeof B);
	if (!mesh)
		return QAWS_STATUS_INVALID_ARGUMENT;
	memset(&D, 0, sizeof D);
	if (desc)
		D = *desc;
	if ((unsigned int)D.method > (unsigned int)QAWS_SMOOTH_BILATERAL || (unsigned int)D.weights > (unsigned int)QAWS_SMOOTH_WEIGHTS_UNIFORM
		|| !(D.step >= 0) || !isfinite((double)D.mu))
		return QAWS_STATUS_INVALID_ARGUMENT;
	iters = D.iterations ? D.iterations : D.method == QAWS_SMOOTH_BILATERAL ? 20 : D.method <= QAWS_SMOOTH_TAUBIN ? 10 : 1;
	step = D.step > 0 ? (double)D.step : D.method == QAWS_SMOOTH_BILATERAL ? 0.2 : D.method <= QAWS_SMOOTH_TAUBIN ? 0.5 : 1;
	st = sm_read(mesh, &x, &t, &nv, &nt);
	if (st != QAWS_STATUS_OK)
		goto done;
	st = QAWS_STATUS_ALLOCATION_FAILURE;
	y = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 3, (size_t)nv + 1));
	mass = (double*)qaws_internal_malloc(sizeof(double) * ((size_t)nv + 1));
	rhs = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 3, (size_t)nv + 1));
	cls = (unsigned char*)qaws_internal_zalloc((size_t)nv + 1, 1);
	move = (unsigned char*)qaws_internal_zalloc((size_t)nv + 1, 1);
	ln = (unsigned int*)qaws_internal_malloc(qaws_internal_mul3(sizeof(unsigned int), 2, (size_t)nv + 1));
	nfe = qaws_mesh_3d_get_feature_edge_count(mesh);
	fe = (unsigned int*)qaws_internal_malloc(qaws_internal_mul3(sizeof(unsigned int), 2, (size_t)nfe + 1));
	if (!y || !mass || !rhs || !cls || !move || !ln || !fe)
		goto done;
	if (nfe && (st = qaws_mesh_3d_get_feature_edges(fe, mesh)) != QAWS_STATUS_OK)
		goto done;
	st = QAWS_STATUS_ALLOCATION_FAILURE;
	{
		/* the lines: feature edges, and the border's with QAWS_SMOOTH_MOVE_BORDER (each edge once); per vertex its line
		   neighbours (deg: how many); class 0 free, 1 on one line, 2 held */
		unsigned long long* he = (unsigned long long*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned long long), 3 * (size_t)nt + nfe + 1));
		unsigned int* deg = (unsigned int*)qaws_internal_zalloc((size_t)nv + 1, sizeof(unsigned int));
		unsigned int n = 0, j;
		unsigned char* border = (unsigned char*)qaws_internal_zalloc((size_t)nv + 1, 1);
		if (!he || !deg || !border)
		{
			qaws_internal_free(he); qaws_internal_free(deg); qaws_internal_free(border);
			goto done;
		}
		/* border edges: used once */
		for (i = 0; i < 3 * nt; i++)
		{
			unsigned long long p = t[i], q = t[3 * (i / 3) + (i % 3 + 1) % 3];
			he[n++] = p < q ? (p << 32 | q) : (q << 32 | p);
		}
		qsort(he, n, sizeof *he, sm_cmp_key);
		for (i = 0, j = 0; i < n; i = j)
		{
			for (j = i + 1; j < n && he[j] == he[i]; j++)
				;
			if (j == i + 1)
			{
				border[(unsigned int)(he[i] >> 32)] = 1;
				border[(unsigned int)he[i]] = 1;
			}
		}
		{
			unsigned int nb = 0;
			for (i = 0, j = 0; i < n; i = j)
			{
				for (j = i + 1; j < n && he[j] == he[i]; j++)
					;
				if (j == i + 1 && (D.flags & QAWS_SMOOTH_MOVE_BORDER))
					he[nb++] = he[i];
			}
			n = nb;
		}
		for (i = 0; i < nfe; i++)
		{
			unsigned long long p = fe[2 * (size_t)i], q = fe[2 * (size_t)i + 1];
			he[n++] = p < q ? (p << 32 | q) : (q << 32 | p);
		}
		qsort(he, n, sizeof *he, sm_cmp_key);
		for (i = 0; i < 2 * nv; i++)
			ln[i] = SM_NONE;
		for (i = 0; i < n; i++)
		{
			unsigned int a = (unsigned int)(he[i] >> 32), b = (unsigned int)he[i];
			if (i && he[i] == he[i - 1])
				continue;
			if (deg[a] < 2)
				ln[2 * (size_t)a + deg[a]] = b;
			if (deg[b] < 2)
				ln[2 * (size_t)b + deg[b]] = a;
			deg[a]++;
			deg[b]++;
		}
		for (i = 0; i < nv; i++)
		{
			int held = (D.held && D.held[i]) || (border[i] && !(D.flags & QAWS_SMOOTH_MOVE_BORDER)) || qaws_mesh_3d_is_feature_corner(mesh, i)
				|| (deg[i] && deg[i] != 2);
			cls[i] = held ? 2 : deg[i] == 2 ? 1 : 0;
			nothing_held &= cls[i] == 0;
		}
		qaws_internal_free(he);
		qaws_internal_free(deg);
		qaws_internal_free(border);
	}
	for (i = 0; i < 3 * nt; i++)
	{
		double const *a = x + 3 * (size_t)t[i], *b = x + 3 * (size_t)t[3 * (i / 3) + (i % 3 + 1) % 3];
		h += sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]));
	}
	h = nt ? h / (3.0 * nt) : 1;
	if (D.method == QAWS_SMOOTH_BILATERAL)
	{
		st = sm_bilateral(x, nv, t, nt, cls, ln, iters, step) ? sm_write(mesh, x, nv) : QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	if (nothing_held && (D.flags & QAWS_SMOOTH_KEEP_VOLUME))
		vol0 = sm_volume(x, t, nt);
	for (it = 0; it < iters; it++)
	{
		unsigned int pass, passes = D.method == QAWS_SMOOTH_TAUBIN ? 2 : 1;
		for (pass = 0; pass < passes; pass++)
		{
			double s = pass ? (D.mu != 0 ? (double)D.mu : -0.53) : step;
			if (D.method <= QAWS_SMOOTH_TAUBIN)
			{
				/* explicit: each moving vertex a fraction s of the way to its neighbours' weighted mean (the lines' along the line) */
				if (!sm_laplacian(&K, mass, x, nv, t, nt, D.weights == QAWS_SMOOTH_WEIGHTS_UNIFORM, 1))
					goto done;
				memcpy(y, x, sizeof(double) * 3 * nv);
				for (i = 0; i < nv; i++)
				{
					double m[3] = { 0, 0, 0 }, ws = 0;
					unsigned int p;
					if (cls[i] == 2)
						continue;
					if (cls[i] == 1)
					{
						for (k = 0; k < 3; k++)
							m[k] = 0.5 * (x[3 * (size_t)ln[2 * (size_t)i] + k] + x[3 * (size_t)ln[2 * (size_t)i + 1] + k]);
						ws = 1;
					}
					else
						for (p = K.rp[i]; p < K.rp[i + 1]; p++)
							if (K.ci[p] != i)
							{
								for (k = 0; k < 3; k++)
									m[k] -= K.v[p] * x[3 * (size_t)K.ci[p] + k];
								ws -= K.v[p];
							}
					if (!(ws > 0))
						continue;
					for (k = 0; k < 3; k++)
						y[3 * (size_t)i + k] = x[3 * (size_t)i + k] + s * (m[k] / ws - x[3 * (size_t)i + k]);
				}
				memcpy(x, y, sizeof(double) * 3 * nv);
				sm_csr_free(&K);
				continue;
			}
			/* implicit: the lines first ((I + s K1) x' = x along them), then the surface with the lines held */
			if (!sm_line_laplacian(&A, ln, cls, nv))
				goto done;
			for (i = 0; i < nv; i++)
			{
				unsigned int p;
				move[i] = cls[i] == 1;
				for (p = A.rp[i]; p < A.rp[i + 1]; p++)
					A.v[p] = (A.ci[p] == i ? 1 : 0) + s * A.v[p];
			}
			memcpy(rhs, x, sizeof(double) * 3 * nv);
			st = sm_solve(&A, x, rhs, move);
			sm_csr_free(&A);
			if (st != QAWS_STATUS_OK)
				goto done;
			st = QAWS_STATUS_ALLOCATION_FAILURE;
			if (!sm_laplacian(&K, mass, x, nv, t, nt, D.weights == QAWS_SMOOTH_WEIGHTS_UNIFORM, 0))
				goto done;
			if (D.method == QAWS_SMOOTH_IMPLICIT)
			{
				double tt = s * (D.weights == QAWS_SMOOTH_WEIGHTS_UNIFORM ? 1 : h * h);
				for (i = 0; i < K.rp[nv]; i++)
					K.v[i] *= tt;
			}
			else
			{
				double tt = s * (D.weights == QAWS_SMOOTH_WEIGHTS_UNIFORM ? 1 : h * h * h * h);
				double* inv = (double*)qaws_internal_malloc(sizeof(double) * ((size_t)nv + 1));
				if (!inv)
					goto done;
				for (i = 0; i < nv; i++)
					inv[i] = 1 / mass[i];
				k = (unsigned int)sm_mul(&B, &K, inv, &K);
				qaws_internal_free(inv);
				sm_csr_free(&K);
				if (!k)
					goto done;
				K = B;
				memset(&B, 0, sizeof B);
				for (i = 0; i < K.rp[nv]; i++)
					K.v[i] *= tt;
			}
			/* + M on the diagonal; the right-hand side M x */
			for (i = 0; i < nv; i++)
			{
				unsigned int p;
				for (p = K.rp[i]; p < K.rp[i + 1]; p++)
					if (K.ci[p] == i)
						K.v[p] += mass[i];
				move[i] = cls[i] == 0;
				for (k = 0; k < 3; k++)
					rhs[3 * (size_t)i + k] = mass[i] * x[3 * (size_t)i + k];
			}
			st = sm_solve(&K, x, rhs, move);
			sm_csr_free(&K);
			if (st != QAWS_STATUS_OK)
				goto done;
		}
		if (vol0 != 0)
		{
			/* back to the volume, about the vertices' mean */
			double c[3] = { 0, 0, 0 }, v = sm_volume(x, t, nt), sc;
			for (i = 0; i < nv; i++)
				for (k = 0; k < 3; k++)
					c[k] += x[3 * (size_t)i + k] / nv;
			sc = v != 0 && vol0 / v > 0 ? cbrt(vol0 / v) : 1;
			for (i = 0; i < nv; i++)
				for (k = 0; k < 3; k++)
					x[3 * (size_t)i + k] = c[k] + sc * (x[3 * (size_t)i + k] - c[k]);
		}
	}
	st = sm_write(mesh, x, nv);
done:
	sm_csr_free(&K);
	sm_csr_free(&A);
	sm_csr_free(&B);
	qaws_internal_free(x); qaws_internal_free(y); qaws_internal_free(t); qaws_internal_free(mass); qaws_internal_free(rhs);
	qaws_internal_free(cls); qaws_internal_free(move); qaws_internal_free(ln); qaws_internal_free(fe);
	return st;
}

static unsigned int sm_find(unsigned int* par, unsigned int i)
{
	while (par[i] != i)
	{
		par[i] = par[par[i]];
		i = par[i];
	}
	return i;
}

qaws_status qaws_mesh_3d_fair(qaws_mesh_3d* mesh, unsigned char const* free_vertices, qaws_fair_desc const* desc)
{
	unsigned int order = desc && desc->order ? desc->order : 2, nv, nt, i, k, pass;
	int uniform = desc && desc->weights == QAWS_SMOOTH_WEIGHTS_UNIFORM;
	double *x = NULL, *mass = NULL, *rhs = NULL, *inv = NULL;
	unsigned int *t = NULL, *par = NULL;
	unsigned char *move = NULL, *anchored = NULL;
	sm_csr K, P, Q;
	qaws_status st;
	memset(&K, 0, sizeof K);
	memset(&P, 0, sizeof P);
	memset(&Q, 0, sizeof Q);
	if (!mesh || !free_vertices || order > 3 || (desc && (unsigned int)desc->weights > (unsigned int)QAWS_SMOOTH_WEIGHTS_UNIFORM))
		return QAWS_STATUS_INVALID_ARGUMENT;
	st = sm_read(mesh, &x, &t, &nv, &nt);
	if (st != QAWS_STATUS_OK)
		goto done;
	st = QAWS_STATUS_ALLOCATION_FAILURE;
	mass = (double*)qaws_internal_malloc(sizeof(double) * ((size_t)nv + 1));
	inv = (double*)qaws_internal_malloc(sizeof(double) * ((size_t)nv + 1));
	rhs = (double*)qaws_internal_zalloc(qaws_internal_mul2(3, (size_t)nv + 1), sizeof(double));
	par = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)nv + 1));
	move = (unsigned char*)qaws_internal_zalloc((size_t)nv + 1, 1);
	anchored = (unsigned char*)qaws_internal_zalloc((size_t)nv + 1, 1);
	if (!mass || !inv || !rhs || !par || !move || !anchored)
		goto done;
	/* the free parts (joined by free-free edges) touching a held vertex move */
	for (i = 0; i < nv; i++)
		par[i] = i;
	for (i = 0; i < 3 * nt; i++)
	{
		unsigned int a = t[i], b = t[3 * (i / 3) + (i % 3 + 1) % 3];
		if (free_vertices[a] && free_vertices[b])
		{
			a = sm_find(par, a);
			b = sm_find(par, b);
			if (a != b)
				par[a < b ? b : a] = a < b ? a : b;
		}
	}
	for (i = 0; i < 3 * nt; i++)
	{
		unsigned int a = t[i], b = t[3 * (i / 3) + (i % 3 + 1) % 3];
		if (free_vertices[a] && !free_vertices[b])
			anchored[sm_find(par, a)] = 1;
		if (free_vertices[b] && !free_vertices[a])
			anchored[sm_find(par, b)] = 1;
	}
	for (i = 0; i < nv; i++)
		move[i] = free_vertices[i] && anchored[sm_find(par, i)];
	/* cotangent weights read the free vertices' own positions, whatever state they are in: solved first with uniform weights,
	   then twice more with the cotangents of the last solution */
	for (pass = uniform ? 2 : 0; pass < 3; pass++)
	{
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		if (!sm_laplacian(&K, mass, x, nv, t, nt, pass == 0 || uniform, 0))
			goto done;
		for (i = 0; i < nv; i++)
			inv[i] = 1 / mass[i];
		/* L^k: K (M^-1 K)^(k - 1), symmetric */
		for (k = 1; k < order; k++)
		{
			if (!sm_mul(&Q, k == 1 ? &K : &P, inv, &K))
				goto done;
			sm_csr_free(&P);
			P = Q;
			memset(&Q, 0, sizeof Q);
		}
		st = sm_solve(order == 1 ? &K : &P, x, rhs, move);
		sm_csr_free(&K);
		sm_csr_free(&P);
		if (st != QAWS_STATUS_OK)
			goto done;
	}
	st = sm_write(mesh, x, nv);
done:
	sm_csr_free(&K);
	sm_csr_free(&P);
	sm_csr_free(&Q);
	qaws_internal_free(x); qaws_internal_free(t); qaws_internal_free(mass); qaws_internal_free(rhs); qaws_internal_free(inv);
	qaws_internal_free(par); qaws_internal_free(move); qaws_internal_free(anchored);
	return st;
}

/* sm_label_cut's flow network: the sites, then the source and the sink; each arc and its reverse side by side */
typedef struct sm_flow
{
	unsigned int *head, *to, *next, *level, *it, *queue, *stack;
	double* cap;
	unsigned int na, nn;
} sm_flow;

static void sm_flow_arc(sm_flow* F, unsigned int a, unsigned int b, double c, double back)
{
	F->to[F->na] = b; F->cap[F->na] = c; F->next[F->na] = F->head[a]; F->head[a] = F->na++;
	F->to[F->na] = a; F->cap[F->na] = back; F->next[F->na] = F->head[b]; F->head[b] = F->na++;
}

/* the most flow from s to t (Dinic: levels by breadth, blocking flows by an explicit stack) */
static void sm_flow_max(sm_flow* F, unsigned int s, unsigned int t)
{
	for (;;)
	{
		unsigned int qh = 0, qt = 0, i;
		for (i = 0; i < F->nn; i++) F->level[i] = 0xFFFFFFFFu;
		F->level[s] = 0;
		F->queue[qt++] = s;
		while (qh < qt)
		{
			unsigned int v = F->queue[qh++], e;
			for (e = F->head[v]; e != 0xFFFFFFFFu; e = F->next[e])
				if (F->cap[e] > 1e-12 && F->level[F->to[e]] == 0xFFFFFFFFu)
				{
					F->level[F->to[e]] = F->level[v] + 1;
					F->queue[qt++] = F->to[e];
				}
		}
		if (F->level[t] == 0xFFFFFFFFu)
			return;
		for (i = 0; i < F->nn; i++) F->it[i] = F->head[i];
		for (;;)
		{
			/* a path s to t along rising levels: its arcs on the stack */
			unsigned int sp = 0, v = s;
			while (v != t)
			{
				unsigned int e = F->it[v];
				while (e != 0xFFFFFFFFu && !(F->cap[e] > 1e-12 && F->level[F->to[e]] == F->level[v] + 1))
					e = F->it[v] = F->next[e];
				if (e == 0xFFFFFFFFu)
				{
					if (!sp)
						break;
					F->level[v] = 0xFFFFFFFFu;   /* (a dead end) */
					v = F->to[F->stack[--sp] ^ 1];
					continue;
				}
				F->stack[sp++] = e;
				v = F->to[e];
			}
			if (v != t)
				break;
			{
				double f = HUGE_VAL;
				for (i = 0; i < sp; i++) f = fmin(f, F->cap[F->stack[i]]);
				for (i = 0; i < sp; i++)
				{
					F->cap[F->stack[i]] -= f;
					F->cap[F->stack[i] ^ 1] += f;
				}
			}
		}
	}
}

static double sm_label_energy(unsigned char const* lab, double const* data, unsigned int n, unsigned int L, unsigned int const* nbr,
	double const* wgt, double const* V)
{
	double e = 0;
	unsigned int i, k;
	for (i = 0; i < n; i++)
	{
		e += data[(size_t)L * i + lab[i]];
		for (k = 0; k < 3; k++)
			if (nbr[3 * (size_t)i + k] != 0xFFFFFFFFu && nbr[3 * (size_t)i + k] > i)
				e += wgt[3 * (size_t)i + k] * V[(size_t)L * lab[i] + lab[nbr[3 * (size_t)i + k]]];
	}
	return e;
}

int sm_label_cut(unsigned char* lab, double const* data, unsigned int n, unsigned int L, unsigned int const* nbr, double const* wgt,
	double const* V, unsigned int rounds)
{
	sm_flow F;
	size_t arcs = 2 * (5 * (size_t)n + 1);
	double *c0 = (double*)qaws_internal_malloc(sizeof(double) * ((size_t)n + 1)), *c1 = (double*)qaws_internal_malloc(sizeof(double) * ((size_t)n + 1));
	unsigned char* trial = (unsigned char*)qaws_internal_malloc((size_t)n + 1);
	unsigned int r, a, i, k, ok = 0;
	double best = sm_label_energy(lab, data, n, L, nbr, wgt, V);
	memset(&F, 0, sizeof F);
	F.nn = n + 2;
	F.head = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * F.nn);
	F.level = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * F.nn);
	F.it = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * F.nn);
	F.queue = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * F.nn);
	F.stack = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * F.nn);
	F.to = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * arcs);
	F.next = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * arcs);
	F.cap = (double*)qaws_internal_malloc(sizeof(double) * arcs);
	if (!c0 || !c1 || !trial || !F.head || !F.level || !F.it || !F.queue || !F.stack || !F.to || !F.next || !F.cap)
		goto done;
	for (r = 0; r < rounds; r++)
	{
		int better = 0;
		for (a = 0; a < L; a++)
		{
			double e;
			/* x_i = 1: site i takes a. The source side keeps its label (cut s -> i: the cost of 1), the sink side takes a (cut
			   i -> t: the cost of 0); a pair E(x_i, x_j) = A + (C - A) x_i - C x_j + (B + C - A) (1 - x_i) x_j */
			for (i = 0; i < n; i++)
			{
				c0[i] = data[(size_t)L * i + lab[i]];
				c1[i] = data[(size_t)L * i + a];
				F.head[i] = 0xFFFFFFFFu;
			}
			F.head[n] = F.head[n + 1] = 0xFFFFFFFFu;
			F.na = 0;
			for (i = 0; i < n; i++)
				for (k = 0; k < 3; k++)
				{
					unsigned int j = nbr[3 * (size_t)i + k];
					double w = wgt[3 * (size_t)i + k], A, B, C;
					if (j == 0xFFFFFFFFu || j < i)
						continue;
					A = w * V[(size_t)L * lab[i] + lab[j]];
					B = w * V[(size_t)L * lab[i] + a];
					C = w * V[(size_t)L * a + lab[j]];
					if (C > A)
						c1[i] += C - A;
					else
						c0[i] += A - C;
					c0[j] += C;
					if (B + C - A > 0)
						sm_flow_arc(&F, i, j, B + C - A, 0);
				}
			for (i = 0; i < n; i++)
			{
				double m = fmin(c0[i], c1[i]);
				if (lab[i] == a)
					continue;   /* (it has a already: either side the same) */
				if (c1[i] - m > 0)
					sm_flow_arc(&F, n, i, c1[i] - m, 0);
				if (c0[i] - m > 0)
					sm_flow_arc(&F, i, n + 1, c0[i] - m, 0);
			}
			sm_flow_max(&F, n, n + 1);
			/* the source's side: reachable from it in the residual network; the rest takes a */
			{
				unsigned int qh = 0, qt = 0;
				for (i = 0; i < F.nn; i++) F.level[i] = 0;
				F.level[n] = 1;
				F.queue[qt++] = n;
				while (qh < qt)
				{
					unsigned int v = F.queue[qh++], g;
					for (g = F.head[v]; g != 0xFFFFFFFFu; g = F.next[g])
						if (F.cap[g] > 1e-12 && !F.level[F.to[g]])
						{
							F.level[F.to[g]] = 1;
							F.queue[qt++] = F.to[g];
						}
				}
			}
			for (i = 0; i < n; i++)
				trial[i] = F.level[i] ? lab[i] : (unsigned char)a;
			e = sm_label_energy(trial, data, n, L, nbr, wgt, V);
			if (e < best - 1e-12 * fabs(best))
			{
				best = e;
				memcpy(lab, trial, n);
				better = 1;
			}
		}
		if (!better)
			break;
	}
	ok = 1;
done:
	qaws_internal_free(c0); qaws_internal_free(c1); qaws_internal_free(trial);
	qaws_internal_free(F.head); qaws_internal_free(F.level); qaws_internal_free(F.it); qaws_internal_free(F.queue); qaws_internal_free(F.stack);
	qaws_internal_free(F.to); qaws_internal_free(F.next); qaws_internal_free(F.cap);
	return (int)ok;
}
