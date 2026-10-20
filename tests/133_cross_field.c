/*
 * Test 133: smoothest cross fields (Knoppel et al. 2013) on a sparse LDL^T
 *
 *   - the factor: a grid Laplacian solved to rounding, against conjugate
 *     gradients; nested dissection's fill against the natural order and
 *     against 12 n log2(sqrt n); a complex Hermitian system in 2 x 2
 *     blocks against dense elimination; a semidefinite Laplacian shifted;
 *   - SMOOTHEST through the public API: indices summing to n chi exactly
 *     (icosphere for n = 2, 4, 6; torus; a disc held by its border; the
 *     octahedron, icosahedron and cubes, whose faces' excesses pass half
 *     a turn), a sphere beside a held grid in one mesh solved as alone, the
 *     can's cones inside its caps, a cube's crease vertices along their
 *     creases to rounding with its cones at the corners; drawn to
 *     curvature on a torus, along the parallels; determinism; refusals;
 *   - INSTANT keeping the best of the fields: the can, the blob and a
 *     torus all quads, their irregular vertices;
 *   - the time of a 50k-vertex field.
 *
 * INSTANT's three meshes (16 s) and the 50k-vertex field run with
 * QAWS_TEST_HEAVY=1 only; test 135 quadrangulates the same can and blob by
 * default, with tighter bounds on their irregular vertices.
 */

#include "test_common.h"
#include "qaws_quadrangulate.h"
#include "qaws_mesh_field.h"
#include "qaws_remesh.h"
#include "internal/qaws_internal_ldlt.h"
#include "internal/qaws_internal_sparse.h"
#include "mesh/qaws_cross_field.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

static unsigned int cft_seed = 12345u;

static double cft_rand(void)
{
	cft_seed = cft_seed * 1664525u + 1013904223u;
	return ((cft_seed >> 8) & 0xFFFFFF) / (double)0x1000000 - 0.5;
}

static qaws_mesh_3d* cft_make(double const* x, unsigned int nv, unsigned int const* f, unsigned int const* fs, unsigned int nf)
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

static double* cft_positions(qaws_mesh_3d const* m)
{
	unsigned int n = qaws_mesh_3d_get_vertex_count(m), i;
	qaws_scalar* p = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * (n + 1));
	double* x = (double*)malloc(sizeof(double) * 3 * (n + 1));
	qaws_mesh_3d_get_positions(p, m);
	for (i = 0; i < 3 * n; i++)
		x[i] = (double)p[i];
	free(p);
	return x;
}

/* ---- the factor -------------------------------------------------------- */

/* the 5-point Laplacian of a k x k grid plus s on the diagonal, full pattern, the points' coordinates */
static unsigned int cft_grid(unsigned int k, double s, unsigned int* rp, unsigned int* col, double* val, double* xyz)
{
	unsigned int i, j, at = 0;
	for (j = 0; j < k; j++)
		for (i = 0; i < k; i++)
		{
			unsigned int r = j * k + i, nb[5], q, nn = 0, d;
			if (j > 0) nb[nn++] = r - k;
			if (i > 0) nb[nn++] = r - 1;
			d = nn;
			nb[nn++] = r;
			if (i + 1 < k) nb[nn++] = r + 1;
			if (j + 1 < k) nb[nn++] = r + k;
			rp[r] = at;
			for (q = 0; q < nn; q++)
			{
				col[at] = nb[q];
				val[at++] = q == d ? s + (nn - 1) : -1;
			}
			xyz[3 * r] = i;
			xyz[3 * r + 1] = j;
			xyz[3 * r + 2] = 0;
		}
	rp[k * k] = at;
	return at;
}

static void cft_matvec(unsigned int n, unsigned int const* rp, unsigned int const* col, double const* val, double const* x, double* y)
{
	unsigned int i, p;
	for (i = 0; i < n; i++)
	{
		y[i] = 0;
		for (p = rp[i]; p < rp[i + 1]; p++)
			y[i] += val[p] * x[col[p]];
	}
}

static double cft_norm(double const* a, unsigned int n)
{
	double s = 0;
	unsigned int i;
	for (i = 0; i < n; i++)
		s += a[i] * a[i];
	return sqrt(s);
}

static void test_ldlt(void)
{
	unsigned int k = 100, n = k * k, i, j, nnz;
	unsigned int* rp = (unsigned int*)malloc(sizeof(unsigned int) * (n + 1));
	unsigned int* col = (unsigned int*)malloc(sizeof(unsigned int) * 5 * n);
	double *val = (double*)malloc(sizeof(double) * 5 * n), *xyz = (double*)malloc(sizeof(double) * 3 * n);
	double *b = (double*)malloc(sizeof(double) * n), *x = (double*)malloc(sizeof(double) * n), *y = (double*)malloc(sizeof(double) * n);
	double res, cg = -1, fill100, fill60nd = 0, fill60nat = 0, bound;
	qaws_internal_ldlt f;
	qaws_status st;
	char msg[400];
	nnz = cft_grid(k, 1e-3, rp, col, val, xyz);
	for (i = 0; i < n; i++)
		b[i] = cft_rand();
	st = qaws_internal_ldlt_factor(&f, rp, col, val, xyz, n, 0);
	if (st != QAWS_STATUS_OK)
	{
		TEST_ASSERT(0, "ldlt: the grid factored");
		goto done;
	}
	fill100 = f.nnz + n;
	qaws_internal_ldlt_solve(x, &f, b);
	cft_matvec(n, rp, col, val, x, y);
	for (i = 0; i < n; i++)
		y[i] -= b[i];
	res = cft_norm(y, n) / cft_norm(b, n);
	qaws_internal_ldlt_free(&f);
#if !QAWS_SCALAR_IS_FLOAT
	{
		qaws_internal_csr A;
		if (qaws_internal_csr_alloc(&A, n, n, nnz) == QAWS_STATUS_OK)
		{
			memcpy(A.row_ptr, rp, sizeof(unsigned int) * (n + 1));
			memcpy(A.col_idx, col, sizeof(unsigned int) * nnz);
			memcpy(A.values, val, sizeof(double) * nnz);
			qaws_internal_cg_solve(y, &A, b, 20000, 1e-12);
			for (i = 0; i < n; i++)
				y[i] -= x[i];
			cg = cft_norm(y, n) / cft_norm(x, n);
			qaws_internal_csr_free(&A);
		}
	}
#else
	cg = 0;
#endif
	/* fill: nested dissection against the natural order on 60 x 60 */
	cft_grid(60, 1e-3, rp, col, val, xyz);
	if (qaws_internal_ldlt_factor(&f, rp, col, val, xyz, 3600, 0) == QAWS_STATUS_OK)
		fill60nd = f.nnz + 3600;
	qaws_internal_ldlt_free(&f);
	if (qaws_internal_ldlt_factor(&f, rp, col, val, NULL, 3600, 0) == QAWS_STATUS_OK)
		fill60nat = f.nnz + 3600;
	qaws_internal_ldlt_free(&f);
	bound = 12.0 * n * log2(sqrt((double)n));
	sprintf(msg, "ldlt: 100 x 100 grid residual %.2g (1e-12), against CG %.2g (1e-8); fill 100 x 100 %.0f (bound %.0f), 60 x 60 nested %.0f, natural %.0f",
		res, cg, fill100, bound, fill60nd, fill60nat);
	printf("  %s\n", msg);
	TEST_ASSERT(res <= 1e-12 && cg >= 0 && cg <= 1e-8 && fill100 <= bound && fill60nd > 0 && fill60nd <= 0.5 * fill60nat, msg);
	/* a complex Hermitian system of 50 unknowns in 2 x 2 blocks (a connection Laplacian, shifted) against dense elimination */
	{
		unsigned int m = 50, N = 2 * m, e, at = 0;
		double *A = (double*)calloc((size_t)N * N, sizeof(double)), *G = (double*)malloc(sizeof(double) * N * (N + 1)), *xg = (double*)malloc(sizeof(double) * N);
		unsigned int* crp = (unsigned int*)malloc(sizeof(unsigned int) * (N + 1));
		unsigned int* ccol = (unsigned int*)malloc(sizeof(unsigned int) * N * N);
		double *cval = (double*)malloc(sizeof(double) * N * N), *cxyz = (double*)malloc(sizeof(double) * 3 * N), err = 0, cres;
		for (i = 0; i < m; i++)
		{
			A[(2 * i) * N + 2 * i] += 0.1;
			A[(2 * i + 1) * N + 2 * i + 1] += 0.1;
		}
		for (e = 0; e < 2 * m; e++)
		{
			unsigned int a = e % m, c = (a + (e < m ? 1 : 7)) % m;
			double w = 1 + cft_rand(), t = 6.2831853 * cft_rand(), zr = -w * cos(t), zi = -w * sin(t);
			/* w |u_c - r u_a|^2: A_aa, A_cc += w; A_ca = -w r, A_ac = -w conj(r) */
			for (j = 0; j < 2; j++)
			{
				A[(2 * a + j) * N + 2 * a + j] += w;
				A[(2 * c + j) * N + 2 * c + j] += w;
			}
			A[(2 * c) * N + 2 * a] += zr; A[(2 * c) * N + 2 * a + 1] -= zi; A[(2 * c + 1) * N + 2 * a] += zi; A[(2 * c + 1) * N + 2 * a + 1] += zr;
			A[(2 * a) * N + 2 * c] += zr; A[(2 * a) * N + 2 * c + 1] += zi; A[(2 * a + 1) * N + 2 * c] -= zi; A[(2 * a + 1) * N + 2 * c + 1] += zr;
		}
		for (i = 0; i < N; i++)
		{
			crp[i] = at;
			for (j = 0; j < N; j++)
				if (A[i * N + j] != 0 || (j ^ 1) == i)
				{
					ccol[at] = j;
					cval[at++] = A[i * N + j];
				}
			cxyz[3 * i] = cos(6.2831853 * (i / 2) / m);
			cxyz[3 * i + 1] = sin(6.2831853 * (i / 2) / m);
			cxyz[3 * i + 2] = 0;
			b[i] = cft_rand();
		}
		crp[N] = at;
		/* (100 rows: dissected) */
		st = qaws_internal_ldlt_factor(&f, crp, ccol, cval, cxyz, N, 0);
		if (st == QAWS_STATUS_OK)
		{
			qaws_internal_ldlt_solve(x, &f, b);
			qaws_internal_ldlt_free(&f);
		}
		/* dense Gaussian elimination with partial pivoting */
		for (i = 0; i < N; i++)
		{
			memcpy(G + i * (N + 1), A + i * N, sizeof(double) * N);
			G[i * (N + 1) + N] = b[i];
		}
		for (i = 0; i < N; i++)
		{
			unsigned int pv = i, q;
			for (j = i + 1; j < N; j++)
				if (fabs(G[j * (N + 1) + i]) > fabs(G[pv * (N + 1) + i]))
					pv = j;
			for (q = 0; q <= N; q++)
			{
				double t = G[i * (N + 1) + q];
				G[i * (N + 1) + q] = G[pv * (N + 1) + q];
				G[pv * (N + 1) + q] = t;
			}
			for (j = i + 1; j < N; j++)
			{
				double l = G[j * (N + 1) + i] / G[i * (N + 1) + i];
				for (q = i; q <= N; q++)
					G[j * (N + 1) + q] -= l * G[i * (N + 1) + q];
			}
		}
		for (i = N; i > 0; i--)
		{
			double s = G[(i - 1) * (N + 1) + N];
			for (j = i; j < N; j++)
				s -= G[(i - 1) * (N + 1) + j] * xg[j];
			xg[i - 1] = s / G[(i - 1) * (N + 1) + i - 1];
		}
		for (i = 0; i < N; i++)
			err = fmax(err, fabs(x[i] - xg[i]));
		cft_matvec(N, crp, ccol, cval, x, y);
		for (i = 0; i < N; i++)
			y[i] -= b[i];
		cres = cft_norm(y, N) / cft_norm(b, N);
		sprintf(msg, "ldlt: complex Hermitian 2 x 2 blocks, 50 unknowns: residual %.2g, against dense elimination %.2g (1e-10)", cres, err);
		printf("  %s\n", msg);
		TEST_ASSERT(st == QAWS_STATUS_OK && cres <= 1e-10 && err <= 1e-10, msg);
		free(A); free(G); free(xg); free(crp); free(ccol); free(cval); free(cxyz);
	}
	/* semidefinite: the plain grid Laplacian, shifted by 1e-8 */
	{
		int finite = 1;
		cft_grid(40, 0, rp, col, val, xyz);
		for (i = 0; i < 1600; i++)
			b[i] = i < 800 ? 1 : -1;
		st = qaws_internal_ldlt_factor(&f, rp, col, val, xyz, 1600, 1e-8);
		if (st == QAWS_STATUS_OK)
		{
			qaws_internal_ldlt_solve(x, &f, b);
			for (i = 0; i < 1600; i++)
				finite &= x[i] == x[i] && fabs(x[i]) < 1e30;
			qaws_internal_ldlt_free(&f);
		}
		sprintf(msg, "ldlt: a semidefinite Laplacian shifted by 1e-8: %s, finite %d", qaws_status_to_string(st), finite);
		TEST_ASSERT(st == QAWS_STATUS_OK && finite, msg);
	}
done:
	free(rp); free(col); free(val); free(xyz); free(b); free(x); free(y);
}

/* ---- surfaces ---------------------------------------------------------- */

/* the icosahedron's faces cut in f along each edge, on the unit sphere (10 f^2 + 2 vertices); points shared by faces
   are computed alike bit for bit (sums of the same products) and welded by a hash of their bits */
static qaws_mesh_3d* cft_sphere(unsigned int f)
{
	static unsigned int const f0[60] = { 0, 11, 5, 0, 5, 1, 0, 1, 7, 0, 7, 10, 0, 10, 11, 1, 5, 9, 5, 11, 4, 11, 10, 2, 10, 7, 6, 7, 1, 8,
		3, 9, 4, 3, 4, 2, 3, 2, 6, 3, 6, 8, 3, 8, 9, 4, 9, 5, 2, 4, 11, 6, 2, 10, 8, 6, 7, 9, 8, 1 };
	double g = (1 + sqrt(5.0)) / 2;
	double const v0[36] = { -1, g, 0, 1, g, 0, -1, -g, 0, 1, -g, 0, 0, -1, g, 0, 1, g, 0, -1, -g, 0, 1, -g, g, 0, -1, g, 0, 1, -g, 0, -1, -g, 0, 1 };
	unsigned int np = 20 * (f + 1) * (f + 2) / 2, cap = 1, nv = 0, nt = 0, t, i, j, k;
	unsigned int *hash, *loc = (unsigned int*)malloc(sizeof(unsigned int) * (f + 1) * (f + 1)), *tri = (unsigned int*)malloc(sizeof(unsigned int) * 60 * f * f);
	double* x = (double*)malloc(sizeof(double) * 3 * np);
	qaws_mesh_3d* m;
	while (cap < 2 * np)
		cap *= 2;
	hash = (unsigned int*)malloc(sizeof(unsigned int) * cap);
	memset(hash, 0xFF, sizeof(unsigned int) * cap);
	for (t = 0; t < 20; t++)
	{
		double const *A = v0 + 3 * f0[3 * t], *B = v0 + 3 * f0[3 * t + 1], *C = v0 + 3 * f0[3 * t + 2];
		for (i = 0; i <= f; i++)
			for (j = 0; i + j <= f; j++)
			{
				double p[3];
				unsigned long long hv = 1469598103934665603ull;
				unsigned int h;
				for (k = 0; k < 3; k++)
				{
					unsigned long long bits;
					p[k] = (A[k] * (f - i - j) + B[k] * i + C[k] * j) / f;
					if (p[k] == 0)
						p[k] = 0;   /* (no negative zero) */
					memcpy(&bits, &p[k], sizeof bits);
					hv = (hv ^ bits) * 1099511628211ull;
				}
				for (h = (unsigned int)(hv ^ (hv >> 29)) & (cap - 1); hash[h] != 0xFFFFFFFFu; h = (h + 1) & (cap - 1))
					if (memcmp(x + 3 * hash[h], p, sizeof p) == 0)
						break;
				if (hash[h] == 0xFFFFFFFFu)
				{
					memcpy(x + 3 * nv, p, sizeof p);
					hash[h] = nv++;
				}
				loc[i * (f + 1) + j] = hash[h];
			}
		for (i = 0; i < f; i++)
			for (j = 0; i + j < f; j++)
			{
				tri[3 * nt] = loc[i * (f + 1) + j]; tri[3 * nt + 1] = loc[(i + 1) * (f + 1) + j]; tri[3 * nt + 2] = loc[i * (f + 1) + j + 1];
				nt++;
				if (i + j + 1 < f)
				{
					tri[3 * nt] = loc[(i + 1) * (f + 1) + j]; tri[3 * nt + 1] = loc[(i + 1) * (f + 1) + j + 1]; tri[3 * nt + 2] = loc[i * (f + 1) + j + 1];
					nt++;
				}
			}
	}
	for (i = 0; i < nv; i++)
	{
		double r = sqrt(x[3 * i] * x[3 * i] + x[3 * i + 1] * x[3 * i + 1] + x[3 * i + 2] * x[3 * i + 2]);
		for (k = 0; k < 3; k++)
			x[3 * i + k] /= r;
	}
	m = cft_make(x, nv, tri, NULL, nt);
	free(hash); free(loc); free(tri); free(x);
	return m;
}

static qaws_mesh_3d* cft_torus(unsigned int nu, unsigned int nv, double R, double r)
{
	double* x = (double*)malloc(sizeof(double) * 3 * nu * nv);
	unsigned int* f = (unsigned int*)malloc(sizeof(unsigned int) * 6 * nu * nv), i, j, k = 0;
	qaws_mesh_3d* m;
	for (j = 0; j < nv; j++)
		for (i = 0; i < nu; i++)
		{
			double u = 2 * M_PI * i / nu, v = 2 * M_PI * (j + 0.5 * (i & 1)) / nv;
			x[3 * (j * nu + i)] = (R + r * cos(v)) * cos(u);
			x[3 * (j * nu + i) + 1] = (R + r * cos(v)) * sin(u);
			x[3 * (j * nu + i) + 2] = r * sin(v);
		}
	for (j = 0; j < nv; j++)
		for (i = 0; i < nu; i++)
		{
			unsigned int a = j * nu + i, b = j * nu + (i + 1) % nu, c = ((j + 1) % nv) * nu + (i + 1) % nu, d = ((j + 1) % nv) * nu + i;
			f[k++] = a; f[k++] = b; f[k++] = c;
			f[k++] = a; f[k++] = c; f[k++] = d;
		}
	m = cft_make(x, nu * nv, f, NULL, 2 * nu * nv);
	free(x);
	free(f);
	return m;
}

/* test 125's can: 8 sides of radius 1, height 1.6 */
static qaws_mesh_3d* cft_can(void)
{
	unsigned int n = 8, i, k = 0;
	double x[3 * 18];
	unsigned int fv[3 * 32];
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
	return cft_make(x, 2 * n + 2, fv, NULL, 4 * n);
}

/* the showcase's blob: an icosahedron split once (edges split, new edges flipped), on a bumped ellipsoid */
static qaws_mesh_3d* cft_blob(void)
{
	static unsigned int const f[60] = { 0, 11, 5, 0, 5, 1, 0, 1, 7, 0, 7, 10, 0, 10, 11, 1, 5, 9, 5, 11, 4, 11, 10, 2, 10, 7, 6, 7, 1, 8,
		3, 9, 4, 3, 4, 2, 3, 2, 6, 3, 6, 8, 3, 8, 9, 4, 9, 5, 2, 4, 11, 6, 2, 10, 8, 6, 7, 9, 8, 1 };
	double g = (1 + sqrt(5.0)) / 2, x[36];
	double const v[36] = { -1, g, 0, 1, g, 0, -1, -g, 0, 1, -g, 0, 0, -1, g, 0, 1, g, 0, -1, -g, 0, 1, -g, g, 0, -1, g, 0, 1, -g, 0, -1, -g, 0, 1 };
	qaws_mesh_3d* m;
	unsigned int i, n, ne, nv0, vv;
	unsigned int *e, *src;
	double* p;
	qaws_scalar* ps;
	for (i = 0; i < 36; i++)
		x[i] = v[i] / sqrt(1 + g * g);
	m = cft_make(x, 12, f, NULL, 20);
	ne = qaws_mesh_3d_get_edge_count(m);
	nv0 = qaws_mesh_3d_get_vertex_count(m);
	e = (unsigned int*)malloc(sizeof(unsigned int) * 2 * ne);
	src = (unsigned int*)malloc(sizeof(unsigned int) * 2 * ne);
	qaws_mesh_3d_get_edges(e, m);
	for (i = 0; i < ne; i++)
		if (qaws_mesh_3d_split_edge(&vv, m, NULL, NULL, e[2 * i], e[2 * i + 1]) == QAWS_STATUS_OK)
		{
			src[2 * (vv - nv0)] = e[2 * i];
			src[2 * (vv - nv0) + 1] = e[2 * i + 1];
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
	p = cft_positions(m);
	ps = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * n);
	for (i = 0; i < n; i++)
	{
		double* q = p + 3 * i, l = sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2]), s;
		q[0] /= l; q[1] /= l; q[2] /= l;
		s = 1 + 0.32 * sin(2.2 * q[0] + 0.4) * cos(1.7 * q[1]) + 0.2 * q[2] * q[2];
		q[0] *= 1.25 * s; q[1] *= s; q[2] *= 0.9 * s;
	}
	for (i = 0; i < 3 * n; i++)
		ps[i] = (qaws_scalar)p[i];
	qaws_mesh_3d_set_positions(m, ps);
	free(ps);
	free(p);
	return m;
}

static qaws_mesh_3d* cft_remesh(qaws_mesh_3d const* m, double L)
{
	qaws_remesh_desc d;
	qaws_mesh_3d* out = NULL;
	memset(&d, 0, sizeof d);
	d.target_length = (qaws_scalar)L;
	qaws_mesh_remesh(&out, NULL, m, &d);
	return out;
}

/* a flat disc of radius 1: a 32-gon fan remeshed */
static qaws_mesh_3d* cft_disc(double L)
{
	double x[3 * 33];
	unsigned int f[3 * 32], i;
	qaws_mesh_3d *m, *out;
	x[0] = x[1] = x[2] = 0;
	for (i = 0; i < 32; i++)
	{
		x[3 * (i + 1)] = cos(2 * M_PI * i / 32);
		x[3 * (i + 1) + 1] = sin(2 * M_PI * i / 32);
		x[3 * (i + 1) + 2] = 0;
		f[3 * i] = 0;
		f[3 * i + 1] = i + 1;
		f[3 * i + 2] = (i + 1) % 32 + 1;
	}
	m = cft_make(x, 33, f, NULL, 32);
	out = cft_remesh(m, L);
	qaws_mesh_3d_destroy(m);
	return out;
}

/* the SMOOTHEST field: its directions (may be NULL) and indices, the counts of +1 and of other non-zero faces */
static qaws_status cft_field(qaws_mesh_3d const* m, unsigned int n, int aligned, qaws_scalar* dir, int* idx, unsigned int* plus, unsigned int* other)
{
	qaws_direction_field_desc d;
	unsigned int f, nf = qaws_mesh_3d_get_face_count(m);
	qaws_status st;
	memset(&d, 0, sizeof d);
	d.method = QAWS_DIRECTION_FIELD_SMOOTHEST;
	d.symmetry = n;
	d.curvature_aligned = aligned;
	st = qaws_mesh_3d_direction_field(dir, idx, m, &d);
	*plus = *other = 0;
	for (f = 0; st == QAWS_STATUS_OK && f < nf; f++)
	{
		*plus += idx[f] == 1;
		*other += idx[f] != 0 && idx[f] != 1;
	}
	return st;
}

static void test_indices(void)
{
	qaws_mesh_3d *s = cft_sphere(8), *t = cft_torus(30, 12, 2, 0.75), *disc = cft_disc(0.08);
	unsigned int nfs = qaws_mesh_3d_get_face_count(s), nft = qaws_mesh_3d_get_face_count(t), nfd = disc ? qaws_mesh_3d_get_face_count(disc) : 0;
	unsigned int p4, o4, p6, o6, p2, o2, pt, ot, pd, od, f;
	int* idx = (int*)malloc(sizeof(int) * (nfs + nft + nfd + 1)), tsum = 0, tbig = 0;
	qaws_status s4, s6, s2, st, sd = QAWS_STATUS_INTERNAL_ERROR;
	char msg[400];
	s4 = cft_field(s, 4, 0, NULL, idx, &p4, &o4);
	s6 = cft_field(s, 6, 0, NULL, idx, &p6, &o6);
	s2 = cft_field(s, 2, 0, NULL, idx, &p2, &o2);
	st = cft_field(t, 4, 0, NULL, idx, &pt, &ot);
	for (f = 0; st == QAWS_STATUS_OK && f < nft; f++)
	{
		tsum += idx[f];
		tbig += idx[f] < -1 || idx[f] > 1;
	}
	if (disc)
		sd = cft_field(disc, 4, 0, NULL, idx, &pd, &od);
	sprintf(msg, "SMOOTHEST indices: icosphere of %u vertices n = 4: %u of +1, %u others (8, 0); n = 6: %u, %u (12, 0); n = 2: %u, %u (4, 0); torus %u of +1, %u others, sum %d (0); disc of %u vertices %u, %u (4, 0)",
		qaws_mesh_3d_get_vertex_count(s), p4, o4, p6, o6, p2, o2, pt, ot, tsum, disc ? qaws_mesh_3d_get_vertex_count(disc) : 0, sd == QAWS_STATUS_OK ? pd : 0,
		sd == QAWS_STATUS_OK ? od : 0);
	printf("  %s\n", msg);
	TEST_ASSERT(s4 == QAWS_STATUS_OK && p4 == 8 && o4 == 0, msg);
	TEST_ASSERT(s6 == QAWS_STATUS_OK && p6 == 12 && o6 == 0, msg);
	TEST_ASSERT(s2 == QAWS_STATUS_OK && p2 == 4 && o2 == 0, msg);
	/* the target was no singular face, but the smoothest field of this fat torus (R = 2, r = 0.75) has pairs (+1 outside,
	   -1 inside): measured 6 and 6 after the 20 steps, 4 and 4 converged (200 steps); asserted: sum 0 in simple
	   singularities; drawn to curvature it has none (below) */
	TEST_ASSERT(st == QAWS_STATUS_OK && tsum == 0 && tbig == 0 && pt <= 8, msg);
	TEST_ASSERT(sd == QAWS_STATUS_OK && pd == 4 && od == 0, msg);
	free(idx);
	qaws_mesh_3d_destroy(s);
	qaws_mesh_3d_destroy(t);
	qaws_mesh_3d_destroy(disc);
}

/* the SMOOTHEST field's index sum over faces [f0, f1) (all when f1 is 0) */
static int cft_sum(qaws_mesh_3d const* m, unsigned int n, unsigned int f0, unsigned int f1, qaws_status* st)
{
	unsigned int nf = qaws_mesh_3d_get_face_count(m), plus, other, f;
	int* idx = (int*)malloc(sizeof(int) * (nf + 1)), sum = 0;
	*st = cft_field(m, n, 0, NULL, idx, &plus, &other);
	for (f = f0; *st == QAWS_STATUS_OK && f < (f1 ? f1 : nf); f++)
		sum += idx[f];
	free(idx);
	return sum;
}

/* coarse closed meshes, where n times a face's rescaled angle excess passes half a turn: the indices still sum to 2 n */
static void test_coarse(void)
{
	static double const ox[18] = { 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1 };
	static unsigned int const of[24] = { 0, 2, 4, 1, 4, 2, 0, 4, 3, 1, 3, 4, 0, 5, 2, 1, 2, 5, 0, 3, 5, 1, 5, 3 };
	static double const cx[24] = { -1, -1, -1, 1, -1, -1, 1, 1, -1, -1, 1, -1, -1, -1, 1, 1, -1, 1, 1, 1, 1, -1, 1, 1 };
	static unsigned int const cf[24] = { 0, 3, 2, 1, 4, 5, 6, 7, 0, 1, 5, 4, 1, 2, 6, 5, 2, 3, 7, 6, 3, 0, 4, 7 };
	static unsigned int const cs[6] = { 4, 4, 4, 4, 4, 4 };
	unsigned int ct[36], i;
	qaws_mesh_3d *oct = cft_make(ox, 6, of, NULL, 8), *ico = cft_sphere(1), *cube, *creased, *quads = cft_make(cx, 8, cf, cs, 6);
	qaws_status st[8];
	int sum[8], ok = 1;
	char msg[400];
	for (i = 0; i < 6; i++)
	{
		ct[6 * i] = cf[4 * i]; ct[6 * i + 1] = cf[4 * i + 1]; ct[6 * i + 2] = cf[4 * i + 2];
		ct[6 * i + 3] = cf[4 * i]; ct[6 * i + 4] = cf[4 * i + 2]; ct[6 * i + 5] = cf[4 * i + 3];
	}
	cube = cft_make(cx, 8, ct, NULL, 12);
	creased = cft_make(cx, 8, ct, NULL, 12);
	qaws_mesh_3d_detect_features(creased, (qaws_scalar)(M_PI / 6), 0, 0);
	sum[0] = cft_sum(oct, 4, 0, 0, &st[0]);
	sum[1] = cft_sum(oct, 6, 0, 0, &st[1]);
	sum[2] = cft_sum(ico, 4, 0, 0, &st[2]);
	sum[3] = cft_sum(ico, 6, 0, 0, &st[3]);
	sum[4] = cft_sum(cube, 4, 0, 0, &st[4]);
	sum[5] = cft_sum(cube, 6, 0, 0, &st[5]);
	sum[6] = cft_sum(creased, 4, 0, 0, &st[6]);
	sum[7] = cft_sum(quads, 4, 0, 0, &st[7]);
	for (i = 0; i < 8; i++)
		ok &= st[i] == QAWS_STATUS_OK && sum[i] == (i == 1 || i == 3 || i == 5 ? 12 : 8);
	sprintf(msg, "coarse index sums: octahedron n = 4: %d (8), n = 6: %d (12); icosahedron %d (8), %d (12); 12-triangle cube %d (8), %d (12), creased %d (8); 6-quad cube %d (8)",
		sum[0], sum[1], sum[2], sum[3], sum[4], sum[5], sum[6], sum[7]);
	printf("  %s\n", msg);
	TEST_ASSERT(ok, msg);
	qaws_mesh_3d_destroy(oct);
	qaws_mesh_3d_destroy(ico);
	qaws_mesh_3d_destroy(cube);
	qaws_mesh_3d_destroy(creased);
	qaws_mesh_3d_destroy(quads);
}

/* the singular faces among [f0, f1): +1 and others */
static void cft_count(int const* idx, unsigned int f0, unsigned int f1, unsigned int* plus, unsigned int* other)
{
	unsigned int f;
	*plus = *other = 0;
	for (f = f0; f < f1; f++)
	{
		*plus += idx[f] == 1;
		*other += idx[f] != 0 && idx[f] != 1;
	}
}

/* two components in one mesh, an icosphere without lines and a flat grid held by its border: the sphere's field its
   own eigenvector (smoothest; drawn to a guide only the grid has; left by the relaxation of the held grid) */
static void test_components(void)
{
	qaws_mesh_3d *s = cft_sphere(8), *m;
	unsigned int ns = qaws_mesh_3d_get_vertex_count(s), nfs = qaws_mesh_3d_get_face_count(s), k = 13, nv = ns + k * k, nf = nfs + 2 * (k - 1) * (k - 1);
	unsigned int i, j, at = 3 * nfs, p[4], o[4];
	double *x = (double*)malloc(sizeof(double) * 3 * nv), *xs = cft_positions(s), *g = (double*)calloc(3 * (size_t)nv, sizeof(double));
	double* q = (double*)malloc(sizeof(double) * 3 * nv);
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 3 * nf);
	int* idx = (int*)calloc(nf, sizeof(int));
	qaws_status st = QAWS_STATUS_INTERNAL_ERROR;
	mf_hier h;
	cf_system* sys = NULL;
	char msg[400];
	memcpy(x, xs, sizeof(double) * 3 * ns);
	qaws_mesh_3d_get_triangles(t, s);
	for (j = 0; j < k; j++)
		for (i = 0; i < k; i++)
		{
			double* y = x + 3 * (ns + j * k + i);
			y[0] = 3 + 0.1 * i; y[1] = 0.1 * j; y[2] = 0;
			g[3 * (ns + j * k + i)] = 1;
		}
	for (j = 0; j + 1 < k; j++)
		for (i = 0; i + 1 < k; i++)
		{
			unsigned int a = ns + j * k + i;
			t[at++] = a; t[at++] = a + 1; t[at++] = a + k + 1;
			t[at++] = a; t[at++] = a + k + 1; t[at++] = a + k;
		}
	m = cft_make(x, nv, t, NULL, nf);
	memset(p, 0, sizeof p);
	memset(o, 0xFF, sizeof o);
	if (cft_field(m, 4, 0, NULL, idx, &i, &j) == QAWS_STATUS_OK)
		cft_count(idx, 0, nfs, &p[0], &o[0]);
	if (mf_build(&h, &m->hm, NULL) == QAWS_STATUS_OK)
	{
		if (cf_create(&m->hm, h.lv[0].nrm, &h.lv[0], 4, &sys) == QAWS_STATUS_OK)
		{
			st = cf_aligned(sys, g, q);
			cf_index(sys, idx);
			cft_count(idx, 0, nfs, &p[1], &o[1]);
			if (st == QAWS_STATUS_OK)
				st = cf_smoothest(sys, q);
			if (st == QAWS_STATUS_OK)
				st = cf_relax(sys, 0.2, 100, q);
			cf_index(sys, idx);
			cft_count(idx, 0, nfs, &p[2], &o[2]);
			cf_free(sys);
		}
		mf_free(&h);
	}
	sprintf(msg, "two components, the sphere's faces: smoothest %u of +1, %u others; drawn to the grid's guide %u, %u; after relaxing the grid %u, %u (8, 0 each)",
		p[0], o[0], p[1], o[1], p[2], o[2]);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && p[0] == 8 && o[0] == 0 && p[1] == 8 && o[1] == 0 && p[2] == 8 && o[2] == 0, msg);
	free(x); free(xs); free(g); free(q); free(t); free(idx);
	qaws_mesh_3d_destroy(s);
	qaws_mesh_3d_destroy(m);
}

static void test_can(void)
{
	qaws_mesh_3d *can = cft_can(), *d = NULL;
	unsigned int nf, nv, f, plus = 0, other = 0, top = 0, bottom = 0, side = 0, *tri = NULL, *e = NULL, ne, same = 1;
	int *idx = NULL, *idx2 = NULL;
	qaws_scalar *dir = NULL, *dir2 = NULL;
	double *x = NULL, mean = 0, nearest = HUGE_VAL;
	qaws_status st;
	char msg[400];
	qaws_mesh_3d_detect_features(can, (qaws_scalar)(50 * M_PI / 180), 0, 0);
	d = cft_remesh(can, 0.06);
	if (!d)
	{
		TEST_ASSERT(0, "can: remeshed");
		goto done;
	}
	nf = qaws_mesh_3d_get_face_count(d);
	nv = qaws_mesh_3d_get_vertex_count(d);
	ne = qaws_mesh_3d_get_edge_count(d);
	idx = (int*)malloc(sizeof(int) * nf);
	idx2 = (int*)malloc(sizeof(int) * nf);
	dir = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	dir2 = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	tri = (unsigned int*)malloc(sizeof(unsigned int) * 3 * nf);
	e = (unsigned int*)malloc(sizeof(unsigned int) * 2 * ne);
	x = cft_positions(d);
	qaws_mesh_3d_get_triangles(tri, d);
	qaws_mesh_3d_get_edges(e, d);
	for (f = 0; f < ne; f++)
	{
		double const *a = x + 3 * e[2 * f], *b = x + 3 * e[2 * f + 1];
		mean += sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2])) / ne;
	}
	st = cft_field(d, 4, 0, dir, idx, &plus, &other);
	for (f = 0; st == QAWS_STATUS_OK && f < nf; f++)
		if (idx[f])
		{
			double c[3] = { 0, 0, 0 };
			unsigned int k, j;
			for (j = 0; j < 3; j++)
				for (k = 0; k < 3; k++)
					c[k] += x[3 * tri[3 * f + j] + k] / 3;
			if (fabs(c[2] - 0.8) < 1e-6)
				top++;
			else if (fabs(c[2] + 0.8) < 1e-6)
				bottom++;
			else
				side++;
			nearest = fmin(nearest, (1 - sqrt(c[0] * c[0] + c[1] * c[1])) / mean);
		}
	/* determinism */
	cft_field(d, 4, 0, dir2, idx2, &f, &f);
	same = memcmp(dir, dir2, sizeof(qaws_scalar) * 3 * nv) == 0 && memcmp(idx, idx2, sizeof(int) * nf) == 0;
	sprintf(msg, "SMOOTHEST on the can of %u vertices: %u of +1, %u others (8, 0); top %u, bottom %u, side %u (4, 4, 0); nearest the rim %.2f edges (2); two runs alike %d",
		nv, plus, other, top, bottom, side, nearest, same);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && plus == 8 && other == 0 && top == 4 && bottom == 4 && side == 0 && nearest >= 2, msg);
	TEST_ASSERT(same, "SMOOTHEST: two runs bitwise alike");
done:
	free(idx); free(idx2); free(dir); free(dir2); free(tri); free(e); free(x);
	qaws_mesh_3d_destroy(can);
	qaws_mesh_3d_destroy(d);
}

static void test_cube(void)
{
	static double const cx[24] = { -1, -1, -1, 1, -1, -1, 1, 1, -1, -1, 1, -1, -1, -1, 1, 1, -1, 1, 1, 1, 1, -1, 1, 1 };
	static unsigned int const cf[24] = { 0, 3, 2, 1, 4, 5, 6, 7, 0, 1, 5, 4, 1, 2, 6, 5, 2, 3, 7, 6, 3, 0, 4, 7 };
	static unsigned int const cs[6] = { 4, 4, 4, 4, 4, 4 };
	qaws_mesh_3d *cube = cft_make(cx, 8, cf, cs, 6), *d;
	unsigned int nv, nf, nfe, f, i, plus, other, sum = 0, away = 0, crease = 0, *tri, *fe;
	int* idx;
	qaws_scalar* dir;
	double *x, err = 0;
	qaws_status st;
	char msg[400];
	qaws_mesh_3d_detect_features(cube, (qaws_scalar)(M_PI / 6), 0, 0);
	d = cft_remesh(cube, 0.2);
	if (!d)
	{
		TEST_ASSERT(0, "cube: remeshed");
		qaws_mesh_3d_destroy(cube);
		return;
	}
	nv = qaws_mesh_3d_get_vertex_count(d);
	nf = qaws_mesh_3d_get_face_count(d);
	nfe = qaws_mesh_3d_get_feature_edge_count(d);
	idx = (int*)malloc(sizeof(int) * nf);
	dir = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	tri = (unsigned int*)malloc(sizeof(unsigned int) * 3 * nf);
	fe = (unsigned int*)malloc(sizeof(unsigned int) * 2 * (nfe + 1));
	x = cft_positions(d);
	qaws_mesh_3d_get_triangles(tri, d);
	qaws_mesh_3d_get_feature_edges(fe, d);
	st = cft_field(d, 4, 0, dir, idx, &plus, &other);
	for (f = 0; st == QAWS_STATUS_OK && f < nf; f++)
	{
		int corner = 0;
		sum += idx[f];
		for (i = 0; i < 3; i++)
		{
			double const* p = x + 3 * tri[3 * f + i];
			corner |= fabs(fabs(p[0]) - 1) < 1e-9 && fabs(fabs(p[1]) - 1) < 1e-9 && fabs(fabs(p[2]) - 1) < 1e-9;
		}
		away += idx[f] != 0 && !corner;
	}
	/* crease vertices off the corners: the direction or its quarter turn along the crease */
	for (i = 0; st == QAWS_STATUS_OK && i < nfe; i++)
	{
		unsigned int j;
		double e[3], l;
		for (j = 0; j < 3; j++)
			e[j] = x[3 * fe[2 * i + 1] + j] - x[3 * fe[2 * i] + j];
		l = sqrt(e[0] * e[0] + e[1] * e[1] + e[2] * e[2]);
		for (j = 0; j < 2; j++)
		{
			double const* p = x + 3 * fe[2 * i + j];
			qaws_scalar const* q = dir + 3 * fe[2 * i + j];
			double dt, cr[3];
			if (fabs(fabs(p[0]) - 1) < 1e-9 && fabs(fabs(p[1]) - 1) < 1e-9 && fabs(fabs(p[2]) - 1) < 1e-9)
				continue;
			dt = (q[0] * e[0] + q[1] * e[1] + q[2] * e[2]) / l;
			cr[0] = (q[1] * e[2] - q[2] * e[1]) / l;
			cr[1] = (q[2] * e[0] - q[0] * e[2]) / l;
			cr[2] = (q[0] * e[1] - q[1] * e[0]) / l;
			err = fmax(err, fmin(fabs(dt), sqrt(cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2])));
			crease++;
		}
	}
	sprintf(msg, "SMOOTHEST on a cube of %u vertices with its creases: index sum %u (8), %u of +1, %u others, %u singular faces off the corners (0); crease directions off by %.2g (1e-9) at %u ends",
		nv, sum, plus, other, away, err, crease);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && sum == 8 && away == 0 && err <= 1e-9 && crease > 0, msg);
	free(idx); free(dir); free(tri); free(fe); free(x);
	qaws_mesh_3d_destroy(cube);
	qaws_mesh_3d_destroy(d);
}

static void test_aligned(void)
{
	qaws_mesh_3d* t = cft_torus(48, 18, 2, 0.75);
	unsigned int nv = qaws_mesh_3d_get_vertex_count(t), nf = qaws_mesh_3d_get_face_count(t), i, good = 0, plus, other;
	int* idx = (int*)malloc(sizeof(int) * nf);
	qaws_scalar* dir = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	double* x = cft_positions(t);
	qaws_status st = cft_field(t, 4, 1, dir, idx, &plus, &other);
	char msg[300];
	for (i = 0; st == QAWS_STATUS_OK && i < nv; i++)
	{
		double const* p = x + 3 * i;
		double rho = sqrt(p[0] * p[0] + p[1] * p[1]), par[3], n[3], mer[3], a, b, l;
		par[0] = -p[1] / rho; par[1] = p[0] / rho; par[2] = 0;
		n[0] = p[0] - 2 * p[0] / rho; n[1] = p[1] - 2 * p[1] / rho; n[2] = p[2];
		l = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
		n[0] /= l; n[1] /= l; n[2] /= l;
		mer[0] = n[1] * par[2] - n[2] * par[1];
		mer[1] = n[2] * par[0] - n[0] * par[2];
		mer[2] = n[0] * par[1] - n[1] * par[0];
		a = fabs(dir[3 * i] * par[0] + dir[3 * i + 1] * par[1] + dir[3 * i + 2] * par[2]);
		b = fabs(dir[3 * i] * mer[0] + dir[3 * i + 1] * mer[1] + dir[3 * i + 2] * mer[2]);
		good += fmax(a, b) >= cos(5 * M_PI / 180);
	}
	sprintf(msg, "SMOOTHEST drawn to curvature on a torus of %u vertices: %u within 5 degrees of a parallel or a meridian (95%%), %u singular faces (0)",
		nv, good, plus + other);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && good >= 0.95 * nv && plus + other == 0, msg);
	free(idx); free(dir); free(x);
	qaws_mesh_3d_destroy(t);
}

static void test_refusals(void)
{
	qaws_mesh_3d* t = cft_torus(12, 6, 2, 0.75);
	qaws_direction_field_desc d;
	qaws_status s1, s2, s3;
	char msg[300];
	memset(&d, 0, sizeof d);
	d.symmetry = 6;
	s1 = qaws_mesh_3d_direction_field(NULL, NULL, t, &d);
	d.method = (qaws_direction_field_method)7;
	s2 = qaws_mesh_3d_direction_field(NULL, NULL, t, &d);
	d.method = QAWS_DIRECTION_FIELD_SMOOTHEST;
	d.symmetry = 3;
	s3 = qaws_mesh_3d_direction_field(NULL, NULL, t, &d);
	sprintf(msg, "refusals: LOCAL 6-symmetric %s, method 7 %s, SMOOTHEST 3-symmetric %s", qaws_status_to_string(s1), qaws_status_to_string(s2),
		qaws_status_to_string(s3));
	TEST_ASSERT(s1 == QAWS_STATUS_UNSUPPORTED_OPERATION && s2 == QAWS_STATUS_INVALID_ARGUMENT && s3 == QAWS_STATUS_UNSUPPORTED_OPERATION, msg);
	qaws_mesh_3d_destroy(t);
}

static unsigned int cft_quads(qaws_mesh_3d const* m, double L, qaws_status* st, double* mean, unsigned int* quads)
{
	qaws_quadrangulate_desc d;
	qaws_mesh_3d* out = NULL;
	qaws_mesh_report rep;
	memset(&d, 0, sizeof d);
	memset(&rep, 0, sizeof rep);
	d.target_length = (qaws_scalar)L;
	d.pure = 1;
	*st = qaws_mesh_quadrangulate(&out, &rep, m, &d);
	qaws_mesh_3d_destroy(out);
	*mean = (double)rep.quality.mean_edge;
	*quads = rep.quality.quads;
	return *st == QAWS_STATUS_OK ? rep.quality.irregular : 0xFFFFFFFFu;
}

static void test_instant(void)
{
	qaws_mesh_3d *can = cft_can(), *blob = cft_blob(), *torus = cft_torus(24, 10, 2, 0.75);
	qaws_status s1, s2, s3;
	unsigned int ic, ib, it, qc, qb, qt;
	double ec, eb, et;
	char msg[400];
	clock_t t0 = clock();
	qaws_mesh_3d_detect_features(can, (qaws_scalar)(50 * M_PI / 180), 0, 0);
	ic = cft_quads(can, 0.14, &s1, &ec, &qc);
	ib = cft_quads(blob, 0.11, &s2, &eb, &qb);
	it = cft_quads(torus, 0.3, &s3, &et, &qt);
	/* the targets: the can 18 or fewer (from 22), the blob no more than 12, the torus none; the mean edges near their
	   targets (the blob's quads grew with the global field, 4% longer) */
	sprintf(msg, "INSTANT, all quads: can %u irregular (18; was 22), %u quads, mean edge %.4f (0.14); blob %u (12), %u quads, mean edge %.4f (0.11); torus %u (0), %u quads, mean edge %.4f (0.3); %.2f s",
		ic, qc, ec, ib, qb, eb, it, qt, et, (double)(clock() - t0) / CLOCKS_PER_SEC);
	printf("  %s\n", msg);
	TEST_ASSERT(s1 == QAWS_STATUS_OK && s2 == QAWS_STATUS_OK && s3 == QAWS_STATUS_OK && ic <= 18 && ib <= 12 && it == 0, msg);
	TEST_ASSERT(fabs(ec - 0.14) < 0.15 * 0.14 && fabs(eb - 0.11) < 0.15 * 0.11 && fabs(et - 0.3) < 0.15 * 0.3, msg);
	qaws_mesh_3d_destroy(can);
	qaws_mesh_3d_destroy(blob);
	qaws_mesh_3d_destroy(torus);
}

static void test_timing(void)
{
	qaws_mesh_3d* s = cft_sphere(71);
	unsigned int nv = qaws_mesh_3d_get_vertex_count(s), nf = qaws_mesh_3d_get_face_count(s), plus, other;
	int* idx = (int*)malloc(sizeof(int) * nf);
	qaws_scalar* dir = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	clock_t t0 = clock();
	qaws_status st = cft_field(s, 4, 0, dir, idx, &plus, &other);
	double sec = (double)(clock() - t0) / CLOCKS_PER_SEC;
	char msg[300];
	/* target 1.5 s in Release; asserted below 10 s (debug builds, loaded machines) */
	sprintf(msg, "SMOOTHEST on an icosphere of %u vertices: %.2f s (1.5), %u of +1, %u others", nv, sec, plus, other);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && sec < 10 && plus == 8 && other == 0, msg);
	free(idx); free(dir);
	qaws_mesh_3d_destroy(s);
}

int test_133_cross_field_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 133: smoothest cross fields\n");
	test_ldlt();
	test_indices();
	test_coarse();
	test_components();
	test_can();
	test_cube();
	test_aligned();
	test_refusals();
	if (test_heavy())
	{
		test_instant();
		test_timing();
	}
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
