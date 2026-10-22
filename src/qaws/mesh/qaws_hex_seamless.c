/*
 * qaws: hexahedra of a frame field with singular edges (an integer-grid map with transitions, CubeCover: Nieser,
 * Reitebuch, Polthier 2011; its extraction as HexEx does it: Lyon, Bommes, Kobbelt 2016)
 */

#include "../internal/qaws_internal_memory.h"
#include "../internal/qaws_internal_ldlt.h"
#include "qaws_mesh_internal.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define HS_ROUNDS 600u      /* the most solves of the rounding */
#define HS_TOL 1e-5         /* a barycentric within this of 0: on that side (the solve's constraints hold to about 1e-7) */

/* the signed permutation P (P[i] = j + 1 or -(j + 1)) taking b's axes onto a's (axis i of a ~ sign axis j of b), a rotation
   nearest b^T a (axes: 3 each, a + 3 i axis i) */
static void hs_match(double const* a, double const* b, int* P)
{
	static int const perm[6][3] = { { 0, 1, 2 }, { 0, 2, 1 }, { 1, 0, 2 }, { 1, 2, 0 }, { 2, 0, 1 }, { 2, 1, 0 } };
	double D[3][3], best = -1;
	int p, i, j, bp = 0;
	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
			D[i][j] = a[3 * i] * b[3 * j] + a[3 * i + 1] * b[3 * j + 1] + a[3 * i + 2] * b[3 * j + 2];
	for (p = 0; p < 6; p++)
	{
		double s = fabs(D[0][perm[p][0]]) + fabs(D[1][perm[p][1]]) + fabs(D[2][perm[p][2]]);
		if (s > best)
		{
			best = s;
			bp = p;
		}
	}
	for (i = 0; i < 3; i++)
		P[i] = D[i][perm[bp][i]] < 0 ? -(perm[bp][i] + 1) : perm[bp][i] + 1;
	{
		int odd = bp == 1 || bp == 2 || bp == 5, neg = (P[0] < 0) + (P[1] < 0) + (P[2] < 0), w = 0;
		if ((odd + neg) & 1)
		{
			for (i = 1; i < 3; i++)
				if (fabs(D[i][perm[bp][i]]) < fabs(D[w][perm[bp][w]]))
					w = i;
			P[w] = -P[w];
		}
	}
}

/* M (3 x 3, row-major, entries -1 0 1) of a signed permutation: M[i][j] = sign when P[i] = +-(j + 1) */
static void hs_mat(int const* P, int* M)
{
	int i;
	memset(M, 0, sizeof(int) * 9);
	for (i = 0; i < 3; i++)
		M[3 * i + (P[i] < 0 ? -P[i] : P[i]) - 1] = P[i] < 0 ? -1 : 1;
}

static void hs_mul(int const* A, int const* B, int* C)
{
	int i, j, k, T[9];
	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
		{
			T[3 * i + j] = 0;
			for (k = 0; k < 3; k++)
				T[3 * i + j] += A[3 * i + k] * B[3 * k + j];
		}
	memcpy(C, T, sizeof T);
}

/* a tetrahedron's volume and its corners' barycentric gradients */
static double hs_grads(qaws_volume_mesh const* t, unsigned int c, double g[4][3])
{
	double const* p[4];
	double M[9], inv[9], det;
	unsigned int i, k;
	for (i = 0; i < 4; i++)
		p[i] = t->x + 3 * (size_t)t->cv[t->c0[c] + i];
	for (i = 0; i < 3; i++)
		for (k = 0; k < 3; k++)
			M[3 * k + i] = p[i + 1][k] - p[0][k];   /* (column i: corner i + 1 less corner 0) */
	inv[0] = M[4] * M[8] - M[5] * M[7]; inv[3] = M[5] * M[6] - M[3] * M[8]; inv[6] = M[3] * M[7] - M[4] * M[6];
	inv[1] = M[7] * M[2] - M[8] * M[1]; inv[4] = M[8] * M[0] - M[6] * M[2]; inv[7] = M[6] * M[1] - M[7] * M[0];
	inv[2] = M[1] * M[5] - M[2] * M[4]; inv[5] = M[2] * M[3] - M[0] * M[5]; inv[8] = M[0] * M[4] - M[1] * M[3];
	det = M[0] * inv[0] + M[1] * inv[3] + M[2] * inv[6];
	for (k = 0; k < 3; k++)
	{
		g[0][k] = 0;
		for (i = 0; i < 3; i++)
		{
			g[i + 1][k] = det != 0 ? inv[3 * i + k] / det : 0;   /* (row i of the inverse) */
			g[0][k] -= g[i + 1][k];
		}
	}
	return det / 6;
}

/* barycentric coordinates of q in a tetrahedron of corners P (3 each); its image's determinant (0: flat) */
static double hs_bary(double const* P, double const* q, double* l)
{
	double M[9], inv[9], det, d[3];
	unsigned int i, k;
	for (i = 0; i < 3; i++)
		for (k = 0; k < 3; k++)
			M[3 * k + i] = P[3 * (i + 1) + k] - P[k];
	inv[0] = M[4] * M[8] - M[5] * M[7]; inv[3] = M[5] * M[6] - M[3] * M[8]; inv[6] = M[3] * M[7] - M[4] * M[6];
	inv[1] = M[7] * M[2] - M[8] * M[1]; inv[4] = M[8] * M[0] - M[6] * M[2]; inv[7] = M[6] * M[1] - M[7] * M[0];
	inv[2] = M[1] * M[5] - M[2] * M[4]; inv[5] = M[2] * M[3] - M[0] * M[5]; inv[8] = M[0] * M[4] - M[1] * M[3];
	det = M[0] * inv[0] + M[1] * inv[3] + M[2] * inv[6];
	for (k = 0; k < 3; k++)
		d[k] = q[k] - P[k];
	l[0] = 1;
	for (i = 0; i < 3; i++)
	{
		l[i + 1] = det != 0 ? (inv[3 * i] * d[0] + inv[3 * i + 1] * d[1] + inv[3 * i + 2] * d[2]) / det : 0;
		l[0] -= l[i + 1];
	}
	return det;
}

static unsigned int hs_uf(unsigned int* uf, unsigned int i)
{
	while (uf[i] != i)
	{
		uf[i] = uf[uf[i]];
		i = uf[i];
	}
	return i;
}

static int hs_cmp_double(void const* a, void const* b)
{
	double x = *(double const*)a, y = *(double const*)b;
	return x < y ? -1 : x > y;
}

static int hs_cmp_u64(void const* a, void const* b)
{
	unsigned long long x = *(unsigned long long const*)a, y = *(unsigned long long const*)b;
	return x < y ? -1 : x > y;
}

/* index of key in the sorted keys, or HM_NONE */
static unsigned int hs_find(unsigned long long const* keys, unsigned int n, unsigned long long key)
{
	unsigned int lo = 0, hi = n;
	while (lo < hi)
	{
		unsigned int mid = (lo + hi) / 2;
		if (keys[mid] < key)
			lo = mid + 1;
		else
			hi = mid;
	}
	return lo < n && keys[lo] == key ? lo : HM_NONE;
}

/* a linear form of the unknowns: sum a_i U[idx_i] (at most 8 terms, like terms merged) */
typedef struct hs_form
{
	unsigned int n, kind;       /* kind (candidates): 0 a singular edge's coordinate, 1 a border plane, 2 a transition */
	unsigned int idx[8];
	double a[8];
} hs_form;

static void hs_term(hs_form* f, unsigned int idx, double a)
{
	unsigned int i;
	for (i = 0; i < f->n; i++)
		if (f->idx[i] == idx)
		{
			f->a[i] += a;
			return;
		}
	f->idx[f->n] = idx;
	f->a[f->n++] = a;
}

static double hs_eval(hs_form const* f, double const* U)
{
	double s = 0;
	unsigned int i;
	for (i = 0; i < f->n; i++)
		s += f->a[i] * U[f->idx[i]];
	return s;
}

typedef struct hs_trip
{
	unsigned long long key;     /* row << 32 | column */
	double v;
} hs_trip;

static int hs_cmp_trip(void const* a, void const* b)
{
	unsigned long long x = ((hs_trip const*)a)->key, y = ((hs_trip const*)b)->key;
	return x < y ? -1 : x > y;
}

typedef struct hs_forms
{
	hs_form* f;
	unsigned int n, cap;
} hs_forms;

static hs_form* hs_add(hs_forms* F, unsigned int kind)
{
	if (F->n == F->cap)
	{
		unsigned int nc = qaws_internal_grow(F->cap, (size_t)F->n + 1, 256);
		hs_form* g = nc ? (hs_form*)qaws_internal_realloc(F->f, qaws_internal_mul2(sizeof(hs_form), nc)) : NULL;
		if (!g)
			return NULL;
		F->f = g;
		F->cap = nc;
	}
	memset(&F->f[F->n], 0, sizeof(hs_form));
	F->f[F->n].kind = kind;
	return &F->f[F->n++];
}

/* the map's energy (each cell's coordinates' gradients against its axes over h) plus w |form - value|^2 over the
   constraints (value 0) and the fixed candidates; U solved; one copy pinned at 0 */
static qaws_status hs_solve(qaws_volume_mesh const* t, double const* R, double h, unsigned int const* copy, unsigned int ncopy, double const* cx,
	hs_forms const* con, hs_forms const* cand, double const* fixv, unsigned char const* isfix, double* U)
{
	unsigned int n = 3 * ncopy, nt = 0, cap, c, q, r, k, i, j, *rp = NULL, *col = NULL;
	hs_trip* T = NULL;
	double *val = NULL, *rhs = NULL, *xyz = NULL, *dg = NULL, *rhs2 = NULL, w, diag = 0;
	qaws_internal_ldlt f;
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	memset(&f, 0, sizeof f);
	cap = 48 * t->nc + 16 * (con->n + cand->n) + n + 3;
	T = (hs_trip*)qaws_internal_malloc(qaws_internal_mul2(sizeof(hs_trip), cap));
	rhs = (double*)qaws_internal_zalloc((size_t)n + 1, sizeof(double));
	dg = (double*)qaws_internal_zalloc((size_t)n + 1, sizeof(double));
	rhs2 = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), (size_t)n + 1));
	xyz = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 3, (size_t)n + 1));
	rp = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (size_t)n + 1));
	if (!T || !rhs || !xyz || !rp || !dg || !rhs2)
		goto done;
	for (c = 0; c < t->nc; c++)
	{
		double g[4][3], V = hs_grads(t, c, g);
		if (!(V > 0))
			continue;
		for (q = 0; q < 4; q++)
			for (k = 0; k < 3; k++)
			{
				double const* ax = R + 9 * (size_t)c + 3 * k;
				rhs[3 * (size_t)copy[4 * c + q] + k] += V * (g[q][0] * ax[0] + g[q][1] * ax[1] + g[q][2] * ax[2]) / h;
				for (r = 0; r < 4; r++)
				{
					double s = V * (g[q][0] * g[r][0] + g[q][1] * g[r][1] + g[q][2] * g[r][2]);
					T[nt].key = ((unsigned long long)(3 * copy[4 * c + q] + k) << 32) | (3 * copy[4 * c + r] + k);
					T[nt++].v = s;
					if (q == r)
						dg[3 * (size_t)copy[4 * c + q] + k] += s;
				}
			}
	}
	/* the weight from the median row (a sliver's huge gradients would drive a mean) */
	memcpy(rhs2, dg, sizeof(double) * n);
	qsort(rhs2, n, sizeof(double), hs_cmp_double);
	diag = n ? rhs2[n / 2] : 1;
	w = 1e7 * diag;
	for (j = 0; j < 2; j++)
	{
		hs_forms const* F = j ? cand : con;
		for (i = 0; i < F->n; i++)
		{
			hs_form const* fm = &F->f[i];
			double b = j ? fixv[i] : 0;
			if (j && !isfix[i])
				continue;   /* (a free candidate) */
			for (q = 0; q < fm->n; q++)
			{
				rhs[fm->idx[q]] += w * fm->a[q] * b;
				for (r = 0; r < fm->n; r++)
				{
					T[nt].key = ((unsigned long long)fm->idx[q] << 32) | fm->idx[r];
					T[nt++].v = w * fm->a[q] * fm->a[r];
				}
			}
		}
	}
	for (k = 0; k < 3; k++)
	{
		T[nt].key = ((unsigned long long)k << 32) | k;
		T[nt++].v = w;   /* (the first copy pinned at 0) */
	}
	/* every unknown pulled to 0 by a hundred-millionth of a typical weight: the parts the cut parts translate freely (their
	   transitions not yet whole) stay put and the matrix positive, factored once */
	for (i = 0; i < n; i++)
	{
		T[nt].key = ((unsigned long long)i << 32) | i;
		T[nt++].v = 1e-8 * diag;
	}
	qsort(T, nt, sizeof *T, hs_cmp_trip);
	for (i = 0, j = 0; i < nt; i++)
		if (j && T[j - 1].key == T[i].key)
			T[j - 1].v += T[i].v;
		else
			T[j++] = T[i];
	nt = j;
	col = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (size_t)nt + 1));
	val = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), (size_t)nt + 1));
	if (!col || !val)
		goto done;
	for (i = 0, r = 0; r <= n; r++)
	{
		while (i < nt && (T[i].key >> 32) < r)
			i++;
		rp[r] = i;
	}
	for (i = 0; i < nt; i++)
	{
		col[i] = (unsigned int)(T[i].key & 0xFFFFFFFFu);
		val[i] = T[i].v;
	}
	for (i = 0; i < n; i++)
		memcpy(xyz + 3 * (size_t)i, cx + 3 * (size_t)(i / 3), 3 * sizeof(double));
	st = qaws_internal_ldlt_factor(&f, rp, col, val, xyz, n, 0);
	if (st == QAWS_STATUS_OK)
		qaws_internal_ldlt_solve(U, &f, rhs);
done:
	qaws_internal_ldlt_free(&f);
	qaws_internal_free(T); qaws_internal_free(rhs); qaws_internal_free(dg); qaws_internal_free(rhs2); qaws_internal_free(xyz); qaws_internal_free(rp); qaws_internal_free(col); qaws_internal_free(val);
	return st;
}

/* a point of the lattice (doubled coordinates: corners even, cube centres odd) seen from a cell */
typedef struct hs_reg
{
	unsigned int c;
	int p[3];
	double l[4];
} hs_reg;

static unsigned long long hs_rkey(unsigned int c, int const* p)
{
	/* 28 bits of cell, 12 bits per coordinate (offset): lattices of up to 2048 steps across */
	return ((unsigned long long)c << 36) | ((unsigned long long)((p[0] + 2048) & 0xFFF) << 24) | ((unsigned long long)((p[1] + 2048) & 0xFFF) << 12) |
		(unsigned long long)((p[2] + 2048) & 0xFFF);
}

/* doubled coordinates p of c's chart into d's across slot s (from c's side): u_d = Pi^T (u_c - t) */
static void hs_across(int const* P, int const* t, int const* p, int scale, int* o)
{
	int i, x[3];
	for (i = 0; i < 3; i++)
		x[i] = p[i] - scale * t[i];
	for (i = 0; i < 3; i++)
	{
		int j = (P[i] < 0 ? -P[i] : P[i]) - 1;
		o[j] = P[i] < 0 ? -x[i] : x[i];
	}
}

static void hs_across_d(int const* P, int const* t, double const* p, double* o)
{
	int i;
	double x[3];
	for (i = 0; i < 3; i++)
		x[i] = p[i] - t[i];
	for (i = 0; i < 3; i++)
	{
		int j = (P[i] < 0 ? -P[i] : P[i]) - 1;
		o[j] = P[i] < 0 ? -x[i] : x[i];
	}
}

/* the lattice's registrations and what finding a cube's corners needs */
typedef struct hs_lat
{
	qaws_volume_mesh const* t;
	double const* par;              /* each cell's corners in the map (12) */
	unsigned char const* flat;
	int const *Ps, *Ts;
	unsigned long long const* rkey;
	unsigned int const *rid, *copy;
	unsigned int* uf;
	unsigned int nreg;
	double const* px;               /* the points' positions */
	unsigned int *bq, *vis, nbq, stamp;
	int* bk;
	double reach2;                  /* how far a corner may lie from where it is sought (squared) */
} hs_lat;

/* the point at doubled coordinates K2 of cell c's chart, sought from fr2 (doubled, in c): walked straight to just short of
   it across the cells and their transitions, else breadth first from c (the coordinates carried across each facet) among
   the points within reach of xyz; HM_NONE when none */
static unsigned int hs_corner(hs_lat* L, unsigned int c, int const* fr2, int const* K2in, double const* xyz)
{
	qaws_volume_mesh const* t = L->t;
	unsigned int cur = c, guard, at = HM_NONE, k, q;
	double from[3], to[3];
	int K2[3];
	for (k = 0; k < 3; k++)
	{
		K2[k] = K2in[k];
		from[k] = 0.5 * fr2[k];
		to[k] = 0.5 * K2[k] + 1e-3 * (from[k] - 0.5 * K2[k]);
	}
	for (guard = 0; guard < 4096; guard++)
	{
		double l0[4], l1[4], sbest = 2;
		unsigned int qb = HM_NONE;
		if (L->flat[cur])
			break;
		hs_bary(L->par + 12 * (size_t)cur, from, l0);
		hs_bary(L->par + 12 * (size_t)cur, to, l1);
		for (q = 0; q < 4; q++)
			if (l1[q] < -1e-12)
			{
				double sq = l0[q] / (l0[q] - l1[q]);
				if (sq < sbest)
				{
					sbest = sq;
					qb = q;
				}
			}
		if (qb == HM_NONE)
		{
			at = hs_find(L->rkey, L->nreg, hs_rkey(cur, K2));
			for (q = 0; q < 4 && at == HM_NONE; q++)
			{
				unsigned int s = t->f0[cur] + q, o = t->fo[s];
				int p2[3];
				if (o == HM_NONE)
					continue;
				hs_across(L->Ps + 3 * (size_t)s, L->Ts + 3 * (size_t)s, K2, 2, p2);
				at = hs_find(L->rkey, L->nreg, hs_rkey(o / 4, p2));
			}
			break;
		}
		{
			unsigned int s = t->f0[cur] + qb, o = t->fo[s];
			double m[3];
			int p2[3];
			if (o == HM_NONE)
				break;
			for (k = 0; k < 3; k++)
				from[k] = from[k] + fmax(0, sbest) * (to[k] - from[k]);
			hs_across_d(L->Ps + 3 * (size_t)s, L->Ts + 3 * (size_t)s, from, m);
			memcpy(from, m, sizeof m);
			hs_across_d(L->Ps + 3 * (size_t)s, L->Ts + 3 * (size_t)s, to, m);
			memcpy(to, m, sizeof m);
			hs_across(L->Ps + 3 * (size_t)s, L->Ts + 3 * (size_t)s, K2, 2, p2);
			memcpy(K2, p2, sizeof p2);
			cur = o / 4;
		}
	}
	if (at != HM_NONE && L->copy[hs_uf(L->uf, L->rid[at])] != HM_NONE)
		return L->copy[hs_uf(L->uf, L->rid[at])];
	{
		unsigned int head = 0, tail = 0, depth = 0, layer;
		L->stamp++;
		L->bq[tail] = c;
		memcpy(L->bk, K2in, 3 * sizeof(int));
		tail++;
		L->vis[c] = L->stamp;
		for (layer = tail; head < tail && depth < 10; depth++, layer = tail)
			for (; head < layer; head++)
			{
				unsigned int cc = L->bq[head], ff, f2 = hs_find(L->rkey, L->nreg, hs_rkey(cc, L->bk + 3 * (size_t)head));
				if (f2 != HM_NONE && L->copy[hs_uf(L->uf, L->rid[f2])] != HM_NONE)
				{
					unsigned int id = L->copy[hs_uf(L->uf, L->rid[f2])];
					double const* pp = L->px + 3 * (size_t)id;
					if ((pp[0] - xyz[0]) * (pp[0] - xyz[0]) + (pp[1] - xyz[1]) * (pp[1] - xyz[1]) + (pp[2] - xyz[2]) * (pp[2] - xyz[2]) < L->reach2)
						return id;
				}
				for (ff = 0; ff < 4 && tail + 1 < L->nbq; ff++)
				{
					unsigned int s = t->f0[cc] + ff, o = t->fo[s];
					if (o == HM_NONE || L->vis[o / 4] == L->stamp)
						continue;
					L->vis[o / 4] = L->stamp;
					L->bq[tail] = o / 4;
					hs_across(L->Ps + 3 * (size_t)s, L->Ts + 3 * (size_t)s, L->bk + 3 * (size_t)head, 2, L->bk + 3 * (size_t)tail);
					tail++;
				}
			}
	}
	return HM_NONE;
}

/* whether point x lies in a cell (3D) reached breadth first from c within 10 facets */
static int hs_inside(hs_lat* L, unsigned int c, double const* x)
{
	qaws_volume_mesh const* t = L->t;
	unsigned int head = 0, tail = 0, depth = 0, layer, q;
	L->stamp++;
	L->bq[tail++] = c;
	L->vis[c] = L->stamp;
	for (layer = tail; head < tail && depth < 10; depth++, layer = tail)
		for (; head < layer; head++)
		{
			unsigned int cc = L->bq[head];
			double P[12], l[4];
			for (q = 0; q < 4; q++)
				memcpy(P + 3 * q, t->x + 3 * (size_t)t->cv[t->c0[cc] + q], 3 * sizeof(double));
			if (hs_bary(P, x, l) > 0 && l[0] >= -1e-9 && l[1] >= -1e-9 && l[2] >= -1e-9 && l[3] >= -1e-9)
				return 1;
			for (q = 0; q < 4 && tail + 1 < L->nbq; q++)
			{
				unsigned int o = t->fo[t->f0[cc] + q];
				if (o == HM_NONE || L->vis[o / 4] == L->stamp)
					continue;
				L->vis[o / 4] = L->stamp;
				L->bq[tail++] = o / 4;
			}
		}
	return 0;
}

/* a table of quadrilateral faces (their sorted points) and how many cubes hold each */
typedef struct hs_faces
{
	unsigned int* key;          /* 4 per slot, HM_NONE first: empty */
	unsigned char* n;
	unsigned int cap, used;
} hs_faces;

static unsigned int hs_face_slot(hs_faces const* F, unsigned int const* s)
{
	unsigned int hsh = (s[0] * 2654435761u) ^ (s[1] * 40503u) ^ (s[2] * 2246822519u) ^ (s[3] * 3266489917u), i = hsh & (F->cap - 1);
	while (F->key[4 * (size_t)i] != HM_NONE && memcmp(F->key + 4 * (size_t)i, s, 4 * sizeof(unsigned int)) != 0)
		i = (i + 1) & (F->cap - 1);
	return i;
}

/* the cube's six faces' counts (out), or with add 1 their counts raised; 0 on allocation failure */
static int hs_cube_faces(hs_faces* F, unsigned int const* v, unsigned int* cnt, int add)
{
	static unsigned char const FQ[6][4] = { { 0, 3, 2, 1 }, { 4, 5, 6, 7 }, { 0, 1, 5, 4 }, { 1, 2, 6, 5 }, { 2, 3, 7, 6 }, { 3, 0, 4, 7 } };
	unsigned int f, a, b;
	if (add && 2 * (F->used + 6) >= F->cap)
	{
		hs_faces G;
		unsigned int i;
		G.cap = F->cap ? 2 * F->cap : 1024;
		while (2 * (F->used + 6) >= G.cap)
			G.cap *= 2;
		G.used = F->used;
		G.key = (unsigned int*)qaws_internal_malloc(qaws_internal_mul3(sizeof(unsigned int), 4, G.cap));
		G.n = (unsigned char*)qaws_internal_zalloc(G.cap, 1);
		if (!G.key || !G.n)
		{
			qaws_internal_free(G.key);
			qaws_internal_free(G.n);
			return 0;
		}
		for (i = 0; i < G.cap; i++)
			G.key[4 * (size_t)i] = HM_NONE;
		for (i = 0; i < F->cap; i++)
			if (F->key[4 * (size_t)i] != HM_NONE)
			{
				unsigned int j = hs_face_slot(&G, F->key + 4 * (size_t)i);
				memcpy(G.key + 4 * (size_t)j, F->key + 4 * (size_t)i, 4 * sizeof(unsigned int));
				G.n[j] = F->n[i];
			}
		qaws_internal_free(F->key);
		qaws_internal_free(F->n);
		*F = G;
	}
	for (f = 0; f < 6; f++)
	{
		unsigned int s[4], j;
		for (a = 0; a < 4; a++)
			s[a] = v[FQ[f][a]];
		for (a = 1; a < 4; a++)
			for (b = a; b > 0 && s[b - 1] > s[b]; b--)
			{
				unsigned int sw = s[b];
				s[b] = s[b - 1];
				s[b - 1] = sw;
			}
		if (!F->cap)
		{
			cnt[f] = 0;
			continue;
		}
		j = hs_face_slot(F, s);
		if (add)
		{
			if (F->key[4 * (size_t)j] == HM_NONE)
			{
				memcpy(F->key + 4 * (size_t)j, s, sizeof s);
				F->used++;
			}
			F->n[j]++;
		}
		cnt[f] = F->key[4 * (size_t)j] == HM_NONE ? 0 : F->n[j];
	}
	return 1;
}

static int hs_cmp_desc(void const* a, void const* b)
{
	double x = *(double const*)a, y = *(double const*)b;
	return x > y ? -1 : x < y;
}

/* a cube's faces against the oriented faces held (F: a face's sorted points; dir: per slot the first holder's points from
   its least, in its turn): 1 when none is held twice or held turned the same way, then its faces added */
static int hs_cube_fits(hs_faces* F, unsigned int* dir, unsigned int const* v)
{
	static unsigned char const FQ[6][4] = { { 0, 3, 2, 1 }, { 4, 5, 6, 7 }, { 0, 1, 5, 4 }, { 1, 2, 6, 5 }, { 2, 3, 7, 6 }, { 3, 0, 4, 7 } };
	unsigned int f, a, b, s[6][4], d[6][4], slot[6];
	for (f = 0; f < 6; f++)
	{
		unsigned int m = 0;
		for (a = 0; a < 4; a++)
		{
			s[f][a] = v[FQ[f][a]];
			if (v[FQ[f][a]] < v[FQ[f][m]])
				m = a;
		}
		for (a = 0; a < 4; a++)
			d[f][a] = v[FQ[f][(m + a) % 4]];
		for (a = 1; a < 4; a++)
			for (b = a; b > 0 && s[f][b - 1] > s[f][b]; b--)
			{
				unsigned int w = s[f][b];
				s[f][b] = s[f][b - 1];
				s[f][b - 1] = w;
			}
		slot[f] = hs_face_slot(F, s[f]);
		if (F->key[4 * (size_t)slot[f]] != HM_NONE && (F->n[slot[f]] >= 2 || memcmp(dir + 4 * (size_t)slot[f], d[f], sizeof d[f]) == 0))
			return 0;
	}
	for (f = 0; f < 6; f++)
	{
		unsigned int j = hs_face_slot(F, s[f]);
		if (F->key[4 * (size_t)j] == HM_NONE)
		{
			memcpy(F->key + 4 * (size_t)j, s[f], sizeof s[f]);
			memcpy(dir + 4 * (size_t)j, d[f], sizeof d[f]);
			F->used++;
		}
		F->n[j]++;
	}
	return 1;
}

static int hs_cmp_cube(void const* a, void const* b)
{
	unsigned int const *x = (unsigned int const*)a + 8, *y = (unsigned int const*)b + 8;
	int k;
	for (k = 0; k < 8; k++)
		if (x[k] != y[k])
			return x[k] < y[k] ? -1 : 1;
	return 0;
}

qaws_status hs_hexes(qaws_volume_mesh const* t, double* R, double h, qaws_volume_mesh** out_mesh, qaws_mesh_report* rep)
{
	unsigned int nturn = 0, nc = t->nc, nv = t->nv, ns = 4 * nc, c, f, i, j, k, q, ne = 0, ncopy = 0, nreg = 0, capreg = 0, round;
	unsigned int *queue = NULL, *uf = NULL, *copy = NULL, *ecell = NULL, *fe = NULL, *ecnt = NULL, *e0 = NULL, *ef = NULL, *stack = NULL;
	unsigned char *tree = NULL, *cut = NULL, *esing = NULL, *ebord = NULL, *seen = NULL, *flat = NULL;
	unsigned long long *ekey = NULL, *rkey = NULL;
	int *Ps = NULL, *Ts = NULL;
	double *U = NULL, *cx = NULL, *fixv = NULL, *par = NULL, *px = NULL;
	unsigned char* isfix = NULL;
	unsigned int *rid = NULL, *cells = NULL, ncell = 0, capc = 0, np = 0, *bq = NULL, *vis = NULL, nbq = 0, cappx = 0;
	int* bk = NULL;
	hs_reg* reg = NULL;
	hs_forms con, cand;
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	memset(&con, 0, sizeof con);
	memset(&cand, 0, sizeof cand);
	queue = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (size_t)ns + 1));
	tree = (unsigned char*)qaws_internal_zalloc((size_t)ns + 1, 1);
	cut = (unsigned char*)qaws_internal_zalloc((size_t)ns + 1, 1);
	seen = (unsigned char*)qaws_internal_zalloc((size_t)nc + 1, 1);
	Ps = (int*)qaws_internal_zalloc(3 * (size_t)ns + 1, sizeof(int));
	Ts = (int*)qaws_internal_zalloc(3 * (size_t)ns + 1, sizeof(int));
	ekey = (unsigned long long*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned long long), 6 * (size_t)nc + 1));
	if (!queue || !tree || !cut || !seen || !Ps || !Ts || !ekey)
		goto done;
	for (c = 0; c < nc; c++)
		if (t->type[c] != QAWS_CELL_TETRAHEDRON)
		{
			st = QAWS_STATUS_UNSUPPORTED_OPERATION;
			goto done;
		}
	/* 0. every frame right-handed (a frame's axes are signless: its third turned where they make a left hand) */
	for (c = 0; c < nc; c++)
	{
		double* A = R + 9 * (size_t)c;
		double d = A[0] * (A[4] * A[8] - A[5] * A[7]) - A[1] * (A[3] * A[8] - A[5] * A[6]) + A[2] * (A[3] * A[7] - A[4] * A[6]);
		if (d < 0)
			for (k = 6; k < 9; k++)
				A[k] = -A[k];
	}
	/* 1. the frames combed along a spanning tree of the cells (its facets then agree without a permutation) */
	for (c = 0; c < nc; c++)
	{
		unsigned int nq = 0, head;
		if (seen[c])
			continue;
		seen[c] = 1;
		queue[nq++] = c;
		for (head = 0; head < nq; head++)
			for (f = 0; f < 4; f++)
			{
				unsigned int s = t->f0[queue[head]] + f, o = t->fo[s], d;
				double old[9];
				int P[3];
				if (o == HM_NONE || seen[o / 4])
					continue;
				d = o / 4;
				hs_match(R + 9 * (size_t)queue[head], R + 9 * (size_t)d, P);
				memcpy(old, R + 9 * (size_t)d, sizeof old);
				for (i = 0; i < 3; i++)
				{
					int jj = (P[i] < 0 ? -P[i] : P[i]) - 1;
					for (k = 0; k < 3; k++)
						R[9 * (size_t)d + 3 * i + k] = (P[i] < 0 ? -1 : 1) * old[3 * jj + k];
				}
				seen[d] = 1;
				tree[s] = tree[o] = 1;
				queue[nq++] = d;
			}
	}
	/* each inner facet's permutation (from its cell's side: u_c = Pi u_d + t); the others cut */
	for (c = 0; c < nc; c++)
		for (f = 0; f < 4; f++)
		{
			unsigned int s = t->f0[c] + f, o = t->fo[s];
			if (o == HM_NONE)
				continue;
			if (tree[s])
			{
				Ps[3 * s] = 1; Ps[3 * s + 1] = 2; Ps[3 * s + 2] = 3;
				continue;
			}
			hs_match(R + 9 * (size_t)c, R + 9 * (size_t)(o / 4), Ps + 3 * (size_t)s);
			cut[s] = 1;
		}
	/* 2. the edges (sorted keys), each cell's six; per edge a cell; singular (the permutations round it compose to a turn)
	   and border edges (their fan opens) */
	for (c = 0; c < nc; c++)
		for (k = 0; k < 6; k++)
		{
			static unsigned char const EA[6] = { 0, 0, 0, 1, 1, 2 }, EB[6] = { 1, 2, 3, 2, 3, 3 };
			unsigned int a = t->cv[t->c0[c] + EA[k]], b = t->cv[t->c0[c] + EB[k]];
			ekey[ne++] = a < b ? ((unsigned long long)a << 32) | b : ((unsigned long long)b << 32) | a;
		}
	qsort(ekey, ne, sizeof *ekey, hs_cmp_u64);
	for (i = 0, j = 0; i < ne; i++)
		if (!j || ekey[j - 1] != ekey[i])
			ekey[j++] = ekey[i];
	ne = j;
	ecell = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (size_t)ne + 1));
	esing = (unsigned char*)qaws_internal_zalloc((size_t)ne + 1, 1);
	ebord = (unsigned char*)qaws_internal_zalloc((size_t)ne + 1, 1);
	ecnt = (unsigned int*)qaws_internal_zalloc((size_t)ne + 1, sizeof(unsigned int));
	fe = (unsigned int*)qaws_internal_malloc(qaws_internal_mul3(sizeof(unsigned int), 3, (size_t)ns + 1));
	if (!ecell || !esing || !ebord || !ecnt || !fe)
		goto done;
	for (c = 0; c < nc; c++)
		for (k = 0; k < 6; k++)
		{
			static unsigned char const EA[6] = { 0, 0, 0, 1, 1, 2 }, EB[6] = { 1, 2, 3, 2, 3, 3 };
			unsigned int a = t->cv[t->c0[c] + EA[k]], b = t->cv[t->c0[c] + EB[k]];
			ecell[hs_find(ekey, ne, a < b ? ((unsigned long long)a << 32) | b : ((unsigned long long)b << 32) | a)] = c;
		}
	/* each facet slot's three edges (edge m: corners m, m + 1 of the facet) */
	for (c = 0; c < nc; c++)
		for (f = 0; f < 4; f++)
		{
			unsigned int n3, s = t->f0[c] + f;
			unsigned char const* loc = vm_facet(QAWS_CELL_TETRAHEDRON, f, &n3);
			for (k = 0; k < 3; k++)
			{
				unsigned int a = t->cv[t->c0[c] + loc[k]], b = t->cv[t->c0[c] + loc[(k + 1) % 3]];
				fe[3 * (size_t)s + k] = hs_find(ekey, ne, a < b ? ((unsigned long long)a << 32) | b : ((unsigned long long)b << 32) | a);
			}
		}
	for (i = 0; i < ne; i++)
	{
		unsigned int a = (unsigned int)(ekey[i] >> 32), b = (unsigned int)(ekey[i] & 0xFFFFFFFFu), cur = ecell[i], prev = HM_NONE, guard, start = ecell[i];
		int M[9] = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
		for (guard = 0; guard < 256; guard++)
		{
			unsigned int nx = HM_NONE, s = HM_NONE;
			for (f = 0; f < 4; f++)
			{
				unsigned int w = t->cv[t->c0[cur] + f], o = t->fo[t->f0[cur] + f];
				if (w == a || w == b)
					continue;   /* (the facet opposite w holds the edge) */
				if (o == HM_NONE)
				{
					ebord[i] = 1;
					continue;
				}
				if (o / 4 == prev)
					continue;
				nx = o / 4;
				s = t->f0[cur] + f;
				break;
			}
			if (nx == HM_NONE || ebord[i])
				break;
			{
				int Pm[9];
				hs_mat(Ps + 3 * (size_t)s, Pm);
				hs_mul(M, Pm, M);
			}
			prev = cur;
			cur = nx;
			if (cur == start)
				break;
		}
		if (!ebord[i] && cur != start)
		{
			st = QAWS_STATUS_UNSUPPORTED_OPERATION;   /* (an edge whose cells do not close round it) */
			goto done;
		}
		esing[i] = !ebord[i] && !(M[0] == 1 && M[4] == 1 && M[8] == 1);
	}
	/* 3. the cut reduced: a cut facet with an edge no other cut facet holds (an inner regular edge, or a border edge where
	   the facet's permutation is none) is glued; the singular edges stay in the cut */
	for (c = 0; c < nc; c++)
		for (f = 0; f < 4; f++)
		{
			unsigned int s = t->f0[c] + f;
			if (cut[s] && t->fo[s] > s)
				for (k = 0; k < 3; k++)
					ecnt[fe[3 * (size_t)s + k]]++;
		}
	{
		/* the cut facets round each edge (CSR) */
		unsigned int tot = 0, nq = 0, head;
		e0 = (unsigned int*)qaws_internal_zalloc((size_t)ne + 2, sizeof(unsigned int));
		if (!e0)
			goto done;
		for (i = 0; i < ne; i++)
		{
			e0[i] = tot;
			tot += ecnt[i];
		}
		e0[ne] = tot;
		ef = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (size_t)tot + 1));
		stack = (unsigned int*)qaws_internal_zalloc((size_t)ne + 1, sizeof(unsigned int));
		qaws_internal_free(queue);
		queue = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), 4 * (size_t)ne + 1));
		if (!ef || !stack || !queue)
			goto done;
		for (c = 0; c < nc; c++)
			for (f = 0; f < 4; f++)
			{
				unsigned int s = t->f0[c] + f;
				if (cut[s] && t->fo[s] > s)
					for (k = 0; k < 3; k++)
					{
						unsigned int e = fe[3 * (size_t)s + k];
						ef[e0[e] + stack[e]++] = s;
					}
			}
		for (i = 0; i < ne; i++)
			if (ecnt[i] == 1 && !esing[i])
				queue[nq++] = i;
		for (head = 0; head < nq; head++)
		{
			unsigned int e = queue[head], s = HM_NONE;
			if (ecnt[e] != 1)
				continue;
			for (k = e0[e]; k < e0[e + 1]; k++)
				if (cut[ef[k]])
					s = ef[k];
			if (s == HM_NONE)
				continue;
			if (ebord[e] && !(Ps[3 * (size_t)s] == 1 && Ps[3 * (size_t)s + 1] == 2 && Ps[3 * (size_t)s + 2] == 3))
				continue;
			cut[s] = cut[t->fo[s]] = 0;
			Ps[3 * (size_t)t->fo[s]] = 1; Ps[3 * (size_t)t->fo[s] + 1] = 2; Ps[3 * (size_t)t->fo[s] + 2] = 3;
			for (k = 0; k < 3; k++)
			{
				unsigned int e2 = fe[3 * (size_t)s + k];
				ecnt[e2]--;
				if (ecnt[e2] == 1 && !esing[e2] && nq < 4 * ne)
					queue[nq++] = e2;
			}
		}
	}
	/* 4. the copies: corners joined across the facets not cut (same vertex) */
	uf = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (size_t)ns + 1));
	copy = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (size_t)ns + 1));
	if (!uf || !copy)
		goto done;
	for (i = 0; i < ns; i++)
		uf[i] = i;
	for (c = 0; c < nc; c++)
		for (f = 0; f < 4; f++)
		{
			unsigned int s = t->f0[c] + f, o = t->fo[s], d;
			if (o == HM_NONE || cut[s] || o < s)
				continue;
			d = o / 4;
			for (q = 0; q < 4; q++)
				for (k = 0; k < 4; k++)
					if (q != f && t->cv[t->c0[c] + q] == t->cv[t->c0[d] + k])
					{
						unsigned int ra = hs_uf(uf, 4 * c + q), rb = hs_uf(uf, 4 * d + k);
						if (ra != rb)
							uf[ra] = rb;
					}
		}
	for (i = 0; i < ns; i++)
		copy[i] = HM_NONE;
	for (i = 0; i < ns; i++)
		if (hs_uf(uf, i) == i)
			copy[i] = ncopy++;
	for (i = 0; i < ns; i++)
		copy[i] = copy[hs_uf(uf, i)];
	cx = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 3, (size_t)ncopy + 1));
	if (!cx)
		goto done;
	/* each copy's point for the factor's ordering: the mean of its cells' centres (the copies of a vertex apart, as their
	   cells are) */
	{
		double* cw = (double*)qaws_internal_zalloc((size_t)ncopy + 1, sizeof(double));
		if (!cw)
			goto done;
		memset(cx, 0, sizeof(double) * 3 * ncopy);
		for (i = 0; i < ns; i++)
		{
			for (k = 0; k < 3; k++)
				for (q = 0; q < 4; q++)
					cx[3 * (size_t)copy[i] + k] += 0.25 * t->x[3 * (size_t)t->cv[t->c0[i / 4] + q] + k];
			cw[copy[i]] += 1;
		}
		for (i = 0; i < 3 * ncopy; i++)
			cx[i] /= cw[i / 3];
		qaws_internal_free(cw);
	}
	/* 5. the constraints and the candidates. A cut facet: its three vertices moved alike (u_c = Pi u_d + t, t the
	   candidate); a border facet: its vertices on one plane of the axis its normal follows (the plane the candidate); a
	   singular edge: its vertices on whole numbers across the axis it turns about */
	for (c = 0; c < nc; c++)
		for (f = 0; f < 4; f++)
		{
			unsigned int s = t->f0[c] + f, o = t->fo[s], n3, qc[3], qd[3];
			unsigned char const* loc = vm_facet(QAWS_CELL_TETRAHEDRON, f, &n3);
			if (o != HM_NONE && (!cut[s] || o < s))
				continue;
			for (k = 0; k < 3; k++)
				qc[k] = loc[k];
			if (o == HM_NONE)
			{
				/* the border: the axis nearest the facet's normal */
				double const *a = t->x + 3 * (size_t)t->cv[t->c0[c] + loc[0]], *b = t->x + 3 * (size_t)t->cv[t->c0[c] + loc[1]],
					*e = t->x + 3 * (size_t)t->cv[t->c0[c] + loc[2]];
				double N[3], dmax = -1;
				unsigned int ax = 0;
				N[0] = (b[1] - a[1]) * (e[2] - a[2]) - (b[2] - a[2]) * (e[1] - a[1]);
				N[1] = (b[2] - a[2]) * (e[0] - a[0]) - (b[0] - a[0]) * (e[2] - a[2]);
				N[2] = (b[0] - a[0]) * (e[1] - a[1]) - (b[1] - a[1]) * (e[0] - a[0]);
				for (k = 0; k < 3; k++)
				{
					double const* A = R + 9 * (size_t)c + 3 * k;
					double dd = fabs(A[0] * N[0] + A[1] * N[1] + A[2] * N[2]);
					if (dd > dmax)
					{
						dmax = dd;
						ax = k;
					}
				}
				for (k = 1; k < 3; k++)
				{
					hs_form* fm = hs_add(&con, 0);
					if (!fm)
						goto done;
					hs_term(fm, 3 * copy[4 * c + qc[k]] + ax, 1);
					hs_term(fm, 3 * copy[4 * c + qc[0]] + ax, -1);
				}
				{
					hs_form* fm = hs_add(&cand, 1);
					if (!fm)
						goto done;
					hs_term(fm, 3 * copy[4 * c + qc[0]] + ax, 1);
				}
				continue;
			}
			/* the cut facet's vertices in the cell across */
			for (k = 0; k < 3; k++)
			{
				unsigned int v = t->cv[t->c0[c] + qc[k]], z;
				qd[k] = HM_NONE;
				for (z = 0; z < 4; z++)
					if (t->cv[t->c0[o / 4] + z] == v)
						qd[k] = z;
			}
			for (i = 0; i < 3; i++)
			{
				int const* P = Ps + 3 * (size_t)s;
				unsigned int jj = (unsigned int)((P[i] < 0 ? -P[i] : P[i]) - 1);
				double sg = P[i] < 0 ? -1 : 1;
				hs_form* fm;
				for (k = 1; k < 3; k++)
				{
					fm = hs_add(&con, 0);
					if (!fm)
						goto done;
					hs_term(fm, 3 * copy[4 * c + qc[k]] + i, 1);
					hs_term(fm, 3 * copy[4 * (o / 4) + qd[k]] + jj, -sg);
					hs_term(fm, 3 * copy[4 * c + qc[0]] + i, -1);
					hs_term(fm, 3 * copy[4 * (o / 4) + qd[0]] + jj, sg);
				}
				fm = hs_add(&cand, 2);
				if (!fm)
					goto done;
				hs_term(fm, 3 * copy[4 * c + qc[0]] + i, 1);
				hs_term(fm, 3 * copy[4 * (o / 4) + qd[0]] + jj, -sg);
			}
		}
	for (i = 0; i < ne; i++)
	{
		unsigned int a = (unsigned int)(ekey[i] >> 32), b = (unsigned int)(ekey[i] & 0xFFFFFFFFu), cur = ecell[i], prev = HM_NONE, guard, z;
		int M[9] = { 1, 0, 0, 0, 1, 0, 0, 0, 1 }, ax = -1;
		if (!esing[i])
			continue;
		for (guard = 0; guard < 256; guard++)
		{
			unsigned int nx = HM_NONE, s = HM_NONE;
			int Pm[9];
			for (f = 0; f < 4; f++)
			{
				unsigned int w = t->cv[t->c0[cur] + f], o = t->fo[t->f0[cur] + f];
				if (w == a || w == b || o == HM_NONE || o / 4 == prev)
					continue;
				nx = o / 4;
				s = t->f0[cur] + f;
				break;
			}
			hs_mat(Ps + 3 * (size_t)s, Pm);
			hs_mul(M, Pm, M);
			prev = cur;
			cur = nx;
			if (cur == ecell[i])
				break;
		}
		/* a quarter turn about one axis of the starting cell's chart: its vertices on whole numbers across it; another turn
		   (a third about a diagonal, where curves meet on an edge) no hexahedra follow: its cells squeezed flat are left out */
		for (k = 0; k < 3; k++)
			if (M[4 * k] == 1)
				ax = (int)k;
		if (ax < 0 || M[0] + M[4] + M[8] != 1)
			continue;
		for (z = 0; z < 4; z++)
		{
			unsigned int v = t->cv[t->c0[ecell[i]] + z];
			if (v != a && v != b)
				continue;
			for (k = 0; k < 3; k++)
				if ((int)k != ax)
				{
					hs_form* fm = hs_add(&cand, 0);
					if (!fm)
						goto done;
					hs_term(fm, 3 * copy[4 * ecell[i] + z] + k, 1);
				}
		}
	}
	/* 6. the rounding (MIQ's greedy order: the singular edges, the border's planes, then the transitions): solved, the
	   candidates already whole fixed there, then the one nearest a whole number of the first kind left (with those within
	   a fiftieth of one), solved again */
	U = (double*)qaws_internal_zalloc(3 * (size_t)ncopy + 1, sizeof(double));
	fixv = (double*)qaws_internal_zalloc((size_t)cand.n + 1, sizeof(double));
	isfix = (unsigned char*)qaws_internal_zalloc((size_t)cand.n + 1, 1);
	if (!U || !fixv || !isfix)
		goto done;
	for (round = 0;; round++)
	{
		unsigned int best = HM_NONE, kind;
		double bd = 2;
		st = hs_solve(t, R, h, copy, ncopy, cx, &con, &cand, fixv, isfix, U);
		if (st != QAWS_STATUS_OK)
			goto done;
		for (kind = 0; kind < 3 && best == HM_NONE; kind++)
			for (i = 0; i < cand.n; i++)
			{
				double v, d;
				if (cand.f[i].kind != kind || isfix[i])
					continue;
				v = hs_eval(&cand.f[i], U);
				d = fabs(v - floor(v + 0.5));
				if (d < 1e-4)
				{
					fixv[i] = floor(v + 0.5);
					isfix[i] = 1;
					continue;
				}
				if (d < bd)
				{
					bd = d;
					best = i;
				}
			}
		if (best == HM_NONE)
			break;
		if (round == HS_ROUNDS)
		{
			st = QAWS_STATUS_NOT_CONVERGED;
			goto done;
		}
		for (i = 0; i < cand.n; i++)
			if (!isfix[i] && cand.f[i].kind == cand.f[best].kind)
			{
				double v = hs_eval(&cand.f[i], U), d = fabs(v - floor(v + 0.5));
				if (i == best || d < 0.02)
				{
					fixv[i] = floor(v + 0.5);
					isfix[i] = 1;
				}
			}
	}
	/* 7. the cut facets' transitions (whole), each cell's chart; every cell turned the right way in the map */
	par = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 12, (size_t)nc + 1));
	flat = (unsigned char*)qaws_internal_zalloc((size_t)nc + 1, 1);
	if (!par || !flat)
		goto done;
	for (c = 0; c < nc; c++)
	{
		for (q = 0; q < 4; q++)
			memcpy(par + 12 * (size_t)c + 3 * q, U + 3 * (size_t)copy[4 * c + q], 3 * sizeof(double));
		{
			double l[4], z[3] = { 0, 0, 0 };
			st = QAWS_STATUS_UNSUPPORTED_OPERATION;
			double g[4][3], V = hs_grads(t, c, g), dm = hs_bary(par + 12 * (size_t)c, z, l) / 6, ref = V / (h * h * h);
			/* a cell flat or turned in the map (squeezed beside a singular edge, or by the rounding) is left out of the lattice,
			   the cubes it hid found from their corners; more than one in a hundred turned is no map */
			flat[c] = !(dm > 1e-6 * ref);
			nturn += dm < -0.05 * ref;
		}
	}
	if (nturn > nc / 100)
		goto done;
	for (c = 0; c < nc; c++)
		for (f = 0; f < 4; f++)
		{
			unsigned int s = t->f0[c] + f, o = t->fo[s], n3, v, z, qd = HM_NONE;
			unsigned char const* loc = vm_facet(QAWS_CELL_TETRAHEDRON, f, &n3);
			if (o == HM_NONE || !cut[s])
				continue;
			v = t->cv[t->c0[c] + loc[0]];
			for (z = 0; z < 4; z++)
				if (t->cv[t->c0[o / 4] + z] == v)
					qd = z;
			for (i = 0; i < 3; i++)
			{
				int const* P = Ps + 3 * (size_t)s;
				unsigned int jj = (unsigned int)((P[i] < 0 ? -P[i] : P[i]) - 1);
				double tv = par[12 * (size_t)c + 3 * loc[0] + i] - (P[i] < 0 ? -1 : 1) * par[12 * (size_t)(o / 4) + 3 * qd + jj];
				Ts[3 * (size_t)s + i] = (int)floor(tv + 0.5);
			}
		}
	/* 8. the lattice's points (doubled coordinates: corners even, cube centres odd) in each cell's closure, registered;
	   one point where a facet holds it seen from both sides */
	for (c = 0; c < nc; c++)
	{
		double lo[3], hi[3];
		double const* P4 = par + 12 * (size_t)c;
		int a0[3], a1[3], p[3];
		if (flat[c])
			continue;
		for (k = 0; k < 3; k++)
		{
			lo[k] = hi[k] = P4[k];
			for (q = 1; q < 4; q++)
			{
				lo[k] = fmin(lo[k], P4[3 * q + k]);
				hi[k] = fmax(hi[k], P4[3 * q + k]);
			}
			a0[k] = (int)ceil(2 * lo[k] - 1e-6);
			a1[k] = (int)floor(2 * hi[k] + 1e-6);
			if (a0[k] < -2000 || a1[k] > 2000)
			{
				st = QAWS_STATUS_UNSUPPORTED_OPERATION;   /* (a lattice too wide for the keys) */
				goto done;
			}
		}
		for (p[0] = a0[0]; p[0] <= a1[0]; p[0]++)
			for (p[1] = a0[1]; p[1] <= a1[1]; p[1]++)
				for (p[2] = a0[2]; p[2] <= a1[2]; p[2]++)
				{
					double u[3], l[4];
					int odd = (p[0] & 1) + (p[1] & 1) + (p[2] & 1);
					if (odd != 0 && odd != 3)
						continue;   /* (edge and face midpoints are not needed) */
					for (k = 0; k < 3; k++)
						u[k] = 0.5 * p[k];
					hs_bary(P4, u, l);
					if (l[0] < -HS_TOL || l[1] < -HS_TOL || l[2] < -HS_TOL || l[3] < -HS_TOL)
						continue;
					if (nreg == capreg)
					{
						unsigned int nn = qaws_internal_grow(capreg, (size_t)nreg + 1, 1024);
						hs_reg* g2 = nn ? (hs_reg*)qaws_internal_realloc(reg, qaws_internal_mul2(sizeof(hs_reg), nn)) : NULL;
						if (!g2)
						{
							st = QAWS_STATUS_ALLOCATION_FAILURE;
							goto done;
						}
						reg = g2;
						capreg = nn;
					}
					reg[nreg].c = c;
					memcpy(reg[nreg].p, p, sizeof p);
					memcpy(reg[nreg].l, l, sizeof l);
					nreg++;
				}
	}
	st = QAWS_STATUS_ALLOCATION_FAILURE;
	rkey = (unsigned long long*)qaws_internal_malloc(qaws_internal_mul3(sizeof(unsigned long long), 2, (size_t)nreg + 1));
	rid = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (size_t)nreg + 1));
	qaws_internal_free(uf);
	uf = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (size_t)nreg + 1));
	if (!rkey || !rid || !uf)
		goto done;
	for (i = 0; i < nreg; i++)
	{
		rkey[2 * (size_t)i] = hs_rkey(reg[i].c, reg[i].p);
		rkey[2 * (size_t)i + 1] = i;
		uf[i] = i;
	}
	qsort(rkey, nreg, 2 * sizeof *rkey, hs_cmp_u64);
	for (i = 0; i < nreg; i++)
		rid[i] = (unsigned int)rkey[2 * (size_t)i + 1];
	for (i = 0; i < nreg; i++)
		rkey[i] = rkey[2 * (size_t)i];   /* (the keys alone, sorted; rid: their registrations) */
	for (i = 0; i < nreg; i++)
		for (q = 0; q < 4; q++)
		{
			unsigned int s = t->f0[reg[i].c] + q, o = t->fo[s], at;
			int p2[3];
			if (reg[i].l[q] > HS_TOL || o == HM_NONE)
				continue;
			hs_across(Ps + 3 * (size_t)s, Ts + 3 * (size_t)s, reg[i].p, 2, p2);
			at = hs_find(rkey, nreg, hs_rkey(o / 4, p2));
			if (at != HM_NONE)
			{
				unsigned int ra = hs_uf(uf, i), rb = hs_uf(uf, rid[at]);
				if (ra != rb)
					uf[ra] = rb;
			}
		}
	/* the corners' points: their positions from a registration */
	qaws_internal_free(copy);
	copy = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (size_t)nreg + 1));   /* (per root its point) */
	px = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 3, (size_t)nreg + 1));
	cappx = nreg + 1;
	if (!copy || !px)
		goto done;
	for (i = 0; i < nreg; i++)
		copy[i] = HM_NONE;
	for (i = 0; i < nreg; i++)
	{
		unsigned int r = hs_uf(uf, i);
		if ((reg[i].p[0] & 1) || copy[r] != HM_NONE)
			continue;
		copy[r] = np;
		for (k = 0; k < 3; k++)
		{
			px[3 * (size_t)np + k] = 0;
			for (q = 0; q < 4; q++)
				px[3 * (size_t)np + k] += reg[i].l[q] * t->x[3 * (size_t)t->cv[t->c0[reg[i].c] + q] + k];
		}
		np++;
	}
	/* 9. the cubes: round each cube centre (one registration per point), then round each point the eight cubes it is a
	   corner of (where a cell left out of the lattice hid the centre); a cube's corners found by hs_corner, kept when they
	   are eight points, the cube positive and (from a point) its middle inside; each cube once (its sorted points) */
	nbq = nc < 4096 ? nc + 1 : 4096;
	bq = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (size_t)nbq + 1));
	bk = (int*)qaws_internal_malloc(qaws_internal_mul3(sizeof(int), 3, (size_t)nbq + 1));
	vis = (unsigned int*)qaws_internal_zalloc((size_t)nc + 1, sizeof(unsigned int));
	qaws_internal_free(seen);
	seen = (unsigned char*)qaws_internal_zalloc((size_t)nreg + 1, 1);
	if (!seen || !bq || !bk || !vis)
		goto done;
	{
		static int const CO[8][3] = { { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 0, 1, 0 }, { 0, 0, 1 }, { 1, 0, 1 }, { 1, 1, 1 }, { 0, 1, 1 } };
		hs_lat L;
		unsigned int pass, oc;
		L.t = t; L.par = par; L.flat = flat; L.Ps = Ps; L.Ts = Ts; L.rkey = rkey; L.rid = rid; L.copy = copy; L.uf = uf; L.nreg = nreg;
		L.px = px; L.bq = bq; L.vis = vis; L.nbq = nbq; L.stamp = 0; L.bk = bk; L.reach2 = 1.5 * 1.5 * 0.75 * h * h * 1.5;
		hs_faces FT;
		memset(&FT, 0, sizeof FT);
		for (pass = 0; pass < 2; pass++)
		{
			if (pass == 1)
			{
				/* the centres' cubes once each, their faces counted */
				qsort(cells, ncell, 16 * sizeof(unsigned int), hs_cmp_cube);
				for (i = 0, k = 0; i < ncell; i++)
					if (!k || hs_cmp_cube(cells + 16 * (size_t)(k - 1), cells + 16 * (size_t)i) != 0)
					{
						memmove(cells + 16 * (size_t)k, cells + 16 * (size_t)i, 16 * sizeof(unsigned int));
						k++;
					}
				ncell = k;
				for (i = 0; i < ncell; i++)
				{
					unsigned int cnt[6];
					if (!hs_cube_faces(&FT, cells + 16 * (size_t)i, cnt, 1))
					{
						qaws_internal_free(FT.key);
						qaws_internal_free(FT.n);
						goto done;
					}
				}
			}
			for (i = 0; i < nreg; i++)
			{
				unsigned int r = hs_uf(uf, i);
				int odd = reg[i].p[0] & 1;
				if (odd != (pass == 0) || seen[r])
					continue;
				seen[r] = 1;
				for (oc = 0; oc < (pass ? 8u : 1u); oc++)
				{
					unsigned int v[16], z, all = 1;
					int C2[3], placed = 0;
					double cxyz[3] = { 0, 0, 0 };
					for (k = 0; k < 3; k++)
					{
						C2[k] = pass ? reg[i].p[k] + 2 * CO[oc][k] - 1 : reg[i].p[k];
						for (q = 0; q < 4; q++)
							cxyz[k] += reg[i].l[q] * t->x[3 * (size_t)t->cv[t->c0[reg[i].c] + q] + k];
					}
					{
						unsigned int miss = 0, mz = 0;
						for (z = 0; z < 8 && miss <= (pass ? 1u : 0u); z++)
						{
							int K2[3];
							for (k = 0; k < 3; k++)
								K2[k] = C2[k] + 2 * CO[z][k] - 1;
							v[z] = hs_corner(&L, reg[i].c, pass ? reg[i].p : C2, K2, cxyz);
							if (v[z] == HM_NONE)
							{
								miss++;
								mz = z;
							}
						}
						all = miss == 0;
						if (pass && miss == 1 && z == 8)
						{
							/* one corner not in the lattice (its cells squeezed flat at the border): placed by the box's rule from the
							   other seven (its edge neighbours, less the face diagonals, plus the opposite corner: exact on a
							   parallelepiped), an existing point within a quarter of a cell taken instead */
							double y[3] = { 0, 0, 0 };
							unsigned int best = HM_NONE, w, bm = (unsigned int)(CO[mz][0] | CO[mz][1] << 1 | CO[mz][2] << 2);
							double bd = 0.0625 * h * h;
							for (w = 0; w < 8; w++)
							{
								unsigned int d = bm ^ (unsigned int)(CO[w][0] | CO[w][1] << 1 | CO[w][2] << 2), hd = (d & 1) + ((d >> 1) & 1) + ((d >> 2) & 1);
								double sg = hd == 2 ? -1 : 1;
								if (!hd)
									continue;
								for (k = 0; k < 3; k++)
									y[k] += sg * px[3 * (size_t)v[w] + k];
							}
							for (w = 0; w < np; w++)
							{
								double d2 = (px[3 * (size_t)w] - y[0]) * (px[3 * (size_t)w] - y[0]) + (px[3 * (size_t)w + 1] - y[1]) * (px[3 * (size_t)w + 1] - y[1]) +
									(px[3 * (size_t)w + 2] - y[2]) * (px[3 * (size_t)w + 2] - y[2]);
								if (d2 < bd)
								{
									bd = d2;
									best = w;
								}
							}
							if (best == HM_NONE)
							{
								if (np == cappx)
								{
									unsigned int nn = qaws_internal_grow(cappx, (size_t)np + 1, 64);
									double* g2 = nn ? (double*)qaws_internal_realloc(px, qaws_internal_mul3(sizeof(double), 3, nn)) : NULL;
									if (!g2)
										goto done;
									px = g2;
									L.px = px;
									cappx = nn;
								}
								memcpy(px + 3 * (size_t)np, y, sizeof y);
								best = np++;
							}
							v[mz] = best;
							all = 1;
							placed = 1;
						}
					}
					for (z = 0; all && z < 8; z++)
						for (k = z + 1; k < 8; k++)
							all &= v[z] != v[k];
					if (!all)
						continue;
					if (pass)
					{
						/* from a point: the cube positive, its middle inside the cells */
						double mid[3] = { 0, 0, 0 };
						unsigned int vv[8];
						for (z = 0; z < 8; z++)
						{
							vv[z] = v[z];
							for (k = 0; k < 3; k++)
								mid[k] += px[3 * (size_t)v[z] + k] / 8;
						}
						if (!(vm_hex_quality(px, vv) > (placed ? 0.2 : -0.5)) || !hs_inside(&L, reg[i].c, mid))
							continue;
						/* a gap filled: no face held twice already, one held once */
						{
							unsigned int cnt[6], ok = 1, touch = 0;
							hs_cube_faces(&FT, vv, cnt, 0);
							for (z = 0; z < 6; z++)
							{
								ok &= cnt[z] <= 1;
								touch |= cnt[z] == 1;
							}
							if (!ok || !touch)
								continue;
							if (!hs_cube_faces(&FT, vv, cnt, 1))
							{
								qaws_internal_free(FT.key);
								qaws_internal_free(FT.n);
								goto done;
							}
						}
					}
					/* the sorted points after the cube's own (the key) */
					for (z = 0; z < 8; z++)
						v[8 + z] = v[z];
					for (z = 1; z < 8; z++)
						for (k = z; k > 0 && v[8 + k - 1] > v[8 + k]; k--)
						{
							unsigned int sw = v[8 + k];
							v[8 + k] = v[8 + k - 1];
							v[8 + k - 1] = sw;
						}
					if (ncell + 1 > capc)
					{
						unsigned int nn = qaws_internal_grow(capc, (size_t)ncell + 1, 256);
						unsigned int* g2 = nn ? (unsigned int*)qaws_internal_realloc(cells, qaws_internal_mul3(sizeof(unsigned int), 16, nn)) : NULL;
						if (!g2)
							goto done;
						cells = g2;
						capc = nn;
					}
					memcpy(cells + 16 * (size_t)ncell++, v, sizeof v);
				}
			}
		}
		qaws_internal_free(FT.key);
		qaws_internal_free(FT.n);
		/* each cube once; then 8 points each */
		qsort(cells, ncell, 16 * sizeof(unsigned int), hs_cmp_cube);
		for (i = 0, k = 0; i < ncell; i++)
			if (!k || hs_cmp_cube(cells + 16 * (size_t)(k - 1), cells + 16 * (size_t)i) != 0)
			{
				memmove(cells + 16 * (size_t)k, cells + 16 * (size_t)i, 16 * sizeof(unsigned int));
				k++;
			}
		ncell = k;
		/* every face held by two cubes at most, turned opposite ways: the cubes taken best first (their scaled Jacobian),
		   one overlapping those taken left out */
		{
			double* qi = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 2, (size_t)ncell + 1));
			unsigned int* keep = (unsigned int*)qaws_internal_malloc(qaws_internal_mul3(sizeof(unsigned int), 16, (size_t)ncell + 1));
			unsigned int* dir = NULL, nk = 0;
			hs_faces G;
			memset(&G, 0, sizeof G);
			G.cap = 1024;
			while (G.cap < 16 * ncell)
				G.cap *= 2;
			G.key = (unsigned int*)qaws_internal_malloc(qaws_internal_mul3(sizeof(unsigned int), 4, G.cap));
			G.n = (unsigned char*)qaws_internal_zalloc(G.cap, 1);
			dir = (unsigned int*)qaws_internal_malloc(qaws_internal_mul3(sizeof(unsigned int), 4, G.cap));
			if (!qi || !keep || !G.key || !G.n || !dir)
			{
				qaws_internal_free(qi); qaws_internal_free(keep); qaws_internal_free(G.key); qaws_internal_free(G.n); qaws_internal_free(dir);
				goto done;
			}
			for (i = 0; i < G.cap; i++)
				G.key[4 * (size_t)i] = HM_NONE;
			for (i = 0; i < ncell; i++)
			{
				qi[2 * (size_t)i] = vm_hex_quality(px, cells + 16 * (size_t)i);
				qi[2 * (size_t)i + 1] = i;
			}
			qsort(qi, ncell, 2 * sizeof(double), hs_cmp_desc);
			for (i = 0; i < ncell; i++)
			{
				unsigned int ci = (unsigned int)qi[2 * (size_t)i + 1];
				if (hs_cube_fits(&G, dir, cells + 16 * (size_t)ci))
					memcpy(keep + 16 * (size_t)nk++, cells + 16 * (size_t)ci, 16 * sizeof(unsigned int));
			}
			memcpy(cells, keep, sizeof(unsigned int) * 16 * nk);
			ncell = nk;
			qaws_internal_free(qi); qaws_internal_free(keep); qaws_internal_free(G.key); qaws_internal_free(G.n); qaws_internal_free(dir);
		}
		for (i = 0; i < ncell; i++)
			memmove(cells + 8 * (size_t)i, cells + 16 * (size_t)i, 8 * sizeof(unsigned int));
	}
	if (!ncell)
	{
		st = QAWS_STATUS_INVALID_ARGUMENT;   /* (no cube fits: a cell size too large) */
		goto done;
	}
	{
		/* the points the hexahedra use, compacted */
		unsigned int* map = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (size_t)np + 1)), nv2 = 0;
		double* x2 = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 3, (size_t)np + 1));
		unsigned char* types = (unsigned char*)qaws_internal_malloc((size_t)ncell + 1);
		if (!map || !x2 || !types)
		{
			qaws_internal_free(map); qaws_internal_free(x2); qaws_internal_free(types);
			goto done;
		}
		for (i = 0; i < np; i++)
			map[i] = HM_NONE;
		for (i = 0; i < 8 * ncell; i++)
		{
			if (map[cells[i]] == HM_NONE)
			{
				memcpy(x2 + 3 * (size_t)nv2, px + 3 * (size_t)cells[i], sizeof(double) * 3);
				map[cells[i]] = nv2++;
			}
			cells[i] = map[cells[i]];
		}
		memset(types, QAWS_CELL_HEXAHEDRON, ncell);
		st = vm_create(x2, nv2, cells, types, ncell, out_mesh);
		qaws_internal_free(map); qaws_internal_free(x2); qaws_internal_free(types);
	}
	if (st == QAWS_STATUS_OK && rep)
	{
		qaws_mesh_3d* bd = NULL;
		if (qaws_volume_mesh_get_boundary(&bd, *out_mesh) == QAWS_STATUS_OK)
			qaws_mesh_3d_measure_quality(&rep->quality, bd);
		qaws_mesh_3d_destroy(bd);
	}
	(void)nv;
done:
	qaws_internal_free(queue); qaws_internal_free(uf); qaws_internal_free(copy); qaws_internal_free(ecell); qaws_internal_free(fe); qaws_internal_free(ecnt);
	qaws_internal_free(e0); qaws_internal_free(ef); qaws_internal_free(stack); qaws_internal_free(tree); qaws_internal_free(cut); qaws_internal_free(esing);
	qaws_internal_free(ebord); qaws_internal_free(seen); qaws_internal_free(ekey); qaws_internal_free(rkey); qaws_internal_free(Ps); qaws_internal_free(Ts);
	qaws_internal_free(U); qaws_internal_free(cx); qaws_internal_free(fixv); qaws_internal_free(par); qaws_internal_free(px); qaws_internal_free(rid);
	qaws_internal_free(cells); qaws_internal_free(reg); qaws_internal_free(con.f); qaws_internal_free(cand.f); qaws_internal_free(isfix); qaws_internal_free(flat);
	qaws_internal_free(bq); qaws_internal_free(bk); qaws_internal_free(vis);
	return st;
}

