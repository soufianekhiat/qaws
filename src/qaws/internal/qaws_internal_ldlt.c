#include "qaws_internal_ldlt.h"
#include "qaws_internal_memory.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define LD_NONE 0xFFFFFFFFu

/* ======================================================================== */
/*  The ordering: geometric nested dissection                               */
/* ======================================================================== */

typedef struct ld_key
{
	double v;
	unsigned int i;
} ld_key;

static int ld_key_cmp(void const* a, void const* b)
{
	ld_key const *p = (ld_key const*)a, *q = (ld_key const*)b;
	if (p->v != q->v)
		return p->v < q->v ? -1 : 1;
	return p->i < q->i ? -1 : p->i > q->i;
}

typedef struct ld_order
{
	unsigned int const *rp, *col;
	double const* xyz;
	unsigned int* stamp;            /* per row: the call whose right half it is in */
	unsigned int* tmp;
	ld_key* key;
	unsigned int gen;
} ld_order;

/* set (n rows) ordered in place: the left half less its separator, the right half, the separator (the left rows with a
   neighbour on the right); both halves dissected again, sets of 64 rows or fewer kept in their order */
static void ld_dissect(ld_order* o, unsigned int* set, unsigned int n)
{
	double lo[3], hi[3];
	unsigned int i, k, a = 0, half = n / 2, nl = 0, ns = 0, g;
	if (n <= 64)
		return;
	for (k = 0; k < 3; k++)
	{
		lo[k] = HUGE_VAL;
		hi[k] = -HUGE_VAL;
	}
	for (i = 0; i < n; i++)
		for (k = 0; k < 3; k++)
		{
			double v = o->xyz[3 * (size_t)set[i] + k];
			lo[k] = v < lo[k] ? v : lo[k];
			hi[k] = v > hi[k] ? v : hi[k];
		}
	for (k = 1; k < 3; k++)
		if (hi[k] - lo[k] > hi[a] - lo[a])
			a = k;
	for (i = 0; i < n; i++)
	{
		o->key[i].v = o->xyz[3 * (size_t)set[i] + a];
		o->key[i].i = set[i];
	}
	qsort(o->key, n, sizeof(ld_key), ld_key_cmp);
	g = ++o->gen;
	for (i = half; i < n; i++)
		o->stamp[o->key[i].i] = g;
	/* the left rows apart from the separator first, the separator at the end */
	for (i = 0; i < half; i++)
	{
		unsigned int r = o->key[i].i, p, sep = 0;
		for (p = o->rp[r]; p < o->rp[r + 1] && !sep; p++)
			sep = o->stamp[o->col[p]] == g;
		if (sep)
			o->tmp[n - 1 - ns++] = r;
		else
			o->tmp[nl++] = r;
	}
	for (i = half; i < n; i++)
		o->tmp[nl + i - half] = o->key[i].i;
	/* (the separator was written backwards) */
	for (i = 0; i < ns / 2; i++)
	{
		unsigned int t = o->tmp[n - 1 - i];
		o->tmp[n - 1 - i] = o->tmp[n - ns + i];
		o->tmp[n - ns + i] = t;
	}
	memcpy(set, o->tmp, sizeof(unsigned int) * n);
	ld_dissect(o, set, nl);
	ld_dissect(o, set + nl, n - half);
}

/* ======================================================================== */
/*  The factor                                                              */
/* ======================================================================== */

void qaws_internal_ldlt_free(qaws_internal_ldlt* f)
{
	qaws_internal_free(f->perm); qaws_internal_free(f->sf); qaws_internal_free(f->rp); qaws_internal_free(f->ri); qaws_internal_free(f->lp); qaws_internal_free(f->lx); qaws_internal_free(f->d); qaws_internal_free(f->w);
	memset(f, 0, sizeof *f);
}

typedef struct ld_work
{
	unsigned int *super, *map, *head, *next, *pos;
	double *acc, *t;
} ld_work;

/* the panels from C's lower triangle by columns (cp, ci, cx) with the shift; 0 on a pivot that is not positive */
static int ld_numeric(qaws_internal_ldlt* f, unsigned int const* cp, unsigned int const* ci, double const* cx, double shift, ld_work* wk)
{
	unsigned int s, c, i, k;
	for (s = 0; s < f->ns; s++)
		wk->head[s] = LD_NONE;
	for (s = 0; s < f->ns; s++)
	{
		unsigned int c0 = f->sf[s], w = f->sf[s + 1] - c0, nr = f->rp[s + 1] - f->rp[s], *R = f->ri + f->rp[s], d, dn;
		double* P = f->lx + f->lp[s];
		memset(P, 0, sizeof(double) * nr * w);
		for (i = 0; i < nr; i++)
			wk->map[R[i]] = i;
		for (c = 0; c < w; c++)
		{
			unsigned int p;
			for (p = cp[c0 + c]; p < cp[c0 + c + 1]; p++)
				P[wk->map[ci[p]] + (size_t)c * nr] += cx[p];
			P[c + (size_t)c * nr] += shift;
		}
		/* the supernodes below whose rows reach these columns: minus L_d (rows from here down) D_d L_d (rows here)^T */
		for (d = wk->head[s], wk->head[s] = LD_NONE; d != LD_NONE; d = dn)
		{
			unsigned int dw = f->sf[d + 1] - f->sf[d], dnr = f->rp[d + 1] - f->rp[d], *Rd = f->ri + f->rp[d], p1 = wk->pos[d], p2, jj, ii;
			double const *Ld = f->lx + f->lp[d], *Dd = f->d + f->sf[d];
			dn = wk->next[d];
			for (p2 = p1; p2 < dnr && Rd[p2] < c0 + w; p2++)
				;
			for (jj = p1; jj < p2; jj++)
			{
				unsigned int len = dnr - jj;
				double* Pc = P + (size_t)(Rd[jj] - c0) * nr;
				for (k = 0; k < dw; k++)
					wk->t[k] = Ld[jj + (size_t)k * dnr] * Dd[k];
				memset(wk->acc, 0, sizeof(double) * len);
				for (k = 0; k < dw; k++)
				{
					double tk = wk->t[k];
					double const* Lk = Ld + (size_t)k * dnr + jj;
					for (ii = 0; ii < len; ii++)
						wk->acc[ii] += Lk[ii] * tk;
				}
				for (ii = 0; ii < len; ii++)
					Pc[wk->map[Rd[jj + ii]]] -= wk->acc[ii];
			}
			wk->pos[d] = p2;
			if (p2 < dnr)
			{
				unsigned int to = wk->super[Rd[p2]];
				wk->next[d] = wk->head[to];
				wk->head[to] = d;
			}
		}
		/* the panel itself, dense, right-looking */
		for (c = 0; c < w; c++)
		{
			double* Pc = P + (size_t)c * nr;
			double dd = Pc[c], inv;
			unsigned int c2;
			if (!(dd > 0))
				return 0;
			f->d[c0 + c] = dd;
			inv = 1 / dd;
			for (i = c + 1; i < nr; i++)
				Pc[i] *= inv;
			for (c2 = c + 1; c2 < w; c2++)
			{
				double tt = Pc[c2] * dd, *P2 = P + (size_t)c2 * nr;
				for (i = c2; i < nr; i++)
					P2[i] -= Pc[i] * tt;
			}
		}
		wk->pos[s] = w;
		if (w < nr)
		{
			unsigned int to = wk->super[R[w]];
			wk->next[s] = wk->head[to];
			wk->head[to] = s;
		}
	}
	return 1;
}

qaws_status qaws_internal_ldlt_factor(qaws_internal_ldlt* f, unsigned int n, unsigned int const* row_ptr, unsigned int const* col,
	double const* val, double const* xyz, double shift)
{
	unsigned int *pinv = NULL, *rp = NULL, *rj = NULL, *cp = NULL, *ci = NULL, *parent = NULL, *flag = NULL, *lnz = NULL, *kids = NULL;
	double *cx = NULL, dmax = 0;
	unsigned int r, p, k, i, s, tries, nnzc, maxr = 1, maxw = 1;
	size_t total = 0;
	ld_work wk;
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	memset(f, 0, sizeof *f);
	memset(&wk, 0, sizeof wk);
	f->n = n;
	f->perm = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * (n + 1));
	f->d = (double*)qaws_internal_malloc(sizeof(double) * (n + 1));
	f->w = (double*)qaws_internal_malloc(sizeof(double) * (2 * (size_t)n + 1));
	f->sf = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * (n + 2));
	pinv = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * (n + 1));
	rp = (unsigned int*)qaws_internal_zalloc((size_t)n + 2, sizeof(unsigned int));
	cp = (unsigned int*)qaws_internal_zalloc((size_t)n + 2, sizeof(unsigned int));
	parent = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * (n + 1));
	flag = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * (n + 1));
	lnz = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * (n + 1));
	kids = (unsigned int*)qaws_internal_zalloc((size_t)n + 1, sizeof(unsigned int));
	wk.super = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * (n + 1));
	wk.map = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * (n + 1));
	if (!f->perm || !f->d || !f->w || !f->sf || !pinv || !rp || !cp || !parent || !flag || !lnz || !kids || !wk.super || !wk.map)
		goto done;
	for (k = 0; k < n; k++)
		f->perm[k] = k;
	if (xyz && n > 64)
	{
		ld_order o;
		o.rp = row_ptr;
		o.col = col;
		o.xyz = xyz;
		o.gen = 0;
		o.stamp = (unsigned int*)qaws_internal_zalloc((size_t)n + 1, sizeof(unsigned int));
		o.tmp = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * (n + 1));
		o.key = (ld_key*)qaws_internal_malloc(sizeof(ld_key) * (n + 1));
		if (o.stamp && o.tmp && o.key)
			ld_dissect(&o, f->perm, n);
		qaws_internal_free(o.stamp); qaws_internal_free(o.tmp); qaws_internal_free(o.key);
		if (o.gen == 0)
			goto done;
	}
	for (k = 0; k < n; k++)
		pinv[f->perm[k]] = k;
	/* C = P A P^T's lower triangle from A's, by rows (rp, rj: the elimination tree's input) and by columns (cp, ci, cx) */
	for (r = 0; r < n; r++)
		for (p = row_ptr[r]; p < row_ptr[r + 1]; p++)
			if (col[p] <= r)
			{
				unsigned int a = pinv[r], b = pinv[col[p]];
				rp[(a > b ? a : b) + 2]++;
				cp[(a < b ? a : b) + 2]++;
			}
	for (k = 0; k < n; k++)
	{
		rp[k + 2] += rp[k + 1];
		cp[k + 2] += cp[k + 1];
	}
	nnzc = cp[n + 1];
	rj = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * (nnzc + 1));
	ci = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * (nnzc + 1));
	cx = (double*)qaws_internal_malloc(sizeof(double) * (nnzc + 1));
	if (!rj || !ci || !cx)
		goto done;
	for (r = 0; r < n; r++)
		for (p = row_ptr[r]; p < row_ptr[r + 1]; p++)
			if (col[p] <= r)
			{
				unsigned int a = pinv[r], b = pinv[col[p]], lo = a < b ? a : b, hi = a < b ? b : a, at = cp[lo + 1]++;
				rj[rp[hi + 1]++] = lo;
				ci[at] = hi;
				cx[at] = val[p];
				if (a == b)
					dmax = fabs(val[p]) > dmax ? fabs(val[p]) : dmax;
			}
	/* the elimination tree and the column counts (marks per row, each path walked once) */
	for (k = 0; k < n; k++)
	{
		parent[k] = LD_NONE;
		flag[k] = k;
		lnz[k] = 0;
		for (p = rp[k]; p < rp[k + 1]; p++)
			for (i = rj[p]; flag[i] != k; i = parent[i])
			{
				if (parent[i] == LD_NONE)
					parent[i] = k;
				lnz[i]++;
				flag[i] = k;
			}
	}
	/* fundamental supernodes: a column joins the one before when it is that column's only child, one row shorter */
	for (k = 0; k < n; k++)
		if (parent[k] != LD_NONE)
			kids[parent[k]]++;
	for (k = 0; k < n; k++)
	{
		if (k == 0 || parent[k - 1] != k || kids[k] != 1 || lnz[k] + 1 != lnz[k - 1])
			f->sf[f->ns++] = k;
		wk.super[k] = f->ns - 1;
	}
	f->sf[f->ns] = n;
	f->rp = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * (f->ns + 1));
	f->lp = (size_t*)qaws_internal_malloc(sizeof(size_t) * (f->ns + 1));
	wk.head = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * (f->ns + 1));
	wk.next = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * (f->ns + 1));
	wk.pos = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * (f->ns + 1));
	if (!f->rp || !f->lp || !wk.head || !wk.next || !wk.pos)
		goto done;
	/* each supernode's rows: its first column and that column's rows below */
	f->rp[0] = 0;
	f->lp[0] = 0;
	for (s = 0; s < f->ns; s++)
	{
		unsigned int w = f->sf[s + 1] - f->sf[s], nr = lnz[f->sf[s]] + 1;
		if (f->rp[s] + nr < f->rp[s])
			goto done;   /* (more than 4G rows) */
		f->rp[s + 1] = f->rp[s] + nr;
		f->lp[s + 1] = f->lp[s] + (size_t)nr * w;
		f->nnz += (double)nr * w - (double)w * (w + 1) / 2;
		maxr = nr > maxr ? nr : maxr;
		maxw = w > maxw ? w : maxw;
	}
	total = f->lp[f->ns];
	f->ri = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)f->rp[f->ns] + 1));
	f->lx = (double*)qaws_internal_malloc(sizeof(double) * (total + 1));
	wk.acc = (double*)qaws_internal_malloc(sizeof(double) * (maxr + 1));
	wk.t = (double*)qaws_internal_malloc(sizeof(double) * (maxw + 1));
	if (!f->ri || !f->lx || !wk.acc || !wk.t)
		goto done;
	/* (gathered row by row along each row's reach, so each list comes sorted) */
	for (s = 0; s < f->ns; s++)
	{
		wk.pos[s] = f->rp[s] + 1;
		f->ri[f->rp[s]] = f->sf[s];
	}
	for (k = 0; k < n; k++)
	{
		flag[k] = k;
		for (p = rp[k]; p < rp[k + 1]; p++)
			for (i = rj[p]; flag[i] != k; i = parent[i])
			{
				flag[i] = k;
				if (f->sf[wk.super[i]] == i)
					f->ri[wk.pos[wk.super[i]]++] = k;
			}
	}
	st = QAWS_STATUS_INVALID_ARGUMENT;
	for (tries = 0; tries < 4; tries++)
	{
		if (ld_numeric(f, cp, ci, cx, shift, &wk))
		{
			st = QAWS_STATUS_OK;
			break;
		}
		shift = shift > 0 ? 10 * shift : 1e-12 * (dmax > 0 ? dmax : 1);
	}
done:
	qaws_internal_free(pinv); qaws_internal_free(rp); qaws_internal_free(rj); qaws_internal_free(cp); qaws_internal_free(ci); qaws_internal_free(cx); qaws_internal_free(parent); qaws_internal_free(flag); qaws_internal_free(lnz); qaws_internal_free(kids);
	qaws_internal_free(wk.super); qaws_internal_free(wk.map); qaws_internal_free(wk.head); qaws_internal_free(wk.next); qaws_internal_free(wk.pos); qaws_internal_free(wk.acc); qaws_internal_free(wk.t);
	if (st != QAWS_STATUS_OK)
		qaws_internal_ldlt_free(f);
	return st;
}

void qaws_internal_ldlt_solve(qaws_internal_ldlt const* f, double const* b, double* x)
{
	unsigned int n = f->n, j, s, c, i;
	double *y = f->w, *t = f->w + n;
	for (j = 0; j < n; j++)
		y[j] = b[f->perm[j]];
	/* per supernode: its own columns' triangle, then the rows below gathered into t (one scatter per row) */
	for (s = 0; s < f->ns; s++)
	{
		unsigned int c0 = f->sf[s], w = f->sf[s + 1] - c0, nr = f->rp[s + 1] - f->rp[s], *R = f->ri + f->rp[s];
		double const* P = f->lx + f->lp[s];
		for (i = w; i < nr; i++)
			t[i] = 0;
		for (c = 0; c < w; c++)
		{
			double const* Pc = P + (size_t)c * nr;
			double yc = y[c0 + c];
			for (i = c + 1; i < w; i++)
				y[c0 + i] -= Pc[i] * yc;
			for (i = w; i < nr; i++)
				t[i] += Pc[i] * yc;
		}
		for (i = w; i < nr; i++)
			y[R[i]] -= t[i];
	}
	for (j = 0; j < n; j++)
		y[j] /= f->d[j];
	for (s = f->ns; s > 0; s--)
	{
		unsigned int c0 = f->sf[s - 1], w = f->sf[s] - c0, nr = f->rp[s] - f->rp[s - 1], *R = f->ri + f->rp[s - 1];
		double const* P = f->lx + f->lp[s - 1];
		for (i = w; i < nr; i++)
			t[i] = y[R[i]];
		for (c = w; c > 0; c--)
		{
			double const* Pc = P + (size_t)(c - 1) * nr;
			double sum = 0;
			for (i = c; i < w; i++)
				sum += Pc[i] * y[c0 + i];
			for (i = w; i < nr; i++)
				sum += Pc[i] * t[i];
			y[c0 + c - 1] -= sum;
		}
	}
	for (j = 0; j < n; j++)
		x[f->perm[j]] = y[j];
}
