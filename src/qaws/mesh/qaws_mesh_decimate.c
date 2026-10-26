/*
 * Decimation (qaws_mesh_decimate.h): quadrics per vertex, a heap of edges
 * by their collapse's error (stale entries skipped by the vertices'
 * stamps), collapses on the corner store without compacting until the end.
 */

#include "../internal/qaws_internal_memory.h"
#include "../qaws_mesh_decimate.h"
#include "../qaws_mesh_reference.h"
#include "qaws_mesh_internal.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* a quadric over the augmented point (x, y, z, the attributes scaled, 1) of n1 coordinates: the symmetric n1 x n1 sum of
   v v^T, packed by rows above the diagonal (without attributes the 10 entries a2 ab ac ad b2 bc bd c2 cd d2); entry (i, j),
   i <= j */
static size_t dc_at(unsigned int n1, unsigned int i, unsigned int j)
{
	return (size_t)i * n1 - (size_t)i * (i ? i - 1 : 0) / 2 + (j - i);
}

/* the plane n . x + d = 0, its squared distance weighted by w */
static void dc_plane(double* Q, unsigned int n1, double const* n, double d, double w)
{
	double p[4] = { n[0], n[1], n[2], d };
	unsigned int at[4] = { 0, 1, 2 }, i, j;
	at[3] = n1 - 1;
	for (i = 0; i < 4; i++)
		for (j = i; j < 4; j++)
			Q[dc_at(n1, at[i], at[j])] += w * p[i] * p[j];
}

/* Q's value at the augmented point z (n1 coordinates, the last 1) */
static double dc_eval_n(double const* Q, double const* z, unsigned int n1)
{
	unsigned int i, j;
	double e = 0;
	for (i = 0; i < n1; i++)
		for (j = i; j < n1; j++)
			e += (i == j ? 1 : 2) * Q[dc_at(n1, i, j)] * z[i] * z[j];
	return e;
}

/* M r' = r for n unknowns (M n x n by rows, overwritten; r the solution), partial pivoting; 0 when near singular */
static int dc_solve(double* M, double* r, unsigned int n)
{
	unsigned int i, j, k;
	double big = 0;
	for (i = 0; i < n * n; i++) big = fmax(big, fabs(M[i]));
	for (k = 0; k < n; k++)
	{
		unsigned int piv = k;
		for (i = k + 1; i < n; i++)
			if (fabs(M[i * n + k]) > fabs(M[piv * n + k]))
				piv = i;
		if (!(fabs(M[piv * n + k]) > 1e-12 * big))
			return 0;
		if (piv != k)
		{
			double t;
			for (j = 0; j < n; j++)
			{
				t = M[k * n + j];
				M[k * n + j] = M[piv * n + j];
				M[piv * n + j] = t;
			}
			t = r[k];
			r[k] = r[piv];
			r[piv] = t;
		}
		for (i = k + 1; i < n; i++)
		{
			double f = M[i * n + k] / M[k * n + k];
			for (j = k; j < n; j++) M[i * n + j] -= f * M[k * n + j];
			r[i] -= f * r[k];
		}
	}
	for (k = n; k-- > 0;)
	{
		for (j = k + 1; j < n; j++) r[k] -= M[k * n + j] * r[j];
		r[k] /= M[k * n + k];
	}
	return 1;
}

static double dc_eval(double const* Q, double const* x)
{
	return Q[0] * x[0] * x[0] + 2 * Q[1] * x[0] * x[1] + 2 * Q[2] * x[0] * x[2] + 2 * Q[3] * x[0] + Q[4] * x[1] * x[1] + 2 * Q[5] * x[1] * x[2]
		+ 2 * Q[6] * x[1] + Q[7] * x[2] * x[2] + 2 * Q[8] * x[2] + Q[9];
}

/* the point where Q is least, by Cramer's rule; 0 when Q's 3 x 3 part is near singular (flat, or along a line) */
static int dc_optimal(double const* Q, double* x)
{
	double a = Q[0], b = Q[1], c = Q[2], e = Q[4], f = Q[5], i = Q[7], r[3] = { -Q[3], -Q[6], -Q[8] };
	double det = a * (e * i - f * f) - b * (b * i - f * c) + c * (b * f - e * c), s = fabs(a) + fabs(e) + fabs(i);
	if (!(fabs(det) > 1e-10 * s * s * s))
		return 0;
	x[0] = (r[0] * (e * i - f * f) - b * (r[1] * i - f * r[2]) + c * (r[1] * f - e * r[2])) / det;
	x[1] = (a * (r[1] * i - f * r[2]) - r[0] * (b * i - f * c) + c * (b * r[2] - r[1] * c)) / det;
	x[2] = (a * (e * r[2] - r[1] * f) - b * (b * r[2] - r[1] * c) + r[0] * (b * f - e * c)) / det;
	return 1;
}

typedef struct dc_entry
{
	double cost;
	unsigned int a, b, sa, sb;
} dc_entry;

typedef struct dc_heap
{
	dc_entry* e;
	unsigned int n, cap;
} dc_heap;

static int dc_push(dc_heap* H, dc_entry const* x)
{
	unsigned int i;
	if (H->n == H->cap)
	{
		unsigned int cap = H->cap ? 2 * H->cap : 1024;
		dc_entry* g = (dc_entry*)qaws_internal_realloc(H->e, qaws_internal_mul2(sizeof(dc_entry), cap));
		if (!g)
			return 0;
		H->e = g;
		H->cap = cap;
	}
	for (i = H->n++; i && H->e[(i - 1) / 2].cost > x->cost; i = (i - 1) / 2)
		H->e[i] = H->e[(i - 1) / 2];
	H->e[i] = *x;
	return 1;
}

static void dc_pop(dc_heap* H, dc_entry* out)
{
	dc_entry last;
	unsigned int i = 0;
	*out = H->e[0];
	last = H->e[--H->n];
	for (;;)
	{
		unsigned int c = 2 * i + 1;
		if (c >= H->n)
			break;
		if (c + 1 < H->n && H->e[c + 1].cost < H->e[c].cost)
			c++;
		if (H->e[c].cost >= last.cost)
			break;
		H->e[i] = H->e[c];
		i = c;
	}
	if (H->n)
		H->e[i] = last;
}

typedef struct dc_work
{
	qaws_hm* m;
	unsigned int n1, K;             /* the augmented coordinates (4 + the attributes' components), a quadric's entries */
	unsigned int *ca, *ck;          /* per attribute component in the quadric: its attribute and component */
	double alpha;                   /* the attributes' scale: a unit of them as this distance */
	double *q, *M, *r, *z;          /* scratch: a quadric, a system, its right side, an augmented point */
	double *Q, *W;
	unsigned int* stamp;
	unsigned char* border;
	unsigned char* locked;          /* neither removed nor moved: a locked border, an importance of 0 or less */
	double cosmax;                  /* the cosine of the largest turn of a face a collapse may make */
	int memless;                    /* QAWS_DECIMATE_MEMORYLESS: costs from the current faces, nothing stored */
	double fw;                      /* the border's and features' weight */
} dc_work;

/* memoryless (Lindstrom and Turk 1998): a and b to one point x placed from the faces round them as they are now. The
   volume under them (each face fanned, its tetrahedra to a) is linear in x: kept, one constraint. Each face moving with
   a or b sweeps a volume n . (x - o) (n its fan triangle's cross product, o the vertex leaving): their squares summed are
   the cost (read as a distance: over the sum of |n|^2), least under the constraint (with the border's and features'
   planes, as now). Lindstrom and Turk's border term: where the collapse touches the border, the area vector of its chain
   there (the sum of v_i x v_i+1 over its edges) kept, the new chain's x x d (d: the far ends of the edges leaving a or b
   less the near ends of those arriving), two more constraints where independent of the volume's; a hair of pull to the
   edge's middle for the directions nothing fixes (their shape term, the ring's squared edges, tried: on an icosphere it
   lowered the smallest angle, in the cost as in the placement). fixed: x is b. */
static void dc_memless(dc_work const* w, unsigned int a, unsigned int b, int fixed, double* x, double* cost)
{
	qaws_hm const* m = w->m;
	double A[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 }, B[3] = { 0, 0, 0 }, C = 0, g[3] = { 0, 0, 0 }, vold = 0, vconst = 0, nn = 0, mid[3], tr, eps, y[3];
	double bd[3] = { 0, 0, 0 }, be[3] = { 0, 0, 0 };   /* the border chain: x x bd = be */
	double const* o = m->x + 3 * (size_t)a;
	unsigned int side, c, k, j, nbd = 0;
	/* the border chain's edges round a and b (relative to a) */
	for (side = 0; side < 2; side++)
	{
		unsigned int v = side ? b : a;
		for (c = m->vc[v]; c != HM_NONE; c = hm_round(m, v, c))
		{
			unsigned int u = hm_to(m, c), pc = hm_prev(m, c), pu = m->cv[pc];
			if (m->co[c] == HM_NONE)
			{
				/* a border edge v -> u */
				double const *pv = m->x + 3 * (size_t)v, *pw = m->x + 3 * (size_t)u;
				double s[3], t[3];
				for (k = 0; k < 3; k++) { s[k] = pv[k] - o[k]; t[k] = pw[k] - o[k]; }
				be[0] += s[1] * t[2] - s[2] * t[1];
				be[1] += s[2] * t[0] - s[0] * t[2];
				be[2] += s[0] * t[1] - s[1] * t[0];
				if (u != a && u != b)
					for (k = 0; k < 3; k++) bd[k] += t[k];
				nbd++;
			}
			if (m->co[pc] == HM_NONE && pu != a && pu != b)
			{
				/* a border edge pu -> v arriving (counted here once: its start is not a or b) */
				double const *pv = m->x + 3 * (size_t)v, *pw = m->x + 3 * (size_t)pu;
				double s[3], t[3];
				for (k = 0; k < 3; k++) { s[k] = pw[k] - o[k]; t[k] = pv[k] - o[k]; }
				be[0] += s[1] * t[2] - s[2] * t[1];
				be[1] += s[2] * t[0] - s[0] * t[2];
				be[2] += s[0] * t[1] - s[1] * t[0];
				for (k = 0; k < 3; k++) bd[k] -= s[k];
				nbd++;
			}
		}
	}
	for (k = 0; k < 3; k++) mid[k] = 0.5 * (m->x[3 * (size_t)b + k] - o[k]);
	for (side = 0; side < 2; side++)
	{
		unsigned int v = side ? b : a;
		for (c = m->vc[v]; c != HM_NONE; c = hm_round(m, v, c))
		{
			unsigned int f = m->cf[c], f0 = m->fc[f], fn = m->fn[f], hasa = 0, hasb = 0, q;
			double P[64][3];
			int mv[64];
			for (q = 0; q < fn; q++)
			{
				hasa |= m->cv[f0 + q] == a;
				hasb |= m->cv[f0 + q] == b;
			}
			if ((side && hasa) || fn > 64)
				continue;   /* (a face round both: counted from a) */
			/* the face relative to a, its moving corners marked */
			for (q = 0; q < fn; q++)
			{
				unsigned int u = m->cv[f0 + q];
				for (k = 0; k < 3; k++) P[q][k] = m->x[3 * (size_t)u + k] - o[k];
				mv[q] = u == a || u == b;
			}
			/* its volume before (fan from its first corner, apex a), and after as constant + g . x */
			for (q = 1; q + 1 < fn; q++)
			{
				double const *p0 = P[0], *p1 = P[q], *p2 = P[q + 1];
				unsigned int nmv = (unsigned int)(mv[0] + mv[q] + mv[q + 1]);
				vold += p0[0] * (p1[1] * p2[2] - p1[2] * p2[1]) - p0[1] * (p1[0] * p2[2] - p1[2] * p2[0]) + p0[2] * (p1[0] * p2[1] - p1[1] * p2[0]);
				if (nmv == 0)
					vconst += p0[0] * (p1[1] * p2[2] - p1[2] * p2[1]) - p0[1] * (p1[0] * p2[2] - p1[2] * p2[0]) + p0[2] * (p1[0] * p2[1] - p1[1] * p2[0]);
				else if (nmv == 1)
				{
					/* det(x, p, q) with the others in turn: x . (p x q) */
					double const *pp = mv[0] ? p1 : mv[q] ? p2 : p0, *qq = mv[0] ? p2 : mv[q] ? p0 : p1;
					g[0] += pp[1] * qq[2] - pp[2] * qq[1];
					g[1] += pp[2] * qq[0] - pp[0] * qq[2];
					g[2] += pp[0] * qq[1] - pp[1] * qq[0];
				}
			}
			/* the swept volumes: a face round one of them only, fanned from the moving corner */
			if (hasa != hasb)
			{
				unsigned int s0 = 0;
				while (!mv[s0]) s0++;
				for (q = 1; q + 1 < fn; q++)
				{
					double const *ov = P[s0], *p1 = P[(s0 + q) % fn], *p2 = P[(s0 + q + 1) % fn];
					double e1[3], e2[3], n[3], d;
					for (k = 0; k < 3; k++) { e1[k] = p1[k] - ov[k]; e2[k] = p2[k] - ov[k]; }
					n[0] = e1[1] * e2[2] - e1[2] * e2[1];
					n[1] = e1[2] * e2[0] - e1[0] * e2[2];
					n[2] = e1[0] * e2[1] - e1[1] * e2[0];
					d = n[0] * ov[0] + n[1] * ov[1] + n[2] * ov[2];
					for (k = 0; k < 3; k++)
					{
						for (j = 0; j < 3; j++) A[3 * k + j] += n[k] * n[j];
						B[k] += n[k] * d;
					}
					C += d * d;
					nn += n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
				}
			}
			/* the border's and features' planes, standing on their edges (weights in the same units: length^4) */
			for (q = 0; q < fn; q++)
			{
				unsigned int cq = f0 + q;
				if (!(mv[q] || mv[(q + 1) % fn]) || !(m->co[cq] == HM_NONE || (m->cflag[cq] & HM_FEATURE)))
					continue;
				{
					double const *pa = P[q], *pb = P[(q + 1) % fn];
					double u[3], fnrm[3] = { 0, 0, 0 }, s[3], l2, ls, d;
					unsigned int r;
					for (r = 1; r + 1 < fn; r++)
					{
						/* (the face's normal) */
						double e1[3], e2[3];
						for (k = 0; k < 3; k++) { e1[k] = P[r][k] - P[0][k]; e2[k] = P[r + 1][k] - P[0][k]; }
						fnrm[0] += e1[1] * e2[2] - e1[2] * e2[1];
						fnrm[1] += e1[2] * e2[0] - e1[0] * e2[2];
						fnrm[2] += e1[0] * e2[1] - e1[1] * e2[0];
					}
					for (k = 0; k < 3; k++) u[k] = pb[k] - pa[k];
					l2 = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
					s[0] = u[1] * fnrm[2] - u[2] * fnrm[1];
					s[1] = u[2] * fnrm[0] - u[0] * fnrm[2];
					s[2] = u[0] * fnrm[1] - u[1] * fnrm[0];
					ls = sqrt(s[0] * s[0] + s[1] * s[1] + s[2] * s[2]);
					if (!(ls > 0))
						continue;
					for (k = 0; k < 3; k++) s[k] /= ls;
					d = s[0] * pa[0] + s[1] * pa[1] + s[2] * pa[2];
					for (k = 0; k < 3; k++)
					{
						for (j = 0; j < 3; j++) A[3 * k + j] += w->fw * l2 * l2 * s[k] * s[j];
						B[k] += w->fw * l2 * l2 * s[k] * d;
					}
					C += w->fw * l2 * l2 * d * d;
				}
			}
		}
	}
	tr = A[0] + A[4] + A[8];
	eps = 1e-9 * (tr > 0 ? tr : 1);
	if (fixed)
		for (k = 0; k < 3; k++) y[k] = m->x[3 * (size_t)b + k] - o[k];
	else
	{
		/* the least of x^T A x - 2 B . x + C (plus eps |x - mid|^2) under the constraints: the volume's
		   g . x = vold - vconst, the border's (d x u_i) . x = u_i . be for u_i across d, each kept when independent of those
		   before (Gram-Schmidt) */
		double rows[3][4], M6[36], r6[6], r4[4], gl = sqrt(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]);
		unsigned int nr = 0, n6, z;
		int ok;
		if (gl > 0)
		{
			for (k = 0; k < 3; k++) rows[0][k] = g[k];
			rows[0][3] = vold - vconst;
			nr = 1;
		}
		if (nbd)
		{
			double dl = sqrt(bd[0] * bd[0] + bd[1] * bd[1] + bd[2] * bd[2]), u[2][3];
			if (dl > 0)
			{
				/* u_0, u_1 across d */
				double ax[3] = { 0, 0, 0 }, l;
				unsigned int im = fabs(bd[0]) <= fabs(bd[1]) && fabs(bd[0]) <= fabs(bd[2]) ? 0 : (fabs(bd[1]) <= fabs(bd[2]) ? 1 : 2);
				ax[im] = 1;
				u[0][0] = bd[1] * ax[2] - bd[2] * ax[1];
				u[0][1] = bd[2] * ax[0] - bd[0] * ax[2];
				u[0][2] = bd[0] * ax[1] - bd[1] * ax[0];
				l = sqrt(u[0][0] * u[0][0] + u[0][1] * u[0][1] + u[0][2] * u[0][2]);
				for (k = 0; k < 3; k++) u[0][k] /= l;
				u[1][0] = (bd[1] * u[0][2] - bd[2] * u[0][1]) / dl;
				u[1][1] = (bd[2] * u[0][0] - bd[0] * u[0][2]) / dl;
				u[1][2] = (bd[0] * u[0][1] - bd[1] * u[0][0]) / dl;
				for (z = 0; z < 2; z++)
				{
					double rr[4], ln;
					unsigned int q;
					/* (d x u_z) . x = u_z . be */
					rr[0] = bd[1] * u[z][2] - bd[2] * u[z][1];
					rr[1] = bd[2] * u[z][0] - bd[0] * u[z][2];
					rr[2] = bd[0] * u[z][1] - bd[1] * u[z][0];
					rr[3] = u[z][0] * be[0] + u[z][1] * be[1] + u[z][2] * be[2];
					ln = sqrt(rr[0] * rr[0] + rr[1] * rr[1] + rr[2] * rr[2]);
					for (q = 0; q < nr; q++)
					{
						double dq = 0, lq = 0, f;
						for (k = 0; k < 3; k++) { dq += rr[k] * rows[q][k]; lq += rows[q][k] * rows[q][k]; }
						f = dq / lq;
						for (k = 0; k < 4; k++) rr[k] -= f * rows[q][k];
					}
					if (sqrt(rr[0] * rr[0] + rr[1] * rr[1] + rr[2] * rr[2]) > 1e-6 * ln)
					{
						memcpy(rows[nr], rr, sizeof rr);
						nr++;
					}
				}
			}
		}
		n6 = 3 + nr;
		memset(M6, 0, sizeof M6);
		for (k = 0; k < 3; k++)
		{
			for (j = 0; j < 3; j++) M6[n6 * k + j] = A[3 * k + j] + (k == j ? eps : 0);
			r6[k] = B[k] + eps * mid[k];
		}
		for (z = 0; z < nr; z++)
		{
			for (k = 0; k < 3; k++)
			{
				M6[n6 * k + 3 + z] = rows[z][k];
				M6[n6 * (3 + z) + k] = rows[z][k];
			}
			r6[3 + z] = rows[z][3];
		}
		ok = nr > 0 && dc_solve(M6, r6, n6);
		for (k = 0; ok && k < 3; k++) r4[k] = r6[k];
		if (!ok)
		{
			double M3[9];
			ok = 0;
			for (k = 0; k < 9; k++) M3[k] = A[k] + (k % 4 == 0 ? eps : 0);
			for (k = 0; k < 3; k++) r4[k] = B[k] + eps * mid[k];
			ok = dc_solve(M3, r4, 3);
		}
		for (k = 0; k < 3; k++) y[k] = ok ? r4[k] : mid[k];
	}
	for (k = 0; k < 3; k++) x[k] = y[k] + o[k];
	{
		double e = C;
		for (k = 0; k < 3; k++)
		{
			e -= 2 * B[k] * y[k];
			for (j = 0; j < 3; j++) e += A[3 * k + j] * y[k] * y[j];
		}
		*cost = (e > 0 ? e : 0) / (nn > 0 ? nn : 1);
	}
}

/* the attributes (scaled) least for Q at position z[0..2], into z[3..n1 - 2]; where Q gives none, a's and b's mean */
static void dc_attributes(dc_work const* w, double const* Q, double* z, unsigned int a, unsigned int b)
{
	unsigned int n1 = w->n1, na = n1 - 4, i, j;
	for (i = 0; i < na; i++)
	{
		for (j = 0; j < na; j++) w->M[i * na + j] = Q[dc_at(n1, 3 + (i < j ? i : j), 3 + (i < j ? j : i))];
		w->r[i] = -Q[dc_at(n1, 3 + i, n1 - 1)];
		for (j = 0; j < 3; j++) w->r[i] -= Q[dc_at(n1, j, 3 + i)] * z[j];
	}
	if (dc_solve(w->M, w->r, na))
		memcpy(z + 3, w->r, sizeof(double) * na);
	else
		for (i = 0; i < na; i++)
			z[3 + i] = 0.5 * w->alpha * (hm_attr_at(w->m, w->ca[i], a)[w->ck[i]] + hm_attr_at(w->m, w->ca[i], b)[w->ck[i]]);
}

/* a merged into b: allowed, and where b goes (and its attributes, scaled, into s when given) and the error there; 0 when
   refused */
static int dc_candidate(dc_work const* w, unsigned int a, unsigned int b, double* x, double* s, double* cost)
{
	qaws_hm const* m = w->m;
	unsigned int c = hm_find_edge(m, a, b), k, n1 = w->n1;
	double *Q = w->q, best = HUGE_VAL, tries[4][3], wsum = w->W[a] + w->W[b];
	int fixed, nt = 0, i;
	if (c == HM_NONE || (m->vflag[a] & HM_CORNER) || (hm_feature_degree(m, a) && !(m->cflag[c] & HM_FEATURE)) || w->locked[a]
		|| !hm_can_collapse(m, c, a))
		return 0;
	for (k = 0; k < w->K; k++)
		Q[k] = w->Q[w->K * (size_t)a + k] + w->Q[w->K * (size_t)b + k];
	/* b stays where it is when it may not leave its place: a corner, a locked border vertex, a feature vertex the edge leaves */
	fixed = (m->vflag[b] & HM_CORNER) || w->locked[b] || (hm_feature_degree(m, b) && !(m->cflag[c] & HM_FEATURE));
	if (w->memless)
	{
		dc_memless(w, a, b, fixed, x, cost);
		return 1;
	}
	memcpy(tries[nt++], m->x + 3 * (size_t)b, 3 * sizeof(double));
	if (!fixed)
	{
		memcpy(tries[nt++], m->x + 3 * (size_t)a, 3 * sizeof(double));
		for (k = 0; k < 3; k++)
			tries[nt][k] = 0.5 * (m->x[3 * (size_t)a + k] + m->x[3 * (size_t)b + k]);
		nt++;
		if (n1 == 4 ? dc_optimal(Q, tries[nt]) : 0)
			nt++;
		else if (n1 > 4)
		{
			/* (with attributes: the position and the attributes least together) */
			unsigned int D = n1 - 1, j;
			for (k = 0; k < D; k++)
			{
				for (j = 0; j < D; j++) w->M[k * D + j] = Q[dc_at(n1, k < j ? k : j, k < j ? j : k)];
				w->r[k] = -Q[dc_at(n1, k, D)];
			}
			if (dc_solve(w->M, w->r, D))
				memcpy(tries[nt++], w->r, 3 * sizeof(double));
		}
	}
	for (i = 0; i < nt; i++)
	{
		double e;
		if (n1 == 4)
			e = dc_eval(Q, tries[i]);
		else
		{
			memcpy(w->z, tries[i], 3 * sizeof(double));
			dc_attributes(w, Q, w->z, a, b);
			w->z[n1 - 1] = 1;
			e = dc_eval_n(Q, w->z, n1);
		}
		if (e < best)
		{
			best = e;
			memcpy(x, tries[i], 3 * sizeof(double));
			if (s && n1 > 4)
				memcpy(s, w->z + 3, sizeof(double) * (n1 - 4));
		}
	}
	*cost = (best > 0 ? best : 0) / (wsum > 0 ? wsum : 1);
	return 1;
}

/* the point of triangle (p0, p1, p2) nearest x: its weights into w, the squared distance returned (Ericson's regions) */
static double dc_closest(double const* x, double const* p0, double const* p1, double const* p2, double* w)
{
	double ab[3], ac[3], ap[3], bp[3], cp[3], d1 = 0, d2 = 0, d3 = 0, d4 = 0, d5 = 0, d6 = 0, va, vb, vc, q[3], r = 0;
	unsigned int k;
	for (k = 0; k < 3; k++)
	{
		ab[k] = p1[k] - p0[k]; ac[k] = p2[k] - p0[k]; ap[k] = x[k] - p0[k]; bp[k] = x[k] - p1[k]; cp[k] = x[k] - p2[k];
		d1 += ab[k] * ap[k]; d2 += ac[k] * ap[k]; d3 += ab[k] * bp[k]; d4 += ac[k] * bp[k]; d5 += ab[k] * cp[k]; d6 += ac[k] * cp[k];
	}
	vc = d1 * d4 - d3 * d2;
	vb = d5 * d2 - d1 * d6;
	va = d3 * d6 - d5 * d4;
	if (d1 <= 0 && d2 <= 0) { w[0] = 1; w[1] = w[2] = 0; }
	else if (d3 >= 0 && d4 <= d3) { w[1] = 1; w[0] = w[2] = 0; }
	else if (d6 >= 0 && d5 <= d6) { w[2] = 1; w[0] = w[1] = 0; }
	else if (vc <= 0 && d1 >= 0 && d3 <= 0) { w[1] = d1 / (d1 - d3); w[0] = 1 - w[1]; w[2] = 0; }
	else if (vb <= 0 && d2 >= 0 && d6 <= 0) { w[2] = d2 / (d2 - d6); w[0] = 1 - w[2]; w[1] = 0; }
	else if (va <= 0 && d4 - d3 >= 0 && d5 - d6 >= 0) { w[2] = (d4 - d3) / ((d4 - d3) + (d5 - d6)); w[1] = 1 - w[2]; w[0] = 0; }
	else
	{
		double den = 1 / (va + vb + vc);
		w[1] = vb * den; w[2] = vc * den; w[0] = 1 - w[1] - w[2];
	}
	for (k = 0; k < 3; k++)
	{
		q[k] = w[0] * p0[k] + w[1] * p1[k] + w[2] * p2[k] - x[k];
		r += q[k] * q[k];
	}
	return r;
}

/* Corner attributes through a collapse of a into b at x (Hoppe's wedges): per corner attribute (of 4 components at most), the
   corners at a and b grouped (equal values at one vertex; across the edge's faces, the corner at a with the one at b), each
   group's value the one at x of its corners' triangles before the collapse, interpolated on the nearest (LINEAR, VECTOR
   and the unit rules; the others take its nearest corner's). Into tab: per
   attribute, pairs (the old value, the new one) for each corner, matched by value after the collapse; *nt their counts.
   shortcut: 64 corners round the edge at most (more: the attribute left as the collapse copies it), upgrade by growing the
   scratch */
#define DC_WEDGE_MAX 64

/* two corner values one wedge's: equal to round-off (tol, relative), q and -q one orientation */
static int dc_wedge_same(double const* p, double const* q, unsigned int n, unsigned int rule, double tol)
{
	unsigned int k;
	int same = 1, opp = rule == QAWS_ATTRIBUTE_QUATERNION;
	for (k = 0; k < n; k++)
	{
		same = same && fabs(p[k] - q[k]) <= tol * (1 + fabs(p[k]));
		opp = opp && fabs(p[k] + q[k]) <= tol * (1 + fabs(p[k]));
	}
	return same || opp;
}

static void dc_wedges_before(qaws_hm const* m, unsigned int a, unsigned int b, double const* x, double* tab, unsigned int* nt, unsigned int stride)
{
	unsigned int ak, off = 0;
	for (ak = 0; ak < m->na; ak++)
	{
		unsigned int n = m->attr[ak].n, cs[DC_WEDGE_MAX], grp[DC_WEDGE_MAX], nc = 0, i, j, v, c;
		double best[DC_WEDGE_MAX], val[DC_WEDGE_MAX * 4];
		nt[ak] = 0;
		int blend = m->attr[ak].rule == QAWS_ATTRIBUTE_LINEAR || m->attr[ak].rule == QAWS_ATTRIBUTE_VECTOR || hm_rule_unit(m->attr[ak].rule);
		if (m->attr[ak].element != QAWS_MESH_CORNER || n > 4)
			continue;
		for (v = 0; v < 2; v++)
			for (c = m->vc[v ? b : a]; c != HM_NONE && nc <= DC_WEDGE_MAX; c = hm_round(m, v ? b : a, c))
				if (nc < DC_WEDGE_MAX)
					cs[nc++] = c;
				else
					nc++;
		if (nc > DC_WEDGE_MAX)
			continue;
		/* the groups: union by the least index */
		for (i = 0; i < nc; i++) grp[i] = i;
		for (i = 0; i < nc; i++)
			for (j = i + 1; j < nc; j++)
			{
				int join = 0;
				if (m->cv[cs[i]] == m->cv[cs[j]])
				{
					join = dc_wedge_same(hm_attr_at(m, ak, cs[i]), hm_attr_at(m, ak, cs[j]), n, m->attr[ak].rule, 1e-12);   /* (as dc_seams reads them) */
				}
				else
					join = m->cf[cs[i]] == m->cf[cs[j]];   /* (a face holding the edge: a and b on the same side of any seam) */
				if (join)
				{
					unsigned int gi = grp[i], gj = grp[j], z;
					while (grp[gi] != gi) gi = grp[gi];
					while (grp[gj] != gj) gj = grp[gj];
					if (gi != gj)
						for (z = 0; z < nc; z++)
							if (grp[z] == (gi > gj ? gi : gj)) grp[z] = gi < gj ? gi : gj;
					grp[gi > gj ? gi : gj] = gi < gj ? gi : gj;
				}
			}
		for (i = 0; i < nc; i++)
		{
			while (grp[grp[i]] != grp[i]) grp[i] = grp[grp[i]];
			best[i] = HUGE_VAL;
		}
		/* each group's value at x: its corners' triangles (the corner and its face's neighbours), the nearest */
		for (i = 0; i < nc; i++)
		{
			unsigned int g = grp[i], cn = hm_next(m, cs[i]), cp = hm_prev(m, cs[i]), k;
			double wt[3], d2 = dc_closest(x, m->x + 3 * (size_t)m->cv[cs[i]], m->x + 3 * (size_t)m->cv[cn], m->x + 3 * (size_t)m->cv[cp], wt);
			if (d2 < best[g])
			{
				best[g] = d2;
				if (blend)
					for (k = 0; k < n; k++)
						val[4 * g + k] = wt[0] * hm_attr_at(m, ak, cs[i])[k] + wt[1] * hm_attr_at(m, ak, cn)[k] + wt[2] * hm_attr_at(m, ak, cp)[k];
				else   /* (other rules: the nearest corner's value, exact at a vertex) */
					memcpy(val + 4 * g, hm_attr_at(m, ak, wt[0] >= wt[1] && wt[0] >= wt[2] ? cs[i] : wt[1] >= wt[2] ? cn : cp), sizeof(double) * n);
			}
		}
		if (hm_rule_unit(m->attr[ak].rule))
			for (i = 0; i < nc; i++)
				if (grp[i] == i)
					hm_space_out(val + 4 * i, n, m->attr[ak].rule);
		for (i = 0; i < nc; i++)
		{
			memcpy(tab + off + 2 * (size_t)n * i, hm_attr_at(m, ak, cs[i]), sizeof(double) * n);
			memcpy(tab + off + 2 * (size_t)n * i + n, val + 4 * grp[i], sizeof(double) * n);
		}
		nt[ak] = nc;
		off += stride;
	}
}

/* after the collapse: b's corners' values replaced by their groups' (matched by their old value) */
static void dc_wedges_after(qaws_hm* m, unsigned int b, double const* tab, unsigned int const* nt, unsigned int stride)
{
	unsigned int ak, off = 0, c, i;
	for (ak = 0; ak < m->na; ak++)
	{
		unsigned int n = m->attr[ak].n;
		if (!nt[ak])
			continue;
		for (c = m->vc[b]; c != HM_NONE; c = hm_round(m, b, c))
		{
			double* y = hm_attr_at(m, ak, c);
			for (i = 0; i < nt[ak]; i++)
				if (dc_wedge_same(y, tab + off + 2 * (size_t)n * i, n, m->attr[ak].rule, 1e-9))
				{
					memcpy(y, tab + off + 2 * (size_t)n * i + n, sizeof(double) * n);
					break;
				}
		}
		off += stride;
	}
}

/* no face round a or b (but the edge's) turned by more than the bound, or flattened, by a and b moving to x */
static int dc_no_fold(qaws_hm const* m, unsigned int a, unsigned int b, double const* x, double cosmax)
{
	unsigned int s, v, c;
	for (s = 0; s < 2; s++)
	{
		v = s ? b : a;
		for (c = m->vc[v]; c != HM_NONE; c = hm_round(m, v, c))
		{
			unsigned int f = m->cf[c], e = m->fc[f] + m->fn[f], q, hasa = 0, hasb = 0;
			double n0[3] = { 0, 0, 0 }, n1[3] = { 0, 0, 0 }, l0, l1;
			for (q = m->fc[f]; q < e; q++)
			{
				hasa |= m->cv[q] == a;
				hasb |= m->cv[q] == b;
			}
			if (hasa && hasb)
				continue;   /* (the edge's own faces go) */
			for (q = m->fc[f]; q < e; q++)
			{
				unsigned int r = q + 1 < e ? q + 1 : m->fc[f];
				double const *p = m->x + 3 * (size_t)m->cv[q], *pr = m->x + 3 * (size_t)m->cv[r];
				double const *P = m->cv[q] == a || m->cv[q] == b ? x : p, *PR = m->cv[r] == a || m->cv[r] == b ? x : pr;
				n0[0] += (p[1] - pr[1]) * (p[2] + pr[2]);
				n0[1] += (p[2] - pr[2]) * (p[0] + pr[0]);
				n0[2] += (p[0] - pr[0]) * (p[1] + pr[1]);
				n1[0] += (P[1] - PR[1]) * (P[2] + PR[2]);
				n1[1] += (P[2] - PR[2]) * (P[0] + PR[0]);
				n1[2] += (P[0] - PR[0]) * (P[1] + PR[1]);
			}
			l0 = sqrt(n0[0] * n0[0] + n0[1] * n0[1] + n0[2] * n0[2]);
			l1 = sqrt(n1[0] * n1[0] + n1[1] * n1[1] + n1[2] * n1[2]);
			if (!(l1 > 1e-6 * l0) || n0[0] * n1[0] + n0[1] * n1[1] + n0[2] * n1[2] < cosmax * l0 * l1)
				return 0;
		}
	}
	return 1;
}

static int dc_push_edges(dc_work const* w, dc_heap* H, unsigned int v, unsigned int** ring, unsigned int* cap)
{
	unsigned int n = hm_ring(w->m, v, *ring, *cap), i;
	if (n > *cap)
	{
		unsigned int* g = (unsigned int*)qaws_internal_realloc(*ring, sizeof(unsigned int) * (size_t)n);
		if (!g)
			return 0;
		*ring = g;
		*cap = n;
		hm_ring(w->m, v, *ring, *cap);
	}
	for (i = 0; i < n; i++)
	{
		unsigned int u = (*ring)[i], s;
		for (s = 0; s < 2; s++)
		{
			dc_entry e;
			double x[3];
			e.a = s ? u : v;
			e.b = s ? v : u;
			if (!dc_candidate(w, e.a, e.b, x, NULL, &e.cost))
				continue;
			e.sa = w->stamp[e.a];
			e.sb = w->stamp[e.b];
			if (!dc_push(H, &e))
				return 0;
		}
	}
	return 1;
}

/* the edges sharper than the angle made features, both corners; the corners recounted (1 or more than 2 feature edges);
   returns how many were made */
static unsigned int dc_creases(qaws_hm* m, double angle)
{
	unsigned int c, v, n = 0;
	double cosd = cos(angle);
	for (c = 0; c < m->nc; c++)
	{
		unsigned int o = m->co[c];
		double n0[3], n1[3], l0, l1;
		if (m->cf[c] == HM_NONE || o == HM_NONE || o < c || (m->cflag[c] & HM_FEATURE))
			continue;
		hm_face_normal(m, m->cf[c], n0);
		hm_face_normal(m, m->cf[o], n1);
		l0 = sqrt(n0[0] * n0[0] + n0[1] * n0[1] + n0[2] * n0[2]);
		l1 = sqrt(n1[0] * n1[0] + n1[1] * n1[1] + n1[2] * n1[2]);
		if (l0 > 0 && l1 > 0 && n0[0] * n1[0] + n0[1] * n1[1] + n0[2] * n1[2] < cosd * l0 * l1)
		{
			m->cflag[c] |= HM_FEATURE;
			m->cflag[o] |= HM_FEATURE;
			n++;
		}
	}
	for (v = 0; v < m->nv; v++)
	{
		unsigned int d = (m->vflag[v] & HM_DELETED) ? 0 : hm_feature_degree(m, v);
		if (d && d != 2)
			m->vflag[v] |= HM_CORNER;
	}
	return n;
}

/* the attributes' seams made features: an inside edge whose two faces' corner attributes differ at either end (UVs cut
   apart, split normals) or whose faces' attributes differ (materials); the corners recounted; returns how many were made */
static unsigned int dc_seams(qaws_hm* m)
{
	unsigned int c, v, n = 0, a;
	for (c = 0; c < m->nc; c++)
	{
		unsigned int o = m->co[c], seam = 0;
		if (m->cf[c] == HM_NONE || o == HM_NONE || o < c || (m->cflag[c] & HM_FEATURE))
			continue;
		for (a = 0; a < m->na && !seam; a++)
		{
			hm_attr const* at = &m->attr[a];
			unsigned int k;
			if (at->element == QAWS_MESH_CORNER)
			{
				/* c runs a -> b in its face, o b -> a in the other: a's values at c and next(o), b's at next(c) and o */
				double const *p0 = hm_attr_at(m, a, c), *p1 = hm_attr_at(m, a, hm_next(m, o)), *q0 = hm_attr_at(m, a, hm_next(m, c)), *q1 = hm_attr_at(m, a, o);
				for (k = 0; k < at->n && !seam; k++)
					seam = fabs(p0[k] - p1[k]) > 1e-12 * (1 + fabs(p0[k])) || fabs(q0[k] - q1[k]) > 1e-12 * (1 + fabs(q0[k]));
			}
			else if (at->element == QAWS_MESH_FACE)
			{
				double const *p0 = hm_attr_at(m, a, m->cf[c]), *p1 = hm_attr_at(m, a, m->cf[o]);
				for (k = 0; k < at->n && !seam; k++)
					seam = fabs(p0[k] - p1[k]) > 1e-12 * (1 + fabs(p0[k]));
			}
		}
		if (seam)
		{
			m->cflag[c] |= HM_FEATURE;
			m->cflag[o] |= HM_FEATURE;
			n++;
		}
	}
	for (v = 0; n && v < m->nv; v++)
	{
		unsigned int d = (m->vflag[v] & HM_DELETED) ? 0 : hm_feature_degree(m, v);
		if (d && d != 2)
			m->vflag[v] |= HM_CORNER;
	}
	return n;
}

/* the view directions: given, and N x N on an octahedral map (cell centres; the hemisphere's map turned 45 degrees) */
static double* dc_views(qaws_decimate_desc const* d, unsigned int* out_n)
{
	unsigned int n = d->octahedral_views, given = d->views ? d->view_count : 0, i, j, k = 0;
	double* v = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 3, (size_t)given + (size_t)n * n + 1));
	if (!v)
		return NULL;
	for (i = 0; i < given; i++, k++)
		for (j = 0; j < 3; j++)
			v[3 * (size_t)k + j] = (double)d->views[3 * (size_t)i + j];
	for (j = 0; j < n; j++)
		for (i = 0; i < n; i++, k++)
		{
			double u = 2 * (i + 0.5) / n - 1, w = 2 * (j + 0.5) / n - 1, x, y, z, l;
			if (d->flags & QAWS_DECIMATE_HEMISPHERE)
			{
				x = 0.5 * (u + w);
				y = 0.5 * (u - w);
				z = 1 - fabs(x) - fabs(y);
			}
			else
			{
				x = u;
				y = w;
				z = 1 - fabs(u) - fabs(w);
				if (z < 0)
				{
					x = (1 - fabs(w)) * (u < 0 ? -1 : 1);
					y = (1 - fabs(u)) * (w < 0 ? -1 : 1);
				}
			}
			l = sqrt(x * x + y * y + z * z);
			v[3 * (size_t)k] = x / l;
			v[3 * (size_t)k + 1] = y / l;
			v[3 * (size_t)k + 2] = z / l;
		}
	*out_n = k;
	return v;
}

/* views s0 to s1's silhouette edges (the first npoints views are eye points, the others directions): with share, each counted
   on its two ends; else the plane through the edge and the view's ray added to its ends' quadrics, weight over the views' count */
static unsigned int dc_silhouettes(dc_work* w, double const* views, unsigned int s0, unsigned int s1, unsigned int npoints, double weight,
	double* share, double* seen)
{
	qaws_hm const* m = w->m;
	unsigned int s, c, k, n = 0, nviews = s1 - s0;
	for (s = s0; s < s1; s++)
	{
		double const* view = views + 3 * (size_t)s;
		int points = s < npoints;
		if (seen)
		{
			/* (the faces turned toward the view: its direction points from the mesh to the eye) */
			unsigned int f;
			for (f = 0; f < m->nf; f++)
				if (!(m->fflag[f] & HM_DELETED))
				{
					double nf[3];
					hm_face_normal(m, f, nf);
					seen[f] += nf[0] * view[0] + nf[1] * view[1] + nf[2] * view[2] > 0;
				}
		}
		for (c = 0; c < m->nc; c++)
		{
			unsigned int o = m->co[c], f, g, q, e;
			double nf[3], ng[3], cf[3] = { 0, 0, 0 }, cg[3] = { 0, 0, 0 }, df[3], dg[3], sf, sg;
			if (m->cf[c] == HM_NONE || o == HM_NONE || o < c)
				continue;
			f = m->cf[c];
			g = m->cf[o];
			hm_face_normal(m, f, nf);
			hm_face_normal(m, g, ng);
			for (q = m->fc[f], e = m->fc[f] + m->fn[f]; q < e; q++)
				for (k = 0; k < 3; k++)
					cf[k] += m->x[3 * (size_t)m->cv[q] + k] / m->fn[f];
			for (q = m->fc[g], e = m->fc[g] + m->fn[g]; q < e; q++)
				for (k = 0; k < 3; k++)
					cg[k] += m->x[3 * (size_t)m->cv[q] + k] / m->fn[g];
			for (k = 0; k < 3; k++)
			{
				df[k] = points ? cf[k] - view[k] : view[k];
				dg[k] = points ? cg[k] - view[k] : view[k];
			}
			sf = nf[0] * df[0] + nf[1] * df[1] + nf[2] * df[2];
			sg = ng[0] * dg[0] + ng[1] * dg[1] + ng[2] * dg[2];
			if ((sf > 0) == (sg > 0))
				continue;
			if (share)
			{
				share[m->cv[c]] += 1;
				share[hm_to(m, c)] += 1;
				n++;
				continue;
			}
			{
				unsigned int a = m->cv[c], b = hm_to(m, c);
				double const *pa = m->x + 3 * (size_t)a, *pb = m->x + 3 * (size_t)b;
				double u[3] = { pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2] }, ray[3], pl[3], l2 = u[0] * u[0] + u[1] * u[1] + u[2] * u[2], ls, d;
				for (k = 0; k < 3; k++)
					ray[k] = points ? 0.5 * (pa[k] + pb[k]) - view[k] : view[k];
				pl[0] = u[1] * ray[2] - u[2] * ray[1];
				pl[1] = u[2] * ray[0] - u[0] * ray[2];
				pl[2] = u[0] * ray[1] - u[1] * ray[0];
				ls = sqrt(pl[0] * pl[0] + pl[1] * pl[1] + pl[2] * pl[2]);
				if (!(ls > 0))
					continue;
				for (k = 0; k < 3; k++)
					pl[k] /= ls;
				d = -(pl[0] * pa[0] + pl[1] * pa[1] + pl[2] * pa[2]);
				dc_plane(w->Q + w->K * (size_t)a, w->n1, pl, d, weight * l2 / nviews);
				dc_plane(w->Q + w->K * (size_t)b, w->n1, pl, d, weight * l2 / nviews);
				n++;
			}
		}
	}
	return n;
}

static unsigned int dc_faces(qaws_hm const* m)
{
	unsigned int f, n = 0;
	for (f = 0; f < m->nf; f++)
		n += !(m->fflag[f] & HM_DELETED);
	return n;
}

qaws_status qaws_mesh_3d_decimate(qaws_decimate_report* report, qaws_mesh_3d* mesh, qaws_decimate_desc const* desc)
{
	qaws_decimate_report local;
	qaws_mesh_reference* ref = NULL;
	qaws_hm* m;
	dc_work w;
	dc_heap H;
	unsigned int *ring = NULL, rcap = 64, f, v, k, c, live, target, na = 0, a0;
	double* sopt = NULL;
	double* wtab = NULL;   /* the corner attributes' wedges round a collapse (dc_wedges_before) */
	unsigned int* wnt = NULL;
	double bound2, *views = NULL, *share = NULL, *seen = NULL, sweight = desc && desc->silhouette_weight > 0 ? (double)desc->silhouette_weight : 100;
	unsigned int nviews = 0;
	int failed = 0;
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	if (!report)
		report = &local;
	memset(report, 0, sizeof *report);
	if (!mesh || (desc && (!(desc->max_error >= 0) || !(desc->crease_angle >= 0) || !(desc->max_normal_angle >= 0) || !(desc->silhouette_weight >= 0)
		|| !(desc->feature_weight >= 0) || !(desc->attribute_weight >= 0) || (desc->view_count && !desc->views)
		|| ((desc->flags & QAWS_DECIMATE_MEMORYLESS) && (desc->view_count || desc->octahedral_views || desc->attribute_weight > 0)))))
		return QAWS_STATUS_INVALID_ARGUMENT;
	m = &mesh->hm;
	memset(&w, 0, sizeof w);
	memset(&H, 0, sizeof H);
	live = dc_faces(m);
	report->faces_before = report->faces_after = live;
	target = desc ? desc->target_faces : live / 2;
	bound2 = desc && desc->max_error > 0 ? (double)desc->max_error * (double)desc->max_error : HUGE_VAL;
	if (desc && (desc->flags & QAWS_DECIMATE_MEASURE))
	{
		/* the input kept as a LINEAR reference for the deviation */
		qaws_mesh_3d* copy = NULL;
		qaws_mesh_reference_desc rd;
		memset(&rd, 0, sizeof rd);
		rd.interpolant = QAWS_INTERPOLANT_LINEAR;
		st = qaws_mesh_3d_clone(&copy, mesh);
		if (st == QAWS_STATUS_OK)
			st = qaws_mesh_reference_create(&ref, copy, &rd);
		qaws_mesh_3d_destroy(copy);
		if (st != QAWS_STATUS_OK)
			return st;
		st = QAWS_STATUS_ALLOCATION_FAILURE;
	}
	w.m = m;
	w.cosmax = desc && desc->max_normal_angle > 0 ? cos((double)desc->max_normal_angle) : 0.2;
	w.memless = desc && (desc->flags & QAWS_DECIMATE_MEMORYLESS);
	w.fw = desc && desc->feature_weight > 0 ? (double)desc->feature_weight : 100;
	if (desc && desc->crease_angle > 0)
		report->creases = dc_creases(m, (double)desc->crease_angle);
	report->seams = dc_seams(m);
	/* the vertex attributes in the quadric (Hoppe 1999): LINEAR and NORMALIZED ones, with attribute_weight */
	w.alpha = desc ? (double)desc->attribute_weight : 0;
	for (a0 = 0; w.alpha > 0 && a0 < m->na; a0++)
		if (m->attr[a0].element == QAWS_MESH_VERTEX && (m->attr[a0].rule == QAWS_ATTRIBUTE_LINEAR || m->attr[a0].rule == QAWS_ATTRIBUTE_VECTOR || hm_rule_unit(m->attr[a0].rule)))
			na += m->attr[a0].n;
	w.n1 = 4 + na;
	w.K = w.n1 * (w.n1 + 1) / 2;
	w.ca = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)na + 1));
	w.ck = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)na + 1));
	w.q = (double*)qaws_internal_malloc(sizeof(double) * ((size_t)w.K + 1));
	w.M = (double*)qaws_internal_malloc(sizeof(double) * ((size_t)w.n1 * w.n1 + 1));
	w.r = (double*)qaws_internal_malloc(sizeof(double) * ((size_t)w.n1 + 1));
	w.z = (double*)qaws_internal_malloc(sizeof(double) * ((size_t)w.n1 + 1));
	sopt = (double*)qaws_internal_malloc(sizeof(double) * ((size_t)na + 1));
	wtab = (double*)qaws_internal_malloc(sizeof(double) * (8 * DC_WEDGE_MAX * (size_t)m->na + 1));
	wnt = (unsigned int*)qaws_internal_zalloc((size_t)m->na + 1, sizeof(unsigned int));
	if (!w.ca || !w.ck || !w.q || !w.M || !w.r || !w.z || !sopt || !wtab || !wnt)
		goto done;
	for (a0 = 0, k = 0; na && a0 < m->na; a0++)
		if (m->attr[a0].element == QAWS_MESH_VERTEX && (m->attr[a0].rule == QAWS_ATTRIBUTE_LINEAR || m->attr[a0].rule == QAWS_ATTRIBUTE_VECTOR || hm_rule_unit(m->attr[a0].rule)))
			for (c = 0; c < m->attr[a0].n; c++, k++)
			{
				w.ca[k] = a0;
				w.ck[k] = c;
			}
	w.Q = (double*)qaws_internal_zalloc(qaws_internal_mul2(w.K, (size_t)m->nv + 1), sizeof(double));
	w.W = (double*)qaws_internal_zalloc((size_t)m->nv + 1, sizeof(double));
	w.stamp = (unsigned int*)qaws_internal_zalloc((size_t)m->nv + 1, sizeof(unsigned int));
	w.border = (unsigned char*)qaws_internal_zalloc((size_t)m->nv + 1, 1);
	w.locked = (unsigned char*)qaws_internal_zalloc((size_t)m->nv + 1, 1);
	ring = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * rcap);
	if (!w.Q || !w.W || !w.stamp || !w.border || !w.locked || !ring)
		goto done;
	for (v = 0; v < m->nv; v++)
	{
		w.border[v] = !(m->vflag[v] & HM_DELETED) && hm_is_border_vertex(m, v);
		w.locked[v] = (w.border[v] && desc && (desc->flags & QAWS_DECIMATE_LOCK_BORDER)) || (desc && desc->importance && !(desc->importance[v] > 0));
	}
	{
		/* the views: given (exact planes, below) and octahedral (each vertex's share of their silhouettes weighs its faces' planes) */
		unsigned int given = desc && desc->views ? desc->view_count : 0;
		if (desc && (given || desc->octahedral_views))
		{
			views = dc_views(desc, &nviews);
			if (!views)
				goto done;
		}
		if (nviews > given)
		{
			share = (double*)qaws_internal_zalloc((size_t)m->nv + 1, sizeof(double));
			if (!share)
				goto done;
			seen = (double*)qaws_internal_zalloc((size_t)m->nf + 1, sizeof(double));
			if (!seen)
				goto done;
			report->silhouette_edges += dc_silhouettes(&w, views, given, nviews, 0, 0, share, seen);
			for (v = 0; v < m->nv; v++)
				share[v] = 1 + sweight * share[v] / (nviews - given);
			/* (a face no view sees keeps a fiftieth of its weight: hidden parts go first, still in one piece) */
			for (f = 0; f < m->nf; f++)
				seen[f] = 0.02 + seen[f] / (nviews - given);
		}
	}
	/* each face's plane on its corners, by area; planes standing on the border's and the features' edges, heavier */
	for (f = 0; f < m->nf; f++)
	{
		double n[3], a, d;
		unsigned int e;
		if (m->fflag[f] & HM_DELETED)
			continue;
		hm_face_normal(m, f, n);
		a = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
		if (!(a > 0))
			continue;
		for (k = 0; k < 3; k++)
			n[k] /= a;
		d = -(n[0] * m->x[3 * (size_t)m->cv[m->fc[f]]] + n[1] * m->x[3 * (size_t)m->cv[m->fc[f]] + 1] + n[2] * m->x[3 * (size_t)m->cv[m->fc[f]] + 2]);
		e = m->fc[f] + m->fn[f];
		for (c = m->fc[f]; c < e; c++)
		{
			dc_plane(w.Q + w.K * (size_t)m->cv[c], w.n1, n, d, a * (desc && desc->importance && desc->importance[m->cv[c]] > 0 ? (double)desc->importance[m->cv[c]] : 1)
				* (share ? share[m->cv[c]] * seen[f] : 1));
			w.W[m->cv[c]] += a;
			if (m->co[c] == HM_NONE || (m->cflag[c] & HM_FEATURE))
			{
				double const *p = m->x + 3 * (size_t)m->cv[c], *q = m->x + 3 * (size_t)hm_to(m, c);
				double u[3] = { q[0] - p[0], q[1] - p[1], q[2] - p[2] }, s[3], l2 = u[0] * u[0] + u[1] * u[1] + u[2] * u[2], ls;
				s[0] = u[1] * n[2] - u[2] * n[1];
				s[1] = u[2] * n[0] - u[0] * n[2];
				s[2] = u[0] * n[1] - u[1] * n[0];
				ls = sqrt(s[0] * s[0] + s[1] * s[1] + s[2] * s[2]);
				if (ls > 0)
				{
					for (k = 0; k < 3; k++)
						s[k] /= ls;
					{
					double ds = -(s[0] * p[0] + s[1] * p[1] + s[2] * p[2]);
					double fw = desc && desc->feature_weight > 0 ? (double)desc->feature_weight : 100;
					dc_plane(w.Q + w.K * (size_t)m->cv[c], w.n1, s, ds, fw * l2);
					dc_plane(w.Q + w.K * (size_t)hm_to(m, c), w.n1, s, ds, fw * l2);
					}
				}
			}
		}
	}
	for (f = 0; na && f < m->nf; f++)
	{
		/* each fan triangle's attributes as linear functions over its plane: (g . x + d - s)^2 by its area on its corners */
		unsigned int c0 = m->fc[f], i;
		double n[3], an;
		if (m->fflag[f] & HM_DELETED)
			continue;
		hm_face_normal(m, f, n);
		an = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
		if (!(an > 0))
			continue;
		for (k = 0; k < 3; k++)
			n[k] /= an;
		for (i = 1; i + 1 < m->fn[f]; i++)
		{
			unsigned int tv[3], j, q;
			double const* P[3];
			double e1[3], e2[3], cr[3], area;
			tv[0] = m->cv[c0];
			tv[1] = m->cv[c0 + i];
			tv[2] = m->cv[c0 + i + 1];
			for (j = 0; j < 3; j++) P[j] = m->x + 3 * (size_t)tv[j];
			for (k = 0; k < 3; k++)
			{
				e1[k] = P[1][k] - P[0][k];
				e2[k] = P[2][k] - P[0][k];
			}
			cr[0] = e1[1] * e2[2] - e1[2] * e2[1];
			cr[1] = e1[2] * e2[0] - e1[0] * e2[2];
			cr[2] = e1[0] * e2[1] - e1[1] * e2[0];
			area = 0.5 * sqrt(cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]);
			if (!(area > 0))
				continue;
			for (q = 0; q < na; q++)
			{
				double M4[16], r4[4];
				for (j = 0; j < 3; j++)
				{
					for (k = 0; k < 3; k++) M4[4 * j + k] = P[j][k];
					M4[4 * j + 3] = 1;
					r4[j] = w.alpha * hm_attr_at(m, w.ca[q], tv[j])[w.ck[q]];
				}
				for (k = 0; k < 3; k++) M4[12 + k] = n[k];
				M4[15] = 0;
				r4[3] = 0;
				if (!dc_solve(M4, r4, 4))
					continue;
				for (k = 0; k < w.n1; k++) w.z[k] = 0;
				for (k = 0; k < 3; k++) w.z[k] = r4[k];
				w.z[3 + q] = -1;
				w.z[w.n1 - 1] = r4[3];
				for (j = 0; j < 3; j++)
				{
					double imp = desc->importance && desc->importance[tv[j]] > 0 ? (double)desc->importance[tv[j]] : 1;
					unsigned int u, t2;
					double* Qv = w.Q + w.K * (size_t)tv[j];
					for (u = 0; u < w.n1; u++)
						for (t2 = u; t2 < w.n1; t2++)
							Qv[dc_at(w.n1, u, t2)] += area * imp * w.z[u] * w.z[t2];
				}
			}
		}
	}
	if (desc && desc->views && desc->view_count)
		report->silhouette_edges += dc_silhouettes(&w, views, 0, desc->view_count, (desc->flags & QAWS_DECIMATE_VIEW_POINTS) ? desc->view_count : 0,
			sweight, NULL, NULL);
	for (v = 0; v < m->nv; v++)
		if (!(m->vflag[v] & HM_DELETED) && !dc_push_edges(&w, &H, v, &ring, &rcap))
			goto done;
	while (H.n && live > target)
	{
		dc_entry e;
		double x[3], cost, t, l2;
		unsigned int a, b, gone = 0;
		dc_pop(&H, &e);
		a = e.a;
		b = e.b;
		if ((m->vflag[a] & HM_DELETED) || (m->vflag[b] & HM_DELETED) || e.sa != w.stamp[a] || e.sb != w.stamp[b])
			continue;   /* (stale) */
		if (e.cost > bound2)
			break;
		if (!dc_candidate(&w, a, b, x, sopt, &cost) || !dc_no_fold(m, a, b, x, w.cosmax))
			continue;
		/* the faces of the edge that go: its triangles */
		for (c = m->vc[a]; c != HM_NONE; c = hm_round(m, a, c))
		{
			unsigned int g = m->cf[c], q, end = m->fc[g] + m->fn[g];
			for (q = m->fc[g]; q < end; q++)
				gone += m->fn[g] == 3 && m->cv[q] == b;
		}
		l2 = 0;
		t = 0;
		for (k = 0; k < 3; k++)
		{
			double ab = m->x[3 * (size_t)b + k] - m->x[3 * (size_t)a + k];
			t += (x[k] - m->x[3 * (size_t)a + k]) * ab;
			l2 += ab * ab;
		}
		t = l2 > 0 ? t / l2 : 1;
		t = t < 0 ? 0 : t > 1 ? 1 : t;
		dc_wedges_before(m, a, b, x, wtab, wnt, 8 * DC_WEDGE_MAX);
		if (!hm_collapse(m, hm_find_edge(m, a, b), a, t, x))
		{
			failed = 1;   /* (out of memory: what was done stays, compacted) */
			break;
		}
		dc_wedges_after(m, b, wtab, wnt, 8 * DC_WEDGE_MAX);
		for (k = 0; k < na; k++)
			hm_attr_at(m, w.ca[k], b)[w.ck[k]] = sopt[k] / w.alpha;
		for (a0 = 0; na && a0 < m->na; a0++)
			if (m->attr[a0].element == QAWS_MESH_VERTEX && hm_rule_unit(m->attr[a0].rule))
				hm_space_out(hm_attr_at(m, a0, b), m->attr[a0].n, m->attr[a0].rule);
		for (k = 0; k < w.K; k++)
			w.Q[w.K * (size_t)b + k] += w.Q[w.K * (size_t)a + k];
		w.W[b] += w.W[a];
		w.border[b] |= w.border[a];
		w.stamp[a]++;
		w.stamp[b]++;
		live -= gone;
		report->collapses++;
		if (sqrt(cost) > (double)report->max_error)
			report->max_error = (qaws_scalar)sqrt(cost);
		if (!dc_push_edges(&w, &H, b, &ring, &rcap))
		{
			failed = 1;
			break;
		}
		if (w.memless)
		{
			/* (memoryless: the faces round b's neighbours changed too, so their edges' costs) */
			unsigned int nr = hm_ring(m, b, ring, rcap), q, *nb2 = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)nr + 1));
			if (!nb2 || nr > rcap)
			{
				qaws_internal_free(nb2);
				failed = !nb2;
				if (failed)
					break;
				continue;
			}
			memcpy(nb2, ring, sizeof(unsigned int) * nr);
			for (q = 0; q < nr && !failed; q++)
			{
				w.stamp[nb2[q]]++;
				failed = !dc_push_edges(&w, &H, nb2[q], &ring, &rcap);
			}
			qaws_internal_free(nb2);
			if (failed)
				break;
		}
	}
	hm_compact(m);
	report->faces_after = dc_faces(m);
	st = failed ? QAWS_STATUS_ALLOCATION_FAILURE : QAWS_STATUS_OK;
	if (ref)
		st = qaws_mesh_3d_measure_deviation(&report->deviation, ref, mesh, NULL);
done:
	qaws_mesh_reference_destroy(ref);
	qaws_internal_free(w.Q); qaws_internal_free(w.W); qaws_internal_free(w.stamp); qaws_internal_free(w.border); qaws_internal_free(w.locked); qaws_internal_free(views); qaws_internal_free(share); qaws_internal_free(seen);
	qaws_internal_free(ring); qaws_internal_free(H.e);
	qaws_internal_free(w.ca); qaws_internal_free(w.ck); qaws_internal_free(w.q); qaws_internal_free(w.M); qaws_internal_free(w.r); qaws_internal_free(w.z);
	qaws_internal_free(sopt);
	qaws_internal_free(wtab);
	qaws_internal_free(wnt);
	return st;
}
