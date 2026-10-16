/*
 * Barnes-Hut evaluation of the tangent-point energies and their gradients
 * (qaws_repulsive.h, R4; Repulsive Curves, Sec. 6.1; Repulsive Surfaces,
 * Sec. 4.1 and 4.2).
 *
 * The energies sum ordered pairs (I, J) of elements (edges, faces) whose
 * term e(I, J) uses I's own tangent or normal. For a target I, the elements
 * J of a tree node far enough away are taken at once, by their total
 * length or area W and weighted centroid X (one-point quadrature):
 *   sum over J in the node of e(I, J) ~ w_I W k(x_I, X, T_I).
 * A node is far for I when max(r_I, r_node) <= theta dist(box_I, box_node)
 * (r the boxes' half diagonals, dist between boxes, so a node holding an
 * element sharing a vertex with I, at distance 0, is never far); the
 * others are opened down to their leaves, where the pairs are the exact
 * discrete terms of the direct evaluator. theta = 0 opens everything: the
 * direct sum.
 *
 * The gradient follows the paper's scheme (Surfaces, Sec. 4.2): each
 * element's vertices get, from its own traversal, the derivative of both
 * orders of each pair, e(I, J) + e(J, I) in I's vertices only: exact for
 * the pairs near I; for a far node, the forward term above and the reverse
 * term, every J's tangent or normal replaced by the node's averaged
 * projector P (sum w (I - T T^T) for edges, sum w N N^T for faces, over W):
 *   sum over J in the node of e(J, I) ~ w_I W (d^T P d)^(alpha/2) / |d|^beta.
 * So the gradient approximates the direct gradient, converging to it as
 * theta falls; it is not the derivative of the approximated energy.
 */

#include "qaws_repulsive.h"
#include "internal/qaws_internal_repulsive.h"
#include "internal/qaws_internal_parallel.h"
#include "internal/qaws_internal_broadphase.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define BH_LEAF 8

/* x^y for x >= 0; when 2 y is a whole number (alpha / 2, beta / 2, p usually
   are) by products and one square root: exact to rounding, and several
   times cheaper than pow, which every pair term calls */
static double bh_pow(double x, double y)
{
	double y2 = y + y, r = 1;
	int n, i;
	if (!(y2 == floor(y2) && fabs(y2) <= 64))
		return pow(x, y);
	n = (int)fabs(y2);
	for (i = 0; i < n / 2; i++) r *= x;
	if (n % 2)
		r *= sqrt(x);
	return y < 0 ? 1 / r : r;
}

/* ------------------------------------------------------------------------ */
/*  A bounding volume hierarchy over items                                  */
/* ------------------------------------------------------------------------ */

typedef struct bh_node
{
	unsigned int begin, end;        /* in perm */
	int child[2];                   /* -1: a leaf */
	double lo[3], hi[3], r;         /* the items' boxes' box, its half diagonal */
	double W, S[3], Q[6];           /* weight; weighted centroid's numerator; weighted projector (xx xy xz yy yz zz) */
} bh_node;

typedef struct bh_tree
{
	bh_node* node;
	unsigned int nn, cap, dim;
	unsigned int* perm;
	double const* ilo;              /* item boxes, dim each */
	double const* ihi;
	double const* ic;               /* item centers */
} bh_tree;

/* perm[begin, end) reordered so that perm[mid] has the median key along
   axis, smaller keys before it, larger after (quickselect; no globals) */
static void bh_select(unsigned int* perm, unsigned int begin, unsigned int end, unsigned int mid, double const* c, unsigned int dim, unsigned int axis)
{
	while (end - begin > 1)
	{
		/* (the lower middle: Hoare's scheme then leaves j < end - 1, so the range shrinks) */
		double pv = c[dim * perm[begin + (end - begin - 1) / 2] + axis];
		unsigned int i = begin, j = end - 1, t;
		for (;;)
		{
			while (c[dim * perm[i] + axis] < pv) i++;
			while (c[dim * perm[j] + axis] > pv) j--;
			if (i >= j)
				break;
			t = perm[i]; perm[i] = perm[j]; perm[j] = t;
			i++; j--;
		}
		if (mid <= j)
			end = j + 1;
		else
			begin = j + 1;
	}
}

static int bh_build(bh_tree* t, unsigned int begin, unsigned int end)
{
	int id;
	unsigned int i, q, axis = 0;
	bh_node* n;
	if (t->nn == t->cap)
	{
		unsigned int nc = t->cap ? 2 * t->cap : 64;
		bh_node* g = (bh_node*)realloc(t->node, sizeof(bh_node) * nc);
		if (!g)
			return -1;
		t->node = g;
		t->cap = nc;
	}
	id = (int)t->nn++;
	n = &t->node[id];
	n->begin = begin; n->end = end; n->child[0] = n->child[1] = -1;
	for (q = 0; q < 3; q++) { n->lo[q] = 0; n->hi[q] = 0; }
	for (q = 0; q < t->dim; q++) { n->lo[q] = 1e300; n->hi[q] = -1e300; }
	for (i = begin; i < end; i++)
		for (q = 0; q < t->dim; q++)
		{
			unsigned int it = t->perm[i];
			n->lo[q] = fmin(n->lo[q], t->ilo[t->dim * it + q]);
			n->hi[q] = fmax(n->hi[q], t->ihi[t->dim * it + q]);
		}
	n->r = 0;
	for (q = 0; q < t->dim; q++)
	{
		n->r += 0.25 * (n->hi[q] - n->lo[q]) * (n->hi[q] - n->lo[q]);
		if (n->hi[q] - n->lo[q] > n->hi[axis] - n->lo[axis]) axis = q;
	}
	n->r = sqrt(n->r);
	if (end - begin > BH_LEAF)
	{
		/* split at the median center along the widest axis */
		unsigned int mid = (begin + end) / 2;
		int c0, c1;
		bh_select(t->perm, begin, end, mid, t->ic, t->dim, axis);
		c0 = bh_build(t, begin, mid);
		if (c0 < 0)
			return -1;
		c1 = bh_build(t, mid, end);
		if (c1 < 0)
			return -1;
		t->node[id].child[0] = c0; t->node[id].child[1] = c1;
	}
	return id;
}

static double bh_boxdist(double const* alo, double const* ahi, double const* blo, double const* bhi, unsigned int dim)
{
	double s = 0;
	unsigned int q;
	for (q = 0; q < dim; q++)
	{
		double g = fmax(0, fmax(alo[q] - bhi[q], blo[q] - ahi[q]));
		s += g * g;
	}
	return sqrt(s);
}

/* the node sums from the items' weights, centers and projectors */
static void bh_sums(bh_tree* t, unsigned int id, double const* w, double const* qp)
{
	bh_node* n = &t->node[id];
	unsigned int i, q;
	n->W = 0;
	for (q = 0; q < 3; q++) n->S[q] = 0;
	for (q = 0; q < 6; q++) n->Q[q] = 0;
	if (n->child[0] >= 0)
	{
		bh_node const *a, *b;
		bh_sums(t, (unsigned int)n->child[0], w, qp);
		bh_sums(t, (unsigned int)n->child[1], w, qp);
		a = &t->node[n->child[0]]; b = &t->node[n->child[1]];
		n->W = a->W + b->W;
		for (q = 0; q < 3; q++) n->S[q] = a->S[q] + b->S[q];
		for (q = 0; q < 6; q++) n->Q[q] = a->Q[q] + b->Q[q];
		return;
	}
	for (i = n->begin; i < n->end; i++)
	{
		unsigned int it = t->perm[i];
		n->W += w[it];
		for (q = 0; q < t->dim; q++) n->S[q] += w[it] * t->ic[t->dim * it + q];
		for (q = 0; q < 6; q++) n->Q[q] += qp[6 * it + q];
	}
}

/* ------------------------------------------------------------------------ */
/*  The traversal                                                           */
/* ------------------------------------------------------------------------ */

typedef struct bh_ctx bh_ctx;

struct bh_ctx
{
	bh_tree t;
	unsigned int nitems, k;         /* items, vertices per item */
	unsigned int const* idx;        /* items' vertices */
	double const* x;
	unsigned int dim;
	double theta, alpha, beta;      /* (meshes: p and 2p) */
	double* grad;                   /* NULL: energy only */
	double* blk;                    /* the current item's gradient block (k dim) */
	qaws_batch_executor const* executor;
	double E;
	unsigned long long near, far;
	int status;
	int skip_shared;                /* (curves: edges sharing a vertex do not interact) */
	/* a far node for I in closed form: its forward term (the energy), and if g the gradient of forward and reverse in I's vertices */
	int (*far_g)(bh_ctx const*, unsigned int I, bh_node const* node, double* E, double* g);
	/* e(A, B) in closed form, and if g its gradient in A's (which 0) or B's (1) vertices */
	int (*near_g)(bh_ctx const*, unsigned int A, unsigned int B, int which, double* E, double* g);
};

static int bh_shares(bh_ctx const* c, unsigned int I, unsigned int J)
{
	unsigned int a, b;
	for (a = 0; a < c->k; a++)
		for (b = 0; b < c->k; b++)
			if (c->idx[c->k * I + a] == c->idx[c->k * J + b])
				return 1;
	return 0;
}

static void bh_visit(bh_ctx* c, unsigned int I, unsigned int id)
{
	bh_node const* n = &c->t.node[id];
	double const* ilo = c->t.ilo + c->dim * I;
	double const* ihi = c->t.ihi + c->dim * I;
	double ri = 0, d;
	unsigned int q, slots = c->grad ? c->k * c->dim : 0, a;
	if (c->status)
		return;
	for (q = 0; q < c->dim; q++) ri += 0.25 * (ihi[q] - ilo[q]) * (ihi[q] - ilo[q]);
	ri = sqrt(ri);
	d = bh_boxdist(ilo, ihi, n->lo, n->hi, c->dim);
	if (c->theta > 0 && fmax(ri, n->r) <= c->theta * d)
	{
		/* far: the node at once */
		double e, gi[9];
		int r;
		if ((r = c->far_g(c, I, n, &e, c->grad ? gi : NULL)) != 0)
		{
			c->status = r;
			return;
		}
		c->E += e;
		c->far++;
		for (a = 0; c->grad && a < slots; a++) c->blk[a] += gi[a];
		return;
	}
	if (n->child[0] >= 0)
	{
		bh_visit(c, I, (unsigned int)n->child[0]);
		bh_visit(c, I, (unsigned int)n->child[1]);
		return;
	}
	for (a = n->begin; a < n->end; a++)
	{
		unsigned int J = c->t.perm[a], s;
		int r;
		double e, gi[9];
		if (J == I || (c->skip_shared && bh_shares(c, I, J)))
			continue;
		/* e(I, J) for the energy; with e(J, I), in I's vertices, for the gradient */
		if ((r = c->near_g(c, I, J, 0, &e, c->grad ? gi : NULL)) != 0)
		{
			c->status = r;
			return;
		}
		c->E += e;
		c->near++;
		if (!c->grad)
			continue;
		for (s = 0; s < slots; s++) c->blk[s] += gi[s];
		if ((r = c->near_g(c, J, I, 1, &e, gi)) != 0)
		{
			c->status = r;
			return;
		}
		for (s = 0; s < slots; s++) c->blk[s] += gi[s];
	}
}

/* ------------------------------------------------------------------------ */
/*  Curves                                                                  */
/* ------------------------------------------------------------------------ */

/* e(A, B) = l_A^(1 - alpha) l_B / 4 sum |e_A x d|^alpha |d|^-beta (the direct
   evaluator's ordered term) in closed form, and if g, its gradient in the
   vertices of A (which 0) or B (which 1), dim each; 0, or
   1 degenerate, 2 coincident */
static int bh_curve_near_g(bh_ctx const* c, unsigned int A, unsigned int B, int which, double* E, double* g)
{
	unsigned int dim = c->dim, i, j, q;
	unsigned int const* P = c->idx + 2 * A;
	unsigned int const* Q = c->idx + 2 * B;
	double e[3] = { 0, 0, 0 }, f[3] = { 0, 0, 0 }, LA = 0, LB = 0, W, S = 0, ge[3] = { 0, 0, 0 };
	double a2 = 0.5 * c->alpha, b2 = 0.5 * c->beta;
	for (q = 0; q < dim; q++)
	{
		e[q] = c->x[dim * P[1] + q] - c->x[dim * P[0] + q]; f[q] = c->x[dim * Q[1] + q] - c->x[dim * Q[0] + q];
		LA += e[q] * e[q]; LB += f[q] * f[q];
	}
	if (!(LA > 0 && LB > 0))
		return 1;
	W = 0.25 * bh_pow(LA, 0.5 * (1 - c->alpha)) * sqrt(LB);
	if (g) for (q = 0; q < 2 * dim; q++) g[q] = 0;
	for (i = 0; i < 2; i++)
		for (j = 0; j < 2; j++)
		{
			double d[3] = { 0, 0, 0 }, r2 = 0, ed = 0, c2, k, kc, kr;
			for (q = 0; q < dim; q++) { d[q] = c->x[dim * P[i] + q] - c->x[dim * Q[j] + q]; r2 += d[q] * d[q]; ed += e[q] * d[q]; }
			if (!(r2 > 0))
				return 2;
			c2 = fmax(LA * r2 - ed * ed, 0);
			k = bh_pow(c2, a2) * bh_pow(r2, -b2);
			S += k;
			if (!g)
				continue;
			/* dk = kc dc2 + kr dr2; dc2 = 2 (r2 e - (e.d) d) . de + 2 (LA d - (e.d) e) . dd; dr2 = 2 d . dd */
			kc = c2 > 0 ? a2 * bh_pow(c2, a2 - 1) * bh_pow(r2, -b2) : (c->alpha == 2 ? bh_pow(r2, -b2) : 0);
			kr = -b2 * k / r2;
			for (q = 0; q < dim; q++)
			{
				double gd = 2 * kc * (LA * d[q] - ed * e[q]) + 2 * kr * d[q];   /* along d */
				if (which == 0)
				{
					double gee = 2 * kc * (r2 * e[q] - ed * d[q]);
					ge[q] += gee;
					g[dim * i + q] += W * gd;
				}
				else
					g[dim * j + q] -= W * gd;
			}
		}
	*E = W * S;
	if (g)
	{
		/* W's own: (1 - alpha) / 2 dLA / LA (A's edge), 1/2 dLB / LB (B's) */
		for (q = 0; q < dim; q++)
		{
			if (which == 0)
			{
				double t = W * ge[q] + (1 - c->alpha) * W * S * e[q] / LA;
				g[dim + q] += t; g[q] -= t;
			}
			else
			{
				double t = W * S * f[q] / LB;
				g[dim + q] += t; g[q] -= t;
			}
		}
	}
	return 0;
}

/* A far node for edge I, in closed form. I's side keeps the direct term's
   trapezoid over its two endpoints; the node's members are taken at their
   centroid X (d_i = x_i - X):
   forward  F = W / 2 l_I^(1 - alpha) sum_i |e_I x d_i|^alpha |d_i|^-beta (the energy);
   reverse  R = W / 2 l_I sum_i (d_i^T P d_i)^(alpha/2) |d_i|^-beta, P the
            node's averaged projector Q / W.
   *E = F; with g, the gradient of F + R in I's two endpoints (g[0..5]).
   0, or 2 when an endpoint sits on X. */
static int bh_curve_far_g(bh_ctx const* c, unsigned int I, bh_node const* node, double* E, double* g)
{
	static unsigned int const ix[3][3] = { { 0, 1, 2 }, { 1, 3, 4 }, { 2, 4, 5 } };
	unsigned int dim = c->dim, i, q, r;
	unsigned int const* P = c->idx + 2 * I;
	double p[2][3] = { { 0, 0, 0 }, { 0, 0, 0 } }, e[3], X[3] = { 0, 0, 0 }, la = 0, l, Wf, Sf = 0, Sr = 0, ge[3] = { 0, 0, 0 };
	double a2 = 0.5 * c->alpha, b2 = 0.5 * c->beta, W = node->W;
	for (i = 0; i < 2; i++)
		for (q = 0; q < dim; q++) p[i][q] = c->x[dim * P[i] + q];
	for (q = 0; q < 3; q++)
	{
		e[q] = p[1][q] - p[0][q];
		la += e[q] * e[q];
		if (q < dim) X[q] = node->S[q] / W;
	}
	if (!(la > 0))
		return 1;
	l = sqrt(la);
	Wf = 0.5 * W * bh_pow(la, 0.5 * (1 - c->alpha));
	if (g)
		for (q = 0; q < 6; q++) g[q] = 0;
	for (i = 0; i < 2; i++)
	{
		double d[3], w[3], r2 = 0, c2, k, inv;
		for (q = 0; q < 3; q++) { d[q] = p[i][q] - X[q]; r2 += d[q] * d[q]; }
		if (!(r2 > 0))
			return 2;
		w[0] = e[1] * d[2] - e[2] * d[1]; w[1] = e[2] * d[0] - e[0] * d[2]; w[2] = e[0] * d[1] - e[1] * d[0];
		c2 = w[0] * w[0] + w[1] * w[1] + w[2] * w[2];
		inv = bh_pow(r2, -b2);
		k = bh_pow(c2, a2) * inv;
		Sf += k;
		if (g)
		{
			/* forward: dk = kc dc2 + kr dr2, dc2 = 2 (d x w) . de + 2 (w x e) . dd */
			double kc = c2 > 0 ? a2 * bh_pow(c2, a2 - 1) * inv : (c->alpha == 2 ? inv : 0), kr = -b2 * k / r2, dw[3], we[3];
			/* reverse: q = d^T P d, kq = q^(a/2) r2^-b2 */
			double Pd[3] = { 0, 0, 0 }, qd = 0, kq, kqq;
			dw[0] = d[1] * w[2] - d[2] * w[1]; dw[1] = d[2] * w[0] - d[0] * w[2]; dw[2] = d[0] * w[1] - d[1] * w[0];
			we[0] = w[1] * e[2] - w[2] * e[1]; we[1] = w[2] * e[0] - w[0] * e[2]; we[2] = w[0] * e[1] - w[1] * e[0];
			for (q = 0; q < 3; q++)
			{
				for (r = 0; r < 3; r++) Pd[q] += node->Q[ix[q][r]] / W * d[r];
				qd += d[q] * Pd[q];
			}
			qd = qd > 0 ? qd : 0;
			kq = bh_pow(qd, a2) * inv;
			kqq = qd > 0 ? a2 * bh_pow(qd, a2 - 1) * inv : (c->alpha == 2 ? inv : 0);
			Sr += kq;
			for (q = 0; q < 3; q++)
			{
				ge[q] += Wf * 2 * kc * dw[q];
				g[3 * i + q] += Wf * (2 * kc * we[q] + 2 * kr * d[q]) + 0.5 * W * l * (2 * kqq * Pd[q] - 2 * b2 * kq * d[q] / r2);
			}
		}
	}
	*E = Wf * Sf;
	if (g)
		for (q = 0; q < 3; q++)
		{
			/* the weights' own: (1 - alpha) e / la for F, e / la for R */
			double t = ge[q] + (1 - c->alpha) * Wf * Sf * e[q] / la + 0.5 * W * l * Sr * e[q] / la;
			g[3 + q] += t; g[q] -= t;
		}
	if (dim == 2 && g)
	{
		/* (2D: back to 2 coordinates per vertex) */
		g[2] = g[3]; g[3] = g[4];
	}
	return 0;
}

/* ------------------------------------------------------------------------ */
/*  Meshes                                                                  */
/* ------------------------------------------------------------------------ */

/* d (u . m) / d corner c of a face with m = e1 x e2 (e1 = x1 - x0, e2 =
   x2 - x0): e2 x u, u x e1, and minus their sum at corner 0 */
static void bh_dm(double const* e1, double const* e2, double const* u, double* out)
{
	double a[3], b[3];
	unsigned int q;
	a[0] = e2[1] * u[2] - e2[2] * u[1]; a[1] = e2[2] * u[0] - e2[0] * u[2]; a[2] = e2[0] * u[1] - e2[1] * u[0];
	b[0] = u[1] * e1[2] - u[2] * e1[1]; b[1] = u[2] * e1[0] - u[0] * e1[2]; b[2] = u[0] * e1[1] - u[1] * e1[0];
	for (q = 0; q < 3; q++) { out[q] = -a[q] - b[q]; out[3 + q] = a[q]; out[6 + q] = b[q]; }
}

/* e(A, B) = a_A a_B (N_A . d)^p |d|^-2p in closed form: 1/4 |m_A|^(1 - p) |m_B| |m_A . d|^p
   |d|^-2p, and if g, its gradient in the corners of A (which 0) or B (1);
   0, or 1 degenerate, 2 coincident */
static int bh_mesh_near_g(bh_ctx const* c, unsigned int A, unsigned int B, int which, double* E, double* g)
{
	double p[2][3][3], e1[2][3], e2[2][3], m[2][3], X[2][3], d[3], a, b, s = 0, r2 = 0, e, al = c->alpha, b2 = 0.5 * c->beta, J[9];
	unsigned int f, k, q;
	for (f = 0; f < 2; f++)
	{
		unsigned int const* t = c->idx + 3 * (f ? B : A);
		for (k = 0; k < 3; k++)
			for (q = 0; q < 3; q++) p[f][k][q] = c->x[3 * t[k] + q];
		for (q = 0; q < 3; q++)
		{
			e1[f][q] = p[f][1][q] - p[f][0][q]; e2[f][q] = p[f][2][q] - p[f][0][q];
			X[f][q] = (p[f][0][q] + p[f][1][q] + p[f][2][q]) / 3;
		}
		m[f][0] = e1[f][1] * e2[f][2] - e1[f][2] * e2[f][1]; m[f][1] = e1[f][2] * e2[f][0] - e1[f][0] * e2[f][2]; m[f][2] = e1[f][0] * e2[f][1] - e1[f][1] * e2[f][0];
	}
	a = sqrt(m[0][0] * m[0][0] + m[0][1] * m[0][1] + m[0][2] * m[0][2]);
	b = sqrt(m[1][0] * m[1][0] + m[1][1] * m[1][1] + m[1][2] * m[1][2]);
	if (!(a > 0 && b > 0))
		return 1;
	for (q = 0; q < 3; q++) { d[q] = X[0][q] - X[1][q]; s += m[0][q] * d[q]; r2 += d[q] * d[q]; }
	if (!(r2 > 0))
		return 2;
	e = 0.25 * bh_pow(a, 1 - al) * b * bh_pow(fabs(s), al) * bh_pow(r2, -b2);
	*E = e;
	if (!g)
		return 0;
	/* de = e ((1 - p) da / a + db / b - b2 dr2 / r2) + 1/4 a^(1 - p) b p |s|^(p - 1) sgn(s) ds r2^-b2 */
	{
		double ks = 0.25 * bh_pow(a, 1 - al) * b * al * bh_pow(fabs(s), al - 1) * (s < 0 ? -1 : 1) * bh_pow(r2, -b2), u[3];
		if (which == 0)
		{
			/* a: u = m_A / a; s: d . dm_A + m_A . dd; d moves by 1/3 of each corner */
			for (q = 0; q < 3; q++) u[q] = e * (1 - al) * m[0][q] / (a * a) + ks * d[q];
			bh_dm(e1[0], e2[0], u, J);
			for (k = 0; k < 3; k++)
				for (q = 0; q < 3; q++) g[3 * k + q] = J[3 * k + q] + (ks * m[0][q] - 2 * b2 * e * d[q] / r2) / 3;
		}
		else
		{
			for (q = 0; q < 3; q++) u[q] = e * m[1][q] / (b * b);
			bh_dm(e1[1], e2[1], u, J);
			for (k = 0; k < 3; k++)
				for (q = 0; q < 3; q++) g[3 * k + q] = J[3 * k + q] - (ks * m[0][q] - 2 * b2 * e * d[q] / r2) / 3;
		}
	}
	return 0;
}

/* A far node for face I, in closed form, d = X_I - X (X the node's centroid):
   forward  F = W / 2 |m_I|^(1 - p) |m_I . d|^p |d|^-2p (the energy);
   reverse  R = W / 2 |m_I| (d^T P d)^(p/2) |d|^-2p, P = Q / W.
   *E = F; with g, the gradient of F + R in I's corners (g[0..8]).
   0, or 1 degenerate, 2 when the barycenters coincide. */
static int bh_mesh_far_g(bh_ctx const* c, unsigned int I, bh_node const* node, double* E, double* g)
{
	static unsigned int const ix[3][3] = { { 0, 1, 2 }, { 1, 3, 4 }, { 2, 4, 5 } };
	unsigned int const* t = c->idx + 3 * I;
	double p[3][3], e1[3], e2[3], m[3], d[3], a, s = 0, r2 = 0, F, inv, al = c->alpha, W = node->W, J[9];
	unsigned int k, q, r;
	for (k = 0; k < 3; k++)
		for (q = 0; q < 3; q++) p[k][q] = c->x[3 * t[k] + q];
	for (q = 0; q < 3; q++)
	{
		e1[q] = p[1][q] - p[0][q]; e2[q] = p[2][q] - p[0][q];
		d[q] = (p[0][q] + p[1][q] + p[2][q]) / 3 - node->S[q] / W;
	}
	m[0] = e1[1] * e2[2] - e1[2] * e2[1]; m[1] = e1[2] * e2[0] - e1[0] * e2[2]; m[2] = e1[0] * e2[1] - e1[1] * e2[0];
	a = sqrt(m[0] * m[0] + m[1] * m[1] + m[2] * m[2]);
	if (!(a > 0))
		return 1;
	for (q = 0; q < 3; q++) { s += m[q] * d[q]; r2 += d[q] * d[q]; }
	if (!(r2 > 0))
		return 2;
	inv = bh_pow(r2, -al);
	F = 0.5 * W * bh_pow(a, 1 - al) * bh_pow(fabs(s), al) * inv;
	*E = F;
	if (!g)
		return 0;
	{
		/* forward: as the near pair with 1/4 |m_B| replaced by W / 2 */
		double ks = 0.5 * W * bh_pow(a, 1 - al) * al * bh_pow(fabs(s), al - 1) * (s < 0 ? -1 : 1) * inv, u[3], Pd[3] = { 0, 0, 0 }, qd = 0, R, kq, gd[3];
		for (q = 0; q < 3; q++)
			for (r = 0; r < 3; r++) Pd[q] += node->Q[ix[q][r]] / W * d[r];
		for (q = 0; q < 3; q++) qd += d[q] * Pd[q];
		qd = qd > 0 ? qd : 0;
		R = 0.5 * W * a * bh_pow(qd, 0.5 * al) * inv;
		/* reverse: dR = R m . dm / a^2 + W/2 a (p q^(p/2 - 1) Pd - 2p q^(p/2) d / r2) r2^-p . dd */
		kq = 0.5 * W * a * (qd > 0 ? al * bh_pow(qd, 0.5 * al - 1) : (al == 2 ? 2 : 0)) * inv;
		for (q = 0; q < 3; q++)
		{
			u[q] = F * (1 - al) * m[q] / (a * a) + ks * d[q] + R * m[q] / (a * a);
			gd[q] = ks * m[q] - 2 * al * F * d[q] / r2 + kq * Pd[q] - 2 * al * R * d[q] / r2;
		}
		bh_dm(e1, e2, u, J);
		for (k = 0; k < 3; k++)
			for (q = 0; q < 3; q++) g[3 * k + q] = J[3 * k + q] + gd[q] / 3;
	}
	return 0;
}

/* ------------------------------------------------------------------------ */
/*  The drivers                                                             */
/* ------------------------------------------------------------------------ */

/* items' boxes, centers, weights and weighted projectors at x: edges (k 2)
   by their midpoints, lengths and l (I - T T^T); faces (k 3) by their
   barycenters, areas and a N N^T */
static int bh_items(bh_ctx* c, double* lo, double* hi, double* ctr, double* w, double* qp)
{
	unsigned int i, a, q, dim = c->dim;
	for (i = 0; i < c->nitems; i++)
	{
		unsigned int const* v = c->idx + c->k * i;
		double u[3] = { 0, 0, 0 }, P[3][3];
		unsigned int r, s;
		for (q = 0; q < dim; q++)
		{
			lo[dim * i + q] = 1e300; hi[dim * i + q] = -1e300; ctr[dim * i + q] = 0;
			for (a = 0; a < c->k; a++)
			{
				double y = c->x[dim * v[a] + q];
				lo[dim * i + q] = fmin(lo[dim * i + q], y);
				hi[dim * i + q] = fmax(hi[dim * i + q], y);
				ctr[dim * i + q] += y / c->k;
			}
		}
		if (c->k == 2)
		{
			for (q = 0; q < dim; q++) u[q] = c->x[dim * v[1] + q] - c->x[dim * v[0] + q];
			w[i] = sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
			if (!(w[i] > 0))
				return 0;
			for (r = 0; r < 3; r++)
				for (s = 0; s < 3; s++) P[r][s] = (r == s && r < dim ? w[i] : 0) - u[r] * u[s] / w[i];
		}
		else
		{
			double e1[3], e2[3], m[3], mm;
			for (q = 0; q < 3; q++) { e1[q] = c->x[3 * v[1] + q] - c->x[3 * v[0] + q]; e2[q] = c->x[3 * v[2] + q] - c->x[3 * v[0] + q]; }
			m[0] = e1[1] * e2[2] - e1[2] * e2[1]; m[1] = e1[2] * e2[0] - e1[0] * e2[2]; m[2] = e1[0] * e2[1] - e1[1] * e2[0];
			mm = sqrt(m[0] * m[0] + m[1] * m[1] + m[2] * m[2]);
			w[i] = 0.5 * mm;
			if (!(w[i] > 0))
				return 0;
			for (r = 0; r < 3; r++)
				for (s = 0; s < 3; s++) P[r][s] = m[r] * m[s] / (2 * mm);   /* a N N^T */
		}
		qp[6 * i] = P[0][0]; qp[6 * i + 1] = P[0][1]; qp[6 * i + 2] = P[0][2];
		qp[6 * i + 3] = P[1][1]; qp[6 * i + 4] = P[1][2]; qp[6 * i + 5] = P[2][2];
	}
	return 1;
}

/* every item's traversal, in chunks through the executor when there is
   one. Each item writes only its own energy share and gradient block, so
   the chunks share nothing; the shares are summed and the blocks gathered
   in item order afterwards: the same result for any executor. */
typedef struct bh_par
{
	bh_ctx const* base;
	bh_ctx* chunk;                  /* one context per chunk */
	double *share, *blocks;
} bh_par;

static qaws_status bh_chunk(void* user, unsigned int chunk, unsigned int begin, unsigned int end)
{
	bh_par* p = (bh_par*)user;
	bh_ctx* c = &p->chunk[chunk];
	unsigned int I, slots = p->base->k * p->base->dim;
	*c = *p->base;
	c->E = 0; c->near = c->far = 0; c->status = 0;
	for (I = begin; I < end && !c->status; I++)
	{
		double e0 = c->E;
		c->blk = p->blocks ? p->blocks + (size_t)slots * I : NULL;
		bh_visit(c, I, 0);
		p->share[I] = c->E - e0;
	}
	return c->status == 0 ? QAWS_STATUS_OK : c->status == 1 ? QAWS_STATUS_DEGENERATE_CURVE : QAWS_STATUS_NUMERICAL_FAILURE;
}

static qaws_status bh_traverse(bh_ctx* c)
{
	enum { GRAIN = 64 };
	unsigned int chunks = qaws_internal_chunk_count(c->nitems, GRAIN), i, a, slots = c->k * c->dim;
	bh_par p;
	qaws_status st;
	p.base = c;
	p.chunk = (bh_ctx*)calloc(chunks + 1, sizeof(bh_ctx));
	p.share = (double*)calloc((size_t)c->nitems + 1, sizeof(double));
	p.blocks = c->grad ? (double*)calloc((size_t)slots * c->nitems + 1, sizeof(double)) : NULL;
	if (!p.chunk || !p.share || (c->grad && !p.blocks))
	{
		free(p.chunk); free(p.share); free(p.blocks);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	st = qaws_internal_parallel(c->nitems, GRAIN, bh_chunk, &p, c->executor);
	if (st == QAWS_STATUS_OK)
	{
		c->E = 0;
		for (i = 0; i < c->nitems; i++) c->E += p.share[i];
		for (i = 0; i < chunks; i++) { c->near += p.chunk[i].near; c->far += p.chunk[i].far; }
		for (i = 0; c->grad && i < c->nitems; i++)
			for (a = 0; a < slots; a++) c->grad[c->dim * c->idx[c->k * i + a / c->dim] + a % c->dim] += p.blocks[(size_t)slots * i + a];
	}
	free(p.chunk); free(p.share); free(p.blocks);
	return st;
}

/* E (and, if grad, grad += its gradient) at x, in double */
static qaws_status bh_run(bh_ctx* c, double const* x, double* grad, double* E, qaws_repulsive_bh_stats* out_stats)
{
	unsigned int i;
	double *lo = NULL, *hi = NULL, *ctr = NULL, *w = NULL, *qp = NULL;
	qaws_status st = QAWS_STATUS_OK;
	lo = (double*)malloc(sizeof(double) * (c->dim * (size_t)c->nitems + 1));
	hi = (double*)malloc(sizeof(double) * (c->dim * (size_t)c->nitems + 1));
	ctr = (double*)malloc(sizeof(double) * (c->dim * (size_t)c->nitems + 1));
	w = (double*)malloc(sizeof(double) * (c->nitems + 1));
	qp = (double*)malloc(sizeof(double) * (6 * (size_t)c->nitems + 1));
	c->t.perm = (unsigned int*)malloc(sizeof(unsigned int) * (c->nitems + 1));
	c->x = x;
	c->grad = grad;
	if (!lo || !hi || !ctr || !w || !qp || !c->t.perm)
	{
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	if (!bh_items(c, lo, hi, ctr, w, qp))
	{
		st = QAWS_STATUS_DEGENERATE_CURVE;
		goto done;
	}
	for (i = 0; i < c->nitems; i++) c->t.perm[i] = i;
	c->t.dim = c->dim; c->t.ilo = lo; c->t.ihi = hi; c->t.ic = ctr;
	if (c->nitems)
	{
		if (bh_build(&c->t, 0, c->nitems) < 0)
		{
			st = QAWS_STATUS_ALLOCATION_FAILURE;
			goto done;
		}
		bh_sums(&c->t, 0, w, qp);
	}
	if ((st = bh_traverse(c)) != QAWS_STATUS_OK)
		goto done;
	*E = c->E;
	if (out_stats) { out_stats->near_pairs = c->near; out_stats->far_pairs = c->far; out_stats->nodes = c->t.nn; }
done:
	free(lo); free(hi); free(ctr); free(w); free(qp); free(c->t.perm); free(c->t.node);
	return st;
}

/* the public form: positions in qaws_scalar, gradient += seed times it */
static qaws_status bh_public(bh_ctx* c, unsigned int nv, qaws_scalar const* positions, qaws_scalar seed, qaws_scalar* gradient, qaws_scalar* out_value,
	qaws_repulsive_bh_stats* out_stats)
{
	unsigned int i;
	double E = 0, *x = (double*)malloc(sizeof(double) * (c->dim * (size_t)nv + 1)), *g = gradient ? (double*)calloc(c->dim * (size_t)nv + 1, sizeof(double)) : NULL;
	qaws_status st;
	if (!x || (gradient && !g))
	{
		free(x); free(g);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (i = 0; i < c->dim * nv; i++) x[i] = (double)positions[i];
	st = bh_run(c, x, g, &E, out_stats);
	if (st == QAWS_STATUS_OK)
	{
		for (i = 0; gradient && i < c->dim * nv; i++) gradient[i] = (qaws_scalar)((double)gradient[i] + (double)seed * g[i]);
		if (out_value) *out_value = (qaws_scalar)E;
	}
	free(x); free(g);
	return st;
}

static void bh_curves_ctx(bh_ctx* c, qaws_repulsive_curves const* curves, double theta)
{
	memset(c, 0, sizeof(*c));
	c->nitems = curves->ne; c->k = 2; c->idx = curves->edges; c->dim = curves->dim; c->theta = theta;
	c->alpha = curves->alpha; c->beta = curves->beta; c->skip_shared = 1; c->executor = curves->executor;
	c->far_g = bh_curve_far_g; c->near_g = bh_curve_near_g;
}

static void bh_mesh_ctx(bh_ctx* c, qaws_repulsive_mesh const* mesh, double theta)
{
	memset(c, 0, sizeof(*c));
	c->nitems = mesh->nf; c->k = 3; c->idx = mesh->tris; c->dim = 3; c->theta = theta;
	c->alpha = mesh->p; c->beta = 2 * mesh->p; c->skip_shared = 0; c->executor = mesh->executor;
	c->far_g = bh_mesh_far_g; c->near_g = bh_mesh_near_g;
}

qaws_status qaws_internal_repulsive_curves_bh(double* E, double* grad, qaws_repulsive_curves const* curves, double const* x, double theta)
{
	bh_ctx c;
	bh_curves_ctx(&c, curves, theta);
	return bh_run(&c, x, grad, E, NULL);
}

qaws_status qaws_internal_repulsive_mesh_bh(double* E, double* grad, qaws_repulsive_mesh const* mesh, double const* x, double theta)
{
	bh_ctx c;
	bh_mesh_ctx(&c, mesh, theta);
	return bh_run(&c, x, grad, E, NULL);
}

qaws_status qaws_repulsive_curves_bh(qaws_scalar* gradient, qaws_scalar* out_value, qaws_repulsive_bh_stats* out_stats, qaws_repulsive_curves const* curves,
	qaws_scalar const* positions, qaws_scalar theta, qaws_scalar seed)
{
	bh_ctx c;
	if (!curves || !positions || !(theta >= 0))
		return QAWS_STATUS_INVALID_ARGUMENT;
	bh_curves_ctx(&c, curves, (double)theta);
	return bh_public(&c, curves->nv, positions, seed, gradient, out_value, out_stats);
}

qaws_status qaws_repulsive_mesh_bh(qaws_scalar* gradient, qaws_scalar* out_value, qaws_repulsive_bh_stats* out_stats, qaws_repulsive_mesh const* mesh,
	qaws_scalar const* positions, qaws_scalar theta, qaws_scalar seed)
{
	bh_ctx c;
	if (!mesh || !positions || !(theta >= 0))
		return QAWS_STATUS_INVALID_ARGUMENT;
	bh_mesh_ctx(&c, mesh, (double)theta);
	return bh_public(&c, mesh->nv, positions, seed, gradient, out_value, out_stats);
}

/* ------------------------------------------------------------------------ */
/*  Hierarchical operators: the metric and the fractional operator         */
/* ------------------------------------------------------------------------ */

/* The metrics' bilinear forms sum, over pairs of items (I, J), a weight
   kB times (g_I - g_J) . (h_I - h_J) and a weight kB0 times (u_I - u_J)
   (v_I - v_J), with g_I = sum_a u(v_a) c_Ia the item's derivative (k
   vertices, gdim components each) and u_I = sum_a u(v_a) / k its average.
   Row I gathers kB (g_I - g_J) and kB0 (u_I - u_J) from every J: the near
   items one by one, a far node's at once through its sums G = sum w_J g_J,
   U = sum w_J u_J with weights per unit w (kB (W g_I - G)). The pairs'
   weights are computed once at x into the near and far lists, so that each
   product is a sweep over them. */

typedef struct rbh_link
{
	unsigned int I, J;              /* J: an item (near) or a node (far) */
	double kB, kB0;
} rbh_link;

struct qaws_internal_rbh_op
{
	bh_tree t;
	unsigned int nitems, k, gdim, nv, dim;
	unsigned int const* idx;        /* items' vertices; NULL: item i is vertex i */
	double *c, *w;                  /* per item: k x gdim coefficients, weight */
	rbh_link *near, *far;
	size_t nnear, nfar, capnear, capfar;
	qaws_internal_rbh_scratch sc;  /* for qaws_internal_rbh_apply */
};

/* the weights of the pair (I, J), or of I and a far node */
typedef int (*rbh_weights)(void const* user, unsigned int I, unsigned int J, bh_node const* far, double* kB, double* kB0);

typedef struct rbh_build
{
	qaws_internal_rbh_op* op;
	double const *ilo, *ihi;
	double theta;
	int skip_shared, status;
	rbh_weights fn;
	void const* user;
	rbh_link *near, *far;            /* this chunk's links, its items in order */
	size_t nnear, nfar, capnear, capfar;
} rbh_build;

static int rbh_push(rbh_link** v, size_t* n, size_t* cap, unsigned int I, unsigned int J, double kB, double kB0)
{
	if (*n == *cap)
	{
		size_t nc = *cap ? 2 * *cap : 1024;
		rbh_link* g = (rbh_link*)realloc(*v, sizeof(rbh_link) * nc);
		if (!g)
			return 0;
		*v = g;
		*cap = nc;
	}
	(*v)[*n].I = I; (*v)[*n].J = J; (*v)[*n].kB = kB; (*v)[*n].kB0 = kB0;
	(*n)++;
	return 1;
}

static unsigned int rbh_vertex(qaws_internal_rbh_op const* op, unsigned int I, unsigned int a)
{
	return op->idx ? op->idx[op->k * I + a] : I;
}

static int rbh_shares(qaws_internal_rbh_op const* op, unsigned int I, unsigned int J)
{
	unsigned int a, b;
	for (a = 0; a < op->k; a++)
		for (b = 0; b < op->k; b++)
			if (rbh_vertex(op, I, a) == rbh_vertex(op, J, b))
				return 1;
	return 0;
}

static void rbh_visit(rbh_build* b, unsigned int I, unsigned int id)
{
	qaws_internal_rbh_op* op = b->op;
	bh_node const* n = &op->t.node[id];
	double const* ilo = b->ilo + op->dim * I;
	double const* ihi = b->ihi + op->dim * I;
	double ri = 0, d, kB, kB0;
	unsigned int q, a;
	if (b->status)
		return;
	for (q = 0; q < op->dim; q++) ri += 0.25 * (ihi[q] - ilo[q]) * (ihi[q] - ilo[q]);
	ri = sqrt(ri);
	d = bh_boxdist(ilo, ihi, n->lo, n->hi, op->dim);
	if (b->theta > 0 && fmax(ri, n->r) <= b->theta * d)
	{
		if (!b->fn(b->user, I, 0, n, &kB, &kB0)) { b->status = 2; return; }
		if (!rbh_push(&b->far, &b->nfar, &b->capfar, I, id, kB, kB0)) b->status = 3;
		return;
	}
	if (n->child[0] >= 0)
	{
		rbh_visit(b, I, (unsigned int)n->child[0]);
		rbh_visit(b, I, (unsigned int)n->child[1]);
		return;
	}
	for (a = n->begin; a < n->end && !b->status; a++)
	{
		unsigned int J = op->t.perm[a];
		if (J == I || (b->skip_shared && rbh_shares(op, I, J)))
			continue;
		if (!b->fn(b->user, I, J, NULL, &kB, &kB0)) { b->status = 2; return; }
		if (!rbh_push(&b->near, &b->nnear, &b->capnear, I, J, kB, kB0)) b->status = 3;
	}
}

static void rbh_scratch_free(qaws_internal_rbh_scratch* sc)
{
	free(sc->gi); free(sc->ui); free(sc->G); free(sc->U); free(sc->accB); free(sc->acc0);
	memset(sc, 0, sizeof(*sc));
}

/* a product's scratch, sized for op (its tree built) */
static int rbh_scratch_make(qaws_internal_rbh_op const* op, qaws_internal_rbh_scratch* sc)
{
	size_t ni = (size_t)op->nitems + 1, nn = (size_t)op->t.nn + 1;
	sc->gi = (double*)malloc(sizeof(double) * (op->gdim * ni + 1));
	sc->ui = (double*)malloc(sizeof(double) * ni);
	sc->accB = (double*)malloc(sizeof(double) * (op->gdim * ni + 1));
	sc->acc0 = (double*)malloc(sizeof(double) * ni);
	sc->G = (double*)malloc(sizeof(double) * (op->gdim * nn + 1));
	sc->U = (double*)malloc(sizeof(double) * nn);
	if (sc->gi && sc->ui && sc->accB && sc->acc0 && sc->G && sc->U)
		return 1;
	rbh_scratch_free(sc);
	return 0;
}

qaws_internal_rbh_scratch* qaws_internal_rbh_scratch_create(qaws_internal_rbh_op const* op)
{
	qaws_internal_rbh_scratch* sc = (qaws_internal_rbh_scratch*)calloc(1, sizeof(*sc));
	if (sc && !rbh_scratch_make(op, sc))
	{
		free(sc);
		return NULL;
	}
	return sc;
}

void qaws_internal_rbh_scratch_destroy(qaws_internal_rbh_scratch* sc)
{
	if (!sc)
		return;
	rbh_scratch_free(sc);
	free(sc);
}

void qaws_internal_rbh_destroy(qaws_internal_rbh_op* op)
{
	if (!op)
		return;
	free(op->t.node); free(op->t.perm); free(op->c); free(op->w); free(op->near); free(op->far);
	rbh_scratch_free(&op->sc);
	free(op);
}

/* the tree and the lists, the items' c and w set; lo, hi, ctr the items'
   boxes and centers, qp their weighted projectors */
/* the lists of items [begin, end), one chunk */
typedef struct rbh_par
{
	rbh_build const* base;
	rbh_build* chunk;
} rbh_par;

static qaws_status rbh_chunk(void* user, unsigned int chunk, unsigned int begin, unsigned int end)
{
	rbh_par* p = (rbh_par*)user;
	rbh_build* b = &p->chunk[chunk];
	unsigned int I;
	*b = *p->base;
	for (I = begin; I < end && !b->status; I++)
		rbh_visit(b, I, 0);
	return b->status == 3 ? QAWS_STATUS_ALLOCATION_FAILURE : b->status ? QAWS_STATUS_NUMERICAL_FAILURE : QAWS_STATUS_OK;
}

static qaws_status rbh_finish(qaws_internal_rbh_op* op, double const* lo, double const* hi, double const* ctr, double const* qp, double theta,
	int skip_shared, rbh_weights fn, void const* user, qaws_batch_executor const* executor)
{
	enum { GRAIN = 64 };
	rbh_build b;
	rbh_par p;
	unsigned int i, chunks = qaws_internal_chunk_count(op->nitems, GRAIN);
	qaws_status st;
	size_t ni = (size_t)op->nitems + 1;
	op->t.perm = (unsigned int*)malloc(sizeof(unsigned int) * ni);
	if (!op->t.perm)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; i < op->nitems; i++) op->t.perm[i] = i;
	op->t.dim = op->dim; op->t.ilo = lo; op->t.ihi = hi; op->t.ic = ctr;
	if (op->nitems)
	{
		if (bh_build(&op->t, 0, op->nitems) < 0)
			return QAWS_STATUS_ALLOCATION_FAILURE;
		bh_sums(&op->t, 0, op->w, qp);
	}
	if (!rbh_scratch_make(op, &op->sc))
		return QAWS_STATUS_ALLOCATION_FAILURE;
	memset(&b, 0, sizeof(b));
	b.op = op; b.ilo = lo; b.ihi = hi; b.theta = theta; b.skip_shared = skip_shared; b.fn = fn; b.user = user;
	/* the items' lists in chunks through the executor, then joined in chunk
	   (so item) order: the same lists for any executor */
	p.base = &b;
	p.chunk = (rbh_build*)calloc(chunks + 1, sizeof(rbh_build));
	if (!p.chunk)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	st = qaws_internal_parallel(op->nitems, GRAIN, rbh_chunk, &p, executor);
	for (i = 0; st == QAWS_STATUS_OK && i < chunks; i++) { op->nnear += p.chunk[i].nnear; op->nfar += p.chunk[i].nfar; }
	if (st == QAWS_STATUS_OK)
	{
		op->near = (rbh_link*)malloc(sizeof(rbh_link) * (op->nnear + 1));
		op->far = (rbh_link*)malloc(sizeof(rbh_link) * (op->nfar + 1));
		op->capnear = op->nnear; op->capfar = op->nfar;
		if (!op->near || !op->far)
			st = QAWS_STATUS_ALLOCATION_FAILURE;
	}
	if (st == QAWS_STATUS_OK)
	{
		size_t an = 0, af = 0;
		for (i = 0; i < chunks; i++)
		{
			if (p.chunk[i].nnear) memcpy(op->near + an, p.chunk[i].near, sizeof(rbh_link) * p.chunk[i].nnear);
			if (p.chunk[i].nfar) memcpy(op->far + af, p.chunk[i].far, sizeof(rbh_link) * p.chunk[i].nfar);
			an += p.chunk[i].nnear; af += p.chunk[i].nfar;
		}
	}
	for (i = 0; i < chunks; i++) { free(p.chunk[i].near); free(p.chunk[i].far); }
	free(p.chunk);
	op->t.ilo = op->t.ihi = op->t.ic = NULL;
	return st;
}

/* y = A u with the caller's scratch: products with one scratch each may run at once */
void qaws_internal_rbh_apply_with(double* y, qaws_internal_rbh_op const* op, qaws_internal_rbh_scratch* sc, double const* u)
{
	unsigned int I, a, q, gd = op->gdim, k = op->k, id;
	size_t l;
	for (I = 0; I < op->nitems; I++)
	{
		double s = 0;
		for (q = 0; q < gd; q++) sc->gi[gd * I + q] = 0;
		for (a = 0; a < k; a++)
		{
			double v = u[rbh_vertex(op, I, a)];
			s += v;
			for (q = 0; q < gd; q++) sc->gi[gd * I + q] += v * op->c[(k * I + a) * gd + q];
		}
		sc->ui[I] = s / k;
	}
	/* the nodes' sums, children (later ids) first */
	for (id = op->t.nn; id-- > 0;)
	{
		bh_node const* n = &op->t.node[id];
		double* G = sc->G + (size_t)gd * id;
		for (q = 0; q < gd; q++) G[q] = 0;
		sc->U[id] = 0;
		if (n->child[0] >= 0)
		{
			unsigned int c0 = (unsigned int)n->child[0], c1 = (unsigned int)n->child[1];
			for (q = 0; q < gd; q++) G[q] = sc->G[(size_t)gd * c0 + q] + sc->G[(size_t)gd * c1 + q];
			sc->U[id] = sc->U[c0] + sc->U[c1];
		}
		else
			for (a = n->begin; a < n->end; a++)
			{
				unsigned int J = op->t.perm[a];
				for (q = 0; q < gd; q++) G[q] += op->w[J] * sc->gi[gd * J + q];
				sc->U[id] += op->w[J] * sc->ui[J];
			}
	}
	memset(sc->accB, 0, sizeof(double) * (gd * (size_t)op->nitems + 1));
	memset(sc->acc0, 0, sizeof(double) * ((size_t)op->nitems + 1));
	for (l = 0; l < op->nnear; l++)
	{
		rbh_link const* e = &op->near[l];
		for (q = 0; q < gd; q++) sc->accB[gd * e->I + q] += e->kB * (sc->gi[gd * e->I + q] - sc->gi[gd * e->J + q]);
		sc->acc0[e->I] += e->kB0 * (sc->ui[e->I] - sc->ui[e->J]);
	}
	for (l = 0; l < op->nfar; l++)
	{
		rbh_link const* e = &op->far[l];
		double W = op->t.node[e->J].W;
		for (q = 0; q < gd; q++) sc->accB[gd * e->I + q] += e->kB * (W * sc->gi[gd * e->I + q] - sc->G[(size_t)gd * e->J + q]);
		sc->acc0[e->I] += e->kB0 * (W * sc->ui[e->I] - sc->U[e->J]);
	}
	memset(y, 0, sizeof(double) * op->nv);
	for (I = 0; I < op->nitems; I++)
		for (a = 0; a < k; a++)
		{
			double s = sc->acc0[I] / k;
			for (q = 0; q < gd; q++) s += op->c[(k * I + a) * gd + q] * sc->accB[gd * I + q];
			y[rbh_vertex(op, I, a)] += s;
		}
}

void qaws_internal_rbh_apply(double* y, qaws_internal_rbh_op* op, double const* u)
{
	qaws_internal_rbh_apply_with(y, op, &op->sc, u);
}

size_t qaws_internal_rbh_links(size_t* out_far, qaws_internal_rbh_op const* op)
{
	if (out_far) *out_far = op->nfar;
	return op->nnear;
}

/* a new operator over nitems items of k vertices, gdim derivative
   components; lo, hi, ctr, qp allocated for the caller to fill */
static qaws_internal_rbh_op* rbh_alloc(unsigned int nitems, unsigned int k, unsigned int gdim, unsigned int dim, unsigned int nv, unsigned int const* idx,
	double** lo, double** hi, double** ctr, double** qp)
{
	qaws_internal_rbh_op* op = (qaws_internal_rbh_op*)calloc(1, sizeof(*op));
	size_t ni = (size_t)nitems + 1;
	*lo = (double*)malloc(sizeof(double) * dim * ni);
	*hi = (double*)malloc(sizeof(double) * dim * ni);
	*ctr = (double*)malloc(sizeof(double) * dim * ni);
	*qp = (double*)calloc(6 * ni, sizeof(double));
	if (op)
	{
		op->nitems = nitems; op->k = k; op->gdim = gdim; op->dim = dim; op->nv = nv; op->idx = idx;
		op->c = (double*)calloc((size_t)k * gdim * ni + 1, sizeof(double));
		op->w = (double*)malloc(sizeof(double) * ni);
	}
	if (!op || !op->c || !op->w || !*lo || !*hi || !*ctr || !*qp)
	{
		qaws_internal_rbh_destroy(op);
		free(*lo); free(*hi); free(*ctr); free(*qp);
		*lo = *hi = *ctr = *qp = NULL;
		return NULL;
	}
	return op;
}

/* ---- curves: Repulsive Curves, Sec. 5.2 ---- */

typedef struct rbh_curves
{
	qaws_repulsive_curves const* c;
	double const* x;
	double ex;                      /* 2 sigma + 1 */
} rbh_curves;

/* the ordered pair (P, Q) as in the dense assembly: w = 1/4 l_P l_Q sum_ij
   r^-ex (symmetric), w0 = 1/4 l_P l_Q sum_ij |T_P x d|^2 / r^4 r^-ex */
static int rbh_curve_pair(rbh_curves const* u, unsigned int const* P, unsigned int const* Q, double* w, double* w0)
{
	unsigned int dim = u->c->dim, i, j, q;
	double eP[3] = { 0, 0, 0 }, lP = 0, lQ = 0;
	for (q = 0; q < dim; q++)
	{
		double f = u->x[dim * Q[1] + q] - u->x[dim * Q[0] + q];
		eP[q] = u->x[dim * P[1] + q] - u->x[dim * P[0] + q];
		lP += eP[q] * eP[q]; lQ += f * f;
	}
	lP = sqrt(lP); lQ = sqrt(lQ);
	*w = *w0 = 0;
	for (i = 0; i < 2; i++)
		for (j = 0; j < 2; j++)
		{
			double r2 = 0, td = 0, inv;
			for (q = 0; q < dim; q++)
			{
				double d = u->x[dim * P[i] + q] - u->x[dim * Q[j] + q];
				r2 += d * d; td += eP[q] * d / lP;
			}
			if (!(r2 > 0))
				return 0;
			inv = bh_pow(r2, -0.5 * u->ex);
			*w += inv;
			*w0 += fmax(r2 - td * td, 0) / (r2 * r2) * inv;
		}
	*w *= 0.25 * lP * lQ;
	*w0 *= 0.25 * lP * lQ;
	return 1;
}

static int rbh_curve_weights(void const* user, unsigned int I, unsigned int J, bh_node const* far, double* kB, double* kB0)
{
	static unsigned int const ix[3][3] = { { 0, 1, 2 }, { 1, 3, 4 }, { 2, 4, 5 } };
	rbh_curves const* u = (rbh_curves const*)user;
	unsigned int dim = u->c->dim, i, q, r, s;
	unsigned int const* P = u->c->edges + 2 * I;
	double e[3] = { 0, 0, 0 }, l = 0, sB = 0, s0 = 0;
	if (!far)
	{
		double w, w0a, w0b;
		if (!rbh_curve_pair(u, P, u->c->edges + 2 * J, &w, &w0a) || !rbh_curve_pair(u, u->c->edges + 2 * J, P, &w, &w0b))
			return 0;
		*kB = 2 * w;
		*kB0 = w0a + w0b;
		return 1;
	}
	/* the node at its centroid X, I's endpoints kept; per unit length:
	   kB = l_I sum_i r_i^-ex, kB0 = l_I / 2 sum_i r_i^-(ex + 4) (|T_I x d_i|^2 + d_i^T Q d_i / W) */
	for (q = 0; q < dim; q++) { e[q] = u->x[dim * P[1] + q] - u->x[dim * P[0] + q]; l += e[q] * e[q]; }
	l = sqrt(l);
	for (i = 0; i < 2; i++)
	{
		double d[3] = { 0, 0, 0 }, r2 = 0, td = 0, dQd = 0, inv;
		for (q = 0; q < dim; q++)
		{
			d[q] = u->x[dim * P[i] + q] - far->S[q] / far->W;
			r2 += d[q] * d[q]; td += e[q] * d[q] / l;
		}
		if (!(r2 > 0))
			return 0;
		for (r = 0; r < dim; r++)
			for (s = 0; s < dim; s++) dQd += d[r] * d[s] * far->Q[ix[r][s]];
		inv = bh_pow(r2, -0.5 * u->ex);
		sB += inv;
		s0 += inv / (r2 * r2) * (fmax(r2 - td * td, 0) + fmax(dQd, 0) / far->W);
	}
	*kB = l * sB;
	*kB0 = 0.5 * l * s0;
	return 1;
}

qaws_status qaws_internal_rbh_curves_metric(qaws_internal_rbh_op** out, qaws_repulsive_curves const* c, double const* x, double theta)
{
	static unsigned int const ix[6][2] = { { 0, 0 }, { 0, 1 }, { 0, 2 }, { 1, 1 }, { 1, 2 }, { 2, 2 } };
	double *lo, *hi, *ctr, *qp;
	unsigned int dim = c->dim, i, a, q, r;
	rbh_curves u;
	qaws_status st = QAWS_STATUS_OK;
	qaws_internal_rbh_op* op = rbh_alloc(c->ne, 2, dim, dim, c->nv, c->edges, &lo, &hi, &ctr, &qp);
	*out = NULL;
	if (!op)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; i < c->ne && st == QAWS_STATUS_OK; i++)
	{
		unsigned int const* v = c->edges + 2 * i;
		double e[3] = { 0, 0, 0 }, l = 0;
		for (q = 0; q < dim; q++)
		{
			double y0 = x[dim * v[0] + q], y1 = x[dim * v[1] + q];
			lo[dim * i + q] = fmin(y0, y1); hi[dim * i + q] = fmax(y0, y1); ctr[dim * i + q] = 0.5 * (y0 + y1);
			e[q] = y1 - y0; l += e[q] * e[q];
		}
		l = sqrt(l);
		if (!(l > 0))
		{
			st = QAWS_STATUS_DEGENERATE_CURVE;
			break;
		}
		op->w[i] = l;
		/* g_I = (u_1 - u_0) / l T: coefficients -T / l, T / l */
		for (a = 0; a < 2; a++)
			for (q = 0; q < dim; q++) op->c[(2 * i + a) * dim + q] = (a ? 1 : -1) * e[q] / (l * l);
		/* l (I - T T^T) */
		for (r = 0; r < 6; r++)
		{
			unsigned int p0 = ix[r][0], p1 = ix[r][1];
			qp[6 * i + r] = p0 < dim && p1 < dim ? (p0 == p1 ? l : 0) - e[p0] * e[p1] / l : 0;
		}
	}
	u.c = c; u.x = x; u.ex = 2 * ((c->beta - 1) / c->alpha - 1) + 1;
	if (st == QAWS_STATUS_OK)
		st = rbh_finish(op, lo, hi, ctr, qp, theta, 1, rbh_curve_weights, &u, c->executor);
	free(lo); free(hi); free(ctr); free(qp);
	if (st != QAWS_STATUS_OK)
	{
		qaws_internal_rbh_destroy(op);
		return st;
	}
	*out = op;
	return QAWS_STATUS_OK;
}

/* ---- meshes: Repulsive Surfaces, Sec. 3.2 (Eqs. 12, 13) ---- */

typedef struct rbh_mesh
{
	double const *a, *N, *X;        /* per face: area, unit normal, barycenter */
	double ex;                      /* 2 sigma + 2 */
} rbh_mesh;

static int rbh_mesh_weights(void const* user, unsigned int S, unsigned int T, bh_node const* far, double* kB, double* kB0)
{
	static unsigned int const ix[3][3] = { { 0, 1, 2 }, { 1, 3, 4 }, { 2, 4, 5 } };
	rbh_mesh const* u = (rbh_mesh const*)user;
	double d[3], r2 = 0, ns = 0, inv, dQd = 0;
	unsigned int q, r;
	for (q = 0; q < 3; q++)
	{
		d[q] = u->X[3 * S + q] - (far ? far->S[q] / far->W : u->X[3 * T + q]);
		r2 += d[q] * d[q]; ns += u->N[3 * S + q] * d[q];
	}
	if (!(r2 > 0))
		return 0;
	inv = bh_pow(r2, -0.5 * u->ex);
	if (!far)
	{
		/* w = a_S a_T r^-ex: kB = 2 w, kB0 = w ((N_S . d)^2 + (N_T . d)^2) / r^4 */
		double nt = 0, w = u->a[S] * u->a[T] * inv;
		for (q = 0; q < 3; q++) nt += u->N[3 * T + q] * d[q];
		*kB = 2 * w;
		*kB0 = w * (ns * ns + nt * nt) / (r2 * r2);
		return 1;
	}
	/* per unit area: kB = 2 a_S r^-ex, kB0 = a_S r^-(ex + 4) ((N_S . d)^2 + d^T Q d / W) */
	for (q = 0; q < 3; q++)
		for (r = 0; r < 3; r++) dQd += d[q] * d[r] * far->Q[ix[q][r]];
	*kB = 2 * u->a[S] * inv;
	*kB0 = u->a[S] * inv / (r2 * r2) * (ns * ns + fmax(dQd, 0) / far->W);
	return 1;
}

qaws_status qaws_internal_rbh_mesh_metric(qaws_internal_rbh_op** out, qaws_repulsive_mesh const* m, double const* x, double theta)
{
	static unsigned int const ix[6][2] = { { 0, 0 }, { 0, 1 }, { 0, 2 }, { 1, 1 }, { 1, 2 }, { 2, 2 } };
	double *lo, *hi, *ctr, *qp, *fa, *fN;
	unsigned int i, c, q, r;
	rbh_mesh u;
	qaws_status st = QAWS_STATUS_OK;
	qaws_internal_rbh_op* op = rbh_alloc(m->nf, 3, 3, 3, m->nv, m->tris, &lo, &hi, &ctr, &qp);
	*out = NULL;
	if (!op)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	fa = (double*)malloc(sizeof(double) * (m->nf + 1));
	fN = (double*)malloc(sizeof(double) * (3 * (size_t)m->nf + 1));
	if (!fa || !fN)
		st = QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; i < m->nf && st == QAWS_STATUS_OK; i++)
	{
		unsigned int const* t = m->tris + 3 * i;
		double e1[3], e2[3], n[3], l;
		for (q = 0; q < 3; q++)
		{
			double y0 = x[3 * t[0] + q], y1 = x[3 * t[1] + q], y2 = x[3 * t[2] + q];
			lo[3 * i + q] = fmin(y0, fmin(y1, y2)); hi[3 * i + q] = fmax(y0, fmax(y1, y2)); ctr[3 * i + q] = (y0 + y1 + y2) / 3;
			e1[q] = y1 - y0; e2[q] = y2 - y0;
		}
		n[0] = e1[1] * e2[2] - e1[2] * e2[1]; n[1] = e1[2] * e2[0] - e1[0] * e2[2]; n[2] = e1[0] * e2[1] - e1[1] * e2[0];
		l = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
		if (!(l > 0))
		{
			st = QAWS_STATUS_DEGENERATE_CURVE;
			break;
		}
		fa[i] = op->w[i] = 0.5 * l;
		for (q = 0; q < 3; q++) fN[3 * i + q] = n[q] / l;
		/* the hat gradients N x e / (2 a), e the opposite edge in turn */
		for (c = 0; c < 3; c++)
		{
			double e[3], g[3];
			for (q = 0; q < 3; q++) e[q] = x[3 * t[(c + 2) % 3] + q] - x[3 * t[(c + 1) % 3] + q];
			g[0] = fN[3 * i + 1] * e[2] - fN[3 * i + 2] * e[1]; g[1] = fN[3 * i + 2] * e[0] - fN[3 * i] * e[2]; g[2] = fN[3 * i] * e[1] - fN[3 * i + 1] * e[0];
			for (q = 0; q < 3; q++) op->c[(3 * i + c) * 3 + q] = g[q] / l;
		}
		/* a N N^T */
		for (r = 0; r < 6; r++) qp[6 * i + r] = fa[i] * fN[3 * i + ix[r][0]] * fN[3 * i + ix[r][1]];
	}
	u.a = fa; u.N = fN; u.X = ctr; u.ex = 2 * (1 - 2 / m->p) + 2;
	if (st == QAWS_STATUS_OK)
		st = rbh_finish(op, lo, hi, ctr, qp, theta, 0, rbh_mesh_weights, &u, m->executor);
	free(lo); free(hi); free(ctr); free(qp); free(fa); free(fN);
	if (st != QAWS_STATUS_OK)
	{
		qaws_internal_rbh_destroy(op);
		return st;
	}
	*out = op;
	return QAWS_STATUS_OK;
}

/* ---- the fractional operator by vertex quadrature ---- */

typedef struct rbh_points
{
	double const *p, *m;
	unsigned int dim;
	double ex;
} rbh_points;

/* 2 m_i m_j r^-ex (the dense form's weight); per unit mass for a far node */
static int rbh_point_weights(void const* user, unsigned int i, unsigned int j, bh_node const* far, double* kB, double* kB0)
{
	rbh_points const* u = (rbh_points const*)user;
	double r2 = 0;
	unsigned int q;
	for (q = 0; q < u->dim; q++)
	{
		double d = u->p[u->dim * i + q] - (far ? far->S[q] / far->W : u->p[u->dim * j + q]);
		r2 += d * d;
	}
	if (!(r2 > 0))
		return 0;
	*kB = 0;
	*kB0 = 2 * u->m[i] * (far ? 1 : u->m[j]) * bh_pow(r2, -0.5 * u->ex);
	return 1;
}

qaws_status qaws_internal_rbh_fractional(qaws_internal_rbh_op** out, double const* p, double const* mass, unsigned int nv, unsigned int dim, double ex,
	double theta, qaws_batch_executor const* executor)
{
	double *lo, *hi, *ctr, *qp;
	unsigned int i, q;
	rbh_points u;
	qaws_status st;
	qaws_internal_rbh_op* op = rbh_alloc(nv, 1, 0, dim, nv, NULL, &lo, &hi, &ctr, &qp);
	*out = NULL;
	if (!op)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; i < nv; i++)
	{
		for (q = 0; q < dim; q++) lo[dim * i + q] = hi[dim * i + q] = ctr[dim * i + q] = p[dim * i + q];
		op->w[i] = mass[i];
	}
	u.p = p; u.m = mass; u.dim = dim; u.ex = ex;
	st = rbh_finish(op, lo, hi, ctr, qp, theta, 0, rbh_point_weights, &u, executor);
	free(lo); free(hi); free(ctr); free(qp);
	if (st != QAWS_STATUS_OK)
	{
		qaws_internal_rbh_destroy(op);
		return st;
	}
	*out = op;
	return QAWS_STATUS_OK;
}

/* ------------------------------------------------------------------------ */
/*  Overlapping boxes: the collision checks' broad phase                    */
/* ------------------------------------------------------------------------ */

/* Every pair i < j of overlapping boxes that accept allows, visited once,
   by a self traversal of the median-split tree above: a pair of nodes
   whose boxes do not overlap is dropped whole. O(n log n + pairs) for
   boxes spread along curves and surfaces, where a uniform grid's cells
   (sized for boxes spread through a volume) are far too coarse. */
typedef struct rbh_pairs
{
	bh_tree t;
	qaws_bp_box const* box;
	int (*accept)(void*, unsigned int, unsigned int);
	qaws_status (*visit)(void*, unsigned int, unsigned int);
	void* user;
	qaws_status st;
} rbh_pairs;

static int rbh_overlap(double const* alo, double const* ahi, double const* blo, double const* bhi, unsigned int dim)
{
	unsigned int q;
	for (q = 0; q < dim; q++)
		if (alo[q] > bhi[q] || blo[q] > ahi[q])
			return 0;
	return 1;
}

static void rbh_pairs_visit(rbh_pairs* P, int a, int b)
{
	bh_node const *na = &P->t.node[a], *nb = &P->t.node[b];
	if (P->st != QAWS_STATUS_OK || !rbh_overlap(na->lo, na->hi, nb->lo, nb->hi, P->t.dim))
		return;
	if (na->child[0] < 0 && nb->child[0] < 0)
	{
		unsigned int i, j;
		for (i = na->begin; i < na->end && P->st == QAWS_STATUS_OK; i++)
			for (j = a == b ? i + 1 : nb->begin; j < nb->end && P->st == QAWS_STATUS_OK; j++)
			{
				unsigned int u = P->t.perm[i], v = P->t.perm[j], lo = u < v ? u : v, hi = u < v ? v : u;
				if (!rbh_overlap(P->box[u].lo, P->box[u].hi, P->box[v].lo, P->box[v].hi, P->t.dim))
					continue;
				if (P->accept && !P->accept(P->user, lo, hi))
					continue;
				P->st = P->visit(P->user, lo, hi);
			}
		return;
	}
	if (a == b)
	{
		rbh_pairs_visit(P, na->child[0], na->child[0]);
		rbh_pairs_visit(P, na->child[1], na->child[1]);
		rbh_pairs_visit(P, na->child[0], na->child[1]);
		return;
	}
	/* (open the larger, or the inner one) */
	if (nb->child[0] < 0 || (na->child[0] >= 0 && na->r >= nb->r))
	{
		rbh_pairs_visit(P, na->child[0], b);
		rbh_pairs_visit(P, na->child[1], b);
	}
	else
	{
		rbh_pairs_visit(P, a, nb->child[0]);
		rbh_pairs_visit(P, a, nb->child[1]);
	}
}

qaws_status qaws_internal_rbh_box_pairs(qaws_bp_box const* box, unsigned int n, unsigned int dim, int (*accept)(void*, unsigned int, unsigned int),
	qaws_status (*visit)(void*, unsigned int, unsigned int), void* user)
{
	rbh_pairs P;
	double *lo = (double*)malloc(sizeof(double) * dim * ((size_t)n + 1)), *hi = (double*)malloc(sizeof(double) * dim * ((size_t)n + 1));
	double* ctr = (double*)malloc(sizeof(double) * dim * ((size_t)n + 1));
	unsigned int i, q;
	memset(&P, 0, sizeof(P));
	P.t.perm = (unsigned int*)malloc(sizeof(unsigned int) * ((size_t)n + 1));
	P.st = QAWS_STATUS_OK;
	if (!lo || !hi || !ctr || !P.t.perm)
		P.st = QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; P.st == QAWS_STATUS_OK && i < n; i++)
	{
		P.t.perm[i] = i;
		for (q = 0; q < dim; q++)
		{
			lo[dim * i + q] = box[i].lo[q]; hi[dim * i + q] = box[i].hi[q];
			ctr[dim * i + q] = 0.5 * (box[i].lo[q] + box[i].hi[q]);
		}
	}
	P.t.dim = dim; P.t.ilo = lo; P.t.ihi = hi; P.t.ic = ctr;
	P.box = box; P.accept = accept; P.visit = visit; P.user = user;
	if (P.st == QAWS_STATUS_OK && n > 1)
	{
		if (bh_build(&P.t, 0, n) < 0)
			P.st = QAWS_STATUS_ALLOCATION_FAILURE;
		else
			rbh_pairs_visit(&P, 0, 0);
	}
	free(lo); free(hi); free(ctr); free(P.t.perm); free(P.t.node);
	return P.st;
}
