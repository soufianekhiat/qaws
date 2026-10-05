#ifndef QAWS_INTERNAL_BIJET_H
#define QAWS_INTERNAL_BIJET_H

#include "../qaws_platform.h"
#include "../core/qaws_dual_core.h"

/*
 * Bivariate jets: all partials d^(a+b)/du^a dv^b up to total order K,
 * stored at index (a+b)(a+b+1)/2 + b (the qaws_surface_jet layout).
 * Every entry is a second-order dual number, so derivatives along a
 * parameter/coordinate tangent ride along with the spatial partials.
 *
 * Products follow the bivariate Leibniz rule
 *   (f g)_{ab} = sum_{i<=a, j<=b} C(a,i) C(b,j) f_{ij} g_{a-i,b-j}
 * and quotients / square roots invert it.
 */

#define QAWS_BIJET_MAX_ORDER 3
#define QAWS_BIJET_COUNT 10

static qaws_scalar const g_bijet_binom[4][4] = {
	{ 1, 0, 0, 0 },
	{ 1, 1, 0, 0 },
	{ 1, 2, 1, 0 },
	{ 1, 3, 3, 1 }
};

static unsigned char const g_bijet_a[QAWS_BIJET_COUNT] = { 0, 1, 0, 2, 1, 0, 3, 2, 1, 0 };
static unsigned char const g_bijet_b[QAWS_BIJET_COUNT] = { 0, 0, 1, 0, 1, 2, 0, 1, 2, 3 };

QAWS_INLINE unsigned int qaws_bijet_index(unsigned int a, unsigned int b)
{
	unsigned int t = a + b;
	return t * (t + 1) / 2 + b;
}

QAWS_INLINE unsigned int qaws_bijet_count(unsigned int order)
{
	return (order + 1) * (order + 2) / 2;
}

/* out = f x g */
QAWS_INLINE void qaws_bijet_cross(qaws_dual3 const* f, qaws_dual3 const* g, unsigned int order, qaws_dual3* out)
{
	unsigned int n = qaws_bijet_count(order), k, i, j;
	for (k = 0; k < n; k++)
	{
		unsigned int a = g_bijet_a[k], b = g_bijet_b[k];
		qaws_dual3 acc = qaws_dual3_const(qaws_v3_zero());
		for (i = 0; i <= a; i++)
			for (j = 0; j <= b; j++)
				acc = qaws_dual3_add(acc, qaws_dual3_mul_const(
					qaws_dual3_cross(f[qaws_bijet_index(i, j)], g[qaws_bijet_index(a - i, b - j)]),
					g_bijet_binom[a][i] * g_bijet_binom[b][j]));
		out[k] = acc;
	}
}

/* out = f . g */
QAWS_INLINE void qaws_bijet_dot(qaws_dual3 const* f, qaws_dual3 const* g, unsigned int order, qaws_dual1* out)
{
	unsigned int n = qaws_bijet_count(order), k, i, j;
	for (k = 0; k < n; k++)
	{
		unsigned int a = g_bijet_a[k], b = g_bijet_b[k];
		qaws_dual1 acc = qaws_dual1_const(QAWS_ZERO);
		for (i = 0; i <= a; i++)
			for (j = 0; j <= b; j++)
			{
				qaws_dual1 term = qaws_dual3_dot(f[qaws_bijet_index(i, j)], g[qaws_bijet_index(a - i, b - j)]);
				qaws_scalar c = g_bijet_binom[a][i] * g_bijet_binom[b][j];
				acc = qaws_dual1_add(acc, qaws_dual1_make(term.v * c, term.t * c, term.tt * c));
			}
		out[k] = acc;
	}
}

/* l = sqrt(q):  2 l_00 l_ab = q_ab - sum_{(i,j) not in {(0,0),(a,b)}} C C l_ij l_{a-i,b-j} */
QAWS_INLINE void qaws_bijet_sqrt(qaws_dual1 const* q, unsigned int order, qaws_dual1* out)
{
	unsigned int n = qaws_bijet_count(order), k, i, j;
	out[0] = qaws_dual1_sqrt(q[0]);
	for (k = 1; k < n; k++)
	{
		unsigned int a = g_bijet_a[k], b = g_bijet_b[k];
		qaws_dual1 acc = q[k];
		for (i = 0; i <= a; i++)
			for (j = 0; j <= b; j++)
			{
				qaws_dual1 term;
				qaws_scalar c;
				if ((i == 0 && j == 0) || (i == a && j == b))
					continue;
				c = g_bijet_binom[a][i] * g_bijet_binom[b][j];
				term = qaws_dual1_mul(out[qaws_bijet_index(i, j)], out[qaws_bijet_index(a - i, b - j)]);
				acc = qaws_dual1_sub(acc, qaws_dual1_make(term.v * c, term.t * c, term.tt * c));
			}
		out[k] = qaws_dual1_div(acc, qaws_dual1_make(QAWS_LITERAL(2.0) * out[0].v,
			QAWS_LITERAL(2.0) * out[0].t, QAWS_LITERAL(2.0) * out[0].tt));
	}
}

/* out = f / l:  l_00 out_ab = f_ab - sum_{(i,j) != (0,0)} C C l_ij out_{a-i,b-j} */
QAWS_INLINE void qaws_bijet_div(qaws_dual3 const* f, qaws_dual1 const* l, unsigned int order, qaws_dual3* out)
{
	unsigned int n = qaws_bijet_count(order), k, i, j;
	for (k = 0; k < n; k++)
	{
		unsigned int a = g_bijet_a[k], b = g_bijet_b[k];
		qaws_dual3 acc = f[k];
		for (i = 0; i <= a; i++)
			for (j = 0; j <= b; j++)
			{
				if (i == 0 && j == 0)
					continue;
				acc = qaws_dual3_sub(acc, qaws_dual3_mul_const(
					qaws_dual3_scale(out[qaws_bijet_index(a - i, b - j)], l[qaws_bijet_index(i, j)]),
					g_bijet_binom[a][i] * g_bijet_binom[b][j]));
			}
		out[k] = qaws_dual3_div(acc, l[0]);
	}
}

/* Shifted views of a jet: (d/du f)_{ab} = f_{a+1,b}, (d/dv f)_{ab} = f_{a,b+1}.
   The input must hold order + 1. */
QAWS_INLINE void qaws_bijet_shift_u(qaws_dual3 const* f, unsigned int order, qaws_dual3* out)
{
	unsigned int n = qaws_bijet_count(order), k;
	for (k = 0; k < n; k++)
		out[k] = f[qaws_bijet_index(g_bijet_a[k] + 1, g_bijet_b[k])];
}

QAWS_INLINE void qaws_bijet_shift_v(qaws_dual3 const* f, unsigned int order, qaws_dual3* out)
{
	unsigned int n = qaws_bijet_count(order), k;
	for (k = 0; k < n; k++)
		out[k] = f[qaws_bijet_index(g_bijet_a[k], g_bijet_b[k] + 1)];
}

#endif /* QAWS_INTERNAL_BIJET_H */
