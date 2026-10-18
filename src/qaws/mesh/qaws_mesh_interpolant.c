/*
 * The interpolants of mesh references as cubic Bezier triangles: PN
 * (Vlachos et al. 2001, creases by their tangents) and LINEAR (the flat
 * triangle degree-elevated), evaluated with their derivatives, and the
 * curvature of a normal field.
 */

#include "qaws_mesh_internal.h"
#include <math.h>
#include <string.h>

/* powers of (b0, b1, b2) and Bernstein factors */
static int const mi_e3[10][3] = { { 3, 0, 0 }, { 0, 3, 0 }, { 0, 0, 3 }, { 2, 1, 0 }, { 1, 2, 0 }, { 0, 2, 1 }, { 0, 1, 2 }, { 1, 0, 2 },
	{ 2, 0, 1 }, { 1, 1, 1 } };
static double const mi_c3[10] = { 1, 1, 1, 3, 3, 3, 3, 3, 3, 6 };
static int const mi_e2[6][3] = { { 2, 0, 0 }, { 0, 2, 0 }, { 0, 0, 2 }, { 1, 1, 0 }, { 0, 1, 1 }, { 1, 0, 1 } };
static double const mi_c2[6] = { 1, 1, 1, 2, 2, 2 };   /* (a quadratic Bezier: the paper writes the mixed terms without the 2, which halves the normals' turning rate at the corners) */

static double mi_dot(double const* a, double const* b)
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static double mi_unit(double* v)
{
	double l = sqrt(mi_dot(v, v));
	if (l > 0)
	{
		v[0] /= l; v[1] /= l; v[2] /= l;
	}
	return l;
}

/* the control point near P[i] on the edge to P[j]: the chord projected on the
   tangent plane (or along the crease tangent T), a third of the way */
static void mi_edge_point(double* out, double const* Pi, double const* Pj, double const* Ni, double const* T, int linear)
{
	double d[3], w;
	int k;
	for (k = 0; k < 3; k++)
		d[k] = Pj[k] - Pi[k];
	if (!linear)
	{
		if (T)
		{
			w = mi_dot(d, T);
			for (k = 0; k < 3; k++)
				d[k] = w * T[k];
		}
		else
		{
			w = mi_dot(d, Ni);
			for (k = 0; k < 3; k++)
				d[k] -= w * Ni[k];
		}
	}
	for (k = 0; k < 3; k++)
		out[k] = Pi[k] + d[k] / 3;
}

void mi_build(mi_tri* t, qaws_interpolant kind, double const* const* P, double const* const* N, double const* const* crease)
{
	/* edge point (index in the net) near corner i toward corner j: 210 (0 to 1), 120 (1 to 0), 021 (1 to 2), 012 (2 to 1), 102 (2 to 0), 201 (0 to 2) */
	static int const near_i[6] = { 0, 1, 1, 2, 2, 0 }, near_j[6] = { 1, 0, 2, 1, 0, 2 }, edge_of[6] = { 0, 0, 1, 1, 2, 2 };
	int linear = kind == QAWS_INTERPOLANT_LINEAR, q, k;
	double E[3] = { 0, 0, 0 }, V[3];
	for (k = 0; k < 3; k++)
	{
		t->p[0][k] = P[0][k];
		t->p[1][k] = P[1][k];
		t->p[2][k] = P[2][k];
		V[k] = (P[0][k] + P[1][k] + P[2][k]) / 3;
	}
	for (q = 0; q < 6; q++)
	{
		int i = near_i[q], j = near_j[q], e = edge_of[q];
		/* the crease tangent at this end: crease[e] holds the edge's start (corner e) then its end */
		double const* T = crease && crease[e] ? crease[e] + (i == e ? 0 : 3) : NULL;
		mi_edge_point(t->p[3 + q], P[i], P[j], N[i], T, linear);
		for (k = 0; k < 3; k++)
			E[k] += t->p[3 + q][k] / 6;
	}
	for (k = 0; k < 3; k++)
		t->p[9][k] = linear ? V[k] : E[k] + (E[k] - V[k]) / 2;
	/* normals: the corners', and mid-edge ones reflected across the edge's mid-plane (PN) or averaged (LINEAR) */
	for (k = 0; k < 3; k++)
	{
		t->n[0][k] = N[0][k];
		t->n[1][k] = N[1][k];
		t->n[2][k] = N[2][k];
	}
	for (q = 0; q < 3; q++)
	{
		int i = q, j = (q + 1) % 3;
		double d[3], s[3], dd, v = 0;
		for (k = 0; k < 3; k++)
		{
			d[k] = P[j][k] - P[i][k];
			s[k] = N[i][k] + N[j][k];
		}
		dd = mi_dot(d, d);
		if (!linear && dd > 0)
			v = 2 * mi_dot(d, s) / dd;
		for (k = 0; k < 3; k++)
			t->n[3 + q][k] = linear ? 0.5 * s[k] : s[k] - v * d[k];   /* (LINEAR: the linear field degree-elevated) */
		if (!linear)
			mi_unit(t->n[3 + q]);
	}
}

static double mi_pow(double x, int e)
{
	return e < 0 ? 0 : e == 0 ? 1 : e == 1 ? x : e == 2 ? x * x : x * x * x;
}

/* b0^i b1^j b2^k, 0 when a power is negative */
static double mi_mono(double const* b, int i, int j, int k)
{
	if (i < 0 || j < 0 || k < 0)
		return 0;
	return mi_pow(b[0], i) * mi_pow(b[1], j) * mi_pow(b[2], k);
}

void mi_eval(mi_tri const* t, double a, double b, double* x, double* xa, double* xb, double* xaa, double* xab, double* xbb,
	double* n, double* na, double* nb)
{
	double bb[3] = { 1 - a - b, a, b }, v[6][3];
	int l, k, o;
	memset(v, 0, sizeof v);
	for (l = 0; l < 10; l++)
	{
		int i = mi_e3[l][0], j = mi_e3[l][1], q = mi_e3[l][2];
		double c = mi_c3[l], w[6];
		w[0] = mi_mono(bb, i, j, q);
		w[1] = -i * mi_mono(bb, i - 1, j, q) + j * mi_mono(bb, i, j - 1, q);
		w[2] = -i * mi_mono(bb, i - 1, j, q) + q * mi_mono(bb, i, j, q - 1);
		w[3] = i * (i - 1) * mi_mono(bb, i - 2, j, q) - 2 * i * j * mi_mono(bb, i - 1, j - 1, q) + j * (j - 1) * mi_mono(bb, i, j - 2, q);
		w[4] = i * (i - 1) * mi_mono(bb, i - 2, j, q) - i * q * mi_mono(bb, i - 1, j, q - 1) - i * j * mi_mono(bb, i - 1, j - 1, q)
			+ j * q * mi_mono(bb, i, j - 1, q - 1);
		w[5] = i * (i - 1) * mi_mono(bb, i - 2, j, q) - 2 * i * q * mi_mono(bb, i - 1, j, q - 1) + q * (q - 1) * mi_mono(bb, i, j, q - 2);
		for (o = 0; o < 6; o++)
			for (k = 0; k < 3; k++)
				v[o][k] += c * w[o] * t->p[l][k];
	}
	for (k = 0; k < 3; k++)
	{
		if (x) x[k] = v[0][k];
		if (xa) xa[k] = v[1][k];
		if (xb) xb[k] = v[2][k];
		if (xaa) xaa[k] = v[3][k];
		if (xab) xab[k] = v[4][k];
		if (xbb) xbb[k] = v[5][k];
	}
	if (n || na || nb)
	{
		double m[3] = { 0, 0, 0 }, ma[3] = { 0, 0, 0 }, mb[3] = { 0, 0, 0 }, len, u[3], ua, ub;
		for (l = 0; l < 6; l++)
		{
			int i = mi_e2[l][0], j = mi_e2[l][1], q = mi_e2[l][2];
			double c = mi_c2[l], w0 = mi_mono(bb, i, j, q), w1 = -i * mi_mono(bb, i - 1, j, q) + j * mi_mono(bb, i, j - 1, q),
				w2 = -i * mi_mono(bb, i - 1, j, q) + q * mi_mono(bb, i, j, q - 1);
			for (k = 0; k < 3; k++)
			{
				m[k] += c * w0 * t->n[l][k];
				ma[k] += c * w1 * t->n[l][k];
				mb[k] += c * w2 * t->n[l][k];
			}
		}
		/* N = m / |m|; its derivative: the part of m' across N, over |m| */
		for (k = 0; k < 3; k++)
			u[k] = m[k];
		len = mi_unit(u);
		ua = mi_dot(u, ma);
		ub = mi_dot(u, mb);
		for (k = 0; k < 3; k++)
		{
			if (n) n[k] = u[k];
			if (na) na[k] = len > 0 ? (ma[k] - ua * u[k]) / len : 0;
			if (nb) nb[k] = len > 0 ? (mb[k] - ub * u[k]) / len : 0;
		}
	}
}

void mi_curvature(double const* xa, double const* xb, double const* n, double const* na, double const* nb, double* k1, double* k2,
	double* d1, double* d2)
{
	double E = mi_dot(xa, xa), F = mi_dot(xa, xb), G = mi_dot(xb, xb);
	double L = mi_dot(na, xa), M = 0.5 * (mi_dot(na, xb) + mi_dot(nb, xa)), N = mi_dot(nb, xb);
	double det = E * G - F * F, H, K, disc, al, be, r1, r2;
	int k;
	if (det <= 0)
	{
		*k1 = *k2 = 0;
		for (k = 0; k < 3; k++)
			d1[k] = d2[k] = 0;
		return;
	}
	H = (L * G + N * E - 2 * M * F) / (2 * det);
	K = (L * N - M * M) / det;
	disc = H * H - K;
	disc = disc > 0 ? sqrt(disc) : 0;
	*k1 = H + disc;
	*k2 = H - disc;
	/* (II - k1 I) (al, be) = 0, from its longer row */
	r1 = fabs(M - *k1 * F) + fabs(*k1 * E - L);
	r2 = fabs(*k1 * G - N) + fabs(M - *k1 * F);
	if (r1 >= r2)
	{
		al = M - *k1 * F;
		be = *k1 * E - L;
	}
	else
	{
		al = *k1 * G - N;
		be = M - *k1 * F;
	}
	if (fabs(al) + fabs(be) < 1e-14 * (fabs(L) + fabs(M) + fabs(N) + 1e-300))
	{
		al = 1;   /* (an umbilic: any direction) */
		be = 0;
	}
	for (k = 0; k < 3; k++)
		d1[k] = al * xa[k] + be * xb[k];
	{
		double w = mi_dot(d1, n);
		for (k = 0; k < 3; k++)
			d1[k] -= w * n[k];
	}
	mi_unit(d1);
	d2[0] = n[1] * d1[2] - n[2] * d1[1];
	d2[1] = n[2] * d1[0] - n[0] * d1[2];
	d2[2] = n[0] * d1[1] - n[1] * d1[0];
	mi_unit(d2);
}
