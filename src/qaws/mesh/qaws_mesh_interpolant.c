/*
 * The interpolants of mesh references as cubic Bezier triangles: PN
 * (Vlachos et al. 2001, creases by their tangents), PHONG (Boubekeur and
 * Alexa 2008) and NAGATA (Nagata 2005), quadratic and degree-elevated,
 * LINEAR (the flat triangle degree-elevated) and G1 (Walton and Meek 1996:
 * PN's boundary curves, a quartic Gregory interior), evaluated with their
 * derivatives, and the curvature of a normal field.
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

#define MI_PHONG_ALPHA 0.75   /* Boubekeur and Alexa's shape factor: the blend of the projections against the flat triangle */

/* Nagata's curvature coefficient of the edge from x0 (normal n0) to x1 (n1): the quadratic x0 + (d - c) t + c t^2 (d = x1 - x0)
   leaves each end in its tangent plane, (d - c) . n0 = 0 and (d + c) . n1 = 0, c taken in the span of n0 and n1; parallel
   normals (no turn to follow) give 0 */
static void mi_nagata_c(double* c, double const* x0, double const* x1, double const* n0, double const* n1)
{
	double d[3], dn0 = 0, dn1 = 0, cs = 0, det, al, be;
	int k;
	for (k = 0; k < 3; k++)
	{
		d[k] = x1[k] - x0[k];
		dn0 += d[k] * n0[k];
		dn1 += d[k] * n1[k];
		cs += n0[k] * n1[k];
	}
	det = 1 - cs * cs;
	if (det < 1e-12)
	{
		c[0] = c[1] = c[2] = 0;
		return;
	}
	/* (al + be cs = d . n0, al cs + be = -d . n1) */
	al = (dn0 + cs * dn1) / det;
	be = (-dn1 - cs * dn0) / det;
	for (k = 0; k < 3; k++)
		c[k] = al * n0[k] + be * n1[k];
}

/* The quadratic interpolants as quadratic Bezier triangles of corners P_i and an edge point Q_ij each, degree-elevated to the
   cubic net. PHONG: p*(w) = (1 - a) p(w) + a sum_i w_i pi_i(p(w)), pi_i the projection on corner i's tangent plane, whose
   edge points are (1 - a) (P_i + P_j) / 2 + a (pi_i(P_j) + pi_j(P_i)) / 2. NAGATA: the edge curve x0 + (d - c) t + c t^2
   (mi_nagata_c), its point (P_i + P_j) / 2 - c / 2. An edge along a crease (tangents given) stays straight. Normals: the
   corners' interpolated linearly. */
static void mi_build_quadratic(mi_tri* t, qaws_interpolant kind, double const* const* P, double const* const* N, double const* const* crease)
{
	static int const near_i[6] = { 0, 1, 1, 2, 2, 0 }, edge_of[6] = { 0, 0, 1, 1, 2, 2 };
	double Q[3][3];
	int e, q, k;
	for (e = 0; e < 3; e++)
	{
		int i = e, j = (e + 1) % 3;
		double pij[3], pji[3], di = 0, dj = 0;
		for (k = 0; k < 3; k++)
		{
			di += (P[j][k] - P[i][k]) * N[i][k];
			dj += (P[i][k] - P[j][k]) * N[j][k];
		}
		double cc[3];
		mi_nagata_c(cc, P[i], P[j], N[i], N[j]);
		for (k = 0; k < 3; k++)
		{
			pij[k] = P[j][k] - di * N[i][k];
			pji[k] = P[i][k] - dj * N[j][k];
			Q[e][k] = crease && crease[e] ? 0.5 * (P[i][k] + P[j][k])
			        : kind == QAWS_INTERPOLANT_NAGATA ? 0.5 * (P[i][k] + P[j][k]) - 0.5 * cc[k]
			                                          : (1 - MI_PHONG_ALPHA) * 0.5 * (P[i][k] + P[j][k]) + MI_PHONG_ALPHA * 0.5 * (pij[k] + pji[k]);
		}
	}
	for (k = 0; k < 3; k++)
	{
		t->p[0][k] = P[0][k];
		t->p[1][k] = P[1][k];
		t->p[2][k] = P[2][k];
		t->p[9][k] = (Q[0][k] + Q[1][k] + Q[2][k]) / 3;
	}
	/* (degree elevation: an edge's cubic points a third and two thirds between its corner and its quadratic point) */
	for (q = 0; q < 6; q++)
		for (k = 0; k < 3; k++)
			t->p[3 + q][k] = (P[near_i[q]][k] + 2 * Q[edge_of[q]][k]) / 3;
	for (k = 0; k < 3; k++)
	{
		t->n[0][k] = N[0][k];
		t->n[1][k] = N[1][k];
		t->n[2][k] = N[2][k];
		for (q = 0; q < 3; q++)
			t->n[3 + q][k] = 0.5 * (N[q][k] + N[(q + 1) % 3][k]);
	}
}

/* PN and LINEAR */
static void mi_build_pn(mi_tri* t, qaws_interpolant kind, double const* const* P, double const* const* N, double const* const* crease)
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

static void mi_cross(double* o, double const* a, double const* b)
{
	o[0] = a[1] * b[2] - a[2] * b[1];
	o[1] = a[2] * b[0] - a[0] * b[2];
	o[2] = a[0] * b[1] - a[1] * b[0];
}

/* the cubic boundary's control point of exponents f */
static double const* mi_c3p(mi_tri const* t, int const* f)
{
	int l;
	for (l = 0; l < 9; l++)
		if (mi_e3[l][0] == f[0] && mi_e3[l][1] == f[1] && mi_e3[l][2] == f[2])
			break;
	return t->p[l];
}

/* the quartic boundary point of exponents f (one of them 0): the cubic degree-elevated, sum_c f_c / 4 C_(f - e_c) */
static void mi_q4(double* o, mi_tri const* t, int const* f)
{
	int c, k;
	o[0] = o[1] = o[2] = 0;
	for (c = 0; c < 3; c++)
		if (f[c])
		{
			int g[3];
			double const* C;
			g[0] = f[0]; g[1] = f[1]; g[2] = f[2];
			g[c]--;
			C = mi_c3p(t, g);
			for (k = 0; k < 3; k++)
				o[k] += 0.25 * f[c] * C[k];
		}
}

/* G1 (Walton and Meek 1996, the cross tangents of Chiyokura and Kimura): PN's cubic boundary curves c, degree-elevated to a
   quartic triangle; along each edge the cross derivative D(t) = k(t) c'(t) + h(t) m(t), k and h linear, m(t) the linear blend
   of N_i x c'(0) and N_j x c'(1), the same from both sides: the tangent planes span c' and m on both, G1. D's ends are the
   corners' (the other edges' tangents), its two middle coefficients set the interior row: two candidates per interior point
   (one per edge at its corner), blended rationally in the evaluation (Gregory). */
static void mi_build_g1(mi_tri* t, double const* const* P, double const* const* N, double const* const* crease)
{
	int e, k, s;
	mi_build_pn(t, QAWS_INTERPOLANT_PN, P, N, crease);
	for (e = 0; e < 3; e++)
	{
		int i = e, j = (e + 1) % 3, o = (e + 2) % 3, f[3];
		double E[4][3], C[4][3], R0[3], R3[3], d[3][3], m0[3], m1[3], D0[3], D3[3], k0, k1, h0, h1, l0, l1;
		for (s = 0; s < 4; s++)
		{
			double const* c;
			f[i] = 4 - s; f[j] = s; f[o] = 0;
			mi_q4(E[s], t, f);
			f[i] = 3 - s;
			c = mi_c3p(t, f);
			for (k = 0; k < 3; k++)
				C[s][k] = c[k];
		}
		f[i] = 3; f[j] = 0; f[o] = 1;
		mi_q4(R0, t, f);
		f[i] = 0; f[j] = 3;
		mi_q4(R3, t, f);
		for (k = 0; k < 3; k++)
		{
			D0[k] = 4 * (R0[k] - E[0][k]);
			D3[k] = 4 * (R3[k] - E[3][k]);
			for (s = 0; s < 3; s++)
				d[s][k] = 3 * (C[s + 1][k] - C[s][k]);
		}
		mi_cross(m0, N[i], d[0]);
		mi_unit(m0);
		mi_cross(m1, N[j], d[2]);
		mi_unit(m1);
		l0 = mi_dot(d[0], d[0]);
		l1 = mi_dot(d[2], d[2]);
		k0 = l0 > 0 ? mi_dot(D0, d[0]) / l0 : 0;
		k1 = l1 > 0 ? mi_dot(D3, d[2]) / l1 : 0;
		h0 = mi_dot(D0, m0);
		h1 = mi_dot(D3, m1);
		for (k = 0; k < 3; k++)
		{
			/* (k c' cubic from linear times quadratic, h m quadratic degree-elevated) */
			double q0 = h0 * m0[k], q1 = 0.5 * (h0 * m1[k] + h1 * m0[k]), q2 = h1 * m1[k];
			double D1 = (2 * k0 * d[1][k] + k1 * d[0][k]) / 3 + (q0 + 2 * q1) / 3, D2 = (k0 * d[2][k] + 2 * k1 * d[1][k]) / 3 + (2 * q1 + q2) / 3;
			t->n[2 * i][k] = E[1][k] + 0.25 * D1;       /* near corner i, from the edge to i + 1 */
			t->n[2 * j + 1][k] = E[2][k] + 0.25 * D2;   /* near corner j, from the edge to j + 2 */
		}
	}
	for (k = 0; k < 3; k++)
		t->p[9][k] = (t->n[0][k] + t->n[1][k] + t->n[2][k] + t->n[3][k] + t->n[4][k] + t->n[5][k]) / 6;
	t->gregory = 1;
}

void mi_build(mi_tri* t, qaws_interpolant kind, double const* const* P, double const* const* N, double const* const* crease)
{
	t->gregory = 0;
	if (kind == QAWS_INTERPOLANT_PHONG || kind == QAWS_INTERPOLANT_NAGATA)
		mi_build_quadratic(t, kind, P, N, crease);
	else if (kind == QAWS_INTERPOLANT_G1)
		mi_build_g1(t, P, N, crease);
	else
		mi_build_pn(t, kind, P, N, crease);
}

static double mi_pow(double x, int e)
{
	return e < 0 ? 0 : e == 0 ? 1 : e == 1 ? x : e == 2 ? x * x : e == 3 ? x * x * x : x * x * x * x;
}

/* b0^i b1^j b2^k, 0 when a power is negative */
static double mi_mono(double const* b, int i, int j, int k)
{
	if (i < 0 || j < 0 || k < 0)
		return 0;
	return mi_pow(b[0], i) * mi_pow(b[1], j) * mi_pow(b[2], k);
}

/* the second derivatives of b0^i b1^j b2^q in (a, b) and its value and first ones: w[6] = (f, fa, fb, faa, fab, fbb) */
static void mi_bern_derivatives(double* w, double const* bb, int i, int j, int q)
{
	w[0] = mi_mono(bb, i, j, q);
	w[1] = -i * mi_mono(bb, i - 1, j, q) + j * mi_mono(bb, i, j - 1, q);
	w[2] = -i * mi_mono(bb, i - 1, j, q) + q * mi_mono(bb, i, j, q - 1);
	w[3] = i * (i - 1) * mi_mono(bb, i - 2, j, q) - 2 * i * j * mi_mono(bb, i - 1, j - 1, q) + j * (j - 1) * mi_mono(bb, i, j - 2, q);
	w[4] = i * (i - 1) * mi_mono(bb, i - 2, j, q) - i * q * mi_mono(bb, i - 1, j, q - 1) - i * j * mi_mono(bb, i - 1, j - 1, q)
		+ j * q * mi_mono(bb, i, j - 1, q - 1);
	w[5] = i * (i - 1) * mi_mono(bb, i - 2, j, q) - 2 * i * q * mi_mono(bb, i - 1, j, q - 1) + q * (q - 1) * mi_mono(bb, i, j, q - 2);
}

static int const mi_e4[15][3] = { { 4, 0, 0 }, { 0, 4, 0 }, { 0, 0, 4 }, { 3, 1, 0 }, { 2, 2, 0 }, { 1, 3, 0 }, { 0, 3, 1 }, { 0, 2, 2 },
	{ 0, 1, 3 }, { 1, 0, 3 }, { 2, 0, 2 }, { 3, 0, 1 }, { 2, 1, 1 }, { 1, 2, 1 }, { 1, 1, 2 } };
static double const mi_c4[15] = { 1, 1, 1, 4, 6, 4, 4, 6, 4, 4, 6, 4, 12, 12, 12 };

/* the quartic Gregory triangle (mi_build_g1): an interior point near corner v is (u G_2v + w G_2v+1) / (u + w), u = b_(v+1),
   w = b_(v+2) (each candidate alone on its edge), differentiated with it; the normal the surface's, xa x xb */
static void mi_eval_gregory(mi_tri const* t, double a, double b, double* x, double* xa, double* xb, double* xaa, double* xab, double* xbb,
	double* n, double* na, double* nb)
{
	static int const db[3][2] = { { -1, -1 }, { 1, 0 }, { 0, 1 } };   /* the barycentrics' derivatives in (a, b) */
	double bb[3] = { 1 - a - b, a, b }, v[6][3];
	int l, k;
	memset(v, 0, sizeof v);
	for (l = 0; l < 15; l++)
	{
		int const* f = mi_e4[l];
		double c = mi_c4[l], w[6], D[6][3];
		memset(D, 0, sizeof D);
		if (f[0] && f[1] && f[2])
		{
			int cv = f[0] == 2 ? 0 : f[1] == 2 ? 1 : 2, a1 = (cv + 1) % 3, a2 = (cv + 2) % 3;
			double u = bb[a1], s = bb[a1] + bb[a2], r = 0.5, ra = 0, rb = 0, raa = 0, rab = 0, rbb = 0;
			double ua = db[a1][0], ub = db[a1][1], sa = db[a1][0] + db[a2][0], sb = db[a1][1] + db[a2][1];
			if (s > 1e-14)
			{
				r = u / s;
				ra = ua / s - u * sa / (s * s);
				rb = ub / s - u * sb / (s * s);
				raa = -2 * ua * sa / (s * s) + 2 * u * sa * sa / (s * s * s);
				rab = -(ua * sb + ub * sa) / (s * s) + 2 * u * sa * sb / (s * s * s);
				rbb = -2 * ub * sb / (s * s) + 2 * u * sb * sb / (s * s * s);
			}
			for (k = 0; k < 3; k++)
			{
				double g = t->n[2 * cv][k] - t->n[2 * cv + 1][k];
				D[0][k] = t->n[2 * cv + 1][k] + r * g;
				D[1][k] = ra * g;
				D[2][k] = rb * g;
				D[3][k] = raa * g;
				D[4][k] = rab * g;
				D[5][k] = rbb * g;
			}
		}
		else
			mi_q4(D[0], t, f);
		mi_bern_derivatives(w, bb, f[0], f[1], f[2]);
		for (k = 0; k < 3; k++)
		{
			v[0][k] += c * w[0] * D[0][k];
			v[1][k] += c * (w[1] * D[0][k] + w[0] * D[1][k]);
			v[2][k] += c * (w[2] * D[0][k] + w[0] * D[2][k]);
			v[3][k] += c * (w[3] * D[0][k] + 2 * w[1] * D[1][k] + w[0] * D[3][k]);
			v[4][k] += c * (w[4] * D[0][k] + w[1] * D[2][k] + w[2] * D[1][k] + w[0] * D[4][k]);
			v[5][k] += c * (w[5] * D[0][k] + 2 * w[2] * D[2][k] + w[0] * D[5][k]);
		}
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
		double m[3], ma[3], mb[3], p1[3], p2[3], u[3], len, ua, ub;
		mi_cross(m, v[1], v[2]);
		mi_cross(p1, v[3], v[2]);
		mi_cross(p2, v[1], v[4]);
		for (k = 0; k < 3; k++)
			ma[k] = p1[k] + p2[k];
		mi_cross(p1, v[4], v[2]);
		mi_cross(p2, v[1], v[5]);
		for (k = 0; k < 3; k++)
		{
			mb[k] = p1[k] + p2[k];
			u[k] = m[k];
		}
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

void mi_eval(mi_tri const* t, double a, double b, double* x, double* xa, double* xb, double* xaa, double* xab, double* xbb,
	double* n, double* na, double* nb)
{
	double bb[3] = { 1 - a - b, a, b }, v[6][3];
	int l, k, o;
	if (t->gregory)
	{
		mi_eval_gregory(t, a, b, x, xa, xb, xaa, xab, xbb, n, na, nb);
		return;
	}
	memset(v, 0, sizeof v);
	for (l = 0; l < 10; l++)
	{
		int i = mi_e3[l][0], j = mi_e3[l][1], q = mi_e3[l][2];
		double c = mi_c3[l], w[6];
		mi_bern_derivatives(w, bb, i, j, q);
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
