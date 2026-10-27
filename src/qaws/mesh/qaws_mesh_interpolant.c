/*
 * The interpolants of mesh references as cubic Bezier triangles: PN
 * (Vlachos et al. 2001, creases by their tangents), PHONG (Boubekeur and
 * Alexa 2008) and NAGATA (Nagata 2005), quadratic and degree-elevated,
 * LINEAR (the flat triangle degree-elevated), G1 (Walton and Meek 1996:
 * PN's boundary curves, a quartic Gregory interior) and SUBDIVISION (the
 * Loop limit surface, creases by Hoppe et al.'s rules), evaluated with their
 * derivatives, and the curvature of a normal field.
 */

#include "qaws_mesh_internal.h"
#include "../internal/qaws_internal_memory.h"
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

static void mi_build_loop(mi_tri* t);

void mi_build(mi_tri* t, qaws_interpolant kind, double const* const* P, double const* const* N, double const* const* crease)
{
	t->gregory = 0;
	if (kind == QAWS_INTERPOLANT_SUBDIVISION)
		mi_build_loop(t);
	else if (kind == QAWS_INTERPOLANT_PHONG || kind == QAWS_INTERPOLANT_NAGATA)
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

static void mi_surface_normal(double v[6][3], double* n, double* na, double* nb);

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
	mi_surface_normal(v, n, na, nb);
}

/* the unit normal xa x xb and its derivatives from v = (x, xa, xb, xaa, xab, xbb); outputs may be NULL */
static void mi_surface_normal(double v[6][3], double* n, double* na, double* nb)
{
	int k;
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

/* ---- SUBDIVISION: the Loop limit surface ----------------------------------
   A triangle whose corners are regular (valence 6, no crease edge at them) is a quartic box spline over its 12 neighbours
   (Stam 1998: mi_box, monomials u^i v^j w^k times 12, in mi_e4's order). Elsewhere the neighbourhood is subdivided, only the
   child holding the point kept with the triangles touching its corners, until that child is regular; at an extraordinary
   vertex or a crease the limit points of a child 2^-52 small, interpolated. Rules: Loop's smooth masks (also at a dart, one
   crease edge), the B-spline's on a crease or the border (edge midpoint, vertex 1/8 3/4 1/8, limit 1/6 4/6 1/6), a corner
   (three or more crease edges, a feature corner on one, a vertex of one face) fixed (Hoppe et al. 1994): those of
   qaws_mesh_3d_subdivide. */

static int const mi_box[12][15] = { { 1, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0 },
	{ 1, 1, 0, 6, 12, 6, 2, 0, 0, 0, 0, 2, 6, 6, 0 }, { 6, 1, 1, 24, 24, 8, 6, 12, 6, 8, 24, 24, 60, 36, 36 },
	{ 1, 0, 1, 2, 0, 0, 0, 0, 2, 6, 12, 6, 6, 0, 6 }, { 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
	{ 1, 6, 1, 8, 24, 24, 24, 24, 8, 6, 12, 6, 36, 60, 36 }, { 1, 1, 6, 6, 12, 6, 8, 24, 24, 24, 24, 8, 36, 36, 60 },
	{ 0, 0, 1, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0 }, { 0, 1, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, 1, 1, 0, 0, 2, 6, 12, 6, 2, 0, 0, 0, 6, 6 }, { 0, 0, 1, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0 } };

typedef struct ml_net
{
	int nv, nt, nc;
	double x[MI_LOOP_MAX][3];
	int t[MI_LOOP_MAX][3];
	int c[MI_LOOP_MAX][2];
	int kp[MI_LOOP_MAX], kq[MI_LOOP_MAX];   /* (while subdividing: each vertex's parent key, a vertex (p, p) or an edge (p < q)) */
	int fx[MI_LOOP_MAX];                    /* 1 a feature corner, 2 in one face only (a vertex's children keep them) */
} ml_net;

static int ml_crease(ml_net const* n, int a, int b)
{
	int i;
	for (i = 0; i < n->nc; i++)
		if ((n->c[i][0] == a && n->c[i][1] == b) || (n->c[i][0] == b && n->c[i][1] == a))
			return 1;
	return 0;
}

/* v's neighbours (nb, nn), the triangles round it, its crease neighbours (cr, ncr) */
static void ml_star(ml_net const* n, int v, int* nb, int* nn, int* ntri, int* cr, int* ncr)
{
	int i, k, q;
	*nn = *ntri = *ncr = 0;
	for (i = 0; i < n->nt; i++)
		if (n->t[i][0] == v || n->t[i][1] == v || n->t[i][2] == v)
		{
			(*ntri)++;
			for (k = 0; k < 3; k++)
			{
				int w = n->t[i][k];
				if (w == v)
					continue;
				for (q = 0; q < *nn && nb[q] != w; q++)
					;
				if (q == *nn && *nn < MI_LOOP_MAX)
					nb[(*nn)++] = w;
			}
		}
	for (i = 0; i < n->nc; i++)
		if ((n->c[i][0] == v || n->c[i][1] == v) && *ncr < MI_LOOP_MAX)
			cr[(*ncr)++] = n->c[i][0] == v ? n->c[i][1] : n->c[i][0];
}

/* v's next position (lim 0) or its limit point (lim 1) */
static void ml_vertex(ml_net const* n, int v, int lim, double* out)
{
	int nb[MI_LOOP_MAX], cr[MI_LOOP_MAX], nn, ntri, ncr, i, k;
	ml_star(n, v, nb, &nn, &ntri, cr, &ncr);
	if (ncr >= 3 || (n->fx[v] & 2) || ((n->fx[v] & 1) && ncr >= 1))
		for (k = 0; k < 3; k++)
			out[k] = n->x[v][k];
	else if (ncr == 2)
		for (k = 0; k < 3; k++)
			out[k] = lim ? (n->x[cr[0]][k] + 4 * n->x[v][k] + n->x[cr[1]][k]) / 6 : 0.125 * (n->x[cr[0]][k] + n->x[cr[1]][k]) + 0.75 * n->x[v][k];
	else
	{
		double c = 0.375 + 0.25 * cos(2 * 3.14159265358979323846 / nn), be = (0.625 - c * c) / nn, w = lim ? 1 / (0.375 / be + nn) : be;
		for (k = 0; k < 3; k++)
		{
			double s = 0;
			for (i = 0; i < nn; i++)
				s += n->x[nb[i]][k];
			out[k] = (1 - nn * w) * n->x[v][k] + w * s;
		}
	}
}

/* the edge (a, b)'s new point: a crease's midpoint, else 3/8 3/8 and 1/8 for each opposite vertex */
static void ml_edge(ml_net const* n, int a, int b, double* out)
{
	int o[2], no = 0, i, k;
	if (!ml_crease(n, a, b))
		for (i = 0; i < n->nt && no < 2; i++)
		{
			int const* t = n->t[i];
			if ((t[0] == a || t[1] == a || t[2] == a) && (t[0] == b || t[1] == b || t[2] == b))
				o[no++] = t[0] != a && t[0] != b ? t[0] : t[1] != a && t[1] != b ? t[1] : t[2];
		}
	for (k = 0; k < 3; k++)
		out[k] = no == 2 ? 0.375 * (n->x[a][k] + n->x[b][k]) + 0.125 * (n->x[o[0]][k] + n->x[o[1]][k]) : 0.5 * (n->x[a][k] + n->x[b][k]);
}

/* the third vertex of a triangle holding a and b but not c, or -1 */
static int ml_opp(ml_net const* n, int a, int b, int c)
{
	int i;
	for (i = 0; i < n->nt; i++)
	{
		int const* t = n->t[i];
		if ((t[0] == a || t[1] == a || t[2] == a) && (t[0] == b || t[1] == b || t[2] == b) && t[0] != c && t[1] != c && t[2] != c)
			return t[0] != a && t[0] != b ? t[0] : t[1] != a && t[1] != b ? t[1] : t[2];
	}
	return -1;
}

/* the key (p, q)'s vertex in m, added (its position from n) when new; -1 when full */
static int ml_key(ml_net* m, ml_net const* n, int p, int q)
{
	int i;
	if (p > q)
	{
		i = p; p = q; q = i;
	}
	for (i = 0; i < m->nv; i++)
		if (m->kp[i] == p && m->kq[i] == q)
			return i;
	if (m->nv == MI_LOOP_MAX)
		return -1;
	m->kp[m->nv] = p;
	m->kq[m->nv] = q;
	m->fx[m->nv] = p == q ? n->fx[p] : 0;
	if (p == q)
		ml_vertex(n, p, 0, m->x[m->nv]);
	else
		ml_edge(n, p, q, m->x[m->nv]);
	return m->nv++;
}

/* a smooth vertex's limit normal from Loop's tangent masks, sum cos (2 pi i / N) p_i and sum sin (2 pi i / N) p_i round its
   ring; 0 for a vertex on a crease, the border or a corner */
static int ml_vertex_normal(ml_net const* n, int v, double* nrm)
{
	int nb[MI_LOOP_MAX], cr[MI_LOOP_MAX], ring[MI_LOOP_MAX], nn, ntri, ncr, i, k, prev, w;
	double t1[3] = { 0, 0, 0 }, t2[3] = { 0, 0, 0 };
	ml_star(n, v, nb, &nn, &ntri, cr, &ncr);
	if (ncr || n->fx[v] || ntri != nn || nn < 3)
		return 0;
	/* (the ring in order: from one neighbour across the triangles round v) */
	prev = -1;
	w = nb[0];
	for (i = 0; i < nn; i++)
	{
		ring[i] = w;
		k = ml_opp(n, v, w, prev < 0 ? -2 : prev);
		if (prev < 0)
		{
			/* (the first step: either way round) */
			int j;
			for (j = 0; j < n->nt; j++)
			{
				int const* t = n->t[j];
				int at = t[0] == v ? 0 : t[1] == v ? 1 : t[2] == v ? 2 : -1;
				if (at >= 0 && t[(at + 1) % 3] == w)
				{
					k = t[(at + 2) % 3];
					break;
				}
			}
		}
		prev = w;
		w = k;
		if (w < 0)
			return 0;
	}
	for (i = 0; i < nn; i++)
		for (k = 0; k < 3; k++)
		{
			t1[k] += cos(2 * 3.14159265358979323846 * i / nn) * n->x[ring[i]][k];
			t2[k] += sin(2 * 3.14159265358979323846 * i / nn) * n->x[ring[i]][k];
		}
	mi_cross(nrm, t1, t2);
	return mi_unit(nrm) > 0;
}

/* the limit point at barycentrics u of n's triangle cur and, regular, its derivatives in (a, b) (J: u's derivatives in them);
   0 when the derivatives are not known (at an extraordinary vertex or a crease), 2 when then the normal is (nrm: at a smooth
   extraordinary vertex) */
static int ml_point(ml_net* n, ml_net* m, int* cur, double* u, double J[3][2], double v[6][3], double* nrm)
{
	int depth, k, l;
	memset(v, 0, sizeof(double) * 18);
	for (depth = 0;; depth++)
	{
		int nb[MI_LOOP_MAX], cr[MI_LOOP_MAX], nn, ntri, ncr, reg = 1, i, j, ch, key[3][2], id[12];
		double u2[3], J2[3][2];
		ml_net* s;
		for (i = 0; i < 3 && reg; i++)
		{
			ml_star(n, cur[i], nb, &nn, &ntri, cr, &ncr);
			reg = ncr == 0 && nn == 6 && ntri == 6;
		}
		if (reg)
		{
			/* the box spline: 4, 7, 8 the corners (u, v, w), the others across the edges */
			id[3] = cur[0]; id[6] = cur[1]; id[7] = cur[2];
			id[2] = ml_opp(n, cur[0], cur[1], cur[2]);
			id[4] = ml_opp(n, cur[0], cur[2], cur[1]);
			id[10] = ml_opp(n, cur[1], cur[2], cur[0]);
			id[0] = id[2] < 0 ? -1 : ml_opp(n, cur[0], id[2], cur[1]);
			id[1] = id[4] < 0 ? -1 : ml_opp(n, cur[0], id[4], cur[2]);
			id[5] = id[2] < 0 ? -1 : ml_opp(n, cur[1], id[2], cur[0]);
			id[9] = id[10] < 0 ? -1 : ml_opp(n, cur[1], id[10], cur[2]);
			id[8] = id[4] < 0 ? -1 : ml_opp(n, cur[2], id[4], cur[0]);
			id[11] = id[10] < 0 ? -1 : ml_opp(n, cur[2], id[10], cur[1]);
			for (i = 0; i < 12 && reg; i++)
				reg = id[i] >= 0;
		}
		if (reg)
		{
			double loc[6][3], aa = J[1][0], ab = J[1][1], ba = J[2][0], bb2 = J[2][1];
			memset(loc, 0, sizeof loc);
			for (l = 0; l < 15; l++)
			{
				double w[6];
				mi_bern_derivatives(w, u, mi_e4[l][0], mi_e4[l][1], mi_e4[l][2]);
				for (i = 0; i < 12; i++)
					if (mi_box[i][l])
						for (j = 0; j < 6; j++)
							for (k = 0; k < 3; k++)
								loc[j][k] += mi_box[i][l] / 12.0 * w[j] * n->x[id[i]][k];
			}
			/* (local (al, be) = (u1, u2), linear in (a, b)) */
			for (k = 0; k < 3; k++)
			{
				v[0][k] = loc[0][k];
				v[1][k] = loc[1][k] * aa + loc[2][k] * ba;
				v[2][k] = loc[1][k] * ab + loc[2][k] * bb2;
				v[3][k] = loc[3][k] * aa * aa + 2 * loc[4][k] * aa * ba + loc[5][k] * ba * ba;
				v[4][k] = loc[3][k] * aa * ab + loc[4][k] * (aa * bb2 + ab * ba) + loc[5][k] * ba * bb2;
				v[5][k] = loc[3][k] * ab * ab + 2 * loc[4][k] * ab * bb2 + loc[5][k] * bb2 * bb2;
			}
			return 1;
		}
		if (depth >= 52 || u[0] > 1 - 1e-15 || u[1] > 1 - 1e-15 || u[2] > 1 - 1e-15)
		{
			for (i = 0; i < 3; i++)
			{
				double y[3];
				ml_vertex(n, cur[i], 1, y);
				for (k = 0; k < 3; k++)
					v[0][k] += u[i] * y[k];
			}
			for (i = 0; i < 3; i++)
				if (u[i] > 1 - 1e-15 && ml_vertex_normal(n, cur[i], nrm))
					return 2;
			return 0;
		}
		/* the child holding the point: a corner's (ci, e(ci, ci+1), e(ci, ci+2)) or the middle (e12, e20, e01) */
		ch = u[0] >= 0.5 ? 0 : u[1] >= 0.5 ? 1 : u[2] >= 0.5 ? 2 : 3;
		for (i = 0; i < 3; i++)
			if (ch < 3)
			{
				int c = cur[ch], o = cur[(ch + i) % 3];
				key[i][0] = c;
				key[i][1] = o;
				u2[i] = i ? 2 * u[(ch + i) % 3] : 2 * u[ch] - 1;
				for (k = 0; k < 2; k++)
					J2[i][k] = 2 * J[(ch + i) % 3][k];
			}
			else
			{
				key[i][0] = cur[(i + 1) % 3];
				key[i][1] = cur[(i + 2) % 3];
				u2[i] = 1 - 2 * u[i];
				for (k = 0; k < 2; k++)
					J2[i][k] = -2 * J[i][k];
			}
		for (i = 0; i < 3; i++)
		{
			if (key[i][0] > key[i][1])
			{
				int w = key[i][0];
				key[i][0] = key[i][1];
				key[i][1] = w;
			}
			u2[i] = u2[i] < 0 ? 0 : u2[i] > 1 ? 1 : u2[i];
		}
		/* the children touching the new triangle's corners */
		m->nv = m->nt = m->nc = 0;
		for (i = 0; i < n->nt; i++)
		{
			int const* t = n->t[i];
			int c4[4][3][2], q;
			for (q = 0; q < 3; q++)
			{
				int p0 = t[q], p1 = t[(q + 1) % 3], p2 = t[(q + 2) % 3];
				c4[q][0][0] = p0; c4[q][0][1] = p0;
				c4[q][1][0] = p0; c4[q][1][1] = p1;
				c4[q][2][0] = p2; c4[q][2][1] = p0;
				c4[3][q][0] = p0; c4[3][q][1] = p1;
			}
			for (q = 0; q < 4; q++)
			{
				int hit = 0, r, z, ids[3];
				for (r = 0; r < 3; r++)
				{
					int a = c4[q][r][0] < c4[q][r][1] ? c4[q][r][0] : c4[q][r][1], b = c4[q][r][0] ^ c4[q][r][1] ^ a;
					c4[q][r][0] = a;
					c4[q][r][1] = b;
					for (z = 0; z < 3; z++)
						hit |= key[z][0] == a && key[z][1] == b;
				}
				if (!hit)
					continue;
				/* (a vertex key is the current triangle's corner: its whole star is here) */
				for (r = 0; r < 3; r++)
					if (c4[q][r][0] == c4[q][r][1] && c4[q][r][0] != cur[0] && c4[q][r][0] != cur[1] && c4[q][r][0] != cur[2])
						hit = 0;
				if (!hit || m->nt == MI_LOOP_MAX)
					return 0;
				for (r = 0; r < 3; r++)
					if ((ids[r] = ml_key(m, n, c4[q][r][0], c4[q][r][1])) < 0)
						return 0;
				for (r = 0; r < 3; r++)
					m->t[m->nt][r] = ids[r];
				m->nt++;
			}
		}
		/* creases: both halves of a crease edge */
		for (i = 0; i < n->nc; i++)
		{
			int h, e[2], a, b, q;
			a = n->c[i][0] < n->c[i][1] ? n->c[i][0] : n->c[i][1];
			b = n->c[i][0] ^ n->c[i][1] ^ a;
			for (h = 0; h < 2; h++)
			{
				int vtx = h ? b : a;
				e[0] = e[1] = -1;
				for (q = 0; q < m->nv; q++)
				{
					if (m->kp[q] == vtx && m->kq[q] == vtx)
						e[0] = q;
					if (m->kp[q] == a && m->kq[q] == b)
						e[1] = q;
				}
				if (e[0] >= 0 && e[1] >= 0 && m->nc < MI_LOOP_MAX)
				{
					m->c[m->nc][0] = e[0];
					m->c[m->nc][1] = e[1];
					m->nc++;
				}
			}
		}
		for (i = 0; i < 3; i++)
		{
			for (j = 0; j < m->nv && !(m->kp[j] == key[i][0] && m->kq[j] == key[i][1]); j++)
				;
			cur[i] = j;
			u[i] = u2[i];
			J[i][0] = J2[i][0];
			J[i][1] = J2[i][1];
		}
		{
			double su = u[0] + u[1] + u[2];
			for (i = 0; i < 3; i++)
				u[i] /= su;
		}
		s = n; n = m; m = s;
	}
}

static int ml_init(ml_net* n, mi_loop const* L)
{
	unsigned int i, k;
	if (L->nv > MI_LOOP_MAX || L->nt > MI_LOOP_MAX || L->nc > MI_LOOP_MAX)
		return 0;
	n->nv = (int)L->nv;
	n->nt = (int)L->nt;
	n->nc = (int)L->nc;
	for (i = 0; i < L->nv; i++)
	{
		for (k = 0; k < 3; k++)
			n->x[i][k] = L->x[3 * (size_t)L->v[i] + k];
		n->fx[i] = i < 3 ? (int)((L->fixed >> i & 1) | (L->fixed >> (3 + i) & 1) << 1) : 0;
	}
	for (i = 0; i < L->nt; i++)
		for (k = 0; k < 3; k++)
			n->t[i][k] = (int)L->t[3 * i + k];
	for (i = 0; i < L->nc; i++)
	{
		n->c[i][0] = (int)L->c[2 * i];
		n->c[i][1] = (int)L->c[2 * i + 1];
	}
	return 1;
}

/* the limit point at (a, b) and its derivatives; where they are not known (an extraordinary vertex, a crease) they are taken
   a hair inside (barycentrics 1e-7 towards the middle), and at a smooth extraordinary vertex its exact normal is nrm (*has) */
static void ml_eval(mi_loop const* L, double a, double b, double v[6][3], double* nrm, int* has)
{
	ml_net *n = (ml_net*)qaws_internal_scratch_malloc(2 * sizeof(ml_net)), *m;
	int pass, k, got;
	memset(v, 0, sizeof(double) * 18);
	*has = 0;
	if (!n)
		return;
	m = n + 1;
	for (pass = 0; pass < 2; pass++)
	{
		double u[3], J[3][2] = { { -1, -1 }, { 1, 0 }, { 0, 1 } }, w[6][3], h = pass ? 1e-7 : 0;
		int cur[3] = { 0, 1, 2 };
		u[0] = (1 - a - b) * (1 - h) + h / 3;
		u[1] = a * (1 - h) + h / 3;
		u[2] = b * (1 - h) + h / 3;
		if (!ml_init(n, L))
			break;
		got = ml_point(n, m, cur, u, J, w, nrm);
		if (!pass && got == 2)
			*has = 1;
		if (got == 1)
		{
			for (k = 0; k < 3; k++)
			{
				int j;
				for (j = pass ? 1 : 0; j < 6; j++)
					v[j][k] = w[j][k];
			}
			break;
		}
		if (!pass)
			for (k = 0; k < 3; k++)
				v[0][k] = w[0][k];
	}
	qaws_internal_free(n);
}

static void mi_eval_loop(mi_tri const* t, double a, double b, double* x, double* xa, double* xb, double* xaa, double* xab, double* xbb,
	double* n, double* na, double* nb)
{
	double v[6][3], nv[3], nn[3];
	int k, has;
	ml_eval(t->loop, a, b, v, nv, &has);
	for (k = 0; k < 3; k++)
	{
		if (x) x[k] = v[0][k];
		if (xa) xa[k] = v[1][k];
		if (xb) xb[k] = v[2][k];
		if (xaa) xaa[k] = v[3][k];
		if (xab) xab[k] = v[4][k];
		if (xbb) xbb[k] = v[5][k];
	}
	mi_surface_normal(v, nn, na, nb);
	/* (the tangent masks' normal, turned to the side of the one a hair inside) */
	if (has && mi_dot(nv, nn) < 0)
		for (k = 0; k < 3; k++)
			nv[k] = -nv[k];
	if (n)
		for (k = 0; k < 3; k++)
			n[k] = has ? nv[k] : nn[k];
}

/* the corners' limit points as p[0..2] (the flat triangle between them in the other slots: starting points, not a net) */
static void mi_build_loop(mi_tri* t)
{
	static double const at[3][2] = { { 0, 0 }, { 1, 0 }, { 0, 1 } };
	int c, k, q;
	for (c = 0; c < 3; c++)
	{
		double v[6][3], nv[3];
		int has;
		ml_eval(t->loop, at[c][0], at[c][1], v, nv, &has);
		for (k = 0; k < 3; k++)
			t->p[c][k] = v[0][k];
	}
	for (q = 3; q < 10; q++)
	{
		double w[3];
		w[0] = mi_e3[q][0] / 3.0;
		w[1] = mi_e3[q][1] / 3.0;
		w[2] = mi_e3[q][2] / 3.0;
		for (k = 0; k < 3; k++)
			t->p[q][k] = w[0] * t->p[0][k] + w[1] * t->p[1][k] + w[2] * t->p[2][k];
	}
	memset(t->n, 0, sizeof t->n);
}

void mi_eval(mi_tri const* t, double a, double b, double* x, double* xa, double* xb, double* xaa, double* xab, double* xbb,
	double* n, double* na, double* nb)
{
	double bb[3] = { 1 - a - b, a, b }, v[6][3];
	int l, k, o;
	if (t->loop)
	{
		mi_eval_loop(t, a, b, x, xa, xb, xaa, xab, xbb, n, na, nb);
		return;
	}
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
