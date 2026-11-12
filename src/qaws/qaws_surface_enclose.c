/*
 * Guaranteed ranges of a surface or a curve and its derivatives over a parameter box (qaws_surface_enclose,
 * qaws_curve_enclose): per knot span the B-spline window - its derivatives' windows for the derivatives - blossomed over
 * the span's part of the box into a Bezier net, in interval arithmetic rounded outward; the nets' hulls bound the
 * homogeneous values, the quotient rule in intervals the rational ones.
 */

#include "internal/qaws_internal_closest.h"
#include "internal/qaws_internal_surface.h"
#include "internal/qaws_internal_types.h"
#include "qaws_inspect.h"
#include "qaws_surface.h"
#include "qaws_platform.h"
#include <math.h>
#include <string.h>

#define EN_MAX QAWS_INTERNAL_NET_MAX

typedef struct en_iv
{
	double lo, hi;
} en_iv;

/* outward by at least an ulp without nextafter (a library call, the cost of every operation here): |x| 2^-51 is two ulps
   or more, so the rounded sum still passes x's neighbour when it crosses into the next binade; 2^-1021 past any
   subnormal's rounding, itself normal (subnormal operands take microcode, a hundred times slower) */
static double en_down(double x) { return isfinite(x) ? x - (fabs(x) * 0x1p-51 + 0x1p-1021) : x; }
static double en_up(double x) { return isfinite(x) ? x + (fabs(x) * 0x1p-51 + 0x1p-1021) : x; }
/* (fmin and fmax are library calls under /fp:precise: the comparisons inline) */
static double en_min(double a, double b) { return a < b ? a : b; }
static double en_max(double a, double b) { return a > b ? a : b; }
static en_iv en_out(double lo, double hi)
{
	en_iv r;
	r.lo = en_down(lo);
	r.hi = en_up(hi);
	return r;
}
static en_iv en_exact(double v)
{
	en_iv r;
	r.lo = r.hi = v;
	return r;
}
/* a point 0 kept exact (a net's zero coordinates would otherwise carry a floor's noise through every product) */
static int en_zero(en_iv a) { return a.lo == 0 && a.hi == 0; }
static en_iv en_add(en_iv a, en_iv b) { return en_zero(a) ? b : en_zero(b) ? a : en_out(a.lo + b.lo, a.hi + b.hi); }
static en_iv en_sub(en_iv a, en_iv b)
{
	en_iv r;
	if (en_zero(b))
		return a;
	if (a.lo == a.hi && b.lo == b.hi && a.lo == b.lo)
	{
		r.lo = r.hi = 0;   /* (x - x, exact: a constant weight's derivative) */
		return r;
	}
	if (en_zero(a))
	{
		r.lo = -b.hi, r.hi = -b.lo;   /* (exact) */
		return r;
	}
	return en_out(a.lo - b.hi, a.hi - b.lo);
}
static en_iv en_mul(en_iv a, en_iv b)
{
	if (en_zero(a) || en_zero(b))
		return en_zero(a) ? a : b;
	double const p0 = a.lo * b.lo, p1 = a.lo * b.hi, p2 = a.hi * b.lo, p3 = a.hi * b.hi;
	return en_out(en_min(en_min(p0, p1), en_min(p2, p3)), en_max(en_max(p0, p1), en_max(p2, p3)));
}
/* b not holding 0 */
static en_iv en_div(en_iv a, en_iv b)
{
	double const p0 = a.lo / b.lo, p1 = a.lo / b.hi, p2 = a.hi / b.lo, p3 = a.hi / b.hi;
	return en_out(en_min(en_min(p0, p1), en_min(p2, p3)), en_max(en_max(p0, p1), en_max(p2, p3)));
}
static en_iv en_hull(en_iv a, en_iv b)
{
	en_iv r;
	r.lo = en_min(a.lo, b.lo);
	r.hi = en_max(a.hi, b.hi);
	return r;
}

/* de Boor's weights for blossoming a window of degree p with its 2p local knots at a or b: per argument (a, b), step r and
   point j, al = (x - k_(j-1)) / (k_(j+p-r) - k_(j-1)) and 1 - al - the same for every row or column of a net */
typedef struct en_weights
{
	unsigned int p;
	en_iv al[2][EN_MAX][EN_MAX], be[2][EN_MAX][EN_MAX];
} en_weights;

static void en_weigh(en_weights* w, double const* uk, unsigned int p, double a, double b)
{
	unsigned int x, r, j;
	w->p = p;
	for (x = 0; x < 2; x++)
		for (r = 1; r <= p; r++)
			for (j = r; j <= p; j++)
			{
				w->al[x][r][j] = en_div(en_sub(en_exact(x ? b : a), en_exact(uk[j - 1])), en_sub(en_exact(uk[j + p - r]), en_exact(uk[j - 1])));
				w->be[x][r][j] = en_sub(en_exact(1), w->al[x][r][j]);
			}
}

/* a window of p + 1 homogeneous points (4 intervals each, `stride` between points; the first nc of them worked)
   blossomed over [a, b] (its weights w): Bezier point i the blossom at (a ^ (p - i), b ^ i), de Boor with those
   arguments */
static void en_blossom(en_iv* pts, unsigned int stride, en_weights const* w, unsigned int nc)
{
	en_iv tmp[EN_MAX * 4], out[EN_MAX * 4];
	unsigned int const p = w->p;
	unsigned int i, r, j, c;
	for (i = 0; i <= p; i++)
	{
		for (j = 0; j <= p; j++)
			memcpy(tmp + j * 4, pts + j * stride, sizeof(en_iv) * 4);
		for (r = 1; r <= p; r++)
		{
			unsigned int const x = r <= p - i ? 0u : 1u;
			for (j = p; j >= r; j--)
			{
				en_iv const al = w->al[x][r][j], be = w->be[x][r][j];
				for (c = 0; c < nc; c++)
					tmp[j * 4 + c] = en_add(en_mul(be, tmp[(j - 1) * 4 + c]), en_mul(al, tmp[j * 4 + c]));
			}
		}
		memcpy(out + i * 4, tmp + p * 4, sizeof(en_iv) * 4);
	}
	for (i = 0; i <= p; i++)
		memcpy(pts + i * stride, out + i * 4, sizeof(en_iv) * 4);
}

/* a window's derivative window: p points p (P_(j+1) - P_j) / (k_(j+p) - k_j), its knots the window's but the ends */
static void en_derive(en_iv* pts, unsigned int stride, double const* uk, unsigned int p, unsigned int nc)
{
	unsigned int j, c;
	for (j = 0; j < p; j++)
	{
		en_iv const f = en_div(en_exact((double)p), en_sub(en_exact(uk[j + p]), en_exact(uk[j])));
		for (c = 0; c < nc; c++)
			pts[j * stride + c] = en_mul(f, en_sub(pts[(j + 1) * stride + c], pts[j * stride + c]));
	}
}

/* the box [a0, a1] x [b0, b1] (inside span (ku, kv)): the hulls of the derivative nets' homogeneous values, d/du^du d/dv^dv,
   into hull (4 intervals); with du = dv = 0 the projected points' hull too. A polynomial net's weights are 1 exactly (its
   derivatives' 0): not worked, so no rounding's noise reaches the quotient rule */
static void en_net(qaws_internal_surface_net const* m, int poly, unsigned int ku, unsigned int kv, double a0, double a1, double b0,
	double b1, unsigned int du, unsigned int dv, en_iv* hull, en_iv* proj, int* weight_ok)
{
	en_iv win[EN_MAX * EN_MAX * 4];
	en_weights w;
	double uk[2 * EN_MAX], vk[2 * EN_MAX];
	unsigned int const pu = m->p[0], pv = m->p[1], s = (pv + 1) * 4, nc = poly ? 3u : 4u;
	unsigned int a, b, c, d, qu, qv;
	if (du > pu || dv > pv)
	{
		/* (past the degree: the derivative is 0) */
		for (c = 0; c < 4; c++)
			hull[c] = en_exact(0);
		return;
	}
	for (a = 0; a <= pu; a++)
		for (b = 0; b <= pv; b++)
		{
			double h[4];
			qaws_internal_surface_net_point(m, ku - pu + a, kv - pv + b, h);
			for (c = 0; c < 4; c++)
				win[a * s + b * 4 + c] = en_exact(h[c]);
		}
	for (a = 0; a < 2 * pu; a++)
		uk[a] = (double)m->k[0][ku - pu + 1 + a];
	for (b = 0; b < 2 * pv; b++)
		vk[b] = (double)m->k[1][kv - pv + 1 + b];
	/* the derivative windows: in u along each column, in v along each row; each time the knots lose their ends */
	qu = pu, qv = pv;
	for (d = 0; d < du; d++, qu--)
		for (b = 0; b <= qv; b++)
			en_derive(win + b * 4, s, uk + d, qu, nc);
	for (d = 0; d < dv; d++, qv--)
		for (a = 0; a <= qu; a++)
			en_derive(win + a * s, 4, vk + d, qv, nc);
	/* blossomed over the box: rows (v) then columns (u) */
	en_weigh(&w, vk + dv, qv, b0, b1);
	for (a = 0; a <= qu; a++)
		en_blossom(win + a * s, 4, &w, nc);
	en_weigh(&w, uk + du, qu, a0, a1);
	for (b = 0; b <= qv; b++)
		en_blossom(win + b * 4, s, &w, nc);
	for (a = 0; a <= qu; a++)
		for (b = 0; b <= qv; b++)
		{
			en_iv* h = win + a * s + b * 4;
			if (poly)
				h[3] = en_exact(du || dv ? 0.0 : 1.0);
			for (c = 0; c < 4; c++)
				hull[c] = a == 0 && b == 0 ? h[c] : en_hull(hull[c], h[c]);
			if (proj)
			{
				if (!(h[3].lo > 0))
				{
					*weight_ok = 0;
					continue;
				}
				for (c = 0; c < 3; c++)
				{
					en_iv const x = poly ? h[c] : en_div(h[c], h[3]);
					proj[c] = a == 0 && b == 0 ? x : en_hull(proj[c], x);
				}
			}
		}
}

qaws_status qaws_surface_enclose(double* out_lo, double* out_hi, qaws_surface const* surface, double u0, double u1, double v0, double v1,
	unsigned int order)
{
	qaws_internal_surface_net m;
	qaws_range ur, vr;
	unsigned int ku, kv, c, q, any = 0;
	int poly;
	en_iv acc[30];
	int weight_ok = 1;
	unsigned int const count = order == 0 ? 3u : order == 1 ? 9u : order == 2 ? 18u : 30u;
	if (!out_lo || !out_hi || !surface || order > 3 || !(u0 <= u1) || !(v0 <= v1) || !isfinite(u0) || !isfinite(u1) || !isfinite(v0)
		|| !isfinite(v1))
		return QAWS_STATUS_INVALID_ARGUMENT;
	ur = qaws_surface_get_u_range(surface);
	vr = qaws_surface_get_v_range(surface);
	if (u0 < (double)ur.min_value || u1 > (double)ur.max_value || v0 < (double)vr.min_value || v1 > (double)vr.max_value)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (!qaws_internal_surface_net_init(surface, &m))
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	poly = m.cp && !m.w;
	for (ku = m.p[0]; ku < m.n[0]; ku++)
	{
		double const ka = (double)m.k[0][ku], kb = (double)m.k[0][ku + 1], a0 = fmax(u0, ka), a1 = fmin(u1, kb);
		if (!(kb > ka) || a0 > a1)
			continue;
		for (kv = m.p[1]; kv < m.n[1]; kv++)
		{
			double const la = (double)m.k[1][kv], lb = (double)m.k[1][kv + 1], b0 = fmax(v0, la), b1 = fmin(v1, lb);
			en_iv base[4], proj[3], X[9][4], S[30], W;
			if (!(lb > la) || b0 > b1)
				continue;
			en_net(&m, poly, ku, kv, a0, a1, b0, b1, 0, 0, base, proj, &weight_ok);
			if (!weight_ok || !(base[3].lo > 0))
				return QAWS_STATUS_NUMERICAL_FAILURE;
			W = base[3];
			for (c = 0; c < 3; c++)
				S[c] = proj[c];
			if (order >= 1)
			{
				/* X = S W: S_u = (X_u - S W_u) / W, S_v likewise */
				en_net(&m, poly, ku, kv, a0, a1, b0, b1, 1, 0, X[0], NULL, NULL);
				en_net(&m, poly, ku, kv, a0, a1, b0, b1, 0, 1, X[1], NULL, NULL);
				for (c = 0; c < 3; c++)
				{
					S[3 + c] = en_div(en_sub(X[0][c], en_mul(S[c], X[0][3])), W);
					S[6 + c] = en_div(en_sub(X[1][c], en_mul(S[c], X[1][3])), W);
				}
			}
			if (order >= 2)
			{
				/* X_uu = S_uu W + 2 S_u W_u + S W_uu; X_uv = S_uv W + S_u W_v + S_v W_u + S W_uv */
				en_iv const two = en_exact(2);
				en_net(&m, poly, ku, kv, a0, a1, b0, b1, 2, 0, X[2], NULL, NULL);
				en_net(&m, poly, ku, kv, a0, a1, b0, b1, 1, 1, X[3], NULL, NULL);
				en_net(&m, poly, ku, kv, a0, a1, b0, b1, 0, 2, X[4], NULL, NULL);
				for (c = 0; c < 3; c++)
				{
					en_iv const su = S[3 + c], sv = S[6 + c], s0 = S[c];
					S[9 + c] = en_div(en_sub(en_sub(X[2][c], en_mul(two, en_mul(su, X[0][3]))), en_mul(s0, X[2][3])), W);
					S[12 + c] = en_div(en_sub(en_sub(en_sub(X[3][c], en_mul(su, X[1][3])), en_mul(sv, X[0][3])), en_mul(s0, X[3][3])), W);
					S[15 + c] = en_div(en_sub(en_sub(X[4][c], en_mul(two, en_mul(sv, X[1][3]))), en_mul(s0, X[4][3])), W);
				}
			}
			if (order >= 3)
			{
				/* X_uuu = S_uuu W + 3 S_uu W_u + 3 S_u W_uu + S W_uuu;
				   X_uuv = S_uuv W + S_uu W_v + 2 S_uv W_u + 2 S_u W_uv + S_v W_uu + S W_uuv; X_uvv and X_vvv alike */
				en_iv const two = en_exact(2), three = en_exact(3);
				en_net(&m, poly, ku, kv, a0, a1, b0, b1, 3, 0, X[5], NULL, NULL);
				en_net(&m, poly, ku, kv, a0, a1, b0, b1, 2, 1, X[6], NULL, NULL);
				en_net(&m, poly, ku, kv, a0, a1, b0, b1, 1, 2, X[7], NULL, NULL);
				en_net(&m, poly, ku, kv, a0, a1, b0, b1, 0, 3, X[8], NULL, NULL);
				for (c = 0; c < 3; c++)
				{
					en_iv const s0 = S[c], su = S[3 + c], sv = S[6 + c], suu = S[9 + c], suv = S[12 + c], svv = S[15 + c];
					en_iv const Wu = X[0][3], Wv = X[1][3], Wuu = X[2][3], Wuv = X[3][3], Wvv = X[4][3];
					en_iv r;
					r = en_sub(en_sub(X[5][c], en_mul(three, en_mul(suu, Wu))), en_mul(three, en_mul(su, Wuu)));
					S[18 + c] = en_div(en_sub(r, en_mul(s0, X[5][3])), W);
					r = en_sub(en_sub(en_sub(X[6][c], en_mul(suu, Wv)), en_mul(two, en_mul(suv, Wu))), en_mul(two, en_mul(su, Wuv)));
					S[21 + c] = en_div(en_sub(en_sub(r, en_mul(sv, Wuu)), en_mul(s0, X[6][3])), W);
					r = en_sub(en_sub(en_sub(X[7][c], en_mul(svv, Wu)), en_mul(two, en_mul(suv, Wv))), en_mul(two, en_mul(sv, Wuv)));
					S[24 + c] = en_div(en_sub(en_sub(r, en_mul(su, Wvv)), en_mul(s0, X[7][3])), W);
					r = en_sub(en_sub(X[8][c], en_mul(three, en_mul(svv, Wv))), en_mul(three, en_mul(sv, Wvv)));
					S[27 + c] = en_div(en_sub(r, en_mul(s0, X[8][3])), W);
				}
			}
			for (q = 0; q < count; q++)
				acc[q] = any ? en_hull(acc[q], S[q]) : S[q];
			any = 1;
		}
	}
	if (!any)
		return QAWS_STATUS_INVALID_ARGUMENT;
	/* the third derivatives over a knot line inside the box across which S is not C2 (multiplicity past the degree less 2):
	   the second derivatives jump there, so no bound holds them over the box - unbounded */
	if (order >= 3)
	{
		unsigned int d, j, i, cut = 0;
		for (d = 0; d < 2 && !cut; d++)
		{
			double const a = d ? v0 : u0, b = d ? v1 : u1;
			for (j = m.p[d] + 1; j < m.n[d] && !cut; j++)
			{
				double const kk = (double)m.k[d][j];
				int mult = 0;
				if (!(kk > a && kk < b))
					continue;
				for (i = 0; i < m.n[d] + m.p[d] + 1; i++)
					mult += (double)m.k[d][i] == kk;
				cut = mult > (int)m.p[d] - 2;
			}
		}
		for (q = 18; cut && q < 30; q++)
			acc[q].lo = -INFINITY, acc[q].hi = INFINITY;
	}
	for (q = 0; q < count; q++)
	{
		out_lo[q] = acc[q].lo;
		out_hi[q] = acc[q].hi;
	}
	return QAWS_STATUS_OK;
}

/* ---- curves ---- */

/* span s of a curve as a window of homogeneous points (4 intervals each, dim + 1 of them used, the weight last) and its
   2p knots: the Bezier and spline families through their B-spline windows; a Hermite span its Bezier points (3 p0,
   3 p0 + m0, 3 p1 - m1, 3 p1) / 3 and a uniform Catmull-Rom span's (6 q1, 6 q1 + q2 - q0, 6 q2 - q3 + q1, 6 q2) / 6, in
   intervals, over the span */
static int en_curve_window(qaws_curve const* c, unsigned int s, unsigned int* p, en_iv* win, double* uk)
{
	unsigned int const dim = (unsigned int)c->dimension, h = dim + 1;
	unsigned int i, j, k;
	if (c->kind == QAWS_CURVE_KIND_HERMITE || c->kind == QAWS_CURVE_KIND_CATMULL_ROM)
	{
		double const ta = (double)c->span_boundaries[s], tb = (double)c->span_boundaries[s + 1];
		en_iv const three = en_exact(3), six = en_exact(6);
		for (k = 0; k < dim; k++)
		{
			en_iv b[4];
			if (c->kind == QAWS_CURVE_KIND_HERMITE)
			{
				qaws_hermite_impl const* m = (qaws_hermite_impl const*)c->impl;
				en_iv const p0 = en_exact((double)m->points[s * dim + k]), p1 = en_exact((double)m->points[(s + 1) * dim + k]);
				en_iv const m0 = en_exact((double)m->tangents[s * dim + k]), m1 = en_exact((double)m->tangents[(s + 1) * dim + k]);
				b[0] = p0, b[3] = p1;
				b[1] = en_div(en_add(en_mul(three, p0), m0), three);
				b[2] = en_div(en_sub(en_mul(three, p1), m1), three);
			}
			else
			{
				qaws_catmull_rom_impl const* m = (qaws_catmull_rom_impl const*)c->impl;
				unsigned int const n = m->control_point_count, b0 = m->closed ? s + n - 1 : s;
				en_iv q[4];
				if (m->parameterization != QAWS_PARAMETERIZATION_UNIFORM || n < 2)
					return 0;
				for (j = 0; j < 4; j++)
					q[j] = en_exact((double)m->control_points[((b0 + j) % n) * dim + k]);
				b[0] = q[1], b[3] = q[2];
				b[1] = en_div(en_sub(en_add(en_mul(six, q[1]), q[2]), q[0]), six);
				b[2] = en_div(en_add(en_sub(en_mul(six, q[2]), q[3]), q[1]), six);
			}
			for (j = 0; j < 4; j++)
				win[j * 4 + k] = b[j];
		}
		for (j = 0; j < 4; j++)
			win[j * 4 + dim] = en_exact(1);
		for (j = 0; j < 6; j++)
			uk[j] = j < 3 ? ta : tb;
		*p = 3;
		return 1;
	}
	{
		double w[EN_MAX * 4], kn[2 * EN_MAX];
		if (c->degree >= EN_MAX || !qaws_internal_curve_window(c, s, p, w, kn))
			return 0;
		for (j = 0; j <= *p; j++)
			for (i = 0; i < h; i++)
				win[j * 4 + i] = en_exact(w[j * h + i]);
		for (j = 0; j < 2 * *p; j++)
			uk[j] = kn[j];
		return 1;
	}
}

qaws_status qaws_curve_enclose(double* out_lo, double* out_hi, qaws_curve const* curve, double t0, double t1, unsigned int order)
{
	qaws_range r;
	unsigned int dim, s, any = 0, q, c, d;
	en_iv acc[9];
	if (!out_lo || !out_hi || !curve || order > 2 || !(t0 <= t1) || !isfinite(t0) || !isfinite(t1))
		return QAWS_STATUS_INVALID_ARGUMENT;
	r = qaws_curve_get_parameter_range(curve);
	if (t0 < (double)r.min_value || t1 > (double)r.max_value)
		return QAWS_STATUS_INVALID_ARGUMENT;
	dim = (unsigned int)curve->dimension;
	if (dim < 1 || dim > 3)
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	for (s = 0; s < curve->span_count; s++)
	{
		double const ta = (double)curve->span_boundaries[s], tb = (double)curve->span_boundaries[s + 1], a0 = fmax(t0, ta), a1 = fmin(t1, tb);
		en_iv win[EN_MAX * 4], X[3][4], C[9], W;
		double uk[2 * EN_MAX];
		unsigned int p, nc;
		int poly = 1;
		if (!(tb > ta) || a0 > a1)
			continue;
		if (!en_curve_window(curve, s, &p, win, uk))
			return QAWS_STATUS_UNSUPPORTED_OPERATION;
		for (q = 0; q <= p; q++)
			poly = poly && win[q * 4 + dim].lo == 1 && win[q * 4 + dim].hi == 1;
		nc = poly ? dim : dim + 1;
		for (d = 0; d <= order; d++)
		{
			en_iv w2[EN_MAX * 4];
			en_weights wt;
			unsigned int qd = p, i;
			if (d > p)
			{
				for (c = 0; c < 4; c++)
					X[d][c] = en_exact(0);
				continue;
			}
			memcpy(w2, win, sizeof(en_iv) * 4 * (p + 1));
			for (i = 0; i < d; i++, qd--)
				en_derive(w2, 4, uk + i, qd, nc);
			en_weigh(&wt, uk + d, qd, a0, a1);
			en_blossom(w2, 4, &wt, nc);
			if (poly)
				for (i = 0; i <= qd; i++)
					w2[i * 4 + dim] = en_exact(d ? 0.0 : 1.0);
			for (i = 0; i <= qd; i++)
				for (c = 0; c <= dim; c++)
					X[d][c] = i == 0 ? w2[c] : en_hull(X[d][c], w2[i * 4 + c]);
			if (d == 0)
				for (i = 0; i <= qd; i++)
				{
					if (!poly && !(w2[i * 4 + dim].lo > 0))
						return QAWS_STATUS_NUMERICAL_FAILURE;
					for (c = 0; c < dim; c++)
					{
						en_iv const x = poly ? w2[i * 4 + c] : en_div(w2[i * 4 + c], w2[i * 4 + dim]);
						C[c] = i == 0 ? x : en_hull(C[c], x);
					}
				}
		}
		W = X[0][dim];
		if (!(W.lo > 0))
			return QAWS_STATUS_NUMERICAL_FAILURE;
		/* X = C W: C' = (X' - C W') / W, C'' = (X'' - 2 C' W' - C W'') / W */
		for (c = 0; order >= 1 && c < dim; c++)
			C[dim + c] = poly ? X[1][c] : en_div(en_sub(X[1][c], en_mul(C[c], X[1][dim])), W);
		for (c = 0; order >= 2 && c < dim; c++)
			C[2 * dim + c] = poly ? X[2][c]
				: en_div(en_sub(en_sub(X[2][c], en_mul(en_exact(2), en_mul(C[dim + c], X[1][dim]))), en_mul(C[c], X[2][dim])), W);
		for (q = 0; q < (order + 1) * dim; q++)
			acc[q] = any ? en_hull(acc[q], C[q]) : C[q];
		any = 1;
	}
	if (!any)
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (q = 0; q < (order + 1) * dim; q++)
	{
		out_lo[q] = acc[q].lo;
		out_hi[q] = acc[q].hi;
	}
	return QAWS_STATUS_OK;
}
