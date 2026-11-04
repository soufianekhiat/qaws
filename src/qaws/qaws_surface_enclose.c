/*
 * Guaranteed ranges of a surface and its derivatives over a parameter box (qaws_surface_enclose): per knot span the
 * B-spline window - its derivatives' windows for the derivatives - blossomed over the span's part of the box into a
 * Bezier net, in interval arithmetic rounded outward; the nets' hulls bound the homogeneous values, the quotient rule in
 * intervals the rational ones.
 */

#include "internal/qaws_internal_closest.h"
#include "internal/qaws_internal_surface.h"
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

static en_iv en_out(double lo, double hi)
{
	en_iv r;
	r.lo = nextafter(lo, -INFINITY);
	r.hi = nextafter(hi, INFINITY);
	return r;
}
static en_iv en_exact(double v)
{
	en_iv r;
	r.lo = r.hi = v;
	return r;
}
static en_iv en_add(en_iv a, en_iv b) { return en_out(a.lo + b.lo, a.hi + b.hi); }
static en_iv en_sub(en_iv a, en_iv b) { return en_out(a.lo - b.hi, a.hi - b.lo); }
static en_iv en_mul(en_iv a, en_iv b)
{
	double const p0 = a.lo * b.lo, p1 = a.lo * b.hi, p2 = a.hi * b.lo, p3 = a.hi * b.hi;
	return en_out(fmin(fmin(p0, p1), fmin(p2, p3)), fmax(fmax(p0, p1), fmax(p2, p3)));
}
/* b not holding 0 */
static en_iv en_div(en_iv a, en_iv b)
{
	double const p0 = a.lo / b.lo, p1 = a.lo / b.hi, p2 = a.hi / b.lo, p3 = a.hi / b.hi;
	return en_out(fmin(fmin(p0, p1), fmin(p2, p3)), fmax(fmax(p0, p1), fmax(p2, p3)));
}
static en_iv en_hull(en_iv a, en_iv b)
{
	en_iv r;
	r.lo = fmin(a.lo, b.lo);
	r.hi = fmax(a.hi, b.hi);
	return r;
}

/* a window of p + 1 homogeneous points (4 intervals each, stride 4 between points of the window, `step` between windows)
   with its 2p local knots, blossomed over [a, b]: Bezier point i the blossom at (a ^ (p - i), b ^ i), de Boor with those
   arguments */
static void en_blossom(en_iv* pts, unsigned int stride, double const* uk, unsigned int p, double a, double b)
{
	en_iv tmp[EN_MAX * 4], out[EN_MAX * 4];
	unsigned int i, r, j, c;
	for (i = 0; i <= p; i++)
	{
		for (j = 0; j <= p; j++)
			memcpy(tmp + j * 4, pts + j * stride, sizeof(en_iv) * 4);
		for (r = 1; r <= p; r++)
		{
			double const x = r <= p - i ? a : b;
			for (j = p; j >= r; j--)
			{
				/* al = (x - k_(j-1)) / (k_(j+p-r) - k_(j-1)) */
				en_iv const al = en_div(en_sub(en_exact(x), en_exact(uk[j - 1])), en_sub(en_exact(uk[j + p - r]), en_exact(uk[j - 1])));
				en_iv const be = en_sub(en_exact(1), al);
				for (c = 0; c < 4; c++)
					tmp[j * 4 + c] = en_add(en_mul(be, tmp[(j - 1) * 4 + c]), en_mul(al, tmp[j * 4 + c]));
			}
		}
		memcpy(out + i * 4, tmp + p * 4, sizeof(en_iv) * 4);
	}
	for (i = 0; i <= p; i++)
		memcpy(pts + i * stride, out + i * 4, sizeof(en_iv) * 4);
}

/* a window's derivative window: p points p (P_(j+1) - P_j) / (k_(j+p) - k_j), its knots the window's but the ends */
static void en_derive(en_iv* pts, unsigned int stride, double const* uk, unsigned int p)
{
	unsigned int j, c;
	for (j = 0; j < p; j++)
	{
		en_iv const f = en_div(en_exact((double)p), en_sub(en_exact(uk[j + p]), en_exact(uk[j])));
		for (c = 0; c < 4; c++)
			pts[j * stride + c] = en_mul(f, en_sub(pts[(j + 1) * stride + c], pts[j * stride + c]));
	}
}

/* the box [a0, a1] x [b0, b1] (inside span (ku, kv)): the hulls of the derivative nets' homogeneous values, d/du^du d/dv^dv,
   into hull (4 intervals); with du = dv = 0 the projected points' hull too */
static void en_net(qaws_internal_surface_net const* m, unsigned int ku, unsigned int kv, double a0, double a1, double b0, double b1,
	unsigned int du, unsigned int dv, en_iv* hull, en_iv* proj, int* weight_ok)
{
	en_iv win[EN_MAX * EN_MAX * 4];
	double uk[2 * EN_MAX], vk[2 * EN_MAX];
	unsigned int const pu = m->p[0], pv = m->p[1], s = (pv + 1) * 4;
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
			en_derive(win + b * 4, s, uk + d, qu);
	for (d = 0; d < dv; d++, qv--)
		for (a = 0; a <= qu; a++)
			en_derive(win + a * s, 4, vk + d, qv);
	/* blossomed over the box: rows (v) then columns (u) */
	for (a = 0; a <= qu; a++)
		en_blossom(win + a * s, 4, vk + dv, qv, b0, b1);
	for (b = 0; b <= qv; b++)
		en_blossom(win + b * 4, s, uk + du, qu, a0, a1);
	for (a = 0; a <= qu; a++)
		for (b = 0; b <= qv; b++)
		{
			en_iv const* h = win + a * s + b * 4;
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
					en_iv const x = en_div(h[c], h[3]);
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
	en_iv acc[18];
	int weight_ok = 1;
	if (!out_lo || !out_hi || !surface || order > 2 || !(u0 <= u1) || !(v0 <= v1) || !isfinite(u0) || !isfinite(u1) || !isfinite(v0)
		|| !isfinite(v1))
		return QAWS_STATUS_INVALID_ARGUMENT;
	ur = qaws_surface_get_u_range(surface);
	vr = qaws_surface_get_v_range(surface);
	if (u0 < (double)ur.min_value || u1 > (double)ur.max_value || v0 < (double)vr.min_value || v1 > (double)vr.max_value)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (!qaws_internal_surface_net_init(surface, &m))
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	for (ku = m.p[0]; ku < m.n[0]; ku++)
	{
		double const ka = (double)m.k[0][ku], kb = (double)m.k[0][ku + 1], a0 = fmax(u0, ka), a1 = fmin(u1, kb);
		if (!(kb > ka) || a0 > a1)
			continue;
		for (kv = m.p[1]; kv < m.n[1]; kv++)
		{
			double const la = (double)m.k[1][kv], lb = (double)m.k[1][kv + 1], b0 = fmax(v0, la), b1 = fmin(v1, lb);
			en_iv base[4], proj[3], X[6][4], S[18], W;
			if (!(lb > la) || b0 > b1)
				continue;
			en_net(&m, ku, kv, a0, a1, b0, b1, 0, 0, base, proj, &weight_ok);
			if (!weight_ok || !(base[3].lo > 0))
				return QAWS_STATUS_NUMERICAL_FAILURE;
			W = base[3];
			for (c = 0; c < 3; c++)
				S[c] = proj[c];
			if (order >= 1)
			{
				/* X = S W: S_u = (X_u - S W_u) / W, S_v likewise */
				en_net(&m, ku, kv, a0, a1, b0, b1, 1, 0, X[0], NULL, NULL);
				en_net(&m, ku, kv, a0, a1, b0, b1, 0, 1, X[1], NULL, NULL);
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
				en_net(&m, ku, kv, a0, a1, b0, b1, 2, 0, X[2], NULL, NULL);
				en_net(&m, ku, kv, a0, a1, b0, b1, 1, 1, X[3], NULL, NULL);
				en_net(&m, ku, kv, a0, a1, b0, b1, 0, 2, X[4], NULL, NULL);
				for (c = 0; c < 3; c++)
				{
					en_iv const su = S[3 + c], sv = S[6 + c], s0 = S[c];
					S[9 + c] = en_div(en_sub(en_sub(X[2][c], en_mul(two, en_mul(su, X[0][3]))), en_mul(s0, X[2][3])), W);
					S[12 + c] = en_div(en_sub(en_sub(en_sub(X[3][c], en_mul(su, X[1][3])), en_mul(sv, X[0][3])), en_mul(s0, X[3][3])), W);
					S[15 + c] = en_div(en_sub(en_sub(X[4][c], en_mul(two, en_mul(sv, X[1][3]))), en_mul(s0, X[4][3])), W);
				}
			}
			for (q = 0; q < (order == 0 ? 3u : order == 1 ? 9u : 18u); q++)
				acc[q] = any ? en_hull(acc[q], S[q]) : S[q];
			any = 1;
		}
	}
	if (!any)
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (q = 0; q < (order == 0 ? 3u : order == 1 ? 9u : 18u); q++)
	{
		out_lo[q] = acc[q].lo;
		out_hi[q] = acc[q].hi;
	}
	return QAWS_STATUS_OK;
}
