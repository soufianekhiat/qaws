/*
 * Validated enclosures of geodesics on polynomial surfaces
 * (qaws_surface_geodesic_enclose, qaws_surface_region.h).
 *
 * Intervals with outward rounding (each bound stepped one ulp out after
 * every operation). The surface is a set of Bezier patches, each with its
 * derivative nets up to the third order, in intervals; a derivative over a
 * box is the patch's Bernstein form over the box, tightened by the mean
 * value form from the next order. The geodesic y = (u, v, u', v') solves
 * y' = f(y), f = (u', v', a), a^k = -G^k_ij w^i w^j; its second derivative
 * y'' = (a, Da) takes the gradient of a in (u, v) (interval duals, through
 * the third derivatives) and in w.
 *
 * A step of length h from the box Y:
 *   - the a priori box B: Y + [0, h] f(B) within B (Picard), found by
 *     inflating Y + [0, h] f(Y); the geodesic over the step lies in B;
 *   - the end: Y + h f(Y) + h^2/2 y''(B), met with Y + h f(B).
 */

#include "qaws_surface_region.h"
#include "qaws_surface.h"
#include "internal/qaws_internal_surface.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------------ */
/*  Intervals                                                               */
/* ------------------------------------------------------------------------ */

typedef struct ge_iv
{
	double lo, hi;
} ge_iv;

/* one rounding (half an ulp, at most |x| 2^-53) stepped over outward; bounds
   under 1e-290 put at it (denormal arithmetic is very slow, and a bound
   moved out stays a bound) */
static double ge_dn(double x) { return fabs(x) < 1e-290 ? -1e-290 : x - fabs(x) * 2.220446049250313e-16; }
static double ge_up(double x) { return fabs(x) < 1e-290 ? 1e-290 : x + fabs(x) * 2.220446049250313e-16; }

static ge_iv ge_pt(double x)
{
	ge_iv r;
	r.lo = r.hi = x;
	return r;
}

static ge_iv ge_mk(double lo, double hi)
{
	ge_iv r;
	r.lo = lo; r.hi = hi;
	return r;
}

static ge_iv ge_add(ge_iv a, ge_iv b) { return ge_mk(ge_dn(a.lo + b.lo), ge_up(a.hi + b.hi)); }
static ge_iv ge_sub(ge_iv a, ge_iv b) { return ge_mk(ge_dn(a.lo - b.hi), ge_up(a.hi - b.lo)); }
static ge_iv ge_neg(ge_iv a) { return ge_mk(-a.hi, -a.lo); }

static ge_iv ge_mul(ge_iv a, ge_iv b)
{
	double p[4], lo, hi;
	unsigned int i;
	p[0] = a.lo * b.lo; p[1] = a.lo * b.hi; p[2] = a.hi * b.lo; p[3] = a.hi * b.hi;
	lo = hi = p[0];
	for (i = 1; i < 4; i++)
	{
		if (p[i] < lo) lo = p[i];
		if (p[i] > hi) hi = p[i];
	}
	return ge_mk(ge_dn(lo), ge_up(hi));
}

static ge_iv ge_scl(ge_iv a, double k) { return ge_mul(a, ge_pt(k)); }

/* a / b; 0 when b holds 0 */
static int ge_div(ge_iv a, ge_iv b, ge_iv* r)
{
	if (!(b.lo > 0 || b.hi < 0))
		return 0;
	*r = ge_mul(a, ge_mk(ge_dn(1 / b.hi), ge_up(1 / b.lo)));
	return 1;
}

static ge_iv ge_hull(ge_iv a, ge_iv b) { return ge_mk(a.lo < b.lo ? a.lo : b.lo, a.hi > b.hi ? a.hi : b.hi); }

/* a met with b (a when they do not meet: both hold the value, so they do) */
static ge_iv ge_meet(ge_iv a, ge_iv b)
{
	ge_iv r = ge_mk(a.lo > b.lo ? a.lo : b.lo, a.hi < b.hi ? a.hi : b.hi);
	return r.lo <= r.hi ? r : a;
}

static double ge_w(ge_iv a) { return a.hi - a.lo; }
static double ge_mid(ge_iv a) { return 0.5 * (a.lo + a.hi); }
static double ge_mag(ge_iv a) { return fabs(a.lo) > fabs(a.hi) ? fabs(a.lo) : fabs(a.hi); }
static int ge_in(ge_iv a, ge_iv b) { return a.lo >= b.lo && a.hi <= b.hi; }

/* t^k for t >= 0 */
static ge_iv ge_powp(ge_iv t, unsigned int k)
{
	ge_iv r = ge_pt(1);
	while (k--)
		r = ge_mul(r, t);
	return r;
}

/* interval duals: a value and its gradient in (u, v) */
typedef struct ge_d
{
	ge_iv v, du, dv;
} ge_d;

static ge_d ge_dadd(ge_d a, ge_d b) { ge_d r; r.v = ge_add(a.v, b.v); r.du = ge_add(a.du, b.du); r.dv = ge_add(a.dv, b.dv); return r; }
static ge_d ge_dsub(ge_d a, ge_d b) { ge_d r; r.v = ge_sub(a.v, b.v); r.du = ge_sub(a.du, b.du); r.dv = ge_sub(a.dv, b.dv); return r; }

static ge_d ge_dmul(ge_d a, ge_d b)
{
	ge_d r;
	r.v = ge_mul(a.v, b.v);
	r.du = ge_add(ge_mul(a.du, b.v), ge_mul(a.v, b.du));
	r.dv = ge_add(ge_mul(a.dv, b.v), ge_mul(a.v, b.dv));
	return r;
}

static int ge_ddiv(ge_d a, ge_d b, ge_d* r)
{
	if (!ge_div(a.v, b.v, &r->v) || !ge_div(ge_sub(a.du, ge_mul(r->v, b.du)), b.v, &r->du) ||
		!ge_div(ge_sub(a.dv, ge_mul(r->v, b.dv)), b.v, &r->dv))
		return 0;
	return 1;
}

static ge_d ge_ddot(ge_d const* a, ge_d const* b)
{
	return ge_dadd(ge_dadd(ge_dmul(a[0], b[0]), ge_dmul(a[1], b[1])), ge_dmul(a[2], b[2]));
}

/* ------------------------------------------------------------------------ */
/*  The surface as Bezier patches                                           */
/* ------------------------------------------------------------------------ */

#define GE_NETS 10     /* derivative orders (a, b), a + b <= 3 */

static unsigned int ge_net_index(unsigned int a, unsigned int b)
{
	static unsigned int const idx[4][4] = { { 0, 2, 5, 9 }, { 1, 4, 8, 99 }, { 3, 7, 99, 99 }, { 6, 99, 99, 99 } };
	return idx[a][b];
}

typedef struct ge_patch
{
	double u0, u1, v0, v1;
	unsigned int m, n;              /* degrees */
	ge_iv* net[GE_NETS];            /* (m - a + 1) x (n - b + 1) points of 3 intervals, u index major */
} ge_patch;

typedef struct ge_surf
{
	ge_patch* p;
	unsigned int np;
	double u0, u1, v0, v1, diag;
} ge_surf;

static void ge_surf_free(ge_surf* s)
{
	unsigned int i, k;
	for (i = 0; i < s->np; i++)
		for (k = 0; k < GE_NETS; k++)
			free(s->p[i].net[k]);
	free(s->p);
	memset(s, 0, sizeof(*s));
}

/* a patch from its control net (intervals, (m + 1) x (n + 1), u major) over
   [u0, u1] x [v0, v1]: its derivative nets */
static qaws_status ge_patch_make(ge_patch* pt, ge_iv const* cp, unsigned int m, unsigned int n, double u0, double u1, double v0, double v1)
{
	unsigned int a, b, i, j, c;
	ge_iv du = ge_sub(ge_pt(u1), ge_pt(u0)), dv = ge_sub(ge_pt(v1), ge_pt(v0));
	memset(pt, 0, sizeof(*pt));
	pt->u0 = u0; pt->u1 = u1; pt->v0 = v0; pt->v1 = v1; pt->m = m; pt->n = n;
	for (a = 0; a <= 3; a++)
		for (b = 0; a + b <= 3; b++)
		{
			unsigned int k = ge_net_index(a, b), rm = m >= a ? m - a + 1 : 0, rn = n >= b ? n - b + 1 : 0;
			ge_iv* q;
			if (!rm || !rn)
			{
				/* a derivative of an order past the degree: zero */
				pt->net[k] = (ge_iv*)calloc(3, sizeof(ge_iv));
				if (!pt->net[k])
					return QAWS_STATUS_ALLOCATION_FAILURE;
				continue;
			}
			q = (ge_iv*)malloc(sizeof(ge_iv) * 3 * rm * rn);
			if (!q)
				return QAWS_STATUS_ALLOCATION_FAILURE;
			pt->net[k] = q;
			if (a == 0 && b == 0)
			{
				memcpy(q, cp, sizeof(ge_iv) * 3 * rm * rn);
				continue;
			}
			{
				/* from the net one order lower: differences times degree over the span */
				int along_u = a > 0;
				unsigned int pa = along_u ? a - 1 : a, pb = along_u ? b : b - 1;
				ge_iv const* s = pt->net[ge_net_index(pa, pb)];
				unsigned int sn = n - pb + 1;
				ge_iv k2;
				if (along_u)
					ge_div(ge_pt((double)(m - pa)), du, &k2);
				else
					ge_div(ge_pt((double)(n - pb)), dv, &k2);
				for (i = 0; i < rm; i++)
					for (j = 0; j < rn; j++)
						for (c = 0; c < 3; c++)
						{
							ge_iv x0 = s[3 * (i * sn + j) + c], x1 = along_u ? s[3 * ((i + 1) * sn + j) + c] : s[3 * (i * sn + j + 1) + c];
							q[3 * (i * rn + j) + c] = ge_mul(ge_sub(x1, x0), k2);
						}
			}
		}
	return QAWS_STATUS_OK;
}

static qaws_status ge_surf_add(ge_surf* s, ge_iv const* cp, unsigned int m, unsigned int n, double u0, double u1, double v0, double v1)
{
	ge_patch* g = (ge_patch*)realloc(s->p, sizeof(ge_patch) * (s->np + 1));
	qaws_status st;
	if (!g)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	s->p = g;
	st = ge_patch_make(&s->p[s->np], cp, m, n, u0, u1, v0, v1);
	s->np++;
	return st;
}

/* knot insertion (Boehm) on rows of 3-interval points: t inserted once
   into the knot vector kn (nk knots, degree d) of cnt points per row,
   for rows rows of stride points; kn and the points grow by one */
static qaws_status ge_insert(ge_iv** pts, unsigned int rows, unsigned int* cnt, double** kn, unsigned int* nk, unsigned int d, double t)
{
	unsigned int c = *cnt, r, i, k, sp = 0;
	ge_iv* q = (ge_iv*)malloc(sizeof(ge_iv) * 3 * rows * (c + 1));
	double* kk = (double*)malloc(sizeof(double) * (*nk + 1));
	if (!q || !kk)
	{
		free(q); free(kk);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	/* the span: kn[sp] <= t < kn[sp + 1] */
	for (i = 0; i + 1 < *nk; i++)
		if ((*kn)[i] <= t && t < (*kn)[i + 1])
			sp = i;
	if (sp < d)
	{
		free(q); free(kk);
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	}
	for (r = 0; r < rows; r++)
	{
		ge_iv const* p = *pts + 3 * r * c;
		ge_iv* o = q + 3 * r * (c + 1);
		for (i = 0; i <= c; i++)
			for (k = 0; k < 3; k++)
			{
				if (i + d <= sp)
					o[3 * i + k] = p[3 * i + k];
				else if (i > sp)
					o[3 * i + k] = p[3 * (i - 1) + k];
				else
				{
					/* alpha = (t - kn[i]) / (kn[i + d] - kn[i]) */
					ge_iv al;
					ge_div(ge_sub(ge_pt(t), ge_pt((*kn)[i])), ge_sub(ge_pt((*kn)[i + d]), ge_pt((*kn)[i])), &al);
					o[3 * i + k] = ge_add(ge_mul(ge_sub(ge_pt(1), al), p[3 * (i - 1) + k]), ge_mul(al, p[3 * i + k]));
				}
			}
	}
	for (i = 0; i <= sp; i++) kk[i] = (*kn)[i];
	kk[sp + 1] = t;
	for (i = sp + 1; i < *nk; i++) kk[i + 1] = (*kn)[i];
	free(*pts); free(*kn);
	*pts = q; *kn = kk;
	(*cnt)++;
	(*nk)++;
	return QAWS_STATUS_OK;
}

/* points (rows x cols, u major) transposed */
static ge_iv* ge_transpose(ge_iv const* p, unsigned int rows, unsigned int cols)
{
	ge_iv* q = (ge_iv*)malloc(sizeof(ge_iv) * 3 * rows * cols);
	unsigned int i, j;
	if (q)
		for (i = 0; i < rows; i++)
			for (j = 0; j < cols; j++)
				memcpy(&q[3 * (j * rows + i)], &p[3 * (i * cols + j)], sizeof(ge_iv) * 3);
	return q;
}

/* a B-spline surface as Bezier patches: every interior knot inserted to the
   degree's multiplicity, in v along the rows, then in u */
static qaws_status ge_from_bspline(ge_surf* s, qaws_surface_bspline_impl const* b, unsigned int m, unsigned int n)
{
	unsigned int cu = b->u_count, cv = b->v_count, nku = b->u_knot_count, nkv = b->v_knot_count, i, j, k, mult;
	ge_iv* pts = (ge_iv*)malloc(sizeof(ge_iv) * 3 * cu * cv);
	double *ku = (double*)malloc(sizeof(double) * nku), *kv = (double*)malloc(sizeof(double) * nkv);
	qaws_status st = QAWS_STATUS_OK;
	if (!pts || !ku || !kv)
	{
		free(pts); free(ku); free(kv);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (i = 0; i < 3 * cu * cv; i++) pts[i] = ge_pt((double)b->control_points[i]);
	for (i = 0; i < nku; i++) ku[i] = (double)b->u_knots[i];
	for (i = 0; i < nkv; i++) kv[i] = (double)b->v_knots[i];
	/* v: rows are u, points along a row v */
	for (k = n; k + n + 1 <= nkv && st == QAWS_STATUS_OK; k++)
	{
		double t = kv[k];
		for (mult = 0, j = 0; j < nkv; j++) if (kv[j] == t) mult++;
		while (mult < n && st == QAWS_STATUS_OK)
		{
			st = ge_insert(&pts, cu, &cv, &kv, &nkv, n, t);
			mult++;
		}
	}
	if (st == QAWS_STATUS_OK)
	{
		ge_iv* tp = ge_transpose(pts, cu, cv);
		free(pts);
		pts = tp;
		if (!pts)
			st = QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (k = m; k + m + 1 <= nku && st == QAWS_STATUS_OK; k++)
	{
		double t = ku[k];
		for (mult = 0, j = 0; j < nku; j++) if (ku[j] == t) mult++;
		while (mult < m && st == QAWS_STATUS_OK)
		{
			st = ge_insert(&pts, cv, &cu, &ku, &nku, m, t);
			mult++;
		}
	}
	if (st == QAWS_STATUS_OK)
	{
		/* pts is now cv rows of cu: back to u major */
		ge_iv* tp = ge_transpose(pts, cv, cu);
		free(pts);
		pts = tp;
		if (!pts)
			st = QAWS_STATUS_ALLOCATION_FAILURE;
	}
	/* the patches: the domain's spans of positive length, every knot of
	   multiplicity m (n) now, their points P[i - m .. i] x P[j - n .. j] */
	for (i = m; i < cu && st == QAWS_STATUS_OK; i++)
	{
		if (!(ku[i] < ku[i + 1]) || ku[i] < s->u0 || ku[i + 1] > s->u1)
			continue;
		for (j = n; j < cv && st == QAWS_STATUS_OK; j++)
		{
			ge_iv* cp;
			unsigned int a, c2, iu = i - m, iv0 = j - n;
			if (!(kv[j] < kv[j + 1]) || kv[j] < s->v0 || kv[j + 1] > s->v1)
				continue;
			cp = (ge_iv*)malloc(sizeof(ge_iv) * 3 * (m + 1) * (n + 1));
			if (!cp)
			{
				st = QAWS_STATUS_ALLOCATION_FAILURE;
				break;
			}
			for (a = 0; a <= m; a++)
				for (c2 = 0; c2 <= n; c2++)
					memcpy(&cp[3 * (a * (n + 1) + c2)], &pts[3 * ((iu + a) * cv + iv0 + c2)], sizeof(ge_iv) * 3);
			st = ge_surf_add(s, cp, m, n, ku[i], ku[i + 1], kv[j], kv[j + 1]);
			free(cp);
		}
	}
	free(pts); free(ku); free(kv);
	return st;
}

static qaws_status ge_surf_init(ge_surf* s, qaws_surface const* surf)
{
	qaws_surface_kind kind = qaws_surface_get_kind(surf);
	unsigned int m = surf->u_degree, n = surf->v_degree, i, j, c;
	qaws_status st = QAWS_STATUS_OK;
	memset(s, 0, sizeof(*s));
	s->u0 = (double)surf->u_range.min_value; s->u1 = (double)surf->u_range.max_value;
	s->v0 = (double)surf->v_range.min_value; s->v1 = (double)surf->v_range.max_value;
	s->diag = hypot(s->u1 - s->u0, s->v1 - s->v0);
	if (kind == QAWS_SURFACE_KIND_BEZIER || kind == QAWS_SURFACE_KIND_BILINEAR || kind == QAWS_SURFACE_KIND_BIQUADRATIC)
	{
		ge_iv* cp;
		if (kind == QAWS_SURFACE_KIND_BEZIER)
		{
			qaws_surface_bezier_impl const* b = (qaws_surface_bezier_impl const*)surf->impl;
			m = b->u_count - 1; n = b->v_count - 1;
			if (b->u_count > 8 || b->v_count > 8)
				return QAWS_STATUS_UNSUPPORTED_OPERATION;   /* (binomials to degree 7) */
		}
		cp = (ge_iv*)malloc(sizeof(ge_iv) * 3 * (m + 1) * (n + 1) + sizeof(ge_iv) * 12);
		if (!cp)
			return QAWS_STATUS_ALLOCATION_FAILURE;
		if (kind == QAWS_SURFACE_KIND_BEZIER)
		{
			qaws_surface_bezier_impl const* b = (qaws_surface_bezier_impl const*)surf->impl;
			for (i = 0; i < 3 * b->u_count * b->v_count; i++) cp[i] = ge_pt((double)b->control_points[i]);
		}
		else if (kind == QAWS_SURFACE_KIND_BILINEAR)
		{
			/* the impl's four corners, as qaws_surface_bilinear_desc lays them out */
			qaws_vec3 const* q = (qaws_vec3 const*)surf->impl;
			qaws_vec3 const* order[4];
			m = n = 1;
			order[0] = &q[0]; order[1] = &q[2]; order[2] = &q[1]; order[3] = &q[3];   /* P00, P01, P10, P11 (u major) */
			for (i = 0; i < 4; i++)
			{
				cp[3 * i] = ge_pt((double)order[i]->x); cp[3 * i + 1] = ge_pt((double)order[i]->y); cp[3 * i + 2] = ge_pt((double)order[i]->z);
			}
		}
		else
		{
			/* cp[j * 3 + i]: i along u, j along v */
			qaws_vec3 const* q = (qaws_vec3 const*)surf->impl;
			m = n = 2;
			for (i = 0; i < 3; i++)
				for (j = 0; j < 3; j++)
				{
					qaws_vec3 const* x = &q[j * 3 + i];
					c = 3 * (i * 3 + j);
					cp[c] = ge_pt((double)x->x); cp[c + 1] = ge_pt((double)x->y); cp[c + 2] = ge_pt((double)x->z);
				}
		}
		st = ge_surf_add(s, cp, m, n, s->u0, s->u1, s->v0, s->v1);
		free(cp);
	}
	else if (kind == QAWS_SURFACE_KIND_BSPLINE)
		st = ge_from_bspline(s, (qaws_surface_bspline_impl const*)surf->impl, m, n);
	else
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	if (st == QAWS_STATUS_OK && s->np == 0)
		st = QAWS_STATUS_UNSUPPORTED_OPERATION;
	return st;
}

/* a net's Bernstein form over the parameter box (t, s within [0, 1]) */
static void ge_bern(ge_iv const* q, unsigned int p, unsigned int r, ge_iv t, ge_iv s, ge_iv* out)
{
	static double const bin[8][8] = { { 1 }, { 1, 1 }, { 1, 2, 1 }, { 1, 3, 3, 1 }, { 1, 4, 6, 4, 1 }, { 1, 5, 10, 10, 5, 1 },
		{ 1, 6, 15, 20, 15, 6, 1 }, { 1, 7, 21, 35, 35, 21, 7, 1 } };
	ge_iv bt[64], bs[64], one = ge_pt(1);
	unsigned int i, j, c;
	for (i = 0; i <= p; i++)
		bt[i] = ge_mul(ge_mul(ge_pt(p < 8 ? bin[p][i] : 1), ge_powp(t, i)), ge_powp(ge_sub(one, t), p - i));
	for (j = 0; j <= r; j++)
		bs[j] = ge_mul(ge_mul(ge_pt(r < 8 ? bin[r][j] : 1), ge_powp(s, j)), ge_powp(ge_sub(one, s), r - j));
	for (c = 0; c < 3; c++) out[c] = ge_pt(0);
	for (i = 0; i <= p; i++)
		for (j = 0; j <= r; j++)
		{
			ge_iv w = ge_mul(bt[i], bs[j]);
			for (c = 0; c < 3; c++)
				out[c] = ge_add(out[c], ge_mul(w, q[3 * (i * (r + 1) + j) + c]));
		}
	/* within the net's hull too (the basis sums to one) */
	for (c = 0; c < 3; c++)
	{
		ge_iv h = q[c];
		for (i = 0; i < (p + 1) * (r + 1); i++)
			h = ge_hull(h, q[3 * i + c]);
		out[c] = ge_meet(out[c], h);
	}
}

/* the derivatives of orders 1 to 3 over the box (within the patch): d[k][c]
   for the nets in ge_net_index order */
static void ge_patch_eval(ge_patch const* pt, ge_iv bu, ge_iv bv, ge_iv (*d)[3])
{
	ge_iv lu, lv, t, s, tc, sc, hu, hv;
	double cu = ge_mid(bu), cv = ge_mid(bv);
	unsigned int a, b;
	lu = ge_sub(ge_pt(pt->u1), ge_pt(pt->u0));
	lv = ge_sub(ge_pt(pt->v1), ge_pt(pt->v0));
	ge_div(ge_sub(bu, ge_pt(pt->u0)), lu, &t);
	ge_div(ge_sub(bv, ge_pt(pt->v0)), lv, &s);
	ge_div(ge_sub(ge_pt(cu), ge_pt(pt->u0)), lu, &tc);
	ge_div(ge_sub(ge_pt(cv), ge_pt(pt->v0)), lv, &sc);
	t = ge_meet(t, ge_mk(0, 1)); s = ge_meet(s, ge_mk(0, 1));
	tc = ge_meet(tc, ge_mk(0, 1)); sc = ge_meet(sc, ge_mk(0, 1));
	hu = ge_sub(bu, ge_pt(cu)); hv = ge_sub(bv, ge_pt(cv));
	/* the third order over the box; lower ones by the mean value form from the next */
	for (a = 0; a <= 3; a++)
	{
		unsigned int k = ge_net_index(a, 3 - a);
		if (pt->m < a || pt->n < 3 - a)
			d[k][0] = d[k][1] = d[k][2] = ge_pt(0);
		else
			ge_bern(pt->net[k], pt->m - a, pt->n - (3 - a), t, s, d[k]);
	}
	for (b = 2; b >= 1; b--)
		for (a = 0; a <= b; a++)
		{
			unsigned int k = ge_net_index(a, b - a), c, ku = ge_net_index(a + 1, b - a), kv = ge_net_index(a, b - a + 1);
			ge_iv at[3], box[3];
			int zero = pt->m < a || pt->n < b - a;
			if (zero)
			{
				d[k][0] = d[k][1] = d[k][2] = ge_pt(0);
				continue;
			}
			ge_bern(pt->net[k], pt->m - a, pt->n - (b - a), tc, sc, at);
			ge_bern(pt->net[k], pt->m - a, pt->n - (b - a), t, s, box);
			for (c = 0; c < 3; c++)
				d[k][c] = ge_meet(ge_add(at[c], ge_add(ge_mul(d[ku][c], hu), ge_mul(d[kv][c], hv))), box[c]);
		}
}

/* the derivatives over the box: the hull over the patches it meets; 0 when
   it leaves the domain */
static int ge_eval(ge_surf const* s, ge_iv bu, ge_iv bv, ge_iv (*d)[3])
{
	unsigned int i, k, c, found = 0;
	if (bu.lo < s->u0 || bu.hi > s->u1 || bv.lo < s->v0 || bv.hi > s->v1)
		return 0;
	for (i = 0; i < s->np; i++)
	{
		ge_patch const* pt = &s->p[i];
		ge_iv pu, pv, e[GE_NETS][3];
		if (bu.hi < pt->u0 || bu.lo > pt->u1 || bv.hi < pt->v0 || bv.lo > pt->v1)
			continue;
		pu = ge_mk(bu.lo > pt->u0 ? bu.lo : pt->u0, bu.hi < pt->u1 ? bu.hi : pt->u1);
		pv = ge_mk(bv.lo > pt->v0 ? bv.lo : pt->v0, bv.hi < pt->v1 ? bv.hi : pt->v1);
		ge_patch_eval(pt, pu, pv, e);
		for (k = 1; k < GE_NETS; k++)
			for (c = 0; c < 3; c++)
				d[k][c] = found ? ge_hull(d[k][c], e[k][c]) : e[k][c];
		found = 1;
	}
	return found;
}

/* ------------------------------------------------------------------------ */
/*  The geodesic equation                                                   */
/* ------------------------------------------------------------------------ */

/* over the box y (u, v, u', v'): the acceleration a (with its gradient in
   (u, v)) and the Christoffel symbols G[k][0..2] = G^k_11, G^k_12, G^k_22 */
static int ge_accel(ge_surf const* s, ge_iv const* y, ge_d* acc, ge_iv (*G)[3])
{
	ge_iv d[GE_NETS][3];
	ge_d Su[3], Sv[3], Suu[3], Suv[3], Svv[3], E, F, Gm, det, c[3][2], Gk[2][3];
	unsigned int i, k;
	if (!ge_eval(s, y[0], y[1], d))
		return 0;
	for (i = 0; i < 3; i++)
	{
		Su[i].v = d[1][i]; Su[i].du = d[3][i]; Su[i].dv = d[4][i];
		Sv[i].v = d[2][i]; Sv[i].du = d[4][i]; Sv[i].dv = d[5][i];
		Suu[i].v = d[3][i]; Suu[i].du = d[6][i]; Suu[i].dv = d[7][i];
		Suv[i].v = d[4][i]; Suv[i].du = d[7][i]; Suv[i].dv = d[8][i];
		Svv[i].v = d[5][i]; Svv[i].du = d[8][i]; Svv[i].dv = d[9][i];
	}
	E = ge_ddot(Su, Su); F = ge_ddot(Su, Sv); Gm = ge_ddot(Sv, Sv);
	det = ge_dsub(ge_dmul(E, Gm), ge_dmul(F, F));
	c[0][0] = ge_ddot(Suu, Su); c[0][1] = ge_ddot(Suu, Sv);
	c[1][0] = ge_ddot(Suv, Su); c[1][1] = ge_ddot(Suv, Sv);
	c[2][0] = ge_ddot(Svv, Su); c[2][1] = ge_ddot(Svv, Sv);
	for (i = 0; i < 3; i++)
	{
		if (!ge_ddiv(ge_dsub(ge_dmul(Gm, c[i][0]), ge_dmul(F, c[i][1])), det, &Gk[0][i]) ||
			!ge_ddiv(ge_dsub(ge_dmul(E, c[i][1]), ge_dmul(F, c[i][0])), det, &Gk[1][i]))
			return 0;
	}
	for (k = 0; k < 2; k++)
	{
		/* a^k = -(G^k_11 w1^2 + 2 G^k_12 w1 w2 + G^k_22 w2^2), w plain */
		ge_iv w11 = ge_mul(y[2], y[2]), w12 = ge_scl(ge_mul(y[2], y[3]), 2), w22 = ge_mul(y[3], y[3]);
		ge_d r;
		r.v = ge_neg(ge_add(ge_add(ge_mul(Gk[k][0].v, w11), ge_mul(Gk[k][1].v, w12)), ge_mul(Gk[k][2].v, w22)));
		r.du = ge_neg(ge_add(ge_add(ge_mul(Gk[k][0].du, w11), ge_mul(Gk[k][1].du, w12)), ge_mul(Gk[k][2].du, w22)));
		r.dv = ge_neg(ge_add(ge_add(ge_mul(Gk[k][0].dv, w11), ge_mul(Gk[k][1].dv, w12)), ge_mul(Gk[k][2].dv, w22)));
		acc[k] = r;
		for (i = 0; i < 3; i++) G[k][i] = Gk[k][i].v;
	}
	return 1;
}

/* f(y) over the box */
static int ge_f(ge_surf const* s, ge_iv const* y, ge_iv* out)
{
	ge_d acc[2];
	ge_iv G[2][3];
	if (!ge_accel(s, y, acc, G))
		return 0;
	out[0] = y[2]; out[1] = y[3]; out[2] = acc[0].v; out[3] = acc[1].v;
	return 1;
}

/* y'' over the box: (a, Da), Da^k = a^k_u w1 + a^k_v w2 + a^k_w1 a1 + a^k_w2 a2 */
static int ge_f2(ge_surf const* s, ge_iv const* y, ge_iv* out)
{
	ge_d acc[2];
	ge_iv G[2][3];
	unsigned int k;
	if (!ge_accel(s, y, acc, G))
		return 0;
	out[0] = acc[0].v; out[1] = acc[1].v;
	for (k = 0; k < 2; k++)
	{
		ge_iv aw1 = ge_neg(ge_scl(ge_add(ge_mul(G[k][0], y[2]), ge_mul(G[k][1], y[3])), 2));
		ge_iv aw2 = ge_neg(ge_scl(ge_add(ge_mul(G[k][1], y[2]), ge_mul(G[k][2], y[3])), 2));
		out[2 + k] = ge_add(ge_add(ge_mul(acc[k].du, y[2]), ge_mul(acc[k].dv, y[3])), ge_add(ge_mul(aw1, acc[0].v), ge_mul(aw2, acc[1].v)));
	}
	return 1;
}

/* Df over the box: rows u', v', a^1, a^2 against (u, v, w1, w2) */
static int ge_jac(ge_surf const* s, ge_iv const* y, ge_iv* J)
{
	ge_d acc[2];
	ge_iv G[2][3];
	unsigned int k, i;
	if (!ge_accel(s, y, acc, G))
		return 0;
	for (i = 0; i < 16; i++) J[i] = ge_pt(0);
	J[2] = ge_pt(1); J[7] = ge_pt(1);
	for (k = 0; k < 2; k++)
	{
		ge_iv* row = J + 4 * (2 + k);
		row[0] = acc[k].du;
		row[1] = acc[k].dv;
		row[2] = ge_neg(ge_scl(ge_add(ge_mul(G[k][0], y[2]), ge_mul(G[k][1], y[3])), 2));
		row[3] = ge_neg(ge_scl(ge_add(ge_mul(G[k][1], y[2]), ge_mul(G[k][2], y[3])), 2));
	}
	return 1;
}

/* ------------------------------------------------------------------------ */
/*  Steps                                                                   */
/* ------------------------------------------------------------------------ */

/* the set y^ + A r (A a point matrix, r intervals) as a box */
static void ge_box(double const* yh, double const* A, ge_iv const* r, ge_iv* Y)
{
	unsigned int i, j;
	for (i = 0; i < 4; i++)
	{
		ge_iv x = ge_pt(yh[i]);
		for (j = 0; j < 4; j++)
			x = ge_add(x, ge_mul(ge_pt(A[4 * i + j]), r[j]));
		Y[i] = x;
	}
}

/* Q of the point matrix M (columns by Gram-Schmidt, the widest of M r
   first, as Lohner orders them) */
static void ge_qr(ge_iv const* M, ge_iv const* r, double* Q)
{
	double c[4][4], key[4];
	unsigned int ord[4] = { 0, 1, 2, 3 }, i, j, k;
	for (j = 0; j < 4; j++)
	{
		double n2 = 0;
		for (i = 0; i < 4; i++)
		{
			c[j][i] = ge_mid(M[4 * i + j]);
			n2 += c[j][i] * c[j][i];
		}
		key[j] = sqrt(n2) * ge_w(r[j]);
	}
	for (i = 0; i < 4; i++)
		for (j = i + 1; j < 4; j++)
			if (key[ord[j]] > key[ord[i]]) { unsigned int t = ord[i]; ord[i] = ord[j]; ord[j] = t; }
	for (k = 0; k < 4; k++)
	{
		double v[4], n = 0;
		for (i = 0; i < 4; i++) v[i] = c[ord[k]][i];
		for (j = 0; j < k; j++)
		{
			double d = 0;
			for (i = 0; i < 4; i++) d += v[i] * Q[4 * i + j];
			for (i = 0; i < 4; i++) v[i] -= d * Q[4 * i + j];
		}
		for (i = 0; i < 4; i++) n += v[i] * v[i];
		n = sqrt(n);
		if (!(n > 1e-300))
		{
			/* a column lost: the unit vector most apart from those before */
			for (i = 0; i < 4; i++) v[i] = i == k ? 1 : 0;
			n = 1;
		}
		for (i = 0; i < 4; i++) Q[4 * i + k] = v[i] / n;
	}
}

/* an enclosure of Q^-1 from C = Q^T: |Q^-1 - C| <= |C| d / (1 - d), d = |I - C Q|
   (row sums); 0 when d >= 1/2 */
static int ge_inv(double const* Q, ge_iv* Inv)
{
	ge_iv E[16];
	double d = 0, cn = 0, eps;
	unsigned int i, j, k;
	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
		{
			ge_iv x = ge_pt(i == j ? 1 : 0);
			for (k = 0; k < 4; k++)
				x = ge_sub(x, ge_mul(ge_pt(Q[4 * k + i]), ge_pt(Q[4 * k + j])));
			E[4 * i + j] = x;
		}
	for (i = 0; i < 4; i++)
	{
		double rs = 0, rc = 0;
		for (j = 0; j < 4; j++)
		{
			rs = ge_up(rs + ge_mag(E[4 * i + j]));
			rc = ge_up(rc + fabs(Q[4 * j + i]));
		}
		if (rs > d) d = rs;
		if (rc > cn) cn = rc;
	}
	if (!(d < 0.5))
		return 0;
	eps = ge_up(ge_up(cn * d) / ge_dn(1 - d));
	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
			Inv[4 * i + j] = ge_add(ge_pt(Q[4 * j + i]), ge_mk(-eps, eps));
	return 1;
}

#define GE_STRIDE 10
struct qaws_geodesic_enclosure
{
	double* steps;          /* GE_STRIDE per step: t0, t1, u lo, hi, v lo, hi, u' lo, hi, v' lo, hi */
	unsigned int n, cap;
	ge_iv end[4];
};

static int ge_push(qaws_geodesic_enclosure* e, double t0, double t1, ge_iv const* b)
{
	if (e->n == e->cap)
	{
		unsigned int nc = e->cap ? 2 * e->cap : 256;
		double* g = (double*)realloc(e->steps, sizeof(double) * GE_STRIDE * nc);
		if (!g)
			return 0;
		e->steps = g;
		e->cap = nc;
	}
	e->steps[GE_STRIDE * e->n] = t0; e->steps[GE_STRIDE * e->n + 1] = t1;
	e->steps[GE_STRIDE * e->n + 2] = b[0].lo; e->steps[GE_STRIDE * e->n + 3] = b[0].hi;
	e->steps[GE_STRIDE * e->n + 4] = b[1].lo; e->steps[GE_STRIDE * e->n + 5] = b[1].hi;
	e->steps[GE_STRIDE * e->n + 6] = b[2].lo; e->steps[GE_STRIDE * e->n + 7] = b[2].hi;
	e->steps[GE_STRIDE * e->n + 8] = b[3].lo; e->steps[GE_STRIDE * e->n + 9] = b[3].hi;
	e->n++;
	return 1;
}

/* the a priori box over [0, h] from Y: Y + [0, h] f(B) within B; 0 when
   the iteration does not close (or leaves the domain) */
static int ge_apriori(ge_surf const* s, ge_iv const* Y, double h, ge_iv* B)
{
	ge_iv f[4], H = ge_mk(0, h), C[4];
	unsigned int i, it;
	if (!ge_f(s, Y, f))
		return 0;
	for (i = 0; i < 4; i++)
	{
		ge_iv x = ge_add(Y[i], ge_mul(H, f[i]));
		double r = 0.2 * ge_w(x) + 1e-15 * (1 + ge_mag(x));
		B[i] = ge_mk(ge_dn(x.lo - r), ge_up(x.hi + r));
	}
	for (it = 0; it < 8; it++)
	{
		int inside = 1;
		if (!ge_f(s, B, f))
			return 0;
		for (i = 0; i < 4; i++)
		{
			C[i] = ge_add(Y[i], ge_mul(H, f[i]));
			if (!ge_in(C[i], B[i]))
				inside = 0;
		}
		if (inside)
		{
			memcpy(B, C, sizeof(C));
			return 1;
		}
		for (i = 0; i < 4; i++)
		{
			ge_iv x = ge_hull(B[i], C[i]);
			double r = 0.5 * ge_w(x) + 1e-15 * (1 + ge_mag(x));
			B[i] = ge_mk(ge_dn(x.lo - r), ge_up(x.hi + r));
		}
	}
	return 0;
}

/* the qaws_scalar on the outer side of x */
static double ge_out(double x, int down)
{
	qaws_scalar f = (qaws_scalar)x;
#if QAWS_SCALAR_IS_FLOAT
	if (down && (double)f > x) f = nextafterf(f, -HUGE_VALF);
	if (!down && (double)f < x) f = nextafterf(f, HUGE_VALF);
#else
	(void)down;
#endif
	return (double)f;
}

qaws_status qaws_surface_geodesic_enclose(qaws_surface const* surface, qaws_vec2 p, qaws_vec2 direction, qaws_scalar length,
	qaws_scalar width, qaws_geodesic_enclosure** out_enclosure)
{
	ge_surf s;
	qaws_geodesic_enclosure* e;
	ge_iv Y[4], d[GE_NETS][3], E, F, G, n2, nrm;
	double L = (double)length, t = 0, h, target, rate;
	qaws_status st;
	unsigned int k;
	if (!surface || !out_enclosure || !(length >= 0))
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_enclosure = NULL;
	if ((st = ge_surf_init(&s, surface)) != QAWS_STATUS_OK)
	{
		ge_surf_free(&s);
		return st;
	}
	e = (qaws_geodesic_enclosure*)calloc(1, sizeof(*e));
	if (!e)
	{
		ge_surf_free(&s);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	/* the start: the point, the direction at unit speed in the metric */
	Y[0] = ge_pt((double)p.x); Y[1] = ge_pt((double)p.y);
	if (!ge_eval(&s, Y[0], Y[1], d))
	{
		st = QAWS_STATUS_OUT_OF_RANGE;
		goto done;
	}
	E = ge_add(ge_add(ge_mul(d[1][0], d[1][0]), ge_mul(d[1][1], d[1][1])), ge_mul(d[1][2], d[1][2]));
	F = ge_add(ge_add(ge_mul(d[1][0], d[2][0]), ge_mul(d[1][1], d[2][1])), ge_mul(d[1][2], d[2][2]));
	G = ge_add(ge_add(ge_mul(d[2][0], d[2][0]), ge_mul(d[2][1], d[2][1])), ge_mul(d[2][2], d[2][2]));
	{
		ge_iv dx = ge_pt((double)direction.x), dy = ge_pt((double)direction.y);
		n2 = ge_add(ge_add(ge_mul(E, ge_mul(dx, dx)), ge_scl(ge_mul(F, ge_mul(dx, dy)), 2)), ge_mul(G, ge_mul(dy, dy)));
		if (!(n2.lo > 0))
		{
			st = QAWS_STATUS_INVALID_ARGUMENT;
			goto done;
		}
		nrm = ge_mk(ge_dn(sqrt(n2.lo)), ge_up(sqrt(n2.hi)));
		ge_div(dx, nrm, &Y[2]);
		ge_div(dy, nrm, &Y[3]);
	}
	target = (width > 0 ? (double)width : 1e-6) * s.diag;
	rate = L > 0 ? target / L : target;
	h = L > 0 ? (L < 0.01 * s.diag ? L : 0.01 * s.diag) : 0;
	/* the set as y^ + A r (Lohner): A orthonormal, turned with the flow, so
	   the boxes do not wrap the sliver the set becomes */
	{
		double yh[4], A[16];
		ge_iv r[4];
		for (k = 0; k < 4; k++)
		{
			yh[k] = ge_mid(Y[k]);
			r[k] = ge_sub(Y[k], ge_pt(yh[k]));
		}
		for (k = 0; k < 16; k++) A[k] = k % 5 == 0 ? 1 : 0;
		while (t < L)
		{
			ge_iv B[4], f2[4], J[16], M[16], Inv[16], P[16], z[4], R[4], Tp[4], yp[4], HH, H2, rn[4];
			double Q[16], yn[4], rem = 0, fac = 1;
			unsigned int i, j, q;
			int ok, edge;
			if (h > L - t) h = L - t;
			ge_box(yh, A, r, Y);
			/* (at the domain's edge: the geodesic leaves it) */
			edge = Y[0].lo - s.u0 < 1e-7 * s.diag || s.u1 - Y[0].hi < 1e-7 * s.diag || Y[1].lo - s.v0 < 1e-7 * s.diag || s.v1 - Y[1].hi < 1e-7 * s.diag;
			if (!(h > 1e-14 * (L + s.diag)))
			{
				st = edge ? QAWS_STATUS_OUT_OF_RANGE : QAWS_STATUS_NUMERICAL_FAILURE;
				break;
			}
			for (k = 0; k < 4; k++) yp[k] = ge_pt(yh[k]);
			memcpy(B, Y, sizeof(B));
			ok = ge_apriori(&s, Y, h, B) && ge_f2(&s, B, f2) && ge_jac(&s, Y, J) && ge_f(&s, yp, Tp);
			if (!ok)
			{
				/* the box did not close or left the domain: out at the edge,
				   else a shorter step */
				if (!edge && h > 1e-9 * s.diag)
				{
					h *= 0.5;
					continue;
				}
				st = edge ? QAWS_STATUS_OUT_OF_RANGE : QAWS_STATUS_NUMERICAL_FAILURE;
				break;
			}
			HH = ge_pt(h);
			H2 = ge_mul(ge_mul(HH, HH), ge_pt(0.5));
			for (k = 0; k < 4; k++)
			{
				/* (a velocity's width moves the point by as much per length left) */
				double wk = ge_w(R[k] = ge_mul(H2, f2[k])) * (k < 2 ? 1 : L - t);
				if (wk > rem) rem = wk;
			}
			/* the step's growth against the rate the end's width allows */
			{
				/* the remainder goes as h^3, the allowance as h: the step scaled by
				   the square root of their ratio (rejected when well past it) */
				double g = rem > 0 ? 0.9 * sqrt(rate * h / rem) : 2;
				fac = g < 0.3 ? 0.3 : g > 2 ? 2 : g;
				if (rem > 2 * rate * h && h > 1e-6 * s.diag)
				{
					h *= fac;
					continue;
				}
			}
			/* y(t + h) = T(y) + h^2/2 y''(B), T(y) = y + h f(y); T(y) - T(y^) in
			   (I + h Df(Y)) (y - y^) = (I + h Df(Y)) A r */
			for (i = 0; i < 4; i++)
			{
				z[i] = ge_add(ge_add(yp[i], ge_mul(HH, Tp[i])), R[i]);
				for (j = 0; j < 4; j++)
				{
					ge_iv m = ge_pt(0);
					for (q = 0; q < 4; q++)
					{
						ge_iv jq = ge_mul(HH, J[4 * i + q]);
						if (i == q) jq = ge_add(jq, ge_pt(1));
						m = ge_add(m, ge_mul(jq, ge_pt(A[4 * q + j])));
					}
					M[4 * i + j] = m;
				}
			}
			for (i = 0; i < 4; i++)
			{
				yn[i] = ge_mid(z[i]);
				z[i] = ge_sub(z[i], ge_pt(yn[i]));
			}
			ge_qr(M, r, Q);
			if (!ge_inv(Q, Inv))
			{
				st = QAWS_STATUS_NUMERICAL_FAILURE;
				break;
			}
			/* r' = (Q^-1 M) r + Q^-1 z */
			for (i = 0; i < 4; i++)
				for (j = 0; j < 4; j++)
				{
					ge_iv m = ge_pt(0);
					for (q = 0; q < 4; q++)
						m = ge_add(m, ge_mul(Inv[4 * i + q], M[4 * q + j]));
					P[4 * i + j] = m;
				}
			for (i = 0; i < 4; i++)
			{
				ge_iv x = ge_pt(0);
				for (j = 0; j < 4; j++)
					x = ge_add(x, ge_add(ge_mul(P[4 * i + j], r[j]), ge_mul(Inv[4 * i + j], z[j])));
				rn[i] = x;
			}
			if (!ge_push(e, t, t + h, B))
			{
				st = QAWS_STATUS_ALLOCATION_FAILURE;
				break;
			}
			memcpy(yh, yn, sizeof(yn));
			memcpy(A, Q, sizeof(Q));
			memcpy(r, rn, sizeof(rn));
			t += h;
			h *= fac;
		}
		ge_box(yh, A, r, Y);
	}
	if (st == QAWS_STATUS_OK && L == 0 && !ge_push(e, 0, 0, Y))
		st = QAWS_STATUS_ALLOCATION_FAILURE;
	memcpy(e->end, Y, sizeof(Y));
done:
	ge_surf_free(&s);
	if (st != QAWS_STATUS_OK)
	{
		qaws_geodesic_enclosure_destroy(e);
		return st;
	}
	*out_enclosure = e;
	return QAWS_STATUS_OK;
}

void qaws_geodesic_enclosure_destroy(qaws_geodesic_enclosure* enclosure)
{
	if (!enclosure)
		return;
	free(enclosure->steps);
	free(enclosure);
}

unsigned int qaws_geodesic_enclosure_get_step_count(qaws_geodesic_enclosure const* enclosure)
{
	return enclosure ? enclosure->n : 0;
}

qaws_status qaws_geodesic_enclosure_get_step(qaws_geodesic_enclosure const* enclosure, unsigned int index, qaws_scalar* out_t0,
	qaws_scalar* out_t1, qaws_vec2* out_lo, qaws_vec2* out_hi)
{
	double const* q;
	if (!enclosure || index >= enclosure->n)
		return QAWS_STATUS_INVALID_ARGUMENT;
	q = enclosure->steps + GE_STRIDE * index;
	if (out_t0) *out_t0 = (qaws_scalar)ge_out(q[0], 1);
	if (out_t1) *out_t1 = (qaws_scalar)ge_out(q[1], 0);
	if (out_lo) { out_lo->x = (qaws_scalar)ge_out(q[2], 1); out_lo->y = (qaws_scalar)ge_out(q[4], 1); }
	if (out_hi) { out_hi->x = (qaws_scalar)ge_out(q[3], 0); out_hi->y = (qaws_scalar)ge_out(q[5], 0); }
	return QAWS_STATUS_OK;
}

qaws_status qaws_geodesic_enclosure_get_end(qaws_geodesic_enclosure const* enclosure, qaws_vec2* out_point_lo, qaws_vec2* out_point_hi,
	qaws_vec2* out_velocity_lo, qaws_vec2* out_velocity_hi)
{
	ge_iv const* y;
	if (!enclosure)
		return QAWS_STATUS_INVALID_ARGUMENT;
	y = enclosure->end;
	if (out_point_lo) { out_point_lo->x = (qaws_scalar)ge_out(y[0].lo, 1); out_point_lo->y = (qaws_scalar)ge_out(y[1].lo, 1); }
	if (out_point_hi) { out_point_hi->x = (qaws_scalar)ge_out(y[0].hi, 0); out_point_hi->y = (qaws_scalar)ge_out(y[1].hi, 0); }
	if (out_velocity_lo) { out_velocity_lo->x = (qaws_scalar)ge_out(y[2].lo, 1); out_velocity_lo->y = (qaws_scalar)ge_out(y[3].lo, 1); }
	if (out_velocity_hi) { out_velocity_hi->x = (qaws_scalar)ge_out(y[2].hi, 0); out_velocity_hi->y = (qaws_scalar)ge_out(y[3].hi, 0); }
	return QAWS_STATUS_OK;
}

/* ------------------------------------------------------------------------ */
/*  Crossings                                                               */
/* ------------------------------------------------------------------------ */

static void ge_step_box(qaws_geodesic_enclosure const* e, unsigned int i, ge_iv* x, ge_iv* w)
{
	double const* q = e->steps + GE_STRIDE * i;
	x[0] = ge_mk(q[2], q[3]); x[1] = ge_mk(q[4], q[5]);
	if (w) { w[0] = ge_mk(q[6], q[7]); w[1] = ge_mk(q[8], q[9]); }
}

static int ge_meets(ge_iv const* a, ge_iv const* b)
{
	return a[0].lo <= b[0].hi && b[0].lo <= a[0].hi && a[1].lo <= b[1].hi && b[1].lo <= a[1].hi;
}

/* cross(f, a - b): its sign when certain (else 0) */
static int ge_side(double const* f, ge_iv const* a, ge_iv const* b)
{
	ge_iv dx = ge_sub(a[0], b[0]), dy = ge_sub(a[1], b[1]);
	ge_iv c = ge_sub(ge_mul(ge_pt(f[0]), dy), ge_mul(ge_pt(f[1]), dx));
	return c.lo > 0 ? 1 : c.hi < 0 ? -1 : 0;
}

/* one stretch where the boxes meet: A's steps [i0, i1], B's [j0, j1], those
   ranges a step or two past the meeting; 1 when exactly one crossing is
   proven there */
static int ge_cross_one(qaws_geodesic_enclosure const* a, qaws_geodesic_enclosure const* b, unsigned int i0, unsigned int i1,
	unsigned int j0, unsigned int j1)
{
	ge_iv WA[2], WB[2], x[2], w[2], y[2], c;
	double e[2], f[2];
	unsigned int i, j, k;
	int s0 = 0, s1 = 0, t0 = 0, t1 = 0;
	/* transversal throughout: no direction of one along any of the other,
	   so no second crossing (the chord between two would lie in both cones) */
	for (i = i0; i <= i1; i++)
	{
		ge_step_box(a, i, x, w);
		for (k = 0; k < 2; k++) WA[k] = i == i0 ? w[k] : ge_hull(WA[k], w[k]);
	}
	for (j = j0; j <= j1; j++)
	{
		ge_step_box(b, j, x, w);
		for (k = 0; k < 2; k++) WB[k] = j == j0 ? w[k] : ge_hull(WB[k], w[k]);
	}
	c = ge_sub(ge_mul(WA[0], WB[1]), ge_mul(WA[1], WB[0]));
	if (!(c.lo > 0 || c.hi < 0))
		return 0;
	/* Miranda on F(s, t) = A(s) - B(t) seen across B (along f) and across A
	   (along e): opposite signs on opposite faces, a zero in between */
	e[0] = ge_mid(WA[0]); e[1] = ge_mid(WA[1]);
	f[0] = ge_mid(WB[0]); f[1] = ge_mid(WB[1]);
	for (j = j0; j <= j1; j++)
	{
		int p, q;
		ge_step_box(b, j, y, NULL);
		ge_step_box(a, i0, x, NULL);
		p = ge_side(f, x, y);
		ge_step_box(a, i1, x, NULL);
		q = ge_side(f, x, y);
		if (!p || !q || (s0 && p != s0) || (s1 && q != s1))
			return 0;
		s0 = p; s1 = q;
	}
	for (i = i0; i <= i1; i++)
	{
		int p, q;
		ge_step_box(a, i, x, NULL);
		ge_step_box(b, j0, y, NULL);
		p = ge_side(e, x, y);
		ge_step_box(b, j1, y, NULL);
		q = ge_side(e, x, y);
		if (!p || !q || (t0 && p != t0) || (t1 && q != t1))
			return 0;
		t0 = p; t1 = q;
	}
	return s0 == -s1 && t0 == -t1;
}

qaws_status qaws_geodesic_enclosure_cross(qaws_geodesic_enclosure const* a, qaws_geodesic_enclosure const* b, unsigned int capacity,
	qaws_vec2* out_lo, qaws_vec2* out_hi, unsigned int* out_count)
{
	enum { CH = 32, EXT = 2, GAP = 2 * EXT + 1 };
	unsigned int nc, i, j, k, count = 0, ci0 = 0, ci1 = 0, cj0 = 0, cj1 = 0, open = 0;
	unsigned int *jlo, *jhi;
	ge_iv *cb, *mh, hx[2];
	qaws_status st = QAWS_STATUS_OK;
	if (!a || !b || !out_count || (capacity && (!out_lo || !out_hi)))
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_count = 0;
	nc = (b->n + CH - 1) / CH;
	cb = (ge_iv*)malloc(sizeof(ge_iv) * 2 * (nc + 1));
	jlo = (unsigned int*)malloc(sizeof(unsigned int) * (a->n + 1));
	jhi = (unsigned int*)malloc(sizeof(unsigned int) * (a->n + 1));
	mh = (ge_iv*)malloc(sizeof(ge_iv) * 2 * (a->n + 1));
	if (!cb || !jlo || !jhi || !mh)
	{
		free(cb); free(jlo); free(jhi); free(mh);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	/* B's steps by chunks, their bounds */
	for (k = 0; k < nc; k++)
		for (j = k * CH; j < b->n && j < (k + 1) * CH; j++)
		{
			ge_iv y[2];
			ge_step_box(b, j, y, NULL);
			cb[2 * k] = j == k * CH ? y[0] : ge_hull(cb[2 * k], y[0]);
			cb[2 * k + 1] = j == k * CH ? y[1] : ge_hull(cb[2 * k + 1], y[1]);
		}
	/* each of A's steps: the range of B's steps its box meets (jlo > jhi:
	   none) and the hull of the common parts */
	for (i = 0; i < a->n; i++)
	{
		ge_iv x[2];
		jlo[i] = b->n; jhi[i] = 0;
		ge_step_box(a, i, x, NULL);
		for (k = 0; k < nc; k++)
		{
			if (!ge_meets(x, &cb[2 * k]))
				continue;
			for (j = k * CH; j < b->n && j < (k + 1) * CH; j++)
			{
				ge_iv y[2], m0, m1;
				ge_step_box(b, j, y, NULL);
				if (!ge_meets(x, y))
					continue;
				m0 = ge_mk(x[0].lo > y[0].lo ? x[0].lo : y[0].lo, x[0].hi < y[0].hi ? x[0].hi : y[0].hi);
				m1 = ge_mk(x[1].lo > y[1].lo ? x[1].lo : y[1].lo, x[1].hi < y[1].hi ? x[1].hi : y[1].hi);
				if (jlo[i] > jhi[i] && jlo[i] == b->n) { mh[2 * i] = m0; mh[2 * i + 1] = m1; }
				else { mh[2 * i] = ge_hull(mh[2 * i], m0); mh[2 * i + 1] = ge_hull(mh[2 * i + 1], m1); }
				if (j < jlo[i]) jlo[i] = j;
				if (j > jhi[i]) jhi[i] = j;
			}
		}
	}
	/* stretches of meeting steps (near in both curves' steps), each certified */
	for (i = 0; i <= a->n && st == QAWS_STATUS_OK; i++)
	{
		int met = i < a->n && jlo[i] <= jhi[i];
		if (met && open && i <= ci1 + GAP && jlo[i] + GAP >= cj0 && jhi[i] <= cj1 + GAP)
		{
			ci1 = i;
			if (jlo[i] < cj0) cj0 = jlo[i];
			if (jhi[i] > cj1) cj1 = jhi[i];
			hx[0] = ge_hull(hx[0], mh[2 * i]); hx[1] = ge_hull(hx[1], mh[2 * i + 1]);
			continue;
		}
		if (open && (i == a->n || met || i > ci1 + GAP))
		{
			unsigned int i0 = ci0 > EXT ? ci0 - EXT : 0, i1 = ci1 + EXT < a->n ? ci1 + EXT : a->n - 1;
			unsigned int j0 = cj0 > EXT ? cj0 - EXT : 0, j1 = cj1 + EXT < b->n ? cj1 + EXT : b->n - 1;
			if (!ge_cross_one(a, b, i0, i1, j0, j1))
				st = QAWS_STATUS_CERTIFICATION_FAILED;
			else
			{
				if (count < capacity)
				{
					out_lo[count].x = (qaws_scalar)ge_out(hx[0].lo, 1); out_lo[count].y = (qaws_scalar)ge_out(hx[1].lo, 1);
					out_hi[count].x = (qaws_scalar)ge_out(hx[0].hi, 0); out_hi[count].y = (qaws_scalar)ge_out(hx[1].hi, 0);
				}
				count++;
			}
			open = 0;
		}
		if (met && !open)
		{
			open = 1;
			ci0 = ci1 = i;
			cj0 = jlo[i]; cj1 = jhi[i];
			hx[0] = mh[2 * i]; hx[1] = mh[2 * i + 1];
		}
	}
	free(cb); free(jlo); free(jhi); free(mh);
	*out_count = count;
	return st;
}
