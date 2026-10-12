/*
 * Mean value coordinates over regions bounded by curves (qaws_mean_value.h).
 *
 * Every boundary curve becomes rational Bezier pieces whose homogeneous
 * control points are combinations of the curve's own control points (the
 * handles): Bezier extraction for B-splines and NURBS, rational quadratic
 * pieces for arcs. On a piece of degree d with Bernstein weights w_j,
 *
 *     handle i's basis = sum_j a_ji B_j(u) / W(u),  W = sum_j w_j B_j,
 *
 * so the integrals J_j = integral B_j k / W du give every handle's
 * coordinate on the piece at once. Lines (degree 1, equal weights) have them
 * in closed form (Floater 2003): J_0 = tan(alpha / 2) / (r_0 w),
 * J_1 = tan(alpha / 2) / (r_1 w), alpha the signed angle the segment
 * subtends at x. Curved pieces are integrated by adaptive Gauss-Kronrod
 * (7-15), on intervals graded away from the piece's point nearest x, where
 * the kernel peaks with a width of the distance.
 */
#include "qaws_mean_value.h"
#include "qaws_inspect.h"
#include "qaws_eval.h"
#include "qaws_curve.h"
#include "internal/qaws_internal_types.h"
#include "internal/qaws_internal_kinds.h"
#include "internal/qaws_internal_basis.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MV_MAXD 12              /* highest piece degree */
#define MV_PI 3.14159265358979323846
#if QAWS_SCALAR_IS_FLOAT
#define MV_MERGE 1e-6           /* points this near (of the extent) are one handle: the float rounding */
#else
#define MV_MERGE 1e-12
#endif

typedef struct mv_piece
{
	unsigned int deg;
	double q[MV_MAXD + 1][2];   /* Bezier control points (Cartesian) */
	double w[MV_MAXD + 1];      /* their weights */
	unsigned int nh;            /* the handles it moves */
	unsigned int h[MV_MAXD + 1];
	double a[(MV_MAXD + 1) * (MV_MAXD + 1)];   /* a[j * nh + i]: handle i's share of Bernstein j */
	int line;                   /* degree 1, equal weights */
	double sign;                /* +1 or -1: its loop turned outside counter-clockwise, holes clockwise */
	unsigned int path, curve;
	int arc;                    /* the curve parameter: t0 + u (t1 - t0), or along an arc */
	double t0, t1;
	double m[2], rho;           /* arc: centre, radius (t = t0 + rho * angle from q[0]) */
} mv_piece;

struct qaws_mean_value
{
	mv_piece* pc;
	unsigned int npc, cappc;
	double* hp;                 /* handles, x y */
	unsigned int nh, caph;
	double ext, tol, eps_on;
	qaws_batch_executor const* executor;
};

static int mv_grow(void** p, unsigned int* cap, unsigned int need, size_t size)
{
	unsigned int n;
	void* q;
	if (need <= *cap)
		return 1;
	n = *cap ? *cap : 16;
	while (n < need) n *= 2;
	q = realloc(*p, size * n);
	if (!q)
		return 0;
	*p = q;
	*cap = n;
	return 1;
}

/* ------------------------------------------------------------------------ */
/*  Pieces                                                                  */
/* ------------------------------------------------------------------------ */

/* Bernstein basis of degree n at u, and its derivative */
static void mv_bernstein(unsigned int n, double u, double* b, double* db)
{
	double t[MV_MAXD + 2];
	unsigned int j, k;
	/* degree n - 1 first (for the derivative), then n */
	for (j = 0; j <= n; j++) t[j] = 0;
	t[0] = 1;
	for (k = 1; k < n; k++)
		for (j = k; j + 1 > 0; j--)
			t[j] = (j ? u * t[j - 1] : 0) + (1 - u) * t[j];
	if (db)
	{
		for (j = 0; j <= n; j++)
			db[j] = n == 0 ? 0 : n * ((j ? t[j - 1] : 0) - (j < n ? t[j] : 0));
	}
	if (n == 0)
	{
		b[0] = 1;
		return;
	}
	for (j = 0; j <= n; j++)
		b[j] = (j ? u * t[j - 1] : 0) + (j < n ? (1 - u) * t[j] : 0);
}

/* the piece's point and derivative at u, and the Bernstein values over W */
static void mv_eval(mv_piece const* p, double u, double* c, double* dc, double* bw)
{
	double b[MV_MAXD + 1], db[MV_MAXD + 1], X[2] = { 0, 0 }, dX[2] = { 0, 0 }, W = 0, dW = 0;
	unsigned int j;
	mv_bernstein(p->deg, u, b, db);
	for (j = 0; j <= p->deg; j++)
	{
		double wq0 = p->w[j] * p->q[j][0], wq1 = p->w[j] * p->q[j][1];
		X[0] += wq0 * b[j]; X[1] += wq1 * b[j];
		dX[0] += wq0 * db[j]; dX[1] += wq1 * db[j];
		W += p->w[j] * b[j];
		dW += p->w[j] * db[j];
	}
	c[0] = X[0] / W;
	c[1] = X[1] / W;
	if (dc)
	{
		dc[0] = (dX[0] - c[0] * dW) / W;
		dc[1] = (dX[1] - c[1] * dW) / W;
	}
	if (bw)
		for (j = 0; j <= p->deg; j++)
			bw[j] = b[j] / W;
}

static unsigned int mv_handle(qaws_mean_value* m, double x, double y, unsigned int merge_with)
{
	if (merge_with != 0xFFFFFFFFu && hypot(m->hp[2 * merge_with] - x, m->hp[2 * merge_with + 1] - y) <= MV_MERGE * (m->ext > 0 ? m->ext : 1) + 1e-300)
		return merge_with;
	if (!mv_grow((void**)&m->hp, &m->caph, m->nh + 1, 2 * sizeof(double)))
		return 0xFFFFFFFFu;
	m->hp[2 * m->nh] = x;
	m->hp[2 * m->nh + 1] = y;
	return m->nh++;
}

static mv_piece* mv_piece_new(qaws_mean_value* m)
{
	mv_piece* p;
	if (!mv_grow((void**)&m->pc, &m->cappc, m->npc + 1, sizeof(mv_piece)))
		return NULL;
	p = &m->pc[m->npc++];
	memset(p, 0, sizeof(*p));
	p->sign = 1;
	return p;
}

/* the last handle of the loop so far (to share with the next curve) */
static unsigned int mv_last(qaws_mean_value const* m, unsigned int first_handle)
{
	return m->nh > first_handle ? m->nh - 1 : 0xFFFFFFFFu;
}

/* a B-spline or NURBS (or Bezier, as one span) of control points cp, weights
   w (NULL: 1), knots kn: one piece per span */
static qaws_status mv_add_spline(qaws_mean_value* m, double const* cp, double const* w, unsigned int n, unsigned int p,
	double const* kn, unsigned int nk, unsigned int path, unsigned int curve, unsigned int loop_first)
{
	unsigned int* hid;
	unsigned int i, s;
	if (p > MV_MAXD || p == 0 || n < p + 1)
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	hid = (unsigned int*)malloc(sizeof(unsigned int) * n);
	if (!hid)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; i < n; i++)
	{
		hid[i] = mv_handle(m, cp[2 * i], cp[2 * i + 1], i == 0 ? mv_last(m, loop_first) : 0xFFFFFFFFu);
		if (hid[i] == 0xFFFFFFFFu)
		{
			free(hid);
			return QAWS_STATUS_ALLOCATION_FAILURE;
		}
	}
	for (s = p; s < n; s++)
	{
		double ka = kn[s], kb = kn[s + 1], M[(MV_MAXD + 1) * (MV_MAXD + 1)], N[(MV_MAXD + 1) * (MV_MAXD + 1)], C[(MV_MAXD + 1) * (MV_MAXD + 1)];
		unsigned int k, j, r, nb = p + 1;
		mv_piece* pc;
		if (!(kb > ka))
			continue;
		/* the span's B-spline basis in Bernstein form: Bernstein values M and
		   B-spline values N at p + 1 points, C = M^-1 N */
		for (k = 0; k < nb; k++)
		{
			double u = (double)k / p, b[MV_MAXD + 1], nv[MV_MAXD + 1];
			qaws_scalar kq[256], nq[MV_MAXD + 1];
			mv_bernstein(p, u, b, NULL);
			for (j = 0; j < nb; j++) M[k * nb + j] = b[j];
			if (nk <= 256)
			{
				for (j = 0; j < nk; j++) kq[j] = (qaws_scalar)kn[j];
				qaws_internal_bspline_basis(kq, nk, p, s, (qaws_scalar)(ka + u * (kb - ka)), nq);
				for (j = 0; j < nb; j++) nv[j] = (double)nq[j];
			}
			else
			{
				/* Cox-de Boor on the double knots */
				double left[MV_MAXD + 2], right[MV_MAXD + 2], t = ka + u * (kb - ka);
				nv[0] = 1;
				for (j = 1; j <= p; j++)
				{
					double saved = 0;
					left[j] = t - kn[s + 1 - j];
					right[j] = kn[s + j] - t;
					for (r = 0; r < j; r++)
					{
						double tmp = nv[r] / (right[r + 1] + left[j - r]);
						nv[r] = saved + right[r + 1] * tmp;
						saved = left[j - r] * tmp;
					}
					nv[j] = saved;
				}
			}
			for (j = 0; j < nb; j++) N[k * nb + j] = nv[j];
		}
		/* Gaussian elimination with partial pivoting: M C = N */
		memcpy(C, N, sizeof(double) * nb * nb);
		for (k = 0; k < nb; k++)
		{
			unsigned int piv = k;
			double d;
			for (r = k + 1; r < nb; r++)
				if (fabs(M[r * nb + k]) > fabs(M[piv * nb + k])) piv = r;
			if (piv != k)
				for (j = 0; j < nb; j++)
				{
					double t1 = M[k * nb + j], t2 = C[k * nb + j];
					M[k * nb + j] = M[piv * nb + j]; M[piv * nb + j] = t1;
					C[k * nb + j] = C[piv * nb + j]; C[piv * nb + j] = t2;
				}
			d = M[k * nb + k];
			for (r = 0; r < nb; r++)
			{
				double f;
				if (r == k) continue;
				f = M[r * nb + k] / d;
				for (j = 0; j < nb; j++)
				{
					M[r * nb + j] -= f * M[k * nb + j];
					C[r * nb + j] -= f * C[k * nb + j];
				}
			}
		}
		for (k = 0; k < nb; k++)
			for (j = 0; j < nb; j++)
				C[k * nb + j] /= M[k * nb + k];
		/* C[j * nb + i]: B-spline function s - p + i's j-th Bernstein coefficient */
		pc = mv_piece_new(m);
		if (!pc)
		{
			free(hid);
			return QAWS_STATUS_ALLOCATION_FAILURE;
		}
		pc->deg = p;
		pc->nh = nb;
		pc->path = path; pc->curve = curve;
		pc->t0 = ka; pc->t1 = kb;
		for (i = 0; i < nb; i++)
			pc->h[i] = hid[s - p + i];
		for (j = 0; j < nb; j++)
		{
			double wsum = 0, x = 0, y = 0;
			for (i = 0; i < nb; i++)
			{
				double wi = w ? w[s - p + i] : 1, aji = C[j * nb + i] * wi;
				if (fabs(aji) < 1e-15) aji = 0;
				pc->a[j * nb + i] = aji;
				wsum += aji;
				x += aji * cp[2 * (s - p + i)];
				y += aji * cp[2 * (s - p + i) + 1];
			}
			pc->w[j] = wsum;
			pc->q[j][0] = wsum != 0 ? x / wsum : 0;
			pc->q[j][1] = wsum != 0 ? y / wsum : 0;
		}
		pc->line = p == 1 && pc->w[0] == pc->w[1];
	}
	free(hid);
	return QAWS_STATUS_OK;
}

static qaws_status mv_add_curve(qaws_mean_value* m, qaws_curve const* c, unsigned int path, unsigned int curve, unsigned int loop_first)
{
	qaws_curve_kind kind = qaws_curve_get_kind(c);
	qaws_status st = QAWS_STATUS_OK;
	if (qaws_curve_get_dimension(c) != QAWS_DIMENSION_2D)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (kind == QAWS_CURVE_KIND_ARC)
	{
		qaws_arc_impl const* impl = (qaws_arc_impl const*)c->impl;
		unsigned int sgi;
		double base = 0;
		for (sgi = 0; sgi < impl->segment_count; sgi++)
		{
			qaws_arc_segment const* g = &impl->segments[sgi];
			double sw = (double)g->angle_end - (double)g->angle_start, r = (double)g->radius, a = (double)g->angle_start;
			double cx = (double)g->center[0], cy = (double)g->center[1], phi;
			unsigned int k, kk;
			if (!(fabs(sw) > 0) || !(r > 0))
				continue;
			k = (unsigned int)ceil(fabs(sw) / (0.5 * MV_PI) - 1e-6);   /* (a float full turn rounds above 2 pi) */
			if (k < 1) k = 1;
			phi = sw / k;
			for (kk = 0; kk < k; kk++)
			{
				double a0 = a + kk * phi, a1 = a0 + phi, am = a0 + 0.5 * phi, cw = cos(0.5 * phi);
				mv_piece* pc = mv_piece_new(m);
				unsigned int h0;
				if (!pc)
					return QAWS_STATUS_ALLOCATION_FAILURE;
				pc->deg = 2;
				pc->q[0][0] = cx + r * cos(a0); pc->q[0][1] = cy + r * sin(a0);
				pc->q[1][0] = cx + r / cw * cos(am); pc->q[1][1] = cy + r / cw * sin(am);
				pc->q[2][0] = cx + r * cos(a1); pc->q[2][1] = cy + r * sin(a1);
				pc->w[0] = 1; pc->w[1] = cw; pc->w[2] = 1;
				h0 = mv_handle(m, pc->q[0][0], pc->q[0][1], mv_last(m, loop_first));
				pc->h[0] = h0;
				pc->h[1] = mv_handle(m, pc->q[1][0], pc->q[1][1], 0xFFFFFFFFu);
				pc->h[2] = mv_handle(m, pc->q[2][0], pc->q[2][1], 0xFFFFFFFFu);
				if (pc->h[0] == 0xFFFFFFFFu || pc->h[1] == 0xFFFFFFFFu || pc->h[2] == 0xFFFFFFFFu)
					return QAWS_STATUS_ALLOCATION_FAILURE;
				pc->nh = 3;
				pc->a[0] = 1; pc->a[4] = cw; pc->a[8] = 1;
				pc->path = path; pc->curve = curve;
				pc->arc = 1;
				pc->t0 = base + r * fabs(phi) * kk;
				pc->m[0] = cx; pc->m[1] = cy; pc->rho = r;
			}
			base += r * fabs(sw);
		}
		return QAWS_STATUS_OK;
	}
	if (kind == QAWS_CURVE_KIND_BEZIER || kind == QAWS_CURVE_KIND_RATIONAL_BEZIER || kind == QAWS_CURVE_KIND_BSPLINE ||
		(kind == QAWS_CURVE_KIND_NURBS && !qaws_curve_is_periodic(c)))
	{
		unsigned int n = 0, nk = 0, nw = 0, i, deg = qaws_curve_get_degree(c);
		qaws_scalar* buf;
		double *cp, *kn, *w = NULL;
		qaws_range rg = qaws_curve_get_parameter_range(c);
		qaws_scalar probe[2];
		qaws_curve_get_control_points(c, probe, 0, &n);       /* the counts */
		if (kind == QAWS_CURVE_KIND_BSPLINE || kind == QAWS_CURVE_KIND_NURBS)
			qaws_bspline_get_knots(c, probe, 0, &nk);
		else
			nk = 2 * (deg + 1);
		if (n == 0 || n > 100000 || nk > 200000)
			return QAWS_STATUS_INVALID_ARGUMENT;
		buf = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (2 * n + nk + n));
		cp = (double*)malloc(sizeof(double) * (2 * n + nk + n));
		if (!buf || !cp)
		{
			free(buf); free(cp);
			return QAWS_STATUS_ALLOCATION_FAILURE;
		}
		kn = cp + 2 * n;
		st = qaws_curve_get_control_points(c, buf, n, &n);
		for (i = 0; i < 2 * n; i++) cp[i] = (double)buf[i];
		if (st == QAWS_STATUS_OK && (kind == QAWS_CURVE_KIND_BSPLINE || kind == QAWS_CURVE_KIND_NURBS))
		{
			st = qaws_bspline_get_knots(c, buf, nk, &nk);
			for (i = 0; i < nk; i++) kn[i] = (double)buf[i];
		}
		else
		{
			/* a Bezier: one span */
			for (i = 0; i <= deg; i++)
			{
				kn[i] = (double)rg.min_value;
				kn[deg + 1 + i] = (double)rg.max_value;
			}
		}
		if (st == QAWS_STATUS_OK && kind == QAWS_CURVE_KIND_NURBS)
		{
			w = kn + nk;
			st = qaws_nurbs_get_weights(c, buf, n, &nw);
			for (i = 0; i < nw && i < n; i++) w[i] = (double)buf[i];
		}
		else if (st == QAWS_STATUS_OK && kind == QAWS_CURVE_KIND_RATIONAL_BEZIER)
		{
			qaws_rational_bezier_impl const* impl = (qaws_rational_bezier_impl const*)c->impl;
			w = kn + nk;
			for (i = 0; i < n; i++) w[i] = (double)impl->weights[i];
		}
		if (st == QAWS_STATUS_OK)
			st = mv_add_spline(m, cp, w, n, deg, kn, nk, path, curve, loop_first);
		free(buf); free(cp);
		return st;
	}
	{
		/* any other kind: a cubic B-spline fitted to it (exact for cubic spans) */
		qaws_curve* fit = NULL;
		qaws_range rg = qaws_curve_get_parameter_range(c);
		st = qaws_internal_curve_fit_cubic(c, rg.min_value, rg.max_value, &fit);
		if (st == QAWS_STATUS_OK)
			st = mv_add_curve(m, fit, path, curve, loop_first);
		qaws_curve_destroy(fit);
		return st;
	}
}

/* ------------------------------------------------------------------------ */
/*  Loops: orientation and nesting                                          */
/* ------------------------------------------------------------------------ */

#define MV_LOOP_SAMPLES 24

/* twice the signed area of pieces [a, b), and the winding number round pt */
static void mv_loop_measure(qaws_mean_value const* m, unsigned int a, unsigned int b, double const* pt, double* area2, double* wind)
{
	unsigned int i, k;
	double prev[2] = { 0, 0 }, first[2] = { 0, 0 }, ang = 0;
	int have = 0;
	*area2 = 0;
	for (i = a; i < b; i++)
		for (k = 0; k < MV_LOOP_SAMPLES; k++)
		{
			double c[2];
			mv_eval(&m->pc[i], (double)k / MV_LOOP_SAMPLES, c, NULL, NULL);
			if (have)
			{
				*area2 += prev[0] * c[1] - c[0] * prev[1];
				if (pt)
					ang += atan2((prev[0] - pt[0]) * (c[1] - pt[1]) - (prev[1] - pt[1]) * (c[0] - pt[0]),
						(prev[0] - pt[0]) * (c[0] - pt[0]) + (prev[1] - pt[1]) * (c[1] - pt[1]));
			}
			else
			{
				first[0] = c[0]; first[1] = c[1];
				have = 1;
			}
			prev[0] = c[0]; prev[1] = c[1];
		}
	if (have)
	{
		*area2 += prev[0] * first[1] - first[0] * prev[1];
		if (pt)
			ang += atan2((prev[0] - pt[0]) * (first[1] - pt[1]) - (prev[1] - pt[1]) * (first[0] - pt[0]),
				(prev[0] - pt[0]) * (first[0] - pt[0]) + (prev[1] - pt[1]) * (first[1] - pt[1]));
	}
	*wind = ang / (2 * MV_PI);
}

/* ------------------------------------------------------------------------ */
/*  Integration                                                             */
/* ------------------------------------------------------------------------ */

static double const mv_xk[8] = { 0.991455371120812639206854697526329, 0.949107912342758524526189684047851,
	0.864864423359769072789712788640926, 0.741531185599394439863864773280788, 0.586087235467691130294144845693013,
	0.405845151377397166906606412076961, 0.207784955007898467600689403773245, 0.0 };
static double const mv_wk[8] = { 0.022935322010529224963732008058970, 0.063092092629978553290700663189204,
	0.104790010322250183839876322541518, 0.140653259715525918745189590510238, 0.169004726639267902826583426598550,
	0.190350578064785409913256402421014, 0.204432940075298892414161999234649, 0.209482141084727828012999174891714 };
static double const mv_wg[4] = { 0.129484966168869693270611432679082, 0.279705391489276667901467771423780,
	0.381830050505118944950369775488975, 0.417959183673469387755102040816327 };

/* what is integrated over a piece: n values at u */
typedef struct mv_job
{
	qaws_mean_value const* m;
	mv_piece const* p;
	double x[2];
	unsigned int n;
	qaws_mean_value_boundary_fn fn;
	void* user;
	unsigned int dim;
	double* tmp;                /* dim values */
	double rtol;                /* the tolerance, no finer than the rounding of c - x allows */
} mv_job;

static double mv_param(mv_piece const* p, double u, double const* c)
{
	if (p->arc)
	{
		double ax = p->q[0][0] - p->m[0], ay = p->q[0][1] - p->m[1], bx = c[0] - p->m[0], by = c[1] - p->m[1];
		return p->t0 + p->rho * fabs(atan2(ax * by - ay * bx, ax * bx + ay * by));
	}
	return p->t0 + u * (p->t1 - p->t0);
}

/* the integrand at u: coordinates (deg + 1 values B_j k / W), or the
   interpolant (dim values f k, then k) */
static void mv_integrand(mv_job* jb, double u, double* out)
{
	mv_piece const* p = jb->p;
	double c[2], dc[2], bw[MV_MAXD + 1], dx, dy, r, k;
	unsigned int j;
	mv_eval(p, u, c, dc, bw);
	dx = c[0] - jb->x[0];
	dy = c[1] - jb->x[1];
	r = hypot(dx, dy);
	k = r > 0 ? p->sign * (dx * dc[1] - dy * dc[0]) / (r * r * r) : 0;
	if (!jb->fn)
	{
		for (j = 0; j <= p->deg; j++)
			out[j] = bw[j] * k;
		return;
	}
	{
		qaws_vec2 pt;
		pt.x = (qaws_scalar)c[0]; pt.y = (qaws_scalar)c[1];
		jb->fn(jb->user, p->path, p->curve, (qaws_scalar)mv_param(p, u, c), pt, (qaws_scalar*)jb->tmp);
		for (j = 0; j < jb->dim; j++)
			out[j] = (double)((qaws_scalar*)jb->tmp)[j] * k;
		out[jb->dim] = k;
	}
}

/* adaptive Gauss-Kronrod on [a, b], added to acc */
static void mv_gk(mv_job* jb, double a, double b, double floor, unsigned int depth, double* acc)
{
	double kr[MV_MAXD + 18], gs[MV_MAXD + 18], f1[MV_MAXD + 18], f2[MV_MAXD + 18];
	double hm = 0.5 * (b - a), cm = 0.5 * (a + b), err = 0, tot = 0;
	unsigned int i, j, n = jb->n;
	for (j = 0; j < n; j++) kr[j] = gs[j] = 0;
	mv_integrand(jb, cm, f1);
	for (j = 0; j < n; j++)
	{
		kr[j] += mv_wk[7] * f1[j];
		gs[j] += mv_wg[3] * f1[j];
	}
	for (i = 0; i < 7; i++)
	{
		mv_integrand(jb, cm - hm * mv_xk[i], f1);
		mv_integrand(jb, cm + hm * mv_xk[i], f2);
		for (j = 0; j < n; j++)
		{
			kr[j] += mv_wk[i] * (f1[j] + f2[j]);
			if (i & 1)
				gs[j] += mv_wg[i >> 1] * (f1[j] + f2[j]);
		}
	}
	for (j = 0; j < n; j++)
	{
		kr[j] *= hm;
		gs[j] *= hm;
		err += fabs(kr[j] - gs[j]);
		tot += fabs(kr[j]);
	}
	if (depth < 40 && err > jb->rtol * tot && err > floor)
	{
		mv_gk(jb, a, cm, 0.5 * floor, depth + 1, acc);
		mv_gk(jb, cm, b, 0.5 * floor, depth + 1, acc);
		return;
	}
	for (j = 0; j < n; j++)
		acc[j] += kr[j];
}

/* the piece's parameter nearest x, and the distance there */
static double mv_nearest(mv_piece const* p, double const* x, double* out_u)
{
	unsigned int k, best = 0;
	double bd = 1e300, lo, hi, c[2];
	if (p->line)
	{
		double vx = p->q[1][0] - p->q[0][0], vy = p->q[1][1] - p->q[0][1], l2 = vx * vx + vy * vy, u;
		u = l2 > 0 ? ((x[0] - p->q[0][0]) * vx + (x[1] - p->q[0][1]) * vy) / l2 : 0;
		if (u < 0) u = 0;
		if (u > 1) u = 1;
		*out_u = u;
		return hypot(p->q[0][0] + u * vx - x[0], p->q[0][1] + u * vy - x[1]);
	}
	for (k = 0; k <= 32; k++)
	{
		double d;
		mv_eval(p, k / 32.0, c, NULL, NULL);
		d = hypot(c[0] - x[0], c[1] - x[1]);
		if (d < bd) { bd = d; best = k; }
	}
	lo = best > 0 ? (best - 1) / 32.0 : 0;
	hi = best < 32 ? (best + 1) / 32.0 : 1;
	/* golden section on the distance */
	for (k = 0; k < 60 && hi - lo > 1e-15; k++)
	{
		double m1 = hi - 0.6180339887498949 * (hi - lo), m2 = lo + 0.6180339887498949 * (hi - lo), d1, d2;
		mv_eval(p, m1, c, NULL, NULL);
		d1 = hypot(c[0] - x[0], c[1] - x[1]);
		mv_eval(p, m2, c, NULL, NULL);
		d2 = hypot(c[0] - x[0], c[1] - x[1]);
		if (d1 < d2) hi = m2; else lo = m1;
	}
	*out_u = 0.5 * (lo + hi);
	mv_eval(p, *out_u, c, NULL, NULL);
	return hypot(c[0] - x[0], c[1] - x[1]);
}

/* the piece's integrals into acc (n values) */
static void mv_piece_integrate(mv_job* jb, double floor, double* acc)
{
	mv_piece const* p = jb->p;
	double u0, d, speed, dc[2], c[2], h;
	if (p->line && !jb->fn)
	{
		double d0x = p->q[0][0] - jb->x[0], d0y = p->q[0][1] - jb->x[1], d1x = p->q[1][0] - jb->x[0], d1y = p->q[1][1] - jb->x[1];
		double r0 = hypot(d0x, d0y), r1 = hypot(d1x, d1y), cr = d0x * d1y - d0y * d1x, dt = d0x * d1x + d0y * d1y, den = r0 * r1 + dt, t;
		if (r0 > 0 && r1 > 0 && den > 0)
		{
			t = p->sign * cr / den;
			acc[0] += t / (r0 * p->w[0]);
			acc[1] += t / (r1 * p->w[1]);
			return;
		}
	}
	d = mv_nearest(p, jb->x, &u0);
	/* at a distance d from the boundary c - x is good to about 1e-16 of the
	   extent over d: so is the kernel */
	jb->rtol = jb->m->tol;
	if (d > 0 && 4e-16 * jb->m->ext / d > jb->rtol)
		jb->rtol = 4e-16 * jb->m->ext / d;
	mv_eval(p, u0, c, dc, NULL);
	speed = hypot(dc[0], dc[1]);
	/* intervals graded away from u0 by doubling from the kernel's width */
	h = speed > 0 ? d / speed : 1;
	if (!(h > 1e-15)) h = 1e-15;
	{
		double a = u0, step = h;
		while (a < 1)
		{
			double b = a + step > 1 || step >= 0.25 ? 1 : a + step;
			mv_gk(jb, a, b, floor, 0, acc);
			a = b;
			step *= 2;
		}
		a = u0; step = h;
		while (a > 0)
		{
			double b = a - step < 0 || step >= 0.25 ? 0 : a - step;
			mv_gk(jb, b, a, floor, 0, acc);
			a = b;
			step *= 2;
		}
	}
}

/* the boundary point nearest x when x is on the boundary: its piece and u */
static int mv_on_boundary(qaws_mean_value const* m, double const* x, unsigned int* out_piece, double* out_u)
{
	unsigned int i;
	for (i = 0; i < m->npc; i++)
	{
		double u;
		mv_piece const* p = &m->pc[i];
		/* a quick reject by the control points' hull box */
		double lo[2] = { 1e300, 1e300 }, hi[2] = { -1e300, -1e300 };
		unsigned int j;
		for (j = 0; j <= p->deg; j++)
		{
			if (p->q[j][0] < lo[0]) lo[0] = p->q[j][0];
			if (p->q[j][0] > hi[0]) hi[0] = p->q[j][0];
			if (p->q[j][1] < lo[1]) lo[1] = p->q[j][1];
			if (p->q[j][1] > hi[1]) hi[1] = p->q[j][1];
		}
		if (x[0] < lo[0] - m->eps_on || x[0] > hi[0] + m->eps_on || x[1] < lo[1] - m->eps_on || x[1] > hi[1] + m->eps_on)
			continue;
		if (mv_nearest(p, x, &u) <= m->eps_on)
		{
			*out_piece = i;
			*out_u = u;
			return 1;
		}
	}
	return 0;
}

static qaws_status mv_coordinates(qaws_mean_value const* m, double const* x, double* w)
{
	unsigned int i, j, piece;
	double acc[MV_MAXD + 1], total = 0, u, floor = m->tol * 2 * MV_PI / (m->ext * (m->npc + 1));
	for (i = 0; i < m->nh; i++) w[i] = 0;
	if (mv_on_boundary(m, x, &piece, &u))
	{
		/* the boundary's own basis there */
		mv_piece const* p = &m->pc[piece];
		double c[2], bw[MV_MAXD + 1];
		mv_eval(p, u, c, NULL, bw);
		for (j = 0; j <= p->deg; j++)
			for (i = 0; i < p->nh; i++)
				w[p->h[i]] += p->a[j * p->nh + i] * bw[j];
		return QAWS_STATUS_OK;
	}
	for (i = 0; i < m->npc; i++)
	{
		mv_piece const* p = &m->pc[i];
		mv_job jb;
		memset(&jb, 0, sizeof(jb));
		jb.m = m; jb.p = p; jb.x[0] = x[0]; jb.x[1] = x[1]; jb.n = p->deg + 1;
		for (j = 0; j <= p->deg; j++) acc[j] = 0;
		mv_piece_integrate(&jb, floor, acc);
		for (j = 0; j <= p->deg; j++)
		{
			unsigned int k;
			for (k = 0; k < p->nh; k++)
			{
				double v = p->a[j * p->nh + k] * acc[j];
				w[p->h[k]] += v;
				total += v;
			}
		}
	}
	if (!(fabs(total) > 1e-13 * 2 * MV_PI / m->ext))
		return QAWS_STATUS_NUMERICAL_FAILURE;
	for (i = 0; i < m->nh; i++)
		w[i] /= total;
	return QAWS_STATUS_OK;
}

/* ------------------------------------------------------------------------ */
/*  The API                                                                 */
/* ------------------------------------------------------------------------ */

qaws_status qaws_mean_value_create(qaws_mean_value_desc const* desc, qaws_mean_value** out_coords)
{
	qaws_mean_value* m;
	unsigned int pi, ci, *loop_piece = NULL, nloops = 0;
	qaws_status st = QAWS_STATUS_OK;
	double lo[2] = { 1e300, 1e300 }, hi[2] = { -1e300, -1e300 };
	if (!desc || !out_coords || !desc->paths || desc->path_count == 0)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_coords = NULL;
	m = (qaws_mean_value*)calloc(1, sizeof(*m));
	loop_piece = (unsigned int*)malloc(sizeof(unsigned int) * (desc->path_count + 1));
	if (!m || !loop_piece)
	{
		free(m); free(loop_piece);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
#if QAWS_SCALAR_IS_FLOAT
	m->tol = desc->tolerance > 0 ? (double)desc->tolerance : 1e-5;
#else
	m->tol = desc->tolerance > 0 ? (double)desc->tolerance : 1e-10;
#endif
	m->executor = desc->executor;
	/* the extent first: handles merge relative to it */
	for (pi = 0; pi < desc->path_count; pi++)
		for (ci = 0; ci < desc->paths[pi].curve_count; ci++)
		{
			qaws_vec2 a, b;
			if (qaws_curve_compute_bounds_2d(desc->paths[pi].curves[ci], &a, &b) == QAWS_STATUS_OK)
			{
				if (a.x < lo[0]) lo[0] = a.x;
				if (a.y < lo[1]) lo[1] = a.y;
				if (b.x > hi[0]) hi[0] = b.x;
				if (b.y > hi[1]) hi[1] = b.y;
			}
		}
	m->ext = hi[0] >= lo[0] ? hypot(hi[0] - lo[0], hi[1] - lo[1]) : 1;
	if (!(m->ext > 0)) m->ext = 1;
#if QAWS_SCALAR_IS_FLOAT
	m->eps_on = 1e-7 * m->ext;
#else
	m->eps_on = 1e-13 * m->ext;
#endif
	for (pi = 0; pi < desc->path_count && st == QAWS_STATUS_OK; pi++)
	{
		qaws_path_2d const* path = &desc->paths[pi];
		unsigned int first_piece = m->npc, first_handle = m->nh;
		if (!path->curves || path->curve_count == 0)
			continue;
		for (ci = 0; ci < path->curve_count && st == QAWS_STATUS_OK; ci++)
			st = mv_add_curve(m, path->curves[ci], pi, ci, first_handle);
		if (st != QAWS_STATUS_OK)
			break;
		if (m->npc == first_piece)
			continue;
		/* the loop's last handle on its first: one handle */
		if (m->nh > first_handle + 1 &&
			hypot(m->hp[2 * (m->nh - 1)] - m->hp[2 * first_handle], m->hp[2 * (m->nh - 1) + 1] - m->hp[2 * first_handle + 1]) <= MV_MERGE * m->ext)
		{
			unsigned int last = m->nh - 1, k, j;
			for (k = first_piece; k < m->npc; k++)
				for (j = 0; j < m->pc[k].nh; j++)
					if (m->pc[k].h[j] == last)
						m->pc[k].h[j] = first_handle;
			m->nh--;
		}
		loop_piece[nloops] = first_piece;
		nloops++;
	}
	loop_piece[nloops] = m->npc;
	if (st == QAWS_STATUS_OK && m->npc == 0)
		st = QAWS_STATUS_INVALID_ARGUMENT;
	if (st == QAWS_STATUS_OK)
	{
		/* outside loops counter-clockwise, holes clockwise: nesting depth by
		   the winding of the others round each loop's first point */
		unsigned int l, o, k;
		for (l = 0; l < nloops; l++)
		{
			double area2, wind, pt[2];
			unsigned int depth = 0;
			mv_eval(&m->pc[loop_piece[l]], 0, pt, NULL, NULL);
			for (o = 0; o < nloops; o++)
			{
				double a2;
				if (o == l)
					continue;
				mv_loop_measure(m, loop_piece[o], loop_piece[o + 1], pt, &a2, &wind);
				if (fabs(wind) > 0.5)
					depth++;
			}
			mv_loop_measure(m, loop_piece[l], loop_piece[l + 1], NULL, &area2, &wind);
			for (k = loop_piece[l]; k < loop_piece[l + 1]; k++)
				m->pc[k].sign = ((depth & 1) ? -1.0 : 1.0) * (area2 >= 0 ? 1.0 : -1.0);
		}
	}
	free(loop_piece);
	if (st != QAWS_STATUS_OK)
	{
		qaws_mean_value_destroy(m);
		return st;
	}
	*out_coords = m;
	return QAWS_STATUS_OK;
}

void qaws_mean_value_destroy(qaws_mean_value* coords)
{
	if (!coords)
		return;
	free(coords->pc);
	free(coords->hp);
	free(coords);
}

unsigned int qaws_mean_value_get_handle_count(qaws_mean_value const* coords)
{
	return coords ? coords->nh : 0;
}

qaws_status qaws_mean_value_get_handle(qaws_mean_value const* coords, unsigned int index, qaws_vec2* out_point)
{
	if (!coords || !out_point || index >= coords->nh)
		return QAWS_STATUS_INVALID_ARGUMENT;
	out_point->x = (qaws_scalar)coords->hp[2 * index];
	out_point->y = (qaws_scalar)coords->hp[2 * index + 1];
	return QAWS_STATUS_OK;
}

qaws_status qaws_mean_value_coordinates(qaws_mean_value const* coords, qaws_vec2 point, qaws_scalar* out_weights)
{
	double x[2], *w;
	unsigned int i;
	qaws_status st;
	if (!coords || !out_weights)
		return QAWS_STATUS_INVALID_ARGUMENT;
	w = (double*)malloc(sizeof(double) * (coords->nh + 1));
	if (!w)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	x[0] = (double)point.x; x[1] = (double)point.y;
	st = mv_coordinates(coords, x, w);
	for (i = 0; i < coords->nh; i++)
		out_weights[i] = (qaws_scalar)w[i];
	free(w);
	return st;
}

typedef struct mv_deform_ctx
{
	qaws_mean_value const* m;
	qaws_vec2 const* moved;
	qaws_vec2 const* points;
	qaws_vec2* out;
	qaws_status st;
} mv_deform_ctx;

static void mv_deform_task(void* ctx_, unsigned int begin, unsigned int end)
{
	mv_deform_ctx* ctx = (mv_deform_ctx*)ctx_;
	double* w = (double*)malloc(sizeof(double) * (ctx->m->nh + 1));
	unsigned int i, k;
	if (!w)
	{
		ctx->st = QAWS_STATUS_ALLOCATION_FAILURE;
		return;
	}
	for (i = begin; i < end; i++)
	{
		double x[2], y[2] = { 0, 0 };
		x[0] = (double)ctx->points[i].x; x[1] = (double)ctx->points[i].y;
		if (mv_coordinates(ctx->m, x, w) != QAWS_STATUS_OK)
		{
			ctx->st = QAWS_STATUS_NUMERICAL_FAILURE;
			ctx->out[i] = ctx->points[i];
			continue;
		}
		for (k = 0; k < ctx->m->nh; k++)
		{
			y[0] += w[k] * (double)ctx->moved[k].x;
			y[1] += w[k] * (double)ctx->moved[k].y;
		}
		ctx->out[i].x = (qaws_scalar)y[0];
		ctx->out[i].y = (qaws_scalar)y[1];
	}
	free(w);
}

qaws_status qaws_mean_value_deform(qaws_mean_value const* coords, qaws_vec2 const* moved_handles,
	qaws_vec2 const* points, unsigned int count, qaws_vec2* out_points)
{
	mv_deform_ctx ctx;
	if (!coords || !moved_handles || (!points && count) || (!out_points && count))
		return QAWS_STATUS_INVALID_ARGUMENT;
	ctx.m = coords; ctx.moved = moved_handles; ctx.points = points; ctx.out = out_points; ctx.st = QAWS_STATUS_OK;
	if (coords->executor && coords->executor->parallel_for && count > 1)
		coords->executor->parallel_for(coords->executor->user, count, mv_deform_task, &ctx);
	else
		mv_deform_task(&ctx, 0, count);
	return ctx.st;
}

qaws_status qaws_mean_value_interpolate(qaws_mean_value const* coords, qaws_vec2 point,
	qaws_mean_value_boundary_fn fn, void* user, unsigned int dim, qaws_scalar* out_values)
{
	double x[2], *acc, total, u, floor;
	qaws_scalar* tmp;
	unsigned int i, j, piece;
	if (!coords || !fn || !out_values || dim == 0 || dim > 16)
		return QAWS_STATUS_INVALID_ARGUMENT;
	x[0] = (double)point.x; x[1] = (double)point.y;
	if (mv_on_boundary(coords, x, &piece, &u))
	{
		mv_piece const* p = &coords->pc[piece];
		double c[2];
		qaws_vec2 pt;
		mv_eval(p, u, c, NULL, NULL);
		pt.x = (qaws_scalar)c[0]; pt.y = (qaws_scalar)c[1];
		fn(user, p->path, p->curve, (qaws_scalar)mv_param(p, u, c), pt, out_values);
		return QAWS_STATUS_OK;
	}
	acc = (double*)calloc(dim + 1, sizeof(double));
	tmp = (qaws_scalar*)malloc(sizeof(qaws_scalar) * dim + sizeof(double) * dim);
	if (!acc || !tmp)
	{
		free(acc); free(tmp);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	floor = coords->tol * 2 * MV_PI / (coords->ext * (coords->npc + 1));
	for (i = 0; i < coords->npc; i++)
	{
		mv_job jb;
		memset(&jb, 0, sizeof(jb));
		jb.m = coords; jb.p = &coords->pc[i]; jb.x[0] = x[0]; jb.x[1] = x[1];
		jb.n = dim + 1; jb.fn = fn; jb.user = user; jb.dim = dim; jb.tmp = (double*)tmp;
		mv_piece_integrate(&jb, floor, acc);
	}
	total = acc[dim];
	if (!(fabs(total) > 1e-13 * 2 * MV_PI / coords->ext))
	{
		free(acc); free(tmp);
		return QAWS_STATUS_NUMERICAL_FAILURE;
	}
	for (j = 0; j < dim; j++)
		out_values[j] = (qaws_scalar)(acc[j] / total);
	free(acc); free(tmp);
	return QAWS_STATUS_OK;
}
