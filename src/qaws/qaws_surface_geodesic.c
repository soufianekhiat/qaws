/*
 * Geodesics and curves on surfaces (qaws_surface_region.h).
 *
 * The geodesic equation, by arc length s, in the state y = (u, v, u', v'):
 *     u'' = -(G1_11 u'^2 + 2 G1_12 u' v' + G1_22 v'^2),
 *     v'' = -(G2_11 u'^2 + 2 G2_12 u' v' + G2_22 v'^2),
 * with G^k_ij from g^-1 (S_u . S_ij, S_v . S_ij). Integrated by the
 * Dormand-Prince 5(4) pair with error control. A path is kept as a cubic
 * Hermite spline in (u, v) through equal arc-length steps (positions and
 * velocities, exact at the nodes), the steps halved until the spline's
 * midpoints are within the tolerance of the solution.
 */
#include "qaws_surface_region.h"
#include "qaws_surface.h"
#include "qaws_curve.h"
#include "qaws_hermite.h"
#include "qaws_eval.h"
#include "qaws_inspect.h"
#include "qaws_operations.h"
#include "internal/qaws_internal_surface_region.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define SG_PI 3.14159265358979323846

/* ------------------------------------------------------------------------ */
/*  The surface's metric                                                    */
/* ------------------------------------------------------------------------ */

typedef qaws_internal_sg_ctx sg_ctx;

struct qaws_surface_curve
{
	qaws_surface const* s;
	qaws_surface_curve_kind kind;
	qaws_curve* uv;
	double length;
	double a[2], b[2];          /* geodesic: its ends */
	double centre[2], radius, angle0, angle1;   /* arc */
	sg_ctx ctx;                 /* the surface's domain and closed directions */
};

static double sg_wrap(double x, double x0, double p)
{
	double t;
	if (p <= 0)
		return x;
	t = fmod(x - x0, p);
	if (t < 0) t += p;
	return x0 + t;
}

static qaws_status sg_eval(sg_ctx const* c, double u, double v, unsigned int flags, qaws_surface_eval_result* r)
{
	return qaws_surface_evaluate(r, c->s, (qaws_scalar)sg_wrap(u, c->u0, c->pu), (qaws_scalar)sg_wrap(v, c->v0, c->pv), flags);
}

static double sg_dot(qaws_vec3 a, qaws_vec3 b)
{
	return (double)a.x * b.x + (double)a.y * b.y + (double)a.z * b.z;
}

/* g = {E, F, G}, gam = {G1_11, G1_12, G1_22, G2_11, G2_12, G2_22} */
static qaws_status sg_metric(sg_ctx const* c, double u, double v, double* g, double* gam)
{
	qaws_surface_eval_result r;
	double E, F, G, det;
	qaws_status st;
	if (c->pu <= 0 && (u < c->u0 - 1e-9 * (c->u1 - c->u0) || u > c->u1 + 1e-9 * (c->u1 - c->u0)))
		return QAWS_STATUS_OUT_OF_RANGE;
	if (c->pv <= 0 && (v < c->v0 - 1e-9 * (c->v1 - c->v0) || v > c->v1 + 1e-9 * (c->v1 - c->v0)))
		return QAWS_STATUS_OUT_OF_RANGE;
	st = sg_eval(c, u, v, QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_DV | (gam ? QAWS_SURFACE_EVAL_DUU | QAWS_SURFACE_EVAL_DUV | QAWS_SURFACE_EVAL_DVV : 0), &r);
	if (st != QAWS_STATUS_OK)
		return st;
	E = sg_dot(r.du, r.du); F = sg_dot(r.du, r.dv); G = sg_dot(r.dv, r.dv);
	det = E * G - F * F;
	if (!(det > 1e-24 * (E + G) * (E + G)) || !(E > 0) || !(G > 0))
		return QAWS_STATUS_NUMERICAL_FAILURE;
	if (g) { g[0] = E; g[1] = F; g[2] = G; }
	if (gam)
	{
		qaws_vec3 const* sij[3];
		unsigned int k;
		sij[0] = &r.duu; sij[1] = &r.duv; sij[2] = &r.dvv;
		for (k = 0; k < 3; k++)
		{
			double a = sg_dot(r.du, *sij[k]), b = sg_dot(r.dv, *sij[k]);
			gam[k] = (G * a - F * b) / det;
			gam[3 + k] = (E * b - F * a) / det;
		}
	}
	return QAWS_STATUS_OK;
}

static qaws_status sg_init(sg_ctx* c, qaws_surface const* s, double tol)
{
	qaws_range ru = qaws_surface_get_u_range(s), rv = qaws_surface_get_v_range(s);
	double size = 0, lo[3] = { 1e300, 1e300, 1e300 }, hi[3] = { -1e300, -1e300, -1e300 };
	unsigned int i, j, k, cu = 1, cv = 1;
	memset(c, 0, sizeof(*c));
	c->s = s;
	c->u0 = (double)ru.min_value; c->u1 = (double)ru.max_value;
	c->v0 = (double)rv.min_value; c->v1 = (double)rv.max_value;
	if (!(c->u1 > c->u0) || !(c->v1 > c->v0))
		return QAWS_STATUS_INVALID_ARGUMENT;
	c->diag = hypot(c->u1 - c->u0, c->v1 - c->v0);
#if QAWS_SCALAR_IS_FLOAT
	c->tol = (tol > 0 ? tol : 1e-5) * c->diag;
#else
	c->tol = (tol > 0 ? tol : 1e-10) * c->diag;
#endif
	/* the size, then closed directions: opposite edges on one another */
	for (i = 0; i <= 6; i++)
		for (j = 0; j <= 6; j++)
		{
			qaws_surface_eval_result r;
			if (qaws_surface_evaluate(&r, s, (qaws_scalar)(c->u0 + (c->u1 - c->u0) * i / 6),
				(qaws_scalar)(c->v0 + (c->v1 - c->v0) * j / 6), QAWS_SURFACE_EVAL_POSITION) != QAWS_STATUS_OK)
				continue;
			if (r.position.x < lo[0]) lo[0] = r.position.x;
			if (r.position.y < lo[1]) lo[1] = r.position.y;
			if (r.position.z < lo[2]) lo[2] = r.position.z;
			if (r.position.x > hi[0]) hi[0] = r.position.x;
			if (r.position.y > hi[1]) hi[1] = r.position.y;
			if (r.position.z > hi[2]) hi[2] = r.position.z;
		}
	for (k = 0; k < 3; k++) size += (hi[k] - lo[k]) * (hi[k] - lo[k]);
	size = sqrt(size);
	for (i = 1; i < 8; i++)
	{
		qaws_surface_eval_result a, b;
		double t = (double)i / 8;
#if QAWS_SCALAR_IS_FLOAT
		double close = 1e-5 * size;
#else
		double close = 1e-9 * size;
#endif
		qaws_surface_evaluate(&a, s, ru.min_value, (qaws_scalar)(c->v0 + (c->v1 - c->v0) * t), QAWS_SURFACE_EVAL_POSITION);
		qaws_surface_evaluate(&b, s, ru.max_value, (qaws_scalar)(c->v0 + (c->v1 - c->v0) * t), QAWS_SURFACE_EVAL_POSITION);
		if (sqrt((a.position.x - b.position.x) * (a.position.x - b.position.x) + (a.position.y - b.position.y) * (a.position.y - b.position.y) +
			(a.position.z - b.position.z) * (a.position.z - b.position.z)) > close)
			cu = 0;
		qaws_surface_evaluate(&a, s, (qaws_scalar)(c->u0 + (c->u1 - c->u0) * t), rv.min_value, QAWS_SURFACE_EVAL_POSITION);
		qaws_surface_evaluate(&b, s, (qaws_scalar)(c->u0 + (c->u1 - c->u0) * t), rv.max_value, QAWS_SURFACE_EVAL_POSITION);
		if (sqrt((a.position.x - b.position.x) * (a.position.x - b.position.x) + (a.position.y - b.position.y) * (a.position.y - b.position.y) +
			(a.position.z - b.position.z) * (a.position.z - b.position.z)) > close)
			cv = 0;
	}
	if (cu) c->pu = c->u1 - c->u0;
	if (cv) c->pv = c->v1 - c->v0;
	return QAWS_STATUS_OK;
}

/* ------------------------------------------------------------------------ */
/*  The geodesic equation                                                   */
/* ------------------------------------------------------------------------ */

static qaws_status sg_rhs(sg_ctx const* c, double const* y, double* dy)
{
	double gam[6];
	qaws_status st = sg_metric(c, y[0], y[1], NULL, gam);
	if (st != QAWS_STATUS_OK)
		return st;
	dy[0] = y[2];
	dy[1] = y[3];
	dy[2] = -(gam[0] * y[2] * y[2] + 2 * gam[1] * y[2] * y[3] + gam[2] * y[3] * y[3]);
	dy[3] = -(gam[3] * y[2] * y[2] + 2 * gam[4] * y[2] * y[3] + gam[5] * y[3] * y[3]);
	return QAWS_STATUS_OK;
}

/* y advanced by the arc length len (may be negative), Dormand-Prince 5(4) */
static qaws_status sg_integrate(sg_ctx const* c, double* y, double len)
{
	static double const a21 = 1.0 / 5;
	static double const a31 = 3.0 / 40, a32 = 9.0 / 40;
	static double const a41 = 44.0 / 45, a42 = -56.0 / 15, a43 = 32.0 / 9;
	static double const a51 = 19372.0 / 6561, a52 = -25360.0 / 2187, a53 = 64448.0 / 6561, a54 = -212.0 / 729;
	static double const a61 = 9017.0 / 3168, a62 = -355.0 / 33, a63 = 46732.0 / 5247, a64 = 49.0 / 176, a65 = -5103.0 / 18656;
	static double const b1 = 35.0 / 384, b3 = 500.0 / 1113, b4 = 125.0 / 192, b5 = -2187.0 / 6784, b6 = 11.0 / 84;
	static double const e1 = 71.0 / 57600, e3 = -71.0 / 16695, e4 = 71.0 / 1920, e5 = -17253.0 / 339200, e6 = 22.0 / 525, e7 = -1.0 / 40;
	double done = 0, total = fabs(len), sgn = len < 0 ? -1 : 1, h, tol = 1e-3 * c->tol;
	unsigned int steps = 0;
	double scale;
	if (total == 0)
		return QAWS_STATUS_OK;
	/* a first step: a few hundredths of the length, scaled by the speed in (u, v) */
	scale = fabs(y[2]) + fabs(y[3]);
	h = total / 16;
	if (scale > 0 && 0.05 * c->diag / scale < h)
		h = 0.05 * c->diag / scale;
	while (done < total)
	{
		double k1[4], k2[4], k3[4], k4[4], k5[4], k6[4], k7[4], t[4], y5[4], err = 0, hs;
		unsigned int i;
		qaws_status st;
		if (++steps > 2000000)
			return QAWS_STATUS_NUMERICAL_FAILURE;
		if (h > total - done) h = total - done;
		hs = sgn * h;
		if ((st = sg_rhs(c, y, k1)) != QAWS_STATUS_OK) return st;
		for (i = 0; i < 4; i++) t[i] = y[i] + hs * a21 * k1[i];
		if ((st = sg_rhs(c, t, k2)) != QAWS_STATUS_OK) return st;
		for (i = 0; i < 4; i++) t[i] = y[i] + hs * (a31 * k1[i] + a32 * k2[i]);
		if ((st = sg_rhs(c, t, k3)) != QAWS_STATUS_OK) return st;
		for (i = 0; i < 4; i++) t[i] = y[i] + hs * (a41 * k1[i] + a42 * k2[i] + a43 * k3[i]);
		if ((st = sg_rhs(c, t, k4)) != QAWS_STATUS_OK) return st;
		for (i = 0; i < 4; i++) t[i] = y[i] + hs * (a51 * k1[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i]);
		if ((st = sg_rhs(c, t, k5)) != QAWS_STATUS_OK) return st;
		for (i = 0; i < 4; i++) t[i] = y[i] + hs * (a61 * k1[i] + a62 * k2[i] + a63 * k3[i] + a64 * k4[i] + a65 * k5[i]);
		if ((st = sg_rhs(c, t, k6)) != QAWS_STATUS_OK) return st;
		for (i = 0; i < 4; i++) y5[i] = y[i] + hs * (b1 * k1[i] + b3 * k3[i] + b4 * k4[i] + b5 * k5[i] + b6 * k6[i]);
		if ((st = sg_rhs(c, y5, k7)) != QAWS_STATUS_OK) return st;
		for (i = 0; i < 4; i++)
		{
			double e = hs * (e1 * k1[i] + e3 * k3[i] + e4 * k4[i] + e5 * k5[i] + e6 * k6[i] + e7 * k7[i]);
			/* positions in (u, v); velocities weighted by the step (their effect on the position) */
			if (i >= 2) e *= h;
			if (fabs(e) > err) err = fabs(e);
		}
		if (err <= tol || h <= 1e-14 * total)
		{
			for (i = 0; i < 4; i++) y[i] = y5[i];
			done += h;
		}
		{
			double f = err > 0 ? 0.9 * pow(tol / err, 0.2) : 5;
			if (f < 0.2) f = 0.2;
			if (f > 5) f = 5;
			h *= f;
		}
	}
	return QAWS_STATUS_OK;
}

/* the unit (u, v) velocity along the direction d at (u, v) */
static qaws_status sg_unit(sg_ctx const* c, double u, double v, double const* d, double* out)
{
	double g[3], n;
	qaws_status st = sg_metric(c, u, v, g, NULL);
	if (st != QAWS_STATUS_OK)
		return st;
	n = sqrt(g[0] * d[0] * d[0] + 2 * g[1] * d[0] * d[1] + g[2] * d[1] * d[1]);
	if (!(n > 0))
		return QAWS_STATUS_INVALID_ARGUMENT;
	out[0] = d[0] / n;
	out[1] = d[1] / n;
	return QAWS_STATUS_OK;
}

/* the unit tangent at angle a from S_u, counter-clockwise about S_u x S_v */
static qaws_status sg_dir(sg_ctx const* c, double u, double v, double a, double* out)
{
	double g[3], det, e1[2], e2[2];
	qaws_status st = sg_metric(c, u, v, g, NULL);
	if (st != QAWS_STATUS_OK)
		return st;
	det = g[0] * g[2] - g[1] * g[1];
	e1[0] = 1 / sqrt(g[0]); e1[1] = 0;
	e2[0] = -g[1] / sqrt(g[0] * det); e2[1] = g[0] / sqrt(g[0] * det);
	out[0] = cos(a) * e1[0] + sin(a) * e2[0];
	out[1] = cos(a) * e1[1] + sin(a) * e2[1];
	return QAWS_STATUS_OK;
}

/* the angle of the (u, v) direction d at (u, v), as sg_dir measures it */
static qaws_status sg_angle(sg_ctx const* c, double u, double v, double const* d, double* out)
{
	double g[3], det, al, be;
	qaws_status st = sg_metric(c, u, v, g, NULL);
	if (st != QAWS_STATUS_OK)
		return st;
	det = g[0] * g[2] - g[1] * g[1];
	/* d = al e1 + be e2: be from the v part, al from the u part */
	be = d[1] * sqrt(g[0] * det) / g[0];
	al = (d[0] + be * g[1] / sqrt(g[0] * det)) * sqrt(g[0]);
	*out = atan2(be, al);
	return QAWS_STATUS_OK;
}

/* ------------------------------------------------------------------------ */
/*  Paths                                                                   */
/* ------------------------------------------------------------------------ */

static qaws_status sg_hermite(double const* pts, double const* der, unsigned int n, qaws_curve** out)
{
	qaws_hermite_desc d;
	qaws_scalar* buf = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 4 * n);
	unsigned int i;
	qaws_status st;
	if (!buf)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; i < 2 * n; i++)
	{
		buf[i] = (qaws_scalar)pts[i];
		buf[2 * n + i] = (qaws_scalar)der[i];
	}
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_2D; d.degree = 3;
	d.points = buf; d.point_count = n;
	d.derivatives = buf + 2 * n; d.derivative_count = n;
	st = qaws_curve_create_hermite(out, &d);
	free(buf);
	return st;
}

/* a fit refined by doubling has stopped improving: its error no longer
   falls (a cubic's would by 16) but sits on the integration's own accuracy,
   which more points cannot beat */
static int sg_stalled(double worst, double prev, unsigned int n)
{
	return n >= 64 && worst > 0.5 * prev;
}

/* the geodesic from the state y0 over the length L as a Hermite spline in
   (u, v) through equal arc-length steps */
static qaws_status sg_path(sg_ctx const* c, double const* y0, double L, qaws_curve** out)
{
	unsigned int n = 8;
	double prev = 1e300;
	*out = NULL;
	while (n <= 1u << 18)
	{
		double* pts = (double*)malloc(sizeof(double) * 4 * (n + 1));
		double *der = pts + 2 * (n + 1), y[4], h = L / n, worst = 0;
		unsigned int k;
		qaws_status st = QAWS_STATUS_OK;
		qaws_curve* cv = NULL;
		if (!pts)
			return QAWS_STATUS_ALLOCATION_FAILURE;
		memcpy(y, y0, sizeof(y));
		for (k = 0; k <= n && st == QAWS_STATUS_OK; k++)
		{
			if (k)
				st = sg_integrate(c, y, h);
			pts[2 * k] = y[0]; pts[2 * k + 1] = y[1];
			der[2 * k] = y[2] * h; der[2 * k + 1] = y[3] * h;
		}
		if (st == QAWS_STATUS_OK)
			st = sg_hermite(pts, der, n + 1, &cv);
		/* the midpoints against the solution */
		for (k = 0; k < n && st == QAWS_STATUS_OK; k++)
		{
			double ym[4];
			qaws_eval_result_2d e;
			ym[0] = pts[2 * k]; ym[1] = pts[2 * k + 1]; ym[2] = der[2 * k] / h; ym[3] = der[2 * k + 1] / h;
			st = sg_integrate(c, ym, 0.5 * h);
			if (st == QAWS_STATUS_OK && qaws_curve_evaluate_2d(&e, cv, (qaws_scalar)(k + 0.5), QAWS_EVAL_FLAG_POSITION) == QAWS_STATUS_OK)
			{
				double d = hypot((double)e.position.x - ym[0], (double)e.position.y - ym[1]);
				if (d > worst) worst = d;
			}
		}
		free(pts);
		if (st != QAWS_STATUS_OK)
		{
			qaws_curve_destroy(cv);
			return st;
		}
		if (worst <= c->tol || n >= 1u << 18 || sg_stalled(worst, prev, n))
		{
			*out = cv;
			return QAWS_STATUS_OK;
		}
		qaws_curve_destroy(cv);
		prev = worst;
		n *= 2;
	}
	return QAWS_STATUS_NUMERICAL_FAILURE;
}

static qaws_surface_curve* sg_new(qaws_surface const* s, qaws_surface_curve_kind kind)
{
	qaws_surface_curve* c = (qaws_surface_curve*)calloc(1, sizeof(qaws_surface_curve));
	if (c)
	{
		c->s = s;
		c->kind = kind;
	}
	return c;
}

/* ------------------------------------------------------------------------ */
/*  Shooting                                                                */
/* ------------------------------------------------------------------------ */

/* the end of the geodesic from a at angle th over length L, and its velocity */
static qaws_status sg_shoot(sg_ctx const* c, double const* a, double th, double L, double* end, double* vel)
{
	double y[4], d[2];
	qaws_status st = sg_dir(c, a[0], a[1], th, d);
	if (st != QAWS_STATUS_OK)
		return st;
	y[0] = a[0]; y[1] = a[1]; y[2] = d[0]; y[3] = d[1];
	st = sg_integrate(c, y, L);
	if (st != QAWS_STATUS_OK)
		return st;
	end[0] = y[0]; end[1] = y[1];
	if (vel) { vel[0] = y[2]; vel[1] = y[3]; }
	return QAWS_STATUS_OK;
}

/* Newton on (angle, length) from (th, L) towards b */
static int sg_newton(sg_ctx const* c, double const* a, double const* b, double* th, double* L)
{
	unsigned int it;
	double end[2], vel[2], r, dth = 1e-6;
	if (sg_shoot(c, a, *th, *L, end, vel) != QAWS_STATUS_OK)
		return 0;
	r = hypot(end[0] - b[0], end[1] - b[1]);
	for (it = 0; it < 40; it++)
	{
		double ep[2], em[2], J[4], det, step[2], fx = end[0] - b[0], fy = end[1] - b[1];
		unsigned int k;
		int ok = 0;
		if (r <= 0.05 * c->tol)
			return 1;
		if (sg_shoot(c, a, *th + dth, *L, ep, NULL) != QAWS_STATUS_OK || sg_shoot(c, a, *th - dth, *L, em, NULL) != QAWS_STATUS_OK)
			return 0;
		J[0] = (ep[0] - em[0]) / (2 * dth); J[1] = vel[0];
		J[2] = (ep[1] - em[1]) / (2 * dth); J[3] = vel[1];
		det = J[0] * J[3] - J[1] * J[2];
		if (!(fabs(det) > 0))
			return 0;
		step[0] = (J[3] * fx - J[1] * fy) / det;
		step[1] = (J[0] * fy - J[2] * fx) / det;
		for (k = 0; k < 10; k++)
		{
			double nth = *th - step[0], nL = *L - step[1], ne[2], nv[2], nr;
			if (nL > 0 && sg_shoot(c, a, nth, nL, ne, nv) == QAWS_STATUS_OK && (nr = hypot(ne[0] - b[0], ne[1] - b[1])) < r)
			{
				*th = nth; *L = nL;
				end[0] = ne[0]; end[1] = ne[1]; vel[0] = nv[0]; vel[1] = nv[1];
				r = nr;
				ok = 1;
				break;
			}
			step[0] *= 0.5; step[1] *= 0.5;
		}
		if (!ok)
			return r <= 0.05 * c->tol;
	}
	return r <= 0.05 * c->tol;
}

/* the metric length of the (u, v) segment a b (Gauss-Legendre, 16 points) */
static double sg_chord_length(sg_ctx const* c, double const* a, double const* b)
{
	static double const x[8] = { 0.0950125098376374, 0.2816035507792589, 0.4580167776572274, 0.6178762444026438,
		0.7554044083550030, 0.8656312023878318, 0.9445750230732326, 0.9894009349916499 };
	static double const w[8] = { 0.1894506104550685, 0.1826034150449236, 0.1691565193950025, 0.1495959888165767,
		0.1246289712555339, 0.0951585116824928, 0.0622535239386479, 0.0271524594117541 };
	double d[2] = { b[0] - a[0], b[1] - a[1] }, s = 0;
	unsigned int i, k;
	for (i = 0; i < 8; i++)
		for (k = 0; k < 2; k++)
		{
			double t = 0.5 + (k ? 0.5 : -0.5) * x[i], g[3];
			if (sg_metric(c, a[0] + t * d[0], a[1] + t * d[1], g, NULL) == QAWS_STATUS_OK)
				s += 0.5 * w[i] * sqrt(g[0] * d[0] * d[0] + 2 * g[1] * d[0] * d[1] + g[2] * d[1] * d[1]);
		}
	return s;
}

/* ------------------------------------------------------------------------ */
/*  The API                                                                 */
/* ------------------------------------------------------------------------ */

qaws_status qaws_surface_get_metric(qaws_scalar* out_g, qaws_scalar* out_christoffel, qaws_surface const* surface,
	qaws_scalar u, qaws_scalar v)
{
	sg_ctx c;
	double g[3], gam[6];
	unsigned int i;
	qaws_status st;
	if (!surface)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if ((st = sg_init(&c, surface, 0)) != QAWS_STATUS_OK)
		return st;
	st = sg_metric(&c, (double)u, (double)v, g, gam);
	if (st != QAWS_STATUS_OK)
		return st;
	if (out_g) for (i = 0; i < 3; i++) out_g[i] = (qaws_scalar)g[i];
	if (out_christoffel) for (i = 0; i < 6; i++) out_christoffel[i] = (qaws_scalar)gam[i];
	return QAWS_STATUS_OK;
}

qaws_status qaws_surface_exp(qaws_vec2* out_point, qaws_vec2* out_velocity, qaws_surface const* surface, qaws_vec2 p,
	qaws_vec2 direction, qaws_scalar length, qaws_scalar tolerance)
{
	sg_ctx c;
	double y[4], d[2], un[2];
	qaws_status st;
	if (!surface || !out_point)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if ((st = sg_init(&c, surface, (double)tolerance)) != QAWS_STATUS_OK)
		return st;
	d[0] = (double)direction.x; d[1] = (double)direction.y;
	if ((st = sg_unit(&c, (double)p.x, (double)p.y, d, un)) != QAWS_STATUS_OK)
		return st;
	y[0] = (double)p.x; y[1] = (double)p.y; y[2] = un[0]; y[3] = un[1];
	if ((st = sg_integrate(&c, y, (double)length)) != QAWS_STATUS_OK)
		return st;
	out_point->x = (qaws_scalar)y[0]; out_point->y = (qaws_scalar)y[1];
	if (out_velocity) { out_velocity->x = (qaws_scalar)y[2]; out_velocity->y = (qaws_scalar)y[3]; }
	return QAWS_STATUS_OK;
}

qaws_status qaws_surface_curve_create_uv(qaws_surface_curve** out_curve, qaws_surface const* surface,
	qaws_curve const* uv_curve)
{
	qaws_surface_curve* c;
	qaws_range rg;
	qaws_status st;
	if (!surface || !uv_curve || !out_curve)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_curve = NULL;
	if (qaws_curve_get_dimension(uv_curve) != QAWS_DIMENSION_2D)
		return QAWS_STATUS_INVALID_DIMENSION;
	c = sg_new(surface, QAWS_SURFACE_CURVE_UV);
	if (!c)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	if ((st = sg_init(&c->ctx, surface, 0)) != QAWS_STATUS_OK)
	{
		free(c);
		return st;
	}
	rg = qaws_curve_get_parameter_range(uv_curve);
	st = qaws_curve_extract(&c->uv, uv_curve, rg.min_value, rg.max_value);
	c->length = -1;
	if (st != QAWS_STATUS_OK)
	{
		qaws_surface_curve_destroy(c);
		return st;
	}
	*out_curve = c;
	return QAWS_STATUS_OK;
}

qaws_status qaws_surface_curve_create_geodesic_ray(qaws_surface_curve** out_curve, qaws_surface const* surface, qaws_vec2 p,
	qaws_vec2 direction, qaws_scalar length, qaws_scalar tolerance)
{
	sg_ctx ctx;
	qaws_surface_curve* c;
	double y[4], d[2], un[2], L = (double)length;
	qaws_status st;
	if (!surface || !out_curve || !(L > 0))
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_curve = NULL;
	if ((st = sg_init(&ctx, surface, (double)tolerance)) != QAWS_STATUS_OK)
		return st;
	d[0] = (double)direction.x; d[1] = (double)direction.y;
	if ((st = sg_unit(&ctx, (double)p.x, (double)p.y, d, un)) != QAWS_STATUS_OK)
		return st;
	c = sg_new(surface, QAWS_SURFACE_CURVE_GEODESIC);
	if (!c)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	c->ctx = ctx;
	y[0] = (double)p.x; y[1] = (double)p.y; y[2] = un[0]; y[3] = un[1];
	st = sg_path(&ctx, y, L, &c->uv);
	if (st != QAWS_STATUS_OK)
	{
		qaws_surface_curve_destroy(c);
		return st;
	}
	c->length = L;
	c->a[0] = y[0]; c->a[1] = y[1];
	{
		qaws_eval_result_2d e;
		qaws_range rg = qaws_curve_get_parameter_range(c->uv);
		qaws_curve_evaluate_2d(&e, c->uv, rg.max_value, QAWS_EVAL_FLAG_POSITION);
		c->b[0] = (double)e.position.x; c->b[1] = (double)e.position.y;
	}
	*out_curve = c;
	return QAWS_STATUS_OK;
}

qaws_status qaws_surface_curve_create_geodesic(qaws_surface_curve** out_curve, qaws_surface const* surface, qaws_vec2 a,
	qaws_vec2 b, qaws_scalar tolerance)
{
	sg_ctx ctx;
	double pa[2], pb[2], d[2], th, L;
	int ok;
	qaws_status st;
	if (!surface || !out_curve)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_curve = NULL;
	if ((st = sg_init(&ctx, surface, (double)tolerance)) != QAWS_STATUS_OK)
		return st;
	pa[0] = (double)a.x; pa[1] = (double)a.y; pb[0] = (double)b.x; pb[1] = (double)b.y;
	d[0] = pb[0] - pa[0]; d[1] = pb[1] - pa[1];
	if (!(hypot(d[0], d[1]) > 0))
		return QAWS_STATUS_INVALID_ARGUMENT;
	if ((st = sg_angle(&ctx, pa[0], pa[1], d, &th)) != QAWS_STATUS_OK)
		return st;
	L = sg_chord_length(&ctx, pa, pb);
	ok = sg_newton(&ctx, pa, pb, &th, &L);
	if (!ok)
	{
		/* continuation along the chord: the targets a + l (b - a) */
		unsigned int m, j;
		for (m = 8; m <= 256 && !ok; m *= 4)
		{
			double t2 = th, L2;
			if (sg_angle(&ctx, pa[0], pa[1], d, &t2) != QAWS_STATUS_OK)
				break;
			ok = 1;
			for (j = 1; j <= m && ok; j++)
			{
				double tg[2], lam = (double)j / m;
				tg[0] = pa[0] + lam * d[0]; tg[1] = pa[1] + lam * d[1];
				L2 = j == 1 ? sg_chord_length(&ctx, pa, tg) : L * lam / ((double)(j - 1) / m);
				ok = sg_newton(&ctx, pa, tg, &t2, &L2);
				th = t2;
				L = L2;
			}
		}
	}
	if (!ok)
		return QAWS_STATUS_NUMERICAL_FAILURE;
	{
		qaws_vec2 dir;
		double dd[2];
		qaws_surface_curve* c = NULL;
		sg_dir(&ctx, pa[0], pa[1], th, dd);
		dir.x = (qaws_scalar)dd[0]; dir.y = (qaws_scalar)dd[1];
		st = qaws_surface_curve_create_geodesic_ray(&c, surface, a, dir, (qaws_scalar)L, tolerance);
		if (st != QAWS_STATUS_OK)
			return st;
		/* the end exactly on b */
		c->b[0] = pb[0]; c->b[1] = pb[1];
		*out_curve = c;
	}
	return QAWS_STATUS_OK;
}

/* polyline self-crossing (non-adjacent segments), closed or not */
static int sg_self_crossing(double const* p, unsigned int n, int closed)
{
	unsigned int i, j, m = closed ? n : n - 1;
	for (i = 0; i < m; i++)
	{
		double const* a = p + 2 * i;
		double const* b = p + 2 * ((i + 1) % n);
		for (j = i + 2; j < m; j++)
		{
			double const* c = p + 2 * j;
			double const* d = p + 2 * ((j + 1) % n);
			double d1, d2, d3, d4;
			if (closed && i == 0 && j == m - 1)
				continue;
			d1 = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
			d2 = (b[0] - a[0]) * (d[1] - a[1]) - (b[1] - a[1]) * (d[0] - a[0]);
			d3 = (d[0] - c[0]) * (a[1] - c[1]) - (d[1] - c[1]) * (a[0] - c[0]);
			d4 = (d[0] - c[0]) * (b[1] - c[1]) - (d[1] - c[1]) * (b[0] - c[0]);
			if (((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0)) && ((d3 > 0 && d4 < 0) || (d3 < 0 && d4 > 0)))
				return 1;
		}
	}
	return 0;
}

/* a circle's centre on a pole (a domain edge where S_u x S_v = 0): 0 when
   not, else +1 at v1, -1 at v0 */
static int sg_pole(sg_ctx const* c, double const* p)
{
	double g[3];
#if QAWS_SCALAR_IS_FLOAT
	double const small = 1e-5;          /* cos(pi / 2) in float is 4e-8, not 0 */
#else
	double const small = 1e-10;
#endif
	if (c->pu <= 0)
		return 0;
	if (sg_metric(c, p[0], p[1], g, NULL) != QAWS_STATUS_NUMERICAL_FAILURE &&
		(qaws_internal_sg_first_form(g, p[0], p[1], c) != QAWS_STATUS_OK || !(g[0] <= small * small * g[2])))
		return 0;
	if (fabs(p[1] - c->v1) <= 1e-9 * (c->v1 - c->v0)) return 1;
	if (fabs(p[1] - c->v0) <= 1e-9 * (c->v1 - c->v0)) return -1;
	return 0;
}

/* the circle's point at angle a, radius r: from the centre, or from a pole
   along the meridian u = u0 + P a / (2 pi) (starting a hair away from it,
   where the metric is regular; meridians of revolutions are geodesics) */
static qaws_status sg_ray(sg_ctx const* c, double const* c0, int pole, double a, double r, double* y)
{
	double d[2];
	qaws_status st;
	if (!pole)
	{
		st = sg_dir(c, c0[0], c0[1], a, d);
		if (st != QAWS_STATUS_OK)
			return st;
		y[0] = c0[0]; y[1] = c0[1]; y[2] = d[0]; y[3] = d[1];
		return sg_integrate(c, y, r);
	}
	{
#if QAWS_SCALAR_IS_FLOAT
		double dv = 1e-3 * (c->v1 - c->v0), g[3], len, away[2];
#else
		double dv = 1e-8 * (c->v1 - c->v0), g[3], len, away[2];
#endif
		/* counter-clockwise seen from outside: at v0 the angle runs against u */
		y[0] = c->u0 + (pole > 0 ? 1 : -1) * c->pu * a / (2 * SG_PI);
		y[1] = pole > 0 ? c->v1 - dv : c->v0 + dv;
		/* the hair's length: |S_v| dv at its middle */
		st = sg_metric(c, y[0], pole > 0 ? c->v1 - 0.5 * dv : c->v0 + 0.5 * dv, g, NULL);
		if (st != QAWS_STATUS_OK)
		{
			double gg[3];
			qaws_surface_eval_result rr;
			st = sg_eval(c, y[0], y[1], QAWS_SURFACE_EVAL_DV, &rr);
			if (st != QAWS_STATUS_OK)
				return st;
			gg[2] = sg_dot(rr.dv, rr.dv);
			g[2] = gg[2];
		}
		len = sqrt(g[2]) * dv;
		away[0] = 0; away[1] = pole > 0 ? -1 : 1;
		st = sg_unit(c, y[0], y[1], away, d);
		if (st != QAWS_STATUS_OK)
			return st;
		y[2] = d[0]; y[3] = d[1];
		return sg_integrate(c, y, r - len);
	}
}

qaws_status qaws_surface_curve_create_geodesic_arc(qaws_surface_curve** out_curve, qaws_surface const* surface, qaws_vec2 centre,
	qaws_scalar radius, qaws_scalar angle0, qaws_scalar angle1, qaws_scalar tolerance)
{
	sg_ctx ctx;
	double c0[2], r = (double)radius, a0 = (double)angle0, a1 = (double)angle1, sweep;
	int whole;
	int pole;
	unsigned int n;
	double prev = 1e300;
	qaws_status st;
	if (!surface || !out_curve || !(r > 0) || !(a1 != a0))
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_curve = NULL;
	if ((st = sg_init(&ctx, surface, (double)tolerance)) != QAWS_STATUS_OK)
		return st;
	c0[0] = (double)centre.x; c0[1] = (double)centre.y;
	sweep = a1 - a0;
	pole = sg_pole(&ctx, c0);
	/* a full turn within the float rounding of 2 pi closes the circle */
	whole = fabs(fabs(sweep) - 2 * SG_PI) <= 1e-6 * 2 * SG_PI;
	if (whole)
		sweep = sweep > 0 ? 2 * SG_PI : -2 * SG_PI;
	else if (fabs(sweep) > 2 * SG_PI)
		return QAWS_STATUS_INVALID_ARGUMENT;
	n = (unsigned int)ceil(16 * fabs(sweep) / (2 * SG_PI));
	if (n < 4) n = 4;
	while (n <= 1u << 16)
	{
		/* points at n + 1 angles, two more each side for the derivatives */
		double* pts = (double*)malloc(sizeof(double) * 2 * (n + 5) + sizeof(double) * 4 * (n + 1));
		double *P = pts + 4, *hp, *hd, dt = sweep / n, worst = 0;
		unsigned int k;
		qaws_curve* cv = NULL;
		if (!pts)
			return QAWS_STATUS_ALLOCATION_FAILURE;
		hp = pts + 2 * (n + 5);
		hd = hp + 2 * (n + 1);
		st = QAWS_STATUS_OK;
		for (k = 0; k < n + 5 && st == QAWS_STATUS_OK; k++)
		{
			double y[4], ang = a0 + ((double)k - 2) * dt;
			st = sg_ray(&ctx, c0, pole, ang, r, y);
			pts[2 * k] = y[0]; pts[2 * k + 1] = y[1];
		}
		if (st == QAWS_STATUS_OK && whole)
		{
			/* a closed curve: its two ends one point (round a pole, one period apart in u) */
			P[2 * n] = P[0] + (pole ? (sweep > 0 ? 1 : -1) * pole * ctx.pu : 0); P[2 * n + 1] = P[1];
		}
		for (k = 0; k <= n && st == QAWS_STATUS_OK; k++)
		{
			double const* q = P + 2 * k;
			int j;                  /* (signed: q[-4 + j]) */
			hp[2 * k] = q[0]; hp[2 * k + 1] = q[1];
			/* d/dk by the five-point stencil (fourth order) */
			for (j = 0; j < 2; j++)
				hd[2 * k + j] = (-q[4 + j] + 8 * q[2 + j] - 8 * q[-2 + j] + q[-4 + j]) / 12;
		}
		if (st == QAWS_STATUS_OK && whole)
		{
			hd[2 * n] = hd[0]; hd[2 * n + 1] = hd[1];
		}
		if (st == QAWS_STATUS_OK)
			st = sg_hermite(hp, hd, n + 1, &cv);
		for (k = 0; k < n && st == QAWS_STATUS_OK; k++)
		{
			double y[4];
			qaws_eval_result_2d e;
			st = sg_ray(&ctx, c0, pole, a0 + (k + 0.5) * dt, r, y);
			if (st == QAWS_STATUS_OK && qaws_curve_evaluate_2d(&e, cv, (qaws_scalar)(k + 0.5), QAWS_EVAL_FLAG_POSITION) == QAWS_STATUS_OK)
			{
				double dd = hypot((double)e.position.x - y[0], (double)e.position.y - y[1]);
				if (dd > worst) worst = dd;
			}
		}
		if (st == QAWS_STATUS_OK && (worst <= ctx.tol || n >= 1u << 16 || sg_stalled(worst, prev, n)))
		{
			qaws_surface_curve* c;
			if (sg_self_crossing(P, whole && !pole ? n : n + 1, whole && !pole))
			{
				free(pts);
				qaws_curve_destroy(cv);
				return QAWS_STATUS_INVALID_ARGUMENT;      /* past the injectivity radius */
			}
			c = sg_new(surface, QAWS_SURFACE_CURVE_GEODESIC_ARC);
			if (!c)
			{
				free(pts);
				qaws_curve_destroy(cv);
				return QAWS_STATUS_ALLOCATION_FAILURE;
			}
			c->ctx = ctx;
			c->uv = cv;
			c->centre[0] = c0[0]; c->centre[1] = c0[1];
			c->radius = r; c->angle0 = a0; c->angle1 = a1;
			c->length = -1;
			c->a[0] = P[0]; c->a[1] = P[1];
			c->b[0] = P[2 * n]; c->b[1] = P[2 * n + 1];
			free(pts);
			*out_curve = c;
			return QAWS_STATUS_OK;
		}
		free(pts);
		qaws_curve_destroy(cv);
		if (st != QAWS_STATUS_OK)
			return st;
		prev = worst;
		n *= 2;
	}
	return QAWS_STATUS_NUMERICAL_FAILURE;
}

void qaws_surface_curve_destroy(qaws_surface_curve* curve)
{
	if (!curve)
		return;
	qaws_curve_destroy(curve->uv);
	free(curve);
}

qaws_surface_curve_kind qaws_surface_curve_get_kind(qaws_surface_curve const* curve)
{
	return curve ? curve->kind : QAWS_SURFACE_CURVE_UV;
}

qaws_surface const* qaws_surface_curve_get_surface(qaws_surface_curve const* curve)
{
	return curve ? curve->s : NULL;
}

qaws_curve const* qaws_surface_curve_get_uv(qaws_surface_curve const* curve)
{
	return curve ? curve->uv : NULL;
}

qaws_status qaws_surface_curve_get_length(qaws_scalar* out_length, qaws_surface_curve const* curve)
{
	if (!curve || !out_length)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (curve->kind == QAWS_SURFACE_CURVE_GEODESIC)
	{
		*out_length = (qaws_scalar)curve->length;
		return QAWS_STATUS_OK;
	}
	{
		/* Gauss-Legendre (8 points) on 64 panels per span of the (u, v) image */
		static double const x[4] = { 0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363 };
		static double const w[4] = { 0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763 };
		qaws_range rg = qaws_curve_get_parameter_range(curve->uv);
		unsigned int spans = qaws_curve_get_span_count(curve->uv), np, p, i, k;
		double s = 0, t0 = (double)rg.min_value, t1 = (double)rg.max_value;
		if (spans < 1) spans = 1;
		np = 16 * spans;
		for (p = 0; p < np; p++)
		{
			double a = t0 + (t1 - t0) * p / np, b = t0 + (t1 - t0) * (p + 1) / np;
			for (i = 0; i < 4; i++)
				for (k = 0; k < 2; k++)
				{
					double t = 0.5 * (a + b) + (k ? 0.5 : -0.5) * (b - a) * x[i], g[3];
					qaws_eval_result_2d e;
					qaws_curve_evaluate_2d(&e, curve->uv, (qaws_scalar)t, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1);
					if (sg_metric(&curve->ctx, (double)e.position.x, (double)e.position.y, g, NULL) != QAWS_STATUS_OK)
						return QAWS_STATUS_NUMERICAL_FAILURE;
					s += 0.5 * (b - a) * w[i] * sqrt(g[0] * e.d1.x * e.d1.x + 2 * g[1] * e.d1.x * e.d1.y + g[2] * e.d1.y * e.d1.y);
				}
		}
		*out_length = (qaws_scalar)s;
		return QAWS_STATUS_OK;
	}
}

qaws_status qaws_surface_curve_evaluate(qaws_vec2* out_uv, qaws_vec3* out_point, qaws_surface_curve const* curve,
	qaws_scalar t)
{
	qaws_eval_result_2d e;
	qaws_status st;
	if (!curve)
		return QAWS_STATUS_INVALID_ARGUMENT;
	st = qaws_curve_evaluate_2d(&e, curve->uv, t, QAWS_EVAL_FLAG_POSITION);
	if (st != QAWS_STATUS_OK)
		return st;
	if (out_uv)
		*out_uv = e.position;
	if (out_point)
	{
		qaws_surface_eval_result r;
		if ((st = sg_eval(&curve->ctx, (double)e.position.x, (double)e.position.y, QAWS_SURFACE_EVAL_POSITION, &r)) != QAWS_STATUS_OK)
			return st;
		*out_point = r.position;
	}
	return QAWS_STATUS_OK;
}

/* ------------------------------------------------------------------------ */
/*  Shared with the region code                                             */
/* ------------------------------------------------------------------------ */

qaws_status qaws_internal_sg_init(qaws_surface const* s, double tol, qaws_internal_sg_ctx* c)
{
	return sg_init(c, s, tol);
}

qaws_status qaws_internal_sg_metric(double* g, double* gam, double u, double v, qaws_internal_sg_ctx const* c)
{
	return sg_metric(c, u, v, g, gam);
}

qaws_status qaws_internal_sg_eval(qaws_surface_eval_result* r, double u, double v, unsigned int flags,
	qaws_internal_sg_ctx const* c)
{
	return sg_eval(c, u, v, flags, r);
}

/* E, F, G without the singular check (a pole's edge has zero length) */
qaws_status qaws_internal_sg_first_form(double* g, double u, double v, qaws_internal_sg_ctx const* c)
{
	qaws_surface_eval_result r;
	qaws_status st = sg_eval(c, u, v, QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_DV, &r);
	if (st != QAWS_STATUS_OK)
		return st;
	g[0] = sg_dot(r.du, r.du);
	g[1] = sg_dot(r.du, r.dv);
	g[2] = sg_dot(r.dv, r.dv);
	return QAWS_STATUS_OK;
}

/* sqrt(E G - F^2), zero at a pole */
qaws_status qaws_internal_sg_area_element(double* out, double u, double v, qaws_internal_sg_ctx const* c)
{
	double g[3], det;
	qaws_status st = qaws_internal_sg_first_form(g, u, v, c);
	if (st != QAWS_STATUS_OK)
		return st;
	det = g[0] * g[2] - g[1] * g[1];
	*out = det > 0 ? sqrt(det) : 0;
	return QAWS_STATUS_OK;
}

/* ------------------------------------------------------------------------ */
/*  Offsets: exp along the conormals                                        */
/* ------------------------------------------------------------------------ */

/* the base's point at t and the unit velocity a quarter turn left of its
   tangent (side > 0) or right of it */
static qaws_status sg_conormal(sg_ctx const* c, qaws_curve const* base, double t, double side, double* y)
{
	qaws_eval_result_2d e;
	double d[2], a;
	qaws_status st = qaws_curve_evaluate_2d(&e, base, (qaws_scalar)t, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1);
	if (st != QAWS_STATUS_OK)
		return st;
	d[0] = (double)e.d1.x; d[1] = (double)e.d1.y;
	y[0] = (double)e.position.x; y[1] = (double)e.position.y;
	if ((st = sg_angle(c, y[0], y[1], d, &a)) != QAWS_STATUS_OK)
		return st;
	return sg_dir(c, y[0], y[1], a + (side > 0 ? 0.5 : -0.5) * SG_PI, y + 2);
}

qaws_status qaws_internal_sg_conormal(double* y, qaws_curve const* base, double t, double side, qaws_internal_sg_ctx const* c)
{
	return sg_conormal(c, base, t, side, y);
}

qaws_status qaws_internal_sg_integrate(double* y, double len, qaws_internal_sg_ctx const* c)
{
	return sg_integrate(c, y, len);
}

qaws_status qaws_internal_sg_offset(qaws_curve** out_offset, qaws_curve** out_ray0, qaws_curve** out_ray1,
	qaws_curve const* base, double dist, qaws_internal_sg_ctx const* c)
{
	qaws_range rg = qaws_curve_get_parameter_range(base);
	double t0 = (double)rg.min_value, t1 = (double)rg.max_value, r = fabs(dist), y[4], q0[2], q1[2];
	unsigned int n = 4 * qaws_curve_get_span_count(base);
	qaws_curve *r0 = NULL, *r1 = NULL, *cv = NULL;
	qaws_status st;
	*out_offset = NULL;
	if (out_ray0) *out_ray0 = NULL;
	if (out_ray1) *out_ray1 = NULL;
	if (!(r > 0) || !(t1 > t0))
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (n < 16) n = 16;
	/* the conormal geodesics at the ends: the offset starts and ends exactly
	   where they do */
	st = sg_conormal(c, base, t0, dist, y);
	if (st == QAWS_STATUS_OK) st = sg_path(c, y, r, &r0);
	if (st == QAWS_STATUS_OK) st = sg_conormal(c, base, t1, dist, y);
	if (st == QAWS_STATUS_OK) st = sg_path(c, y, r, &r1);
	if (st == QAWS_STATUS_OK)
	{
		qaws_eval_result_2d e;
		qaws_range ra = qaws_curve_get_parameter_range(r0), rb = qaws_curve_get_parameter_range(r1);
		qaws_curve_evaluate_2d(&e, r0, ra.max_value, QAWS_EVAL_FLAG_POSITION);
		q0[0] = (double)e.position.x; q0[1] = (double)e.position.y;
		qaws_curve_evaluate_2d(&e, r1, rb.max_value, QAWS_EVAL_FLAG_POSITION);
		q1[0] = (double)e.position.x; q1[1] = (double)e.position.y;
	}
	while (st == QAWS_STATUS_OK)
	{
		double* pts = (double*)malloc(sizeof(double) * 4 * (n + 1));
		double *der, worst = 0;
		unsigned int k;
		if (!pts)
		{
			st = QAWS_STATUS_ALLOCATION_FAILURE;
			break;
		}
		der = pts + 2 * (n + 1);
		for (k = 0; k <= n && st == QAWS_STATUS_OK; k++)
		{
			if (k == 0) { pts[0] = q0[0]; pts[1] = q0[1]; continue; }
			if (k == n) { pts[2 * n] = q1[0]; pts[2 * n + 1] = q1[1]; continue; }
			st = sg_conormal(c, base, t0 + (t1 - t0) * k / n, dist, y);
			if (st == QAWS_STATUS_OK) st = sg_integrate(c, y, r);
			pts[2 * k] = y[0]; pts[2 * k + 1] = y[1];
		}
		/* d/dk by five-point stencils, one-sided at the ends */
		for (k = 0; k <= n && st == QAWS_STATUS_OK; k++)
		{
			double const* p = pts;
			unsigned int j;
			for (j = 0; j < 2; j++)
			{
				double v;
				if (k >= 2 && k + 2 <= n)
					v = (-p[2 * (k + 2) + j] + 8 * p[2 * (k + 1) + j] - 8 * p[2 * (k - 1) + j] + p[2 * (k - 2) + j]) / 12;
				else if (k == 0)
					v = (-25 * p[j] + 48 * p[2 + j] - 36 * p[4 + j] + 16 * p[6 + j] - 3 * p[8 + j]) / 12;
				else if (k == 1)
					v = (-3 * p[j] - 10 * p[2 + j] + 18 * p[4 + j] - 6 * p[6 + j] + p[8 + j]) / 12;
				else if (k == n)
					v = (25 * p[2 * n + j] - 48 * p[2 * (n - 1) + j] + 36 * p[2 * (n - 2) + j] - 16 * p[2 * (n - 3) + j] + 3 * p[2 * (n - 4) + j]) / 12;
				else
					v = (3 * p[2 * n + j] + 10 * p[2 * (n - 1) + j] - 18 * p[2 * (n - 2) + j] + 6 * p[2 * (n - 3) + j] - p[2 * (n - 4) + j]) / 12;
				der[2 * k + j] = v;
			}
		}
		if (st == QAWS_STATUS_OK)
			st = sg_hermite(pts, der, n + 1, &cv);
		for (k = 0; k < n && st == QAWS_STATUS_OK; k++)
		{
			qaws_eval_result_2d e;
			st = sg_conormal(c, base, t0 + (t1 - t0) * (k + 0.5) / n, dist, y);
			if (st == QAWS_STATUS_OK) st = sg_integrate(c, y, r);
			if (st == QAWS_STATUS_OK && qaws_curve_evaluate_2d(&e, cv, (qaws_scalar)(k + 0.5), QAWS_EVAL_FLAG_POSITION) == QAWS_STATUS_OK)
			{
				double dd = hypot((double)e.position.x - y[0], (double)e.position.y - y[1]);
				if (dd > worst) worst = dd;
			}
		}
		free(pts);
		/* (past a cusp of the offset, where the base curves more than 1/r, the
		   fit stops refining: that run lies inside the swept strip) */
		if (st != QAWS_STATUS_OK || worst <= c->tol || n >= 1u << 12)
			break;
		qaws_curve_destroy(cv);
		cv = NULL;
		n *= 2;
	}
	if (st != QAWS_STATUS_OK)
	{
		qaws_curve_destroy(cv); qaws_curve_destroy(r0); qaws_curve_destroy(r1);
		return st;
	}
	*out_offset = cv;
	if (out_ray0) *out_ray0 = r0; else qaws_curve_destroy(r0);
	if (out_ray1) *out_ray1 = r1; else qaws_curve_destroy(r1);
	return QAWS_STATUS_OK;
}

qaws_status qaws_surface_curve_create_offset(qaws_surface_curve** out_curve, qaws_surface_curve const* base, qaws_scalar distance,
	qaws_scalar tolerance)
{
	qaws_surface_curve* c;
	qaws_status st;
	if (!base || !out_curve || !(distance != 0))
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_curve = NULL;
	c = sg_new(base->s, QAWS_SURFACE_CURVE_GEODESIC_OFFSET);
	if (!c)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	if ((st = sg_init(&c->ctx, base->s, (double)tolerance)) == QAWS_STATUS_OK)
		st = qaws_internal_sg_offset(&c->uv, NULL, NULL, base->uv, (double)distance, &c->ctx);
	if (st != QAWS_STATUS_OK)
	{
		qaws_surface_curve_destroy(c);
		return st;
	}
	c->length = -1;
	*out_curve = c;
	return QAWS_STATUS_OK;
}

/* ------------------------------------------------------------------------ */
/*  Distances to the boundary                                               */
/* ------------------------------------------------------------------------ */

/* x moved by whole periods to the copy nearest e */
static void sg_near_copy(sg_ctx const* c, double const* e, double const* x, double* out)
{
	out[0] = x[0]; out[1] = x[1];
	if (c->pu > 0) out[0] += c->pu * floor((e[0] - x[0]) / c->pu + 0.5);
	if (c->pv > 0) out[1] += c->pv * floor((e[1] - x[1]) / c->pv + 0.5);
}

/* the conormal geodesic from base(s) on the side, run t: its end and velocity */
static qaws_status sg_foot_shot(sg_ctx const* c, qaws_curve const* base, double side, double s, double t, double* end, double* vel)
{
	double y[4];
	qaws_status st = sg_conormal(c, base, s, side, y);
	if (st == QAWS_STATUS_OK)
		st = sg_integrate(c, y, t);
	end[0] = y[0]; end[1] = y[1];
	if (vel) { vel[0] = y[2]; vel[1] = y[3]; }
	return st;
}

int qaws_internal_sg_foot(double* s, double* t, qaws_curve const* base, double const* x, double side,
	double lo, double hi, qaws_internal_sg_ctx const* c)
{
	unsigned int it;
	double end[2], vel[2], tx[2], r, ds = 1e-6 * (hi - lo);
	if (*s < lo) *s = lo;
	if (*s > hi) *s = hi;
	if (sg_foot_shot(c, base, side, *s, *t, end, vel) != QAWS_STATUS_OK)
		return 0;
	sg_near_copy(c, end, x, tx);
	r = hypot(end[0] - tx[0], end[1] - tx[1]);
	for (it = 0; it < 40; it++)
	{
		double ep[2], em[2], J[4], det, step[2], fx = end[0] - tx[0], fy = end[1] - tx[1], sp, sm;
		unsigned int k;
		int ok = 0;
		if (r <= 0.05 * c->tol)
			return 1;
		sp = *s + ds; sm = *s - ds;
		if (sp > hi) sp = hi;
		if (sm < lo) sm = lo;
		if (sg_foot_shot(c, base, side, sp, *t, ep, NULL) != QAWS_STATUS_OK || sg_foot_shot(c, base, side, sm, *t, em, NULL) != QAWS_STATUS_OK)
			return 0;
		J[0] = (ep[0] - em[0]) / (sp - sm); J[1] = vel[0];
		J[2] = (ep[1] - em[1]) / (sp - sm); J[3] = vel[1];
		det = J[0] * J[3] - J[1] * J[2];
		if (!(fabs(det) > 0))
			return 0;
		step[0] = (J[3] * fx - J[1] * fy) / det;
		step[1] = (J[0] * fy - J[2] * fx) / det;
		for (k = 0; k < 12; k++)
		{
			double ns = *s - step[0], nt = *t - step[1], ne[2], nv[2], nr, nx[2];
			if (ns < lo) ns = lo;
			if (ns > hi) ns = hi;
			if (nt > 0 && sg_foot_shot(c, base, side, ns, nt, ne, nv) == QAWS_STATUS_OK)
			{
				sg_near_copy(c, ne, x, nx);
				nr = hypot(ne[0] - nx[0], ne[1] - nx[1]);
				if (nr < r)
				{
					*s = ns; *t = nt;
					end[0] = ne[0]; end[1] = ne[1]; vel[0] = nv[0]; vel[1] = nv[1];
					tx[0] = nx[0]; tx[1] = nx[1];
					r = nr;
					ok = 1;
					break;
				}
			}
			step[0] *= 0.5; step[1] *= 0.5;
		}
		if (!ok)
			return r <= 0.05 * c->tol;
	}
	return r <= 0.05 * c->tol;
}

int qaws_internal_sg_reach(double* th, double* L, double const* a, double const* x, qaws_internal_sg_ctx const* c)
{
	double end[2], tx[2];
	if (sg_shoot(c, a, *th, *L, end, NULL) != QAWS_STATUS_OK)
		return 0;
	sg_near_copy(c, end, x, tx);
	return sg_newton(c, a, tx, th, L);
}

qaws_status qaws_internal_sg_angle(double* out, double const* d, double u, double v, qaws_internal_sg_ctx const* c)
{
	return sg_angle(c, u, v, d, out);
}

qaws_status qaws_internal_sg_shoot(double* end, double* vel, double const* a, double th, double L, qaws_internal_sg_ctx const* c)
{
	return sg_shoot(c, a, th, L, end, vel);
}
