#include "qaws_diff_functionals.h"
#include "qaws_diff.h"
#include "internal/qaws_internal_types.h"
#include "internal/qaws_internal_curve.h"
#include "internal/qaws_internal_surface.h"
#include "internal/qaws_internal_diff.h"
#include "core/qaws_dual_core.h"
#include "internal/qaws_internal_basis.h"
#include <string.h>

/* ================================================================== */
/*  Gauss-Legendre rules on [-1, 1]                                   */
/* ================================================================== */

static void gauss_rule(unsigned int n, double* x, double* w)
{
	static double const x2[2] = { -0.5773502691896258, 0.5773502691896258 };
	static double const w2[2] = { 1.0, 1.0 };
	static double const x3[3] = { -0.7745966692414834, 0.0, 0.7745966692414834 };
	static double const w3[3] = { 0.5555555555555556, 0.8888888888888889, 0.5555555555555556 };
	static double const x4[4] = { -0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526 };
	static double const w4[4] = { 0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538 };
	static double const x5[5] = { -0.9061798459386640, -0.5384693101056831, 0.0, 0.5384693101056831, 0.9061798459386640 };
	static double const w5[5] = { 0.2369268850561891, 0.4786286704993665, 0.5688888888888889, 0.4786286704993665, 0.2369268850561891 };
	static double const x6[6] = { -0.9324695142031521, -0.6612093864662645, -0.2386191860831969,
		0.2386191860831969, 0.6612093864662645, 0.9324695142031521 };
	static double const w6[6] = { 0.1713244923791704, 0.3607615730481386, 0.4679139345726910,
		0.4679139345726910, 0.3607615730481386, 0.1713244923791704 };
	static double const x7[7] = { -0.9491079123427585, -0.7415311855993945, -0.4058451513773972, 0.0,
		0.4058451513773972, 0.7415311855993945, 0.9491079123427585 };
	static double const w7[7] = { 0.1294849661688697, 0.2797053914892766, 0.3818300505051189, 0.4179591836734694,
		0.3818300505051189, 0.2797053914892766, 0.1294849661688697 };
	static double const x8[8] = { -0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498,
		0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363 };
	static double const w8[8] = { 0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620,
		0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763 };
	double const* xs[9] = { NULL, NULL, x2, x3, x4, x5, x6, x7, x8 };
	double const* ws[9] = { NULL, NULL, w2, w3, w4, w5, w6, w7, w8 };
	memcpy(x, xs[n], sizeof(double) * n);
	memcpy(w, ws[n], sizeof(double) * n);
}

static unsigned int clamp_rule(unsigned int n, unsigned int fallback)
{
	if (n == 0)
		n = fallback;
	if (n < 2)
		n = 2;
	if (n > 8)
		n = 8;
	return n;
}

/* ================================================================== */
/*  Integrands on dual jets                                           */
/* ================================================================== */

typedef qaws_dual1 (*integrand_fn)(qaws_dual3 const* y);

static qaws_dual1 integrand_length(qaws_dual3 const* y)
{
	return qaws_dual3_length(y[1]);
}

static qaws_dual1 integrand_bending(qaws_dual3 const* y)
{
	return qaws_dual3_dot(y[2], y[2]);
}

/* kappa^2 |C'| = |C' x C''|^2 / |C'|^5 */
static qaws_dual1 integrand_curvature_squared(qaws_dual3 const* y)
{
	qaws_dual3 b = qaws_dual3_cross(y[1], y[2]);
	qaws_dual1 s = qaws_dual3_length(y[1]);
	qaws_dual1 s2, s5;
	if (!(s.v > QAWS_EPSILON))
		return qaws_dual1_const(QAWS_ZERO);
	s2 = qaws_dual1_mul(s, s);
	s5 = qaws_dual1_mul(qaws_dual1_mul(s2, s2), s);
	return qaws_dual1_div(qaws_dual3_dot(b, b), s5);
}

static qaws_dual1 integrand_area(qaws_dual3 const* y)
{
	return qaws_dual3_length(qaws_dual3_cross(y[1], y[2]));
}

static qaws_dual1 integrand_thin_plate(qaws_dual3 const* y)
{
	qaws_dual1 a = qaws_dual3_dot(y[3], y[3]);
	qaws_dual1 b = qaws_dual3_dot(y[4], y[4]);
	qaws_dual1 c = qaws_dual3_dot(y[5], y[5]);
	return qaws_dual1_add(qaws_dual1_add(a, qaws_dual1_make(2 * b.v, 2 * b.t, 2 * b.tt)), c);
}

/* H^2 |Su x Sv| */
static qaws_dual1 integrand_willmore(qaws_dual3 const* y)
{
	qaws_dual3 n_raw = qaws_dual3_cross(y[1], y[2]);
	qaws_dual1 nl = qaws_dual3_length(n_raw);
	qaws_dual3 N;
	qaws_dual1 E, F, G, L, M, Nn, det, num, H;
	if (!(nl.v > QAWS_EPSILON))
		return qaws_dual1_const(QAWS_ZERO);
	N = qaws_dual3_div(n_raw, nl);
	E = qaws_dual3_dot(y[1], y[1]);
	F = qaws_dual3_dot(y[1], y[2]);
	G = qaws_dual3_dot(y[2], y[2]);
	L = qaws_dual3_dot(y[3], N);
	M = qaws_dual3_dot(y[4], N);
	Nn = qaws_dual3_dot(y[5], N);
	det = qaws_dual1_sub(qaws_dual1_mul(E, G), qaws_dual1_mul(F, F));
	num = qaws_dual1_sub(qaws_dual1_add(qaws_dual1_mul(E, Nn), qaws_dual1_mul(G, L)),
		qaws_dual1_mul(qaws_dual1_make(2 * F.v, 2 * F.t, 2 * F.tt), M));
	H = qaws_dual1_div(num, qaws_dual1_make(2 * det.v, 2 * det.t, 2 * det.tt));
	return qaws_dual1_mul(qaws_dual1_mul(H, H), nl);
}

static void set_comp(qaws_vec3* v, unsigned int c, qaws_scalar x)
{
	if (c == 0) v->x = x;
	else if (c == 1) v->y = x;
	else v->z = x;
}

static qaws_scalar get_comp(qaws_vec3 const* v, unsigned int c)
{
	return c == 0 ? v->x : (c == 1 ? v->y : v->z);
}

/* q(d) = d^T H d of the integrand at y. */
static qaws_scalar integrand_quadratic(integrand_fn f, qaws_vec3 const* y, qaws_vec3 const* d, unsigned int n)
{
	qaws_dual3 yy[6];
	unsigned int k;
	for (k = 0; k < n; k++)
		yy[k] = qaws_dual3_make(y[k], d[k], qaws_v3_zero());
	return f(yy).tt;
}

/* Gradient of the integrand with respect to the jet: one seed per component. */
static void integrand_gradient(integrand_fn f, qaws_vec3 const* y, unsigned int n, qaws_vec3* g)
{
	qaws_dual3 yy[6];
	unsigned int k, i, c;
	for (i = 0; i < n; i++)
		for (c = 0; c < 3; c++)
		{
			for (k = 0; k < n; k++)
				yy[k] = qaws_dual3_make(y[k], qaws_v3_zero(), qaws_v3_zero());
			set_comp(&yy[i].t, c, QAWS_ONE);
			set_comp(&g[i], c, f(yy).t);
		}
}

/* Hessian of the integrand times ydot, by polarization of q(d). */
static void integrand_hess_vec(integrand_fn f, qaws_vec3 const* y, qaws_vec3 const* ydot, unsigned int n, qaws_vec3* hv)
{
	qaws_vec3 d[6];
	qaws_scalar q_ydot = integrand_quadratic(f, y, ydot, n);
	unsigned int i, c, k;
	for (i = 0; i < n; i++)
		for (c = 0; c < 3; c++)
		{
			qaws_scalar q_e, q_sum;
			for (k = 0; k < n; k++)
				d[k] = qaws_v3_zero();
			set_comp(&d[i], c, QAWS_ONE);
			q_e = integrand_quadratic(f, y, d, n);
			for (k = 0; k < n; k++)
				d[k] = ydot[k];
			set_comp(&d[i], c, get_comp(&ydot[i], c) + QAWS_ONE);
			q_sum = integrand_quadratic(f, y, d, n);
			set_comp(&hv[i], c, QAWS_LITERAL(0.5) * (q_sum - q_e - q_ydot));
		}
}

/* ================================================================== */
/*  Curves                                                            */
/* ================================================================== */

static integrand_fn curve_integrand(qaws_curve_functional f)
{
	switch (f)
	{
	case QAWS_FUNCTIONAL_LENGTH: return integrand_length;
	case QAWS_FUNCTIONAL_BENDING: return integrand_bending;
	case QAWS_FUNCTIONAL_CURVATURE_SQUARED: return integrand_curvature_squared;
	default: return NULL;
	}
}

#define CURVE_CHANNELS (QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2)

/*
 * Composite rule: every span is split into CURVE_PIECES equal pieces of n
 * Gauss points (square roots of polynomials such as the speed are then
 * integrated to near machine precision).
 *
 * Knots as parameters: span boundaries are knot values, so moving a knot
 * also moves the quadrature nodes and weights of the spans it bounds.
 * With a = knot[left], b = knot[left + 1], a node at fraction s of the span,
 * t = a (1 - s) + b s, and weight W = (b - a) c, the quadrature
 * differentiates exactly as
 *   dQ/da += -c f + W f_t (1 - s),   dQ/db += c f + W f_t s
 * on top of the fixed-node knot derivative of f.
 */
#define CURVE_PIECES 8
typedef struct curve_job
{
	qaws_diff_context const* ctx;
	qaws_curve const* curve;
	integrand_fn f;
	qaws_diff_views const* direction;
	qaws_diff_views* sink;
	double value, tangent, tangent2;
	/* knot terms (NULL when knots are not among the parameters) */
	qaws_field_view const* knot_in;
	qaws_field_view* knot_out;
	qaws_scalar* knots;
	unsigned int knot_count, cp_count;
	/* current node: fraction s of its span, weight W = (b - a) c */
	unsigned int left_knot;
	double s, c;
} curve_job;

typedef qaws_status (*curve_point_fn)(curve_job* job, qaws_scalar t, qaws_scalar weight);

static qaws_status curve_quadrature(curve_job* job, unsigned int n, curve_point_fn fn)
{
	qaws_curve const* curve = job->curve;
	double x[8], w[8];
	unsigned int s, q;
	gauss_rule(n, x, w);
	for (s = 0; s < curve->span_count; s++)
	{
		double a = curve->span_boundaries[s], b = curve->span_boundaries[s + 1];
		if (job->knots)
			job->left_knot = qaws_internal_find_knot_span(job->knots, job->knot_count, curve->degree, job->cp_count,
				(qaws_scalar)(0.5 * (a + b)));
		for (q = 0; q < n * CURVE_PIECES; q++)
		{
			qaws_status st;
			job->s = ((q / n) + 0.5 * (1 + x[q % n])) / CURVE_PIECES;
			job->c = 0.5 * w[q % n] / CURVE_PIECES;
			st = fn(job, (qaws_scalar)(a + (b - a) * job->s), (qaws_scalar)((b - a) * job->c));
			if (st != QAWS_STATUS_OK)
				return st;
		}
	}
	return QAWS_STATUS_OK;
}

/* Sets up knot terms when `views` carries a knot field the curve can
   differentiate. */
static qaws_status curve_job_knots(curve_job* job, qaws_diff_views const* views)
{
	qaws_field_desc fields[8];
	qaws_field_view const* kv;
	unsigned int n = 0, i, got = 0;
	int ok = 0;
	if (!views)
		return QAWS_STATUS_OK;
	kv = qaws_diff_views_find(views, QAWS_FIELD_KNOTS);
	if (!kv || !kv->data)
		return QAWS_STATUS_OK;
	if (qaws_curve_describe_fields(job->curve, fields, 8, &n) != QAWS_STATUS_OK)
		return QAWS_STATUS_OK;
	for (i = 0; i < n && i < 8; i++)
	{
		if (fields[i].field == QAWS_FIELD_KNOTS && (fields[i].capabilities & QAWS_CAP_TANGENT))
		{
			ok = 1;
			job->knot_count = fields[i].count;
		}
		if (fields[i].field == QAWS_FIELD_CONTROL_POINTS)
			job->cp_count = fields[i].count;
	}
	if (!ok)
		return QAWS_STATUS_OK;
	job->knots = (qaws_scalar*)qaws_internal_alloc(NULL, (unsigned long)(sizeof(qaws_scalar) * (job->knot_count + 1)));
	if (!job->knots)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	if (qaws_curve_read_field(job->curve, QAWS_FIELD_KNOTS, job->knots, job->knot_count, &got) != QAWS_STATUS_OK)
	{
		qaws_internal_dealloc(NULL, job->knots);
		job->knots = NULL;
		return QAWS_STATUS_OK;
	}
	job->knot_in = kv;
	return QAWS_STATUS_OK;
}

/* Rates of the node and weight along the knot direction. */
static void curve_node_rates(curve_job const* job, double* t_dot, double* w_dot)
{
	qaws_scalar da = QAWS_ZERO, db = QAWS_ZERO;
	*t_dot = 0;
	*w_dot = 0;
	if (!job->knots || !job->knot_in)
		return;
	qaws_internal_view_read(job->knot_in, job->left_knot, 1, &da);
	qaws_internal_view_read(job->knot_in, job->left_knot + 1, 1, &db);
	*t_dot = (1 - job->s) * da + job->s * db;
	*w_dot = job->c * ((double)db - da);
}

static void curve_jet_values(qaws_curve_jet_3d const* j, qaws_vec3* y)
{
	y[0] = j->d[0];
	y[1] = j->d[1];
	y[2] = j->d[2];
}

static qaws_status curve_eval_point(curve_job* job, qaws_scalar t, qaws_scalar weight)
{
	qaws_curve_jet_3d p, tg, tt;
	qaws_dual3 y[3];
	qaws_dual1 v;
	unsigned int k;
	double t_dot, w_dot;
	qaws_status st;
	curve_node_rates(job, &t_dot, &w_dot);
	st = qaws_internal_curve_tangent_any(job->ctx, job->curve, t, (qaws_scalar)t_dot, CURVE_CHANNELS,
		job->direction, &p, &tg, job->direction ? &tt : NULL);
	if (st != QAWS_STATUS_OK)
		return st;
	for (k = 0; k < 3; k++)
		y[k] = qaws_dual3_make(p.d[k], job->direction ? tg.d[k] : qaws_v3_zero(),
			job->direction ? tt.d[k] : qaws_v3_zero());
	v = job->f(y);
	/* (W f)' = W f' + W' f,  (W f)'' = W f'' + 2 W' f'  (W is linear in the knots) */
	job->value += (double)weight * v.v;
	job->tangent += (double)weight * v.t + w_dot * v.v;
	job->tangent2 += (double)weight * v.tt + 2.0 * w_dot * v.t;
	return QAWS_STATUS_OK;
}

static qaws_status curve_gradient_point(curve_job* job, qaws_scalar t, qaws_scalar weight)
{
	qaws_curve_jet_3d p, tg, bar;
	qaws_vec3 y[3], g[3];
	qaws_scalar fv, t_adj = QAWS_ZERO;
	unsigned int k;
	qaws_status st = qaws_internal_curve_tangent_any(job->ctx, job->curve, t, QAWS_ZERO, CURVE_CHANNELS, NULL, &p, &tg, NULL);
	if (st != QAWS_STATUS_OK)
		return st;
	curve_jet_values(&p, y);
	{
		qaws_dual3 yy[3];
		for (k = 0; k < 3; k++)
			yy[k] = qaws_dual3_const(y[k]);
		fv = job->f(yy).v;
		job->value += (double)weight * fv;
	}
	integrand_gradient(job->f, y, 3, g);
	memset(&bar, 0, sizeof(bar));
	for (k = 0; k < 3; k++)
		bar.d[k] = qaws_v3_scale(g[k], weight);
	bar.channels = CURVE_CHANNELS;
	st = qaws_internal_curve_adjoint_any(job->ctx, job->curve, t, CURVE_CHANNELS, &bar, job->sink,
		job->knot_out ? &t_adj : NULL);
	if (st != QAWS_STATUS_OK || !job->knot_out)
		return st;
	{
		/* moving span boundaries: t_adj = W f_t */
		qaws_scalar ga = (qaws_scalar)(-job->c * fv + (1 - job->s) * t_adj);
		qaws_scalar gb = (qaws_scalar)(job->c * fv + job->s * t_adj);
		qaws_internal_view_add(job->knot_out, job->left_knot, 1, &ga);
		qaws_internal_view_add(job->knot_out, job->left_knot + 1, 1, &gb);
	}
	return QAWS_STATUS_OK;
}

static qaws_status curve_hvp_point(curve_job* job, qaws_scalar t, qaws_scalar weight)
{
	qaws_curve_jet_3d p, tg, bar;
	qaws_vec3 y[3], yd[3], hv[3];
	unsigned int k;
	qaws_status st = qaws_internal_curve_tangent_any(job->ctx, job->curve, t, QAWS_ZERO, CURVE_CHANNELS,
		job->direction, &p, &tg, NULL);
	if (st != QAWS_STATUS_OK)
		return st;
	curve_jet_values(&p, y);
	curve_jet_values(&tg, yd);
	integrand_hess_vec(job->f, y, yd, 3, hv);
	memset(&bar, 0, sizeof(bar));
	for (k = 0; k < 3; k++)
		bar.d[k] = qaws_v3_scale(hv[k], weight);
	bar.channels = CURVE_CHANNELS;
	return qaws_internal_curve_adjoint_any(job->ctx, job->curve, t, CURVE_CHANNELS, &bar, job->sink, NULL);
}

static qaws_status curve_job_init(curve_job* job, qaws_diff_context const* ctx, qaws_curve const* curve,
	qaws_curve_functional functional)
{
	memset(job, 0, sizeof(*job));
	if (!curve)
		return QAWS_STATUS_INVALID_ARGUMENT;
	job->ctx = ctx;
	job->curve = curve;
	job->f = curve_integrand(functional);
	return job->f ? QAWS_STATUS_OK : QAWS_STATUS_INVALID_ARGUMENT;
}

qaws_status qaws_curve_functional_eval(
	qaws_diff_context const* ctx, qaws_curve const* curve, qaws_curve_functional functional,
	unsigned int quadrature, qaws_diff_views const* direction,
	qaws_scalar* out_value, qaws_scalar* out_tangent, qaws_scalar* out_tangent2)
{
	curve_job job;
	qaws_status st = curve_job_init(&job, ctx, curve, functional);
	if (st != QAWS_STATUS_OK)
		return st;
	job.direction = direction;
	st = curve_job_knots(&job, direction);
	if (st == QAWS_STATUS_OK)
		st = curve_quadrature(&job, clamp_rule(quadrature, 6), curve_eval_point);
	qaws_internal_dealloc(NULL, job.knots);
	if (st != QAWS_STATUS_OK)
		return st;
	if (out_value) *out_value = (qaws_scalar)job.value;
	if (out_tangent) *out_tangent = (qaws_scalar)job.tangent;
	if (out_tangent2) *out_tangent2 = (qaws_scalar)job.tangent2;
	return QAWS_STATUS_OK;
}

qaws_status qaws_curve_functional_gradient(
	qaws_diff_context const* ctx, qaws_curve const* curve, qaws_curve_functional functional,
	unsigned int quadrature, qaws_diff_views* gradient, qaws_scalar* out_value)
{
	curve_job job;
	qaws_status st = curve_job_init(&job, ctx, curve, functional);
	if (st != QAWS_STATUS_OK)
		return st;
	job.sink = gradient;
	st = curve_job_knots(&job, gradient);
	job.knot_out = (qaws_field_view*)job.knot_in;
	if (st == QAWS_STATUS_OK)
		st = curve_quadrature(&job, clamp_rule(quadrature, 6), curve_gradient_point);
	qaws_internal_dealloc(NULL, job.knots);
	if (st == QAWS_STATUS_OK && out_value)
		*out_value = (qaws_scalar)job.value;
	return st;
}

qaws_status qaws_curve_functional_hvp(
	qaws_diff_context const* ctx, qaws_curve const* curve, qaws_curve_functional functional,
	unsigned int quadrature, qaws_diff_views const* direction, qaws_diff_views* out_hv)
{
	curve_job job;
	qaws_status st = curve_job_init(&job, ctx, curve, functional);
	if (st != QAWS_STATUS_OK)
		return st;
	/* Direct HVP needs the jet to be linear in the parameters. */
	if (!(qaws_curve_get_diff_capabilities(curve) & QAWS_CAP_LINEAR))
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	/* Knots enter non-linearly. */
	if ((direction && qaws_diff_views_find(direction, QAWS_FIELD_KNOTS)) ||
	    (out_hv && qaws_diff_views_find(out_hv, QAWS_FIELD_KNOTS)))
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	job.direction = direction;
	job.sink = out_hv;
	return curve_quadrature(&job, clamp_rule(quadrature, 6), curve_hvp_point);
}

/* ================================================================== */
/*  Surfaces                                                          */
/* ================================================================== */

static integrand_fn surface_integrand(qaws_surface_functional f)
{
	switch (f)
	{
	case QAWS_FUNCTIONAL_AREA: return integrand_area;
	case QAWS_FUNCTIONAL_THIN_PLATE: return integrand_thin_plate;
	case QAWS_FUNCTIONAL_WILLMORE: return integrand_willmore;
	default: return NULL;
	}
}

typedef qaws_status (*surface_point_fn)(void* user, qaws_scalar u, qaws_scalar v, qaws_scalar weight);

static qaws_status surface_quadrature(qaws_surface const* surface, unsigned int cells, surface_point_fn fn, void* user)
{
	double x[4], w[4];
	double u0 = surface->u_range.min_value, u1 = surface->u_range.max_value;
	double v0 = surface->v_range.min_value, v1 = surface->v_range.max_value;
	double du = (u1 - u0) / cells, dv = (v1 - v0) / cells;
	unsigned int i, j, a, b;
	gauss_rule(4, x, w);
	for (i = 0; i < cells; i++)
		for (j = 0; j < cells; j++)
			for (a = 0; a < 4; a++)
				for (b = 0; b < 4; b++)
				{
					double u = u0 + du * (i + 0.5 + 0.5 * x[a]);
					double v = v0 + dv * (j + 0.5 + 0.5 * x[b]);
					qaws_status st = fn(user, (qaws_scalar)u, (qaws_scalar)v, (qaws_scalar)(0.25 * du * dv * w[a] * w[b]));
					if (st != QAWS_STATUS_OK)
						return st;
				}
	return QAWS_STATUS_OK;
}

/*
 * Knots as surface parameters. Interior knots change the integrand at fixed
 * nodes (the surface jet rules carry them). The domain is
 * [u_knots[p], u_knots[n_u]] x [v_knots[q], v_knots[n_v]]: moving those end
 * knots moves every node, u = u0 + a (u1 - u0) at fraction a, and scales
 * the weights W = (u1 - u0)(v1 - v0) c, so
 *   (W f)'  = W f' + W' f,  (W f)'' = W f'' + 2 W' f' + W'' f,
 *   W' / W = (u1' - u0') / (u1 - u0) + (v1' - v0') / (v1 - v0),
 *   W'' / W = 2 (u1' - u0')(v1' - v0') / ((u1 - u0)(v1 - v0)),
 * with the node rates (1 - a) u0' + a u1' passed as coordinate tangents.
 */
typedef struct surface_job
{
	qaws_diff_context const* ctx;
	qaws_surface const* surface;
	integrand_fn f;
	qaws_diff_views const* direction;
	qaws_diff_views* sink;
	double value, tangent, tangent2;
	/* domain knots (indices of the end knots, their rates or adjoint views) */
	int knots;
	unsigned int u_lo, u_hi, v_lo, v_hi;
	double du0, du1, dv0, dv1;
	qaws_field_view* uk_out;
	qaws_field_view* vk_out;
} surface_job;

/* Finds the domain knots of the U / V knot views (forward rates when
   `rates`, adjoint views otherwise). */
static void surface_job_knots(surface_job* job, qaws_diff_views const* views, int rates)
{
	qaws_field_view const* uv = views ? qaws_diff_views_find(views, QAWS_FIELD_U_KNOTS) : NULL;
	qaws_field_view const* vv = views ? qaws_diff_views_find(views, QAWS_FIELD_V_KNOTS) : NULL;
	qaws_field_desc fields[8];
	unsigned int n = 0, i, nu = 0, nv = 0;
	if (!uv && !vv)
		return;
	if (qaws_surface_describe_fields(job->surface, fields, 8, &n) != QAWS_STATUS_OK)
		return;
	for (i = 0; i < n && i < 8; i++)
	{
		if (fields[i].field == QAWS_FIELD_U_KNOTS) nu = fields[i].count;
		if (fields[i].field == QAWS_FIELD_V_KNOTS) nv = fields[i].count;
	}
	job->knots = 1;
	job->u_lo = job->surface->u_degree;
	job->u_hi = nu > job->surface->u_degree ? nu - job->surface->u_degree - 1 : 0;
	job->v_lo = job->surface->v_degree;
	job->v_hi = nv > job->surface->v_degree ? nv - job->surface->v_degree - 1 : 0;
	if (rates)
	{
		qaws_scalar r = QAWS_ZERO;
		if (uv && nu)
		{
			qaws_internal_view_read(uv, job->u_lo, 1, &r); job->du0 = r;
			qaws_internal_view_read(uv, job->u_hi, 1, &r); job->du1 = r;
		}
		if (vv && nv)
		{
			qaws_internal_view_read(vv, job->v_lo, 1, &r); job->dv0 = r;
			qaws_internal_view_read(vv, job->v_hi, 1, &r); job->dv1 = r;
		}
	}
	else
	{
		job->uk_out = nu ? (qaws_field_view*)uv : NULL;
		job->vk_out = nv ? (qaws_field_view*)vv : NULL;
	}
}

static void surface_fractions(surface_job const* job, qaws_scalar u, qaws_scalar v, double* au, double* av, double* lu, double* lv)
{
	*lu = (double)job->surface->u_range.max_value - job->surface->u_range.min_value;
	*lv = (double)job->surface->v_range.max_value - job->surface->v_range.min_value;
	*au = ((double)u - job->surface->u_range.min_value) / *lu;
	*av = ((double)v - job->surface->v_range.min_value) / *lv;
}

static qaws_status surface_eval_point(void* user, qaws_scalar u, qaws_scalar v, qaws_scalar weight)
{
	surface_job* job = (surface_job*)user;
	qaws_surface_jet p, tg, tt;
	qaws_dual3 y[6];
	qaws_dual1 r;
	qaws_scalar ut = QAWS_ZERO, vt = QAWS_ZERO;
	double wd = 0, wdd = 0;
	unsigned int k;
	qaws_status st;
	if (job->knots)
	{
		double au, av, lu, lv;
		surface_fractions(job, u, v, &au, &av, &lu, &lv);
		ut = (qaws_scalar)((1 - au) * job->du0 + au * job->du1);
		vt = (qaws_scalar)((1 - av) * job->dv0 + av * job->dv1);
		wd = (double)weight * ((job->du1 - job->du0) / lu + (job->dv1 - job->dv0) / lv);
		wdd = (double)weight * 2 * (job->du1 - job->du0) * (job->dv1 - job->dv0) / (lu * lv);
	}
	if (job->direction)
		st = qaws_surface_eval_batch_tangent2(job->ctx, job->surface, &u, &v, &ut, &vt, 1, QAWS_SJET_ORDER2,
			job->direction, &p, &tg, &tt);
	else
		st = qaws_surface_eval_jet(job->surface, u, v, QAWS_SJET_ORDER2, &p);
	if (st != QAWS_STATUS_OK)
		return st;
	for (k = 0; k < 6; k++)
		y[k] = qaws_dual3_make(p.d[k], job->direction ? tg.d[k] : qaws_v3_zero(), job->direction ? tt.d[k] : qaws_v3_zero());
	r = job->f(y);
	job->value += (double)weight * r.v;
	job->tangent += (double)weight * r.t + wd * r.v;
	job->tangent2 += (double)weight * r.tt + 2 * wd * r.t + wdd * r.v;
	return QAWS_STATUS_OK;
}

static qaws_status surface_gradient_point(void* user, qaws_scalar u, qaws_scalar v, qaws_scalar weight)
{
	surface_job* job = (surface_job*)user;
	qaws_surface_jet p, bar;
	qaws_vec3 g[6];
	qaws_scalar ua = QAWS_ZERO, va = QAWS_ZERO;
	double fv;
	unsigned int k;
	qaws_status st = qaws_surface_eval_jet(job->surface, u, v, QAWS_SJET_ORDER2, &p);
	if (st != QAWS_STATUS_OK)
		return st;
	{
		qaws_dual3 yy[6];
		for (k = 0; k < 6; k++)
			yy[k] = qaws_dual3_const(p.d[k]);
		fv = job->f(yy).v;
		job->value += (double)weight * fv;
	}
	integrand_gradient(job->f, p.d, 6, g);
	memset(&bar, 0, sizeof(bar));
	for (k = 0; k < 6; k++)
		bar.d[k] = qaws_v3_scale(g[k], weight);
	bar.channels = QAWS_SJET_ORDER2;
	st = qaws_surface_eval_adjoint(job->ctx, job->surface, u, v, QAWS_SJET_ORDER2, &bar, job->sink,
		job->knots ? &ua : NULL, job->knots ? &va : NULL);
	if (st != QAWS_STATUS_OK || !job->knots)
		return st;
	{
		/* moving domain: node motion (W f_u, W f_v) and weight scaling */
		double au, av, lu, lv, wf = (double)weight * fv;
		qaws_scalar g0, g1;
		surface_fractions(job, u, v, &au, &av, &lu, &lv);
		if (job->uk_out)
		{
			g0 = (qaws_scalar)(-wf / lu + (1 - au) * ua);
			g1 = (qaws_scalar)(wf / lu + au * ua);
			qaws_internal_view_add(job->uk_out, job->u_lo, 1, &g0);
			qaws_internal_view_add(job->uk_out, job->u_hi, 1, &g1);
		}
		if (job->vk_out)
		{
			g0 = (qaws_scalar)(-wf / lv + (1 - av) * va);
			g1 = (qaws_scalar)(wf / lv + av * va);
			qaws_internal_view_add(job->vk_out, job->v_lo, 1, &g0);
			qaws_internal_view_add(job->vk_out, job->v_hi, 1, &g1);
		}
	}
	return QAWS_STATUS_OK;
}

static qaws_status surface_hvp_point(void* user, qaws_scalar u, qaws_scalar v, qaws_scalar weight)
{
	surface_job* job = (surface_job*)user;
	qaws_surface_jet p, tg, bar;
	qaws_vec3 hv[6];
	unsigned int k;
	qaws_status st = qaws_surface_eval_tangent(job->ctx, job->surface, u, v, QAWS_ZERO, QAWS_ZERO, QAWS_SJET_ORDER2,
		job->direction, &p, &tg);
	if (st != QAWS_STATUS_OK)
		return st;
	integrand_hess_vec(job->f, p.d, tg.d, 6, hv);
	memset(&bar, 0, sizeof(bar));
	for (k = 0; k < 6; k++)
		bar.d[k] = qaws_v3_scale(hv[k], weight);
	bar.channels = QAWS_SJET_ORDER2;
	return qaws_surface_eval_adjoint(job->ctx, job->surface, u, v, QAWS_SJET_ORDER2, &bar, job->sink, NULL, NULL);
}

static qaws_status surface_job_init(surface_job* job, qaws_diff_context const* ctx, qaws_surface const* surface,
	qaws_surface_functional functional)
{
	memset(job, 0, sizeof(*job));
	if (!surface)
		return QAWS_STATUS_INVALID_ARGUMENT;
	job->ctx = ctx;
	job->surface = surface;
	job->f = surface_integrand(functional);
	return job->f ? QAWS_STATUS_OK : QAWS_STATUS_INVALID_ARGUMENT;
}

/* Knots enter the surface functionals non-linearly: no direct HVP for
   knot views (as for curves). */
static int surface_knot_terms(qaws_diff_views const* views)
{
	return views && (qaws_diff_views_find(views, QAWS_FIELD_U_KNOTS) || qaws_diff_views_find(views, QAWS_FIELD_V_KNOTS));
}

qaws_status qaws_surface_functional_eval(
	qaws_diff_context const* ctx, qaws_surface const* surface, qaws_surface_functional functional,
	unsigned int quadrature, qaws_diff_views const* direction,
	qaws_scalar* out_value, qaws_scalar* out_tangent, qaws_scalar* out_tangent2)
{
	surface_job job;
	qaws_status st = surface_job_init(&job, ctx, surface, functional);
	if (st != QAWS_STATUS_OK)
		return st;
	surface_job_knots(&job, direction, 1);
	job.direction = direction;
	st = surface_quadrature(surface, quadrature ? quadrature : 8u, surface_eval_point, &job);
	if (st != QAWS_STATUS_OK)
		return st;
	if (out_value) *out_value = (qaws_scalar)job.value;
	if (out_tangent) *out_tangent = (qaws_scalar)job.tangent;
	if (out_tangent2) *out_tangent2 = (qaws_scalar)job.tangent2;
	return QAWS_STATUS_OK;
}

qaws_status qaws_surface_functional_gradient(
	qaws_diff_context const* ctx, qaws_surface const* surface, qaws_surface_functional functional,
	unsigned int quadrature, qaws_diff_views* gradient, qaws_scalar* out_value)
{
	surface_job job;
	qaws_status st = surface_job_init(&job, ctx, surface, functional);
	if (st != QAWS_STATUS_OK)
		return st;
	job.sink = gradient;
	surface_job_knots(&job, gradient, 0);
	st = surface_quadrature(surface, quadrature ? quadrature : 8u, surface_gradient_point, &job);
	if (st == QAWS_STATUS_OK && out_value)
		*out_value = (qaws_scalar)job.value;
	return st;
}

qaws_status qaws_surface_functional_hvp(
	qaws_diff_context const* ctx, qaws_surface const* surface, qaws_surface_functional functional,
	unsigned int quadrature, qaws_diff_views const* direction, qaws_diff_views* out_hv)
{
	surface_job job;
	qaws_status st = surface_job_init(&job, ctx, surface, functional);
	if (st != QAWS_STATUS_OK)
		return st;
	if (surface_knot_terms(direction) || surface_knot_terms(out_hv))
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	if (!(qaws_surface_get_diff_capabilities(surface) & QAWS_CAP_LINEAR))
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	job.direction = direction;
	job.sink = out_hv;
	return surface_quadrature(surface, quadrature ? quadrature : 8u, surface_hvp_point, &job);
}
