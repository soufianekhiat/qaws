#include "qaws_diff_sampling.h"
#include "qaws_diff.h"
#include "internal/qaws_internal_types.h"
#include "internal/qaws_internal_curve.h"
#include "internal/qaws_internal_surface.h"
#include "internal/qaws_internal_diff.h"
#include "core/qaws_dual_core.h"
#include "internal/qaws_internal_basis.h"
#include <math.h>
#include <string.h>

/* ================================================================== */
/*  Quadrature                                                        */
/* ================================================================== */

static void sampling_gauss_rule(unsigned int n, double* x, double* w)
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

/* Composite rule: every interval is split into CDF_PIECES equal pieces of
   n Gauss points, so that the measure (a square root of polynomials, or a
   user density along the curve) is integrated to near machine precision. */
#define CDF_PIECES 8

#define JET_CHANNELS (QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2)

static double v3_dot(qaws_vec3 a, qaws_vec3 b)
{
	return (double)a.x * b.x + (double)a.y * b.y + (double)a.z * b.z;
}

static qaws_vec3 v3_make(double x, double y, double z)
{
	qaws_vec3 r;
	r.x = (qaws_scalar)x;
	r.y = (qaws_scalar)y;
	r.z = (qaws_scalar)z;
	return r;
}

static qaws_vec3 v3_axpy(double a, qaws_vec3 x, qaws_vec3 y)
{
	return v3_make(a * x.x + y.x, a * x.y + y.y, a * x.z + y.z);
}

static void v3_set(qaws_vec3* v, unsigned int c, qaws_scalar x)
{
	if (c == 0) v->x = x;
	else if (c == 1) v->y = x;
	else v->z = x;
}

static qaws_scalar v3_get(qaws_vec3 const* v, unsigned int c)
{
	return c == 0 ? v->x : (c == 1 ? v->y : v->z);
}

/* ------------------------------------------------------------------ */
/*  HVP by polarization of a second order forward pass                */
/* ------------------------------------------------------------------ */

/* q(x) = x^T H x: the second directional derivative of the scalar
   objective along the parameter direction x. */
typedef qaws_status (*hvp_quadratic_fn)(void const* user, qaws_diff_views const* x, double* out);

/*
 * e_j^T H d = (q(d + e_j) - q(d) - q(e_j)) / 2 for every active component j
 * of out_hv (added): exact, for any parameters the forward pass
 * differentiates to second order (knots, weights), at about two passes per
 * parameter. Views with children are refused.
 */
static qaws_status hvp_polarize(void const* user, hvp_quadratic_fn q, qaws_diff_views const* direction, qaws_diff_views* out_hv)
{
	qaws_field_view* fields;
	qaws_diff_views dv;
	qaws_scalar *data, *keep;
	unsigned int f, total = 0, off, e, c;
	double q_d = 0;
	qaws_status st;
	if (out_hv->child_count || direction->child_count)
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	for (f = 0; f < out_hv->field_count; f++)
		total += out_hv->fields[f].count * out_hv->fields[f].components;
	fields = (qaws_field_view*)qaws_internal_alloc(NULL, (unsigned long)(sizeof(qaws_field_view) * (out_hv->field_count + 1)));
	data = (qaws_scalar*)qaws_internal_alloc(NULL, (unsigned long)(sizeof(qaws_scalar) * 2 * (total + 1)));
	if (!fields || !data)
	{
		qaws_internal_dealloc(NULL, fields);
		qaws_internal_dealloc(NULL, data);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	keep = data + total + 1;
	/* the direction laid out like out_hv */
	off = 0;
	for (f = 0; f < out_hv->field_count; f++)
	{
		qaws_field_view const* src = qaws_diff_views_find(direction, out_hv->fields[f].field);
		unsigned int comps = out_hv->fields[f].components;
		fields[f] = qaws_field_view_make(out_hv->fields[f].field, data + off, out_hv->fields[f].count, comps);
		for (e = 0; e < out_hv->fields[f].count; e++)
		{
			qaws_scalar tmp[3] = { 0, 0, 0 };
			if (src)
				qaws_internal_view_read(src, e, comps, tmp);
			for (c = 0; c < comps; c++)
				data[off + e * comps + c] = tmp[c];
		}
		off += out_hv->fields[f].count * comps;
	}
	memcpy(keep, data, sizeof(qaws_scalar) * total);
	dv.fields = fields;
	dv.field_count = out_hv->field_count;
	dv.children = NULL;
	dv.child_count = 0;
	st = q(user, &dv, &q_d);
	off = 0;
	for (f = 0; f < out_hv->field_count && st == QAWS_STATUS_OK; f++)
	{
		qaws_field_view* ov = &out_hv->fields[f];
		unsigned int comps = ov->components;
		for (e = 0; e < ov->count && st == QAWS_STATUS_OK; e++)
		{
			if (!qaws_internal_view_element_active(ov, e))
				continue;
			for (c = 0; c < comps && st == QAWS_STATUS_OK; c++)
			{
				double q_sum = 0, q_e = 0;
				qaws_scalar g[3] = { 0, 0, 0 };
				unsigned int k = off + e * comps + c;
				if (!qaws_internal_view_component_active(ov, c))
					continue;
				/* q(d + e_j) */
				data[k] = keep[k] + QAWS_ONE;
				st = q(user, &dv, &q_sum);
				/* q(e_j) */
				memset(data, 0, sizeof(qaws_scalar) * total);
				data[k] = QAWS_ONE;
				if (st == QAWS_STATUS_OK)
					st = q(user, &dv, &q_e);
				memcpy(data, keep, sizeof(qaws_scalar) * total);
				g[c] = (qaws_scalar)(0.5 * (q_sum - q_d - q_e));
				qaws_internal_view_add(ov, e, comps, g);
			}
		}
		off += ov->count * comps;
	}
	qaws_internal_dealloc(NULL, fields);
	qaws_internal_dealloc(NULL, data);
	return st;
}

/* Shared state of one call. */
typedef struct cdf_job
{
	qaws_diff_context const* ctx;
	qaws_curve const* curve;
	qaws_sample_measure_desc measure;
	qaws_diff_views const* direction;
	unsigned int n, channels;
	double x[8], w[8];
	unsigned int spans;
	double* cum;      /* cum[k]: measure of spans [0, k) */
	double* cum1;     /* first directional derivative of cum */
	double* cum2;     /* second directional derivative of cum */
	/* knot terms (NULL when knots are not among the parameters) */
	qaws_field_view const* knot_in;
	qaws_field_view* knot_out;
	unsigned int* left;   /* left[k]: knot index of span k's start */
} cdf_job;

/* ------------------------------------------------------------------ */
/*  The measure integrand m(C, C', C'') on dual jets                  */
/* ------------------------------------------------------------------ */

/* rho(C) on a dual point: value, rho_g . C', C'^T H C' + rho_g . C'' */
static qaws_dual1 density_dual(qaws_sample_measure_desc const* measure, qaws_dual3 c)
{
	qaws_vec3 g;
	qaws_scalar H[6], v;
	double hq, gt, gtt;
	memset(H, 0, sizeof(H));
	g = qaws_v3_zero();
	v = measure->density(c.v, measure->density_user_data, &g, H);
	gt = v3_dot(g, c.t);
	gtt = v3_dot(g, c.tt);
	hq = (double)H[0] * c.t.x * c.t.x + (double)H[3] * c.t.y * c.t.y + (double)H[5] * c.t.z * c.t.z +
		2.0 * ((double)H[1] * c.t.x * c.t.y + (double)H[2] * c.t.x * c.t.z + (double)H[4] * c.t.y * c.t.z);
	return qaws_dual1_make(v, (qaws_scalar)gt, (qaws_scalar)(hq + gtt));
}

static qaws_dual1 measure_eval(cdf_job const* job, qaws_dual3 const* y)
{
	qaws_dual1 s = qaws_dual3_length(y[1]);
	switch (job->measure.kind)
	{
	case QAWS_MEASURE_CURVATURE:
	{
		/* sqrt(floor^2 |C'|^2 + |C' x C''|^2 / |C'|^4): (kappa |C'|)^2 = |C' x C''|^2 / |C'|^4 */
		qaws_dual3 x = qaws_dual3_cross(y[1], y[2]);
		qaws_dual1 xx = qaws_dual3_dot(x, x), s2 = qaws_dual1_mul(s, s), s4 = qaws_dual1_mul(s2, s2);
		qaws_scalar f2 = job->measure.curvature_floor * job->measure.curvature_floor;
		qaws_dual1 a = qaws_dual1_make(f2 * s2.v, f2 * s2.t, f2 * s2.tt);
		return qaws_dual1_sqrt(qaws_dual1_add(a, qaws_dual1_div(xx, s4)));
	}
	case QAWS_MEASURE_DENSITY:
		return qaws_dual1_mul(density_dual(&job->measure, y[0]), s);
	default:
		return s;
	}
}

/* Dual jets from a primal jet and its first / second tangents (NULL: zero). */
static void jets_to_dual(qaws_curve_jet_3d const* p, qaws_curve_jet_3d const* tg, qaws_curve_jet_3d const* tt, qaws_dual3* y)
{
	unsigned int k;
	for (k = 0; k < 3; k++)
		y[k] = qaws_dual3_make(p->d[k], tg ? tg->d[k] : qaws_v3_zero(), tt ? tt->d[k] : qaws_v3_zero());
}

/* Value of m at plain jets. */
static double measure_value(cdf_job const* job, qaws_curve_jet_3d const* p)
{
	qaws_dual3 y[3];
	jets_to_dual(p, NULL, NULL, y);
	return measure_eval(job, y).v;
}

/* t-derivative of m: the jets move along t (C', C'', C'''). */
static double measure_rate_t(cdf_job const* job, qaws_curve_jet_3d const* p)
{
	qaws_dual3 y[3];
	unsigned int k;
	for (k = 0; k < 3; k++)
		y[k] = qaws_dual3_make(p->d[k], p->d[k + 1], qaws_v3_zero());
	return measure_eval(job, y).t;
}

/* Gradient of m with respect to the jet values (C, C', C''). */
static void measure_gradient(cdf_job const* job, qaws_curve_jet_3d const* p, qaws_vec3* g)
{
	qaws_dual3 y[3];
	unsigned int i, c, k;
	for (i = 0; i < 3; i++)
	{
		g[i] = qaws_v3_zero();
		if (i == 0 && job->measure.kind != QAWS_MEASURE_DENSITY)
			continue;
		if (i == 2 && job->measure.kind != QAWS_MEASURE_CURVATURE)
			continue;
		for (c = 0; c < 3; c++)
		{
			for (k = 0; k < 3; k++)
				y[k] = qaws_dual3_const(p->d[k]);
			v3_set(&y[i].t, c, QAWS_ONE);
			v3_set(&g[i], c, measure_eval(job, y).t);
		}
	}
}

/* d^T H d of m at p. */
static double measure_quadratic(cdf_job const* job, qaws_curve_jet_3d const* p, qaws_vec3 const* d)
{
	qaws_dual3 y[3];
	unsigned int k;
	for (k = 0; k < 3; k++)
		y[k] = qaws_dual3_make(p->d[k], d[k], qaws_v3_zero());
	return measure_eval(job, y).tt;
}

/* H ydot of m at p, by polarization of the quadratic form. */
static void measure_hess_vec(cdf_job const* job, qaws_curve_jet_3d const* p, qaws_vec3 const* ydot, qaws_vec3* hv)
{
	qaws_vec3 d[3];
	double q_ydot = measure_quadratic(job, p, ydot);
	unsigned int i, c, k;
	for (i = 0; i < 3; i++)
	{
		hv[i] = qaws_v3_zero();
		for (c = 0; c < 3; c++)
		{
			double q_e, q_sum;
			for (k = 0; k < 3; k++)
				d[k] = qaws_v3_zero();
			v3_set(&d[i], c, QAWS_ONE);
			q_e = measure_quadratic(job, p, d);
			for (k = 0; k < 3; k++)
				d[k] = ydot[k];
			v3_set(&d[i], c, v3_get(&ydot[i], c) + QAWS_ONE);
			q_sum = measure_quadratic(job, p, d);
			v3_set(&hv[i], c, (qaws_scalar)(0.5 * (q_sum - q_e - q_ydot)));
		}
	}
}

/* ------------------------------------------------------------------ */
/*  Job                                                               */
/* ------------------------------------------------------------------ */

/* Node q (0 .. n * CDF_PIECES - 1) of [a, b] and its weight. */
/* Node q of [a, b]: t = a + s (b - a), W = (b - a) c (s, c may be NULL). */
static void cdf_node(cdf_job const* job, double a, double b, unsigned int q, double* t, double* W, double* s, double* c)
{
	unsigned int piece = q / job->n, k = q % job->n;
	double h = (b - a) / CDF_PIECES, lo = a + piece * h;
	*t = lo + 0.5 * h * (1 + job->x[k]);
	*W = 0.5 * h * job->w[k];
	if (s) *s = (piece + 0.5 * (1 + job->x[k])) / CDF_PIECES;
	if (c) *c = 0.5 * job->w[k] / CDF_PIECES;
}

/* Rates of the start and end knots of span k along the knot direction. */
static void cdf_span_rates(cdf_job const* job, unsigned int k, double* da, double* db)
{
	qaws_scalar ra = QAWS_ZERO, rb = QAWS_ZERO;
	if (job->left && job->knot_in)
	{
		qaws_internal_view_read(job->knot_in, job->left[k], 1, &ra);
		qaws_internal_view_read(job->knot_in, job->left[k] + 1, 1, &rb);
	}
	*da = ra;
	*db = rb;
}

/* Jets at t moving at rate t_dot, along direction; tangent2 may be NULL. */
static qaws_status cdf_jets(cdf_job const* job, double t, double t_dot, qaws_diff_views const* direction,
	qaws_curve_jet_3d* p, qaws_curve_jet_3d* tg, qaws_curve_jet_3d* tt)
{
	return qaws_internal_curve_tangent_any(job->ctx, job->curve, (qaws_scalar)t, (qaws_scalar)t_dot, job->channels,
		direction, p, tg, tt);
}

/* Measure of [a, b] and its first and second directional derivatives along
   the job direction (d1, d2 may be NULL); da, db are the rates of a and b
   when they are knots: the nodes then move at (1 - s) da + s db and the
   weights at c (db - da), so (W m)' = W m' + W' m, (W m)'' = W m'' + 2 W' m'. */
static qaws_status cdf_integrate(cdf_job const* job, double a, double b, double da, double db,
	double* v, double* d1, double* d2)
{
	unsigned int q;
	*v = 0;
	if (d1) *d1 = 0;
	if (d2) *d2 = 0;
	if (!(b > a))
		return QAWS_STATUS_OK;
	for (q = 0; q < job->n * CDF_PIECES; q++)
	{
		double t, W, s, c, wd;
		qaws_curve_jet_3d p, tg, tt;
		qaws_dual3 y[3];
		qaws_dual1 m;
		qaws_status st;
		cdf_node(job, a, b, q, &t, &W, &s, &c);
		wd = d1 ? c * (db - da) : 0;
		st = cdf_jets(job, t, d1 ? (1 - s) * da + s * db : 0, job->direction, &p, &tg, job->direction ? &tt : NULL);
		if (st != QAWS_STATUS_OK)
			return st;
		jets_to_dual(&p, job->direction ? &tg : NULL, job->direction ? &tt : NULL, y);
		m = measure_eval(job, y);
		*v += W * m.v;
		if (d1) *d1 += W * m.t + wd * m.v;
		if (d2) *d2 += W * m.tt + 2 * wd * m.t;
	}
	return QAWS_STATUS_OK;
}

static void cdf_job_free(cdf_job* job)
{
	qaws_internal_dealloc(NULL, job->cum);
	qaws_internal_dealloc(NULL, job->left);
	job->cum = NULL;
	job->left = NULL;
}

/* Finds the left knot of every span when `views` carries a knot field the
   curve differentiates (other curves ignore the view). */
static qaws_status cdf_job_knots(cdf_job* job, qaws_diff_views const* views)
{
	qaws_field_desc fields[8];
	qaws_field_view const* kv = views ? qaws_diff_views_find(views, QAWS_FIELD_KNOTS) : NULL;
	qaws_scalar* knots;
	unsigned int n = 0, i, k, got = 0, knot_count = 0, cp_count = 0;
	int ok = 0;
	if (!kv || !kv->data)
		return QAWS_STATUS_OK;
	if (qaws_curve_describe_fields(job->curve, fields, 8, &n) != QAWS_STATUS_OK)
		return QAWS_STATUS_OK;
	for (i = 0; i < n && i < 8; i++)
	{
		if (fields[i].field == QAWS_FIELD_KNOTS && (fields[i].capabilities & QAWS_CAP_TANGENT))
		{
			ok = 1;
			knot_count = fields[i].count;
		}
		if (fields[i].field == QAWS_FIELD_CONTROL_POINTS)
			cp_count = fields[i].count;
	}
	if (!ok || knot_count < 2)
		return QAWS_STATUS_OK;
	knots = (qaws_scalar*)qaws_internal_alloc(NULL, (unsigned long)(sizeof(qaws_scalar) * knot_count));
	job->left = (unsigned int*)qaws_internal_alloc(NULL, (unsigned long)(sizeof(unsigned int) * job->spans));
	if (!knots || !job->left)
	{
		qaws_internal_dealloc(NULL, knots);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	if (qaws_curve_read_field(job->curve, QAWS_FIELD_KNOTS, knots, knot_count, &got) != QAWS_STATUS_OK)
	{
		qaws_internal_dealloc(NULL, knots);
		qaws_internal_dealloc(NULL, job->left);
		job->left = NULL;
		return QAWS_STATUS_OK;
	}
	for (k = 0; k < job->spans; k++)
		job->left[k] = qaws_internal_find_knot_span(knots, knot_count, job->curve->degree, cp_count,
			(qaws_scalar)(0.5 * ((double)job->curve->span_boundaries[k] + job->curve->span_boundaries[k + 1])));
	qaws_internal_dealloc(NULL, knots);
	return QAWS_STATUS_OK;
}

/* Validates the measure, finds the knot terms (knot views move the span
   boundaries; `knots` = 0 rejects them) and builds the cumulative span
   measures. */
static qaws_status cdf_job_init(cdf_job* job, qaws_diff_context const* ctx, qaws_curve const* curve,
	qaws_sample_measure_desc const* measure, unsigned int quadrature, qaws_diff_views const* direction,
	qaws_diff_views* other, int knots)
{
	unsigned int k;
	qaws_status st;
	memset(job, 0, sizeof(*job));
	if (!curve || curve->span_count == 0 || !curve->span_boundaries)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (measure)
	{
		job->measure = *measure;
		if ((measure->kind == QAWS_MEASURE_CURVATURE && !(measure->curvature_floor > 0)) ||
		    (measure->kind == QAWS_MEASURE_DENSITY && !measure->density) ||
		    (unsigned int)measure->kind > (unsigned int)QAWS_MEASURE_DENSITY)
			return QAWS_STATUS_INVALID_ARGUMENT;
	}
	if (!knots && ((direction && qaws_diff_views_find(direction, QAWS_FIELD_KNOTS)) ||
	    (other && qaws_diff_views_find(other, QAWS_FIELD_KNOTS))))
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	job->ctx = ctx;
	job->curve = curve;
	job->direction = direction;
	/* the curvature measure reads C'', whose t-rate needs C''' */
	job->channels = JET_CHANNELS | (job->measure.kind == QAWS_MEASURE_CURVATURE ? QAWS_EVAL_FLAG_D3 : 0u);
	job->n = quadrature == 0 ? 8 : (quadrature < 2 ? 2 : (quadrature > 8 ? 8 : quadrature));
	sampling_gauss_rule(job->n, job->x, job->w);
	job->spans = curve->span_count;
	if (knots)
	{
		st = cdf_job_knots(job, direction ? direction : other);
		if (st != QAWS_STATUS_OK)
			return st;
		if (job->left && direction)
			job->knot_in = qaws_diff_views_find(direction, QAWS_FIELD_KNOTS);
		else if (job->left)
			job->knot_out = (qaws_field_view*)qaws_diff_views_find(other, QAWS_FIELD_KNOTS);
	}
	job->cum = (double*)qaws_internal_alloc(NULL, (unsigned long)(sizeof(double) * 3 * (job->spans + 1)));
	if (!job->cum)
	{
		cdf_job_free(job);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	job->cum1 = job->cum + job->spans + 1;
	job->cum2 = job->cum1 + job->spans + 1;
	job->cum[0] = job->cum1[0] = job->cum2[0] = 0;
	for (k = 0; k < job->spans; k++)
	{
		double v, d1, d2, da, db;
		cdf_span_rates(job, k, &da, &db);
		st = cdf_integrate(job, curve->span_boundaries[k], curve->span_boundaries[k + 1], da, db, &v, &d1, &d2);
		if (st != QAWS_STATUS_OK)
		{
			cdf_job_free(job);
			return st;
		}
		job->cum[k + 1] = job->cum[k] + v;
		job->cum1[k + 1] = job->cum1[k] + d1;
		job->cum2[k + 1] = job->cum2[k] + d2;
	}
	return QAWS_STATUS_OK;
}

static double cdf_total(cdf_job const* job)
{
	return job->cum[job->spans];
}

/* Parameter and span of the measure sigma (clamped to the curve):
   safeguarded Newton on the span holding sigma. */
static qaws_status cdf_solve(cdf_job const* job, double sigma, double* out_t, unsigned int* out_span)
{
	double a, b, lo, hi, t, total = cdf_total(job);
	unsigned int k = 0, it;
	qaws_scalar const* bounds = job->curve->span_boundaries;
	if (sigma < 0) sigma = 0;
	if (sigma > total) sigma = total;
	while (k + 1 < job->spans && job->cum[k + 1] <= sigma)
		k++;
	a = bounds[k];
	b = bounds[k + 1];
	lo = a;
	hi = b;
	{
		double len = job->cum[k + 1] - job->cum[k];
		t = len > 0 ? a + (b - a) * (sigma - job->cum[k]) / len : a;
	}
	for (it = 0; it < 60; it++)
	{
		qaws_curve_jet_3d p, tg;
		double v, f, m;
		qaws_status st = cdf_integrate(job, a, t, 0, 0, &v, NULL, NULL);
		if (st != QAWS_STATUS_OK)
			return st;
		f = job->cum[k] + v - sigma;
		if (f > 0) hi = t; else lo = t;
		if (fabs(f) <= 1e-15 * (total > 1 ? total : 1))
			break;
		st = cdf_jets(job, t, 0, NULL, &p, &tg, NULL);
		if (st != QAWS_STATUS_OK)
			return st;
		m = measure_value(job, &p);
		t = m > 0 ? t - f / m : 0.5 * (lo + hi);
		if (!(t > lo && t < hi))
			t = 0.5 * (lo + hi);
		if (hi - lo <= 1e-15 * (fabs(b - a) + 1e-300))
			break;
	}
	*out_t = t;
	*out_span = k;
	return QAWS_STATUS_OK;
}

static double target_sigma(qaws_cdf_target const* tg, double total)
{
	return (double)tg->distance + (double)tg->fraction * total;
}

/* ================================================================== */
/*  Forward                                                           */
/* ================================================================== */

qaws_status qaws_curve_cdf_sample_tangent(
	qaws_diff_context const* ctx,
	qaws_curve const* curve,
	qaws_sample_measure_desc const* measure,
	qaws_cdf_target const* targets,
	qaws_scalar const* distance_tangent,
	qaws_scalar const* distance_tangent2,
	unsigned int count,
	unsigned int quadrature,
	qaws_diff_views const* param_tangent,
	qaws_cdf_sample* out_value,
	qaws_cdf_sample* out_tangent,
	qaws_cdf_sample* out_tangent2,
	qaws_scalar* out_total)
{
	cdf_job job;
	unsigned int i;
	qaws_status st;
	if (!targets && count)
		return QAWS_STATUS_INVALID_ARGUMENT;
	st = cdf_job_init(&job, ctx, curve, measure, quadrature, param_tangent, NULL, 1);
	if (st != QAWS_STATUS_OK)
		return st;
	if (out_total)
		*out_total = (qaws_scalar)cdf_total(&job);
	for (i = 0; i < count && st == QAWS_STATUS_OK; i++)
	{
		double total = cdf_total(&job), t, m, me, mt, sig1, sig2, M1, M2, v, td, tdd;
		unsigned int k;
		qaws_curve_jet_3d p, tg, tt;
		qaws_dual3 y[3];
		st = cdf_solve(&job, target_sigma(&targets[i], total), &t, &k);
		if (st != QAWS_STATUS_OK)
			break;
		st = cdf_jets(&job, t, 0, param_tangent, &p, &tg, NULL);
		if (st != QAWS_STATUS_OK)
			break;
		if (out_value)
		{
			out_value[i].t = (qaws_scalar)t;
			out_value[i].position = p.d[0];
		}
		if (!out_tangent && !out_tangent2)
			continue;
		jets_to_dual(&p, param_tangent ? &tg : NULL, NULL, y);
		{
			qaws_dual1 md = measure_eval(&job, y);
			m = md.v;
			me = md.t;
		}
		if (!(m > QAWS_EPSILON))
		{
			qaws_internal_diff_report_note(ctx, QAWS_DIFF_SMOOTH, QAWS_DIFF_ILL_CONDITIONED, 0, i);
			if (out_tangent) memset(&out_tangent[i], 0, sizeof(out_tangent[i]));
			if (out_tangent2) memset(&out_tangent2[i], 0, sizeof(out_tangent2[i]));
			continue;
		}
		qaws_internal_diff_report_note(ctx, QAWS_DIFF_SMOOTH, QAWS_DIFF_VALID, 0, i);
		{
			/* [a_k, t] at fixed t: only its start knot moves */
			double da, db;
			cdf_span_rates(&job, k, &da, &db);
			st = cdf_integrate(&job, curve->span_boundaries[k], t, da, 0, &v, &M1, &M2);
		}
		if (st != QAWS_STATUS_OK)
			break;
		M1 += job.cum1[k];
		M2 += job.cum2[k];
		sig1 = (distance_tangent ? (double)distance_tangent[i] : 0) + (double)targets[i].fraction * job.cum1[job.spans];
		sig2 = (distance_tangent2 ? (double)distance_tangent2[i] : 0) + (double)targets[i].fraction * job.cum2[job.spans];
		td = (sig1 - M1) / m;
		mt = out_tangent2 ? measure_rate_t(&job, &p) : 0;
		tdd = (sig2 - M2 - 2 * me * td - mt * td * td) / m;
		/* the point along the curve moving at t' (straight path in (t, theta)),
		   plus C' t'' for the curvature of the path in t */
		st = cdf_jets(&job, t, td, param_tangent, &p, &tg, out_tangent2 ? &tt : NULL);
		if (st != QAWS_STATUS_OK)
			break;
		if (out_tangent)
		{
			out_tangent[i].t = (qaws_scalar)td;
			out_tangent[i].position = tg.d[0];
		}
		if (out_tangent2)
		{
			out_tangent2[i].t = (qaws_scalar)tdd;
			out_tangent2[i].position = v3_axpy(tdd, p.d[1], tt.d[0]);
		}
	}
	cdf_job_free(&job);
	return st;
}

/* ================================================================== */
/*  Backward                                                          */
/* ================================================================== */

/* Per-sample state shared by the adjoint and the HVP. */
typedef struct cdf_sample
{
	double t, m, lambda, lambda_dot, t_dot;
	unsigned int span;
	qaws_vec3 g[3];       /* gradient of m with respect to (C, C', C'') at t */
} cdf_sample;

static qaws_status cdf_add_jet(cdf_job const* job, double t, qaws_vec3 const* bar3, qaws_diff_views* sink,
	qaws_scalar* t_adj)
{
	qaws_curve_jet_3d bar;
	memset(&bar, 0, sizeof(bar));
	bar.d[0] = bar3[0];
	bar.d[1] = bar3[1];
	bar.d[2] = bar3[2];
	bar.channels = JET_CHANNELS;
	return qaws_internal_curve_adjoint_any(job->ctx, job->curve, (qaws_scalar)t, JET_CHANNELS, &bar, sink, t_adj);
}

/* Quadrature node of span `span` at fraction s with W = (b - a) c; `full`
   when its end b is a knot too (a partial span ends at a sample). */
typedef struct cdf_knot_node
{
	unsigned int span;
	double s, c;
	int full;
} cdf_knot_node;

/* Pullback of c1 * dm/dy + c2 * H_m ydot at the node t, weight W; with knot
   adjoints, the moving node and weight add
     a_bar += c1 (-c m + (1 - s) W m_t),  b_bar += c1 (c m + s W m_t). */
static qaws_status cdf_node_pullback(cdf_job const* job, double t, double W, double c1, double c2,
	cdf_knot_node const* node, qaws_diff_views* sink)
{
	qaws_curve_jet_3d p, tg;
	qaws_vec3 g[3], hv[3], bar[3];
	unsigned int k;
	double m;
	qaws_scalar t_adj = QAWS_ZERO;
	int knots = job->knot_out && node && c2 == 0;
	qaws_status st = cdf_jets(job, t, 0, c2 != 0 ? job->direction : NULL, &p, &tg, NULL);
	if (st != QAWS_STATUS_OK)
		return st;
	m = measure_value(job, &p);
	if (!(m > QAWS_EPSILON))
		return QAWS_STATUS_OK;
	measure_gradient(job, &p, g);
	for (k = 0; k < 3; k++)
		bar[k] = v3_make(W * c1 * g[k].x, W * c1 * g[k].y, W * c1 * g[k].z);
	if (c2 != 0)
	{
		qaws_vec3 yd[3];
		for (k = 0; k < 3; k++)
			yd[k] = tg.d[k];
		measure_hess_vec(job, &p, yd, hv);
		for (k = 0; k < 3; k++)
			bar[k] = v3_axpy(W * c2, hv[k], bar[k]);
	}
	st = cdf_add_jet(job, t, bar, sink, knots ? &t_adj : NULL);
	if (st != QAWS_STATUS_OK || !knots)
		return st;
	{
		qaws_scalar ga = (qaws_scalar)(-c1 * node->c * m + (1 - node->s) * t_adj);
		qaws_scalar gb = (qaws_scalar)(c1 * node->c * m + node->s * t_adj);
		qaws_internal_view_add(job->knot_out, job->left[node->span], 1, &ga);
		if (node->full)
			qaws_internal_view_add(job->knot_out, job->left[node->span] + 1, 1, &gb);
	}
	return QAWS_STATUS_OK;
}

/*
 * Pulls back w1(tau) dm/dy + w2(tau) H_m ydot over the measure integrals:
 * every full span k gets
 *   w = F - sum of the samples whose span lies beyond k
 * and every sample subtracts its own weight over [a_k, t_i]. c1 / c2 are the
 * per-sample weights of the two terms (c2 may be NULL).
 */
static qaws_status cdf_pullback_measure(cdf_job const* job, cdf_sample const* smp, unsigned int count,
	double const* c1, double const* c2, double F1, double F2, qaws_diff_views* sink)
{
	unsigned int S = job->spans, k, i, q;
	double* beyond = (double*)qaws_internal_alloc(NULL, (unsigned long)(sizeof(double) * 2 * (S + 1)));
	qaws_status st = QAWS_STATUS_OK;
	if (!beyond)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	memset(beyond, 0, sizeof(double) * 2 * (S + 1));
	for (i = 0; i < count; i++)
		if (smp[i].m > QAWS_EPSILON && smp[i].span > 0)
		{
			beyond[smp[i].span - 1] += c1[i];
			if (c2) beyond[S + 1 + smp[i].span - 1] += c2[i];
		}
	for (k = S; k-- > 0;)
		if (k + 1 < S)
		{
			beyond[k] += beyond[k + 1];
			beyond[S + 1 + k] += beyond[S + 1 + k + 1];
		}
	for (k = 0; k < S && st == QAWS_STATUS_OK; k++)
	{
		double w1 = F1 - beyond[k], w2 = c2 ? F2 - beyond[S + 1 + k] : 0;
		double a = job->curve->span_boundaries[k], b = job->curve->span_boundaries[k + 1];
		if (w1 == 0 && w2 == 0)
			continue;
		for (q = 0; q < job->n * CDF_PIECES && st == QAWS_STATUS_OK; q++)
		{
			double t, W;
			cdf_knot_node node;
			node.span = k;
			node.full = 1;
			cdf_node(job, a, b, q, &t, &W, &node.s, &node.c);
			st = cdf_node_pullback(job, t, W, w1, w2, &node, sink);
		}
	}
	for (i = 0; i < count && st == QAWS_STATUS_OK; i++)
	{
		double a = job->curve->span_boundaries[smp[i].span], b = smp[i].t;
		if (!(smp[i].m > QAWS_EPSILON) || !(b > a))
			continue;
		for (q = 0; q < job->n * CDF_PIECES && st == QAWS_STATUS_OK; q++)
		{
			double t, W;
			cdf_knot_node node;
			node.span = smp[i].span;
			node.full = 0;
			cdf_node(job, a, b, q, &t, &W, &node.s, &node.c);
			st = cdf_node_pullback(job, t, W, -c1[i], c2 ? -c2[i] : 0, &node, sink);
		}
	}
	qaws_internal_dealloc(NULL, beyond);
	return st;
}

/* Solves the samples and their lambda = (t_bar + p_bar . C') / m; the
   position pullback J^T p_bar is added into sink when it is given. */
static qaws_status cdf_samples(cdf_job const* job, qaws_cdf_target const* targets, unsigned int count,
	qaws_cdf_sample const* adjoint, qaws_diff_views* sink, cdf_sample* smp)
{
	unsigned int i;
	double total = cdf_total(job);
	for (i = 0; i < count; i++)
	{
		qaws_curve_jet_3d p, tg;
		qaws_status st = cdf_solve(job, target_sigma(&targets[i], total), &smp[i].t, &smp[i].span);
		qaws_scalar t_adj = QAWS_ZERO;
		if (st != QAWS_STATUS_OK)
			return st;
		st = cdf_jets(job, smp[i].t, 0, NULL, &p, &tg, NULL);
		if (st != QAWS_STATUS_OK)
			return st;
		smp[i].m = measure_value(job, &p);
		smp[i].lambda = 0;
		smp[i].lambda_dot = 0;
		smp[i].t_dot = 0;
		if (!(smp[i].m > QAWS_EPSILON))
		{
			qaws_internal_diff_report_note(job->ctx, QAWS_DIFF_SMOOTH, QAWS_DIFF_ILL_CONDITIONED, 0, i);
			continue;
		}
		qaws_internal_diff_report_note(job->ctx, QAWS_DIFF_SMOOTH, QAWS_DIFF_VALID, 0, i);
		measure_gradient(job, &p, smp[i].g);
		if (sink)
		{
			qaws_curve_jet_3d bar;
			memset(&bar, 0, sizeof(bar));
			bar.d[0] = adjoint[i].position;
			bar.channels = QAWS_EVAL_FLAG_POSITION;
			st = qaws_internal_curve_adjoint_any(job->ctx, job->curve, (qaws_scalar)smp[i].t, QAWS_EVAL_FLAG_POSITION,
				&bar, sink, &t_adj);
			if (st != QAWS_STATUS_OK)
				return st;
		}
		else
			t_adj = (qaws_scalar)v3_dot(adjoint[i].position, p.d[1]);
		smp[i].lambda = ((double)adjoint[i].t + t_adj) / smp[i].m;
	}
	return QAWS_STATUS_OK;
}

qaws_status qaws_curve_cdf_sample_adjoint(
	qaws_diff_context const* ctx,
	qaws_curve const* curve,
	qaws_sample_measure_desc const* measure,
	qaws_cdf_target const* targets,
	unsigned int count,
	unsigned int quadrature,
	qaws_cdf_sample const* adjoint,
	qaws_diff_views* param_adjoint,
	qaws_scalar* distance_adjoint)
{
	cdf_job job;
	cdf_sample* smp;
	double* c1;
	double F = 0;
	unsigned int i;
	qaws_status st;
	if ((!targets || !adjoint) && count)
		return QAWS_STATUS_INVALID_ARGUMENT;
	st = cdf_job_init(&job, ctx, curve, measure, quadrature, NULL, param_adjoint, 1);
	if (st != QAWS_STATUS_OK)
		return st;
	smp = (cdf_sample*)qaws_internal_alloc(NULL, (unsigned long)((sizeof(cdf_sample) + sizeof(double)) * (count + 1)));
	if (!smp)
	{
		cdf_job_free(&job);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	c1 = (double*)(smp + count + 1);
	st = cdf_samples(&job, targets, count, adjoint, param_adjoint, smp);
	for (i = 0; i < count && st == QAWS_STATUS_OK; i++)
	{
		c1[i] = smp[i].lambda;
		F += smp[i].lambda * (double)targets[i].fraction;
		if (distance_adjoint)
			distance_adjoint[i] += (qaws_scalar)smp[i].lambda;
	}
	if (st == QAWS_STATUS_OK && param_adjoint)
		st = cdf_pullback_measure(&job, smp, count, c1, NULL, F, 0, param_adjoint);
	qaws_internal_dealloc(NULL, smp);
	cdf_job_free(&job);
	return st;
}

/*
 * HVP, forward over reverse. With y = (C, C', C'') the jet, the gradient is
 *   G = sum_i J(t_i)^T p_bar_i + lambda_i (f_i grad M_total - grad M(t_i)),
 *   grad M(t) = integral of J_y^T dm/dy,
 * and, for a family linear in its fields (J, J_y independent of the
 * parameters), its derivative along the direction is
 *   sum_i  J_D1(t_i)^T p_bar_i t'_i - lambda_i J_y(t_i)^T dm/dy(t_i) t'_i
 *        + lambda'_i (f_i grad M_total - grad M(t_i))
 *        + lambda_i  (f_i grad' M_total - grad' M(t_i))
 * with grad' M the integral of J_y^T H_m ydot, lambda' = (p_bar . dC'(t, t')
 * - lambda m') / m and m' the rate of m along (t', direction).
 */
/* The quadratic form q(x) = sum adjoint . sample''(x) of curve sampling. */
typedef struct cdf_quadratic
{
	qaws_diff_context const* ctx;
	qaws_curve const* curve;
	qaws_sample_measure_desc const* measure;
	qaws_cdf_target const* targets;
	unsigned int count, quadrature;
	qaws_cdf_sample const* adjoint;
	qaws_cdf_sample* scratch;
} cdf_quadratic;

static qaws_status cdf_second(void const* user, qaws_diff_views const* dir, double* out)
{
	cdf_quadratic const* a = (cdf_quadratic const*)user;
	unsigned int i;
	qaws_status st = qaws_curve_cdf_sample_tangent(a->ctx, a->curve, a->measure, a->targets, NULL, NULL, a->count, a->quadrature, dir,
		NULL, NULL, a->scratch, NULL);
	*out = 0;
	if (st != QAWS_STATUS_OK)
		return st;
	for (i = 0; i < a->count; i++)
		*out += (double)a->adjoint[i].t * a->scratch[i].t + v3_dot(a->adjoint[i].position, a->scratch[i].position);
	return QAWS_STATUS_OK;
}

qaws_status qaws_curve_cdf_sample_hvp(
	qaws_diff_context const* ctx,
	qaws_curve const* curve,
	qaws_sample_measure_desc const* measure,
	qaws_cdf_target const* targets,
	unsigned int count,
	unsigned int quadrature,
	qaws_cdf_sample const* adjoint,
	qaws_diff_views const* direction,
	qaws_diff_views* out_hv)
{
	cdf_job job;
	cdf_sample* smp;
	double *c1, *c2;
	double F1 = 0, F2 = 0;
	unsigned int i, k;
	qaws_status st;
	if ((!targets || !adjoint) && count)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (!curve)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (!direction || !out_hv)
		return QAWS_STATUS_OK;
	if (!(qaws_curve_get_diff_capabilities(curve) & QAWS_CAP_LINEAR) || qaws_diff_views_find(direction, QAWS_FIELD_KNOTS) ||
	    qaws_diff_views_find(out_hv, QAWS_FIELD_KNOTS))
	{
		/* knots move the quadrature spans and rational families depend on
		   their weights nonlinearly: polarize the second order forward pass */
		cdf_quadratic a;
		a.ctx = ctx;
		a.curve = curve;
		a.measure = measure;
		a.targets = targets;
		a.count = count;
		a.quadrature = quadrature;
		a.adjoint = adjoint;
		a.scratch = (qaws_cdf_sample*)qaws_internal_alloc(NULL, (unsigned long)(sizeof(qaws_cdf_sample) * (count + 1)));
		if (!a.scratch)
			return QAWS_STATUS_ALLOCATION_FAILURE;
		st = hvp_polarize(&a, cdf_second, direction, out_hv);
		qaws_internal_dealloc(NULL, a.scratch);
		return st;
	}
	st = cdf_job_init(&job, ctx, curve, measure, quadrature, direction, out_hv, 0);
	if (st != QAWS_STATUS_OK)
		return st;
	smp = (cdf_sample*)qaws_internal_alloc(NULL, (unsigned long)((sizeof(cdf_sample) + 2 * sizeof(double)) * (count + 1)));
	if (!smp)
	{
		cdf_job_free(&job);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	c1 = (double*)(smp + count + 1);
	c2 = c1 + count + 1;
	st = cdf_samples(&job, targets, count, adjoint, NULL, smp);
	for (i = 0; i < count && st == QAWS_STATUS_OK; i++)
	{
		qaws_curve_jet_3d p, tg;
		qaws_dual3 y[3];
		double v, M1, mdot;
		qaws_vec3 bar[3];
		cdf_sample* a = &smp[i];
		c1[i] = c2[i] = 0;
		if (!(a->m > QAWS_EPSILON))
			continue;
		st = cdf_integrate(&job, curve->span_boundaries[a->span], a->t, 0, 0, &v, &M1, NULL);
		if (st != QAWS_STATUS_OK)
			break;
		M1 += job.cum1[a->span];
		a->t_dot = ((double)targets[i].fraction * job.cum1[job.spans] - M1) / a->m;
		/* jets moving along (t', direction): tg.d[1] = dC'(t, t') */
		st = cdf_jets(&job, a->t, a->t_dot, direction, &p, &tg, NULL);
		if (st != QAWS_STATUS_OK)
			break;
		jets_to_dual(&p, &tg, NULL, y);
		mdot = measure_eval(&job, y).t;
		a->lambda_dot = (v3_dot(adjoint[i].position, tg.d[1]) - a->lambda * mdot) / a->m;
		c1[i] = a->lambda_dot;
		c2[i] = a->lambda;
		F1 += a->lambda_dot * (double)targets[i].fraction;
		F2 += a->lambda * (double)targets[i].fraction;
		/* moving evaluation point: J_D1^T p_bar t' - lambda J_y^T dm/dy t' */
		for (k = 0; k < 3; k++)
			bar[k] = v3_make(-a->lambda * a->t_dot * a->g[k].x, -a->lambda * a->t_dot * a->g[k].y,
				-a->lambda * a->t_dot * a->g[k].z);
		bar[1] = v3_axpy(a->t_dot, adjoint[i].position, bar[1]);
		st = cdf_add_jet(&job, a->t, bar, out_hv, NULL);
	}
	if (st == QAWS_STATUS_OK)
		st = cdf_pullback_measure(&job, smp, count, c1, c2, F1, F2, out_hv);
	qaws_internal_dealloc(NULL, smp);
	cdf_job_free(&job);
	return st;
}

/* ================================================================== */
/*  Surfaces                                                          */
/* ================================================================== */

typedef struct scdf_job
{
	qaws_diff_context const* ctx;
	qaws_surface const* surface;
	qaws_sample_measure_desc measure;
	qaws_diff_views const* direction;   /* parameter direction of the forward pass (may be NULL) */
	unsigned int cells, n;
	double x[8], w[8];
	double u0, u1, v0, v1, hu, hv;
	/* rates of the domain end knots along the direction: the grid of
	   cells and every quadrature node move with them */
	double du0, du1, dv0, dv1;
	qaws_dual1* cumA;                   /* cells + 1: marginal measure of the full u cells, along direction */
} scdf_job;

static qaws_dual1 dual_c(double v)
{
	return qaws_dual1_const((qaws_scalar)v);
}

/* a x + b */
static qaws_dual1 dual_axpy(double a, qaws_dual1 x, double b)
{
	return qaws_dual1_make((qaws_scalar)(a * x.v + b), (qaws_scalar)(a * x.t), (qaws_scalar)(a * x.tt));
}

/* a x + y */
static qaws_dual1 dual_lin(double a, qaws_dual1 x, qaws_dual1 y)
{
	return qaws_dual1_make((qaws_scalar)(a * x.v + y.v), (qaws_scalar)(a * x.t + y.t), (qaws_scalar)(a * x.tt + y.tt));
}

/* Start and cell size of the u (or v) grid, moving with the end knots when
   `on` (knots are linear in the direction: no second rates). */
static void scdf_grid(scdf_job const* job, int v_dir, int on, qaws_dual1* start, qaws_dual1* step)
{
	double a = v_dir ? job->v0 : job->u0, h = v_dir ? job->hv : job->hu;
	double da = on ? (v_dir ? job->dv0 : job->du0) : 0, db = on ? (v_dir ? job->dv1 : job->du1) : 0;
	*start = qaws_dual1_make((qaws_scalar)a, (qaws_scalar)da, 0);
	*step = qaws_dual1_make((qaws_scalar)h, (qaws_scalar)((db - da) / job->cells), 0);
}

/* The end of the v domain, moving with its knot when `on`. */
static qaws_dual1 scdf_v1(scdf_job const* job, int on)
{
	return qaws_dual1_make((qaws_scalar)job->v1, (qaws_scalar)(on ? job->dv1 : 0), 0);
}

/* Jet entries the measure reads: (S, S_u, S_v), plus (S_uu, S_uv, S_vv)
   for the curvature measure. */
#define SCDF_Y 6

static unsigned int scdf_y_count(scdf_job const* job)
{
	return job->measure.kind == QAWS_MEASURE_CURVATURE ? 6u : 3u;
}

/* Dual jets y at the dual point (u, v) along dir: the first and second
   rates of u and v enter as the straight-path tangents plus the second
   order corrections y_u u'' + y_v v''. */
static qaws_status scdf_jet(scdf_job const* job, qaws_dual1 u, qaws_dual1 v, qaws_diff_views const* dir, qaws_dual3* y)
{
	/* jet index of the u and v derivatives of entry k */
	static unsigned char const du[6] = { 1, 3, 4, 6, 7, 8 };
	static unsigned char const dv[6] = { 2, 4, 5, 7, 8, 9 };
	qaws_surface_jet p, tg, tt;
	qaws_scalar uu = u.v, vv = v.v, ut = u.t, vt = v.t;
	unsigned int k, n = scdf_y_count(job);
	qaws_status st = qaws_surface_eval_batch_tangent2(job->ctx, job->surface, &uu, &vv, &ut, &vt, 1,
		n > 3 ? QAWS_SJET_ORDER3 : QAWS_SJET_ORDER2, dir, &p, &tg, &tt);
	if (st != QAWS_STATUS_OK)
		return st;
	for (k = 0; k < n; k++)
		y[k] = qaws_dual3_make(p.d[k], tg.d[k], v3_axpy(v.tt, p.d[dv[k]], v3_axpy(u.tt, p.d[du[k]], tt.d[k])));
	return QAWS_STATUS_OK;
}

/*
 * w = rho |S_u x S_v|. The curvature measure takes rho = sqrt(floor^2 +
 * k1^2 + k2^2) with k1^2 + k2^2 = 4 H^2 - 2 K; with D = E G - F^2 and the
 * second form on the unnormalized normal n = S_u x S_v (L = S_uu . n, ...):
 *   4 H^2 = (E N - 2 F M + G L)^2 / D^3,  K = (L N - M^2) / D^2.
 */
static qaws_dual1 scdf_w(scdf_job const* job, qaws_dual3 const* y)
{
	qaws_dual3 n = qaws_dual3_cross(y[1], y[2]);
	qaws_dual1 a = qaws_dual3_length(n);
	if (job->measure.kind == QAWS_MEASURE_DENSITY)
		a = qaws_dual1_mul(density_dual(&job->measure, y[0]), a);
	else if (job->measure.kind == QAWS_MEASURE_CURVATURE)
	{
		qaws_dual1 E = qaws_dual3_dot(y[1], y[1]), F = qaws_dual3_dot(y[1], y[2]), G = qaws_dual3_dot(y[2], y[2]);
		qaws_dual1 L = qaws_dual3_dot(y[3], n), M = qaws_dual3_dot(y[4], n), N = qaws_dual3_dot(y[5], n);
		qaws_dual1 D = qaws_dual1_sub(qaws_dual1_mul(E, G), qaws_dual1_mul(F, F));
		qaws_dual1 D2 = qaws_dual1_mul(D, D);
		qaws_dual1 h = qaws_dual1_add(qaws_dual1_sub(qaws_dual1_mul(E, N), dual_axpy(2, qaws_dual1_mul(F, M), 0)),
			qaws_dual1_mul(G, L));
		qaws_dual1 k = qaws_dual1_sub(qaws_dual1_mul(L, N), qaws_dual1_mul(M, M));
		qaws_dual1 c2 = qaws_dual1_sub(qaws_dual1_div(qaws_dual1_mul(h, h), qaws_dual1_mul(D2, D)), dual_axpy(2, qaws_dual1_div(k, D2), 0));
		double f = job->measure.curvature_floor;
		a = qaws_dual1_mul(qaws_dual1_sqrt(dual_axpy(1, c2, f * f)), a);
	}
	return a;
}

/* B(b; u) = int_{v0}^{b} w(u, .): the full v cells, then the partial cell
   whose nodes and weights move with the dual upper limit b. */
static qaws_status scdf_B(scdf_job const* job, qaws_dual1 u, qaws_dual1 b, qaws_diff_views const* dir, qaws_dual1* out)
{
	unsigned int c, q, k = (unsigned int)((b.v - job->v0) / job->hv);
	qaws_dual1 sum = dual_c(0), len, v0, hv, vk;
	qaws_dual3 y[SCDF_Y];
	qaws_status st;
	if (k >= job->cells)
		k = job->cells - 1;
	scdf_grid(job, 1, dir != NULL, &v0, &hv);
	for (c = 0; c < k; c++)
		for (q = 0; q < job->n; q++)
		{
			st = scdf_jet(job, u, dual_lin(c + 0.5 + 0.5 * job->x[q], hv, v0), dir, y);
			if (st != QAWS_STATUS_OK)
				return st;
			sum = qaws_dual1_add(sum, qaws_dual1_mul(dual_axpy(0.5 * job->w[q], hv, 0), scdf_w(job, y)));
		}
	vk = dual_lin(k, hv, v0);
	len = qaws_dual1_sub(b, vk);
	for (q = 0; q < job->n; q++)
	{
		st = scdf_jet(job, u, dual_lin(0.5 * (1 + job->x[q]), len, vk), dir, y);
		if (st != QAWS_STATUS_OK)
			return st;
		sum = qaws_dual1_add(sum, qaws_dual1_mul(dual_axpy(0.5 * job->w[q], len, 0), scdf_w(job, y)));
	}
	*out = sum;
	return QAWS_STATUS_OK;
}

/* A(u) = int_{u0}^{u} B(v1; .): the full u cells (job->cumA, along the job
   direction when use_dir) and the partial cell moving with the dual u. */
static qaws_status scdf_A(scdf_job const* job, qaws_dual1 u, int use_dir, qaws_dual1* out)
{
	unsigned int q, k = (unsigned int)((u.v - job->u0) / job->hu);
	qaws_dual1 sum, len, u0, hu, uk, v1 = scdf_v1(job, use_dir);
	if (k >= job->cells)
		k = job->cells - 1;
	scdf_grid(job, 0, use_dir, &u0, &hu);
	uk = dual_lin(k, hu, u0);
	sum = use_dir ? job->cumA[k] : dual_c(job->cumA[k].v);
	len = qaws_dual1_sub(u, uk);
	for (q = 0; q < job->n; q++)
	{
		qaws_dual1 b1;
		qaws_status st = scdf_B(job, dual_lin(0.5 * (1 + job->x[q]), len, uk), v1, use_dir ? job->direction : NULL, &b1);
		if (st != QAWS_STATUS_OK)
			return st;
		sum = qaws_dual1_add(sum, qaws_dual1_mul(dual_axpy(0.5 * job->w[q], len, 0), b1));
	}
	*out = sum;
	return QAWS_STATUS_OK;
}

/* F(u) = A(u) - xi_u A_total */
static qaws_status scdf_F(scdf_job const* job, qaws_dual1 u, qaws_dual1 xi, int use_dir, qaws_dual1* out)
{
	qaws_dual1 a, total = use_dir ? job->cumA[job->cells] : dual_c(job->cumA[job->cells].v);
	qaws_status st = scdf_A(job, u, use_dir, &a);
	if (st != QAWS_STATUS_OK)
		return st;
	*out = qaws_dual1_sub(a, qaws_dual1_mul(xi, total));
	return QAWS_STATUS_OK;
}

/* G(v; u) = B(v; u) - xi_v B(v1; u) */
static qaws_status scdf_G(scdf_job const* job, qaws_dual1 u, qaws_dual1 v, qaws_dual1 xi, int use_dir, qaws_dual1* out)
{
	qaws_dual1 b = dual_c(0), b1 = dual_c(0);
	qaws_diff_views const* dir = use_dir ? job->direction : NULL;
	qaws_status st = scdf_B(job, u, v, dir, &b);
	if (st == QAWS_STATUS_OK)
		st = scdf_B(job, u, scdf_v1(job, use_dir), dir, &b1);
	if (st != QAWS_STATUS_OK)
		return st;
	*out = qaws_dual1_sub(b, qaws_dual1_mul(xi, b1));
	return QAWS_STATUS_OK;
}

static void scdf_job_free(scdf_job* job)
{
	qaws_internal_dealloc(NULL, job->cumA);
	job->cumA = NULL;
}

static qaws_status scdf_job_init(scdf_job* job, qaws_diff_context const* ctx, qaws_surface const* surface,
	qaws_sample_measure_desc const* measure, unsigned int cells, unsigned int quadrature, qaws_diff_views const* direction,
	qaws_diff_views const* other)
{
	unsigned int c, q;
	memset(job, 0, sizeof(*job));
	if (!surface)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (measure)
	{
		job->measure = *measure;
		if ((measure->kind == QAWS_MEASURE_CURVATURE && !(measure->curvature_floor > 0)) ||
		    (measure->kind == QAWS_MEASURE_DENSITY && !measure->density) ||
		    (unsigned int)measure->kind > (unsigned int)QAWS_MEASURE_DENSITY)
			return QAWS_STATUS_INVALID_ARGUMENT;
	}
	if (!(qaws_surface_get_diff_capabilities(surface) & QAWS_CAP_TANGENT2))
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	/* knot sinks are filled by the callers (forward passes per knot) */
	if (other && (qaws_diff_views_find(other, QAWS_FIELD_U_KNOTS) || qaws_diff_views_find(other, QAWS_FIELD_V_KNOTS)))
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	if (direction)
	{
		/* rates of the domain end knots: u_knots[p] and u_knots[n_u] */
		qaws_field_view const* uv = qaws_diff_views_find(direction, QAWS_FIELD_U_KNOTS);
		qaws_field_view const* vv = qaws_diff_views_find(direction, QAWS_FIELD_V_KNOTS);
		qaws_field_desc fields[8];
		unsigned int nf = 0, f, nu = 0, nv = 0;
		if ((uv || vv) && qaws_surface_describe_fields(surface, fields, 8, &nf) == QAWS_STATUS_OK)
		{
			qaws_scalar r = QAWS_ZERO;
			for (f = 0; f < nf && f < 8; f++)
			{
				if (fields[f].field == QAWS_FIELD_U_KNOTS) nu = fields[f].count;
				if (fields[f].field == QAWS_FIELD_V_KNOTS) nv = fields[f].count;
			}
			if (uv && nu > surface->u_degree + 1)
			{
				qaws_internal_view_read(uv, surface->u_degree, 1, &r); job->du0 = r;
				qaws_internal_view_read(uv, nu - surface->u_degree - 1, 1, &r); job->du1 = r;
			}
			if (vv && nv > surface->v_degree + 1)
			{
				qaws_internal_view_read(vv, surface->v_degree, 1, &r); job->dv0 = r;
				qaws_internal_view_read(vv, nv - surface->v_degree - 1, 1, &r); job->dv1 = r;
			}
		}
	}
	job->ctx = ctx;
	job->surface = surface;
	job->direction = direction;
	job->cells = cells == 0 ? 6 : cells;
	job->n = quadrature == 0 ? 8 : (quadrature < 2 ? 2 : (quadrature > 8 ? 8 : quadrature));
	sampling_gauss_rule(job->n, job->x, job->w);
	job->u0 = surface->u_range.min_value;
	job->u1 = surface->u_range.max_value;
	job->v0 = surface->v_range.min_value;
	job->v1 = surface->v_range.max_value;
	job->hu = (job->u1 - job->u0) / job->cells;
	job->hv = (job->v1 - job->v0) / job->cells;
	job->cumA = (qaws_dual1*)qaws_internal_alloc(NULL, (unsigned long)(sizeof(qaws_dual1) * (job->cells + 1)));
	if (!job->cumA)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	job->cumA[0] = dual_c(0);
	{
		qaws_dual1 u0, hu, v1 = scdf_v1(job, direction != NULL);
		scdf_grid(job, 0, direction != NULL, &u0, &hu);
		for (c = 0; c < job->cells; c++)
		{
			qaws_dual1 sum = job->cumA[c];
			for (q = 0; q < job->n; q++)
			{
				qaws_dual1 b1;
				qaws_status st = scdf_B(job, dual_lin(c + 0.5 + 0.5 * job->x[q], hu, u0), v1, direction, &b1);
				if (st != QAWS_STATUS_OK)
				{
					scdf_job_free(job);
					return st;
				}
				sum = qaws_dual1_add(sum, qaws_dual1_mul(dual_axpy(0.5 * job->w[q], hu, 0), b1));
			}
			job->cumA[c + 1] = sum;
		}
	}
	return QAWS_STATUS_OK;
}

/* Safeguarded Newton for an increasing f on [lo, hi]; fn returns f and its slope. */
typedef qaws_status (*scdf_fn)(scdf_job const* job, void const* user, double x, double* f, double* slope);

static qaws_status scdf_newton(scdf_job const* job, void const* user, scdf_fn fn, double lo, double hi, double scale, double* out)
{
	double x = 0.5 * (lo + hi), width = hi - lo;
	unsigned int it;
	for (it = 0; it < 60; it++)
	{
		double f, s;
		qaws_status st = fn(job, user, x, &f, &s);
		if (st != QAWS_STATUS_OK)
			return st;
		if (f > 0) hi = x; else lo = x;
		if (fabs(f) <= 1e-15 * (scale > 1 ? scale : 1))
			break;
		x = s > 0 ? x - f / s : 0.5 * (lo + hi);
		if (!(x > lo && x < hi))
			x = 0.5 * (lo + hi);
		if (hi - lo <= 1e-15 * width)
			break;
	}
	*out = x;
	return QAWS_STATUS_OK;
}

typedef struct scdf_target
{
	double xi, u;
} scdf_target;

static qaws_status scdf_fu(scdf_job const* job, void const* user, double u, double* f, double* slope)
{
	scdf_target const* tg = (scdf_target const*)user;
	qaws_dual1 F;
	qaws_status st = scdf_F(job, qaws_dual1_make((qaws_scalar)u, 1, 0), dual_c(tg->xi), 0, &F);
	*f = F.v;
	*slope = F.t;
	return st;
}

static qaws_status scdf_gv(scdf_job const* job, void const* user, double v, double* f, double* slope)
{
	scdf_target const* tg = (scdf_target const*)user;
	qaws_dual1 G;
	qaws_status st = scdf_G(job, dual_c(tg->u), qaws_dual1_make((qaws_scalar)v, 1, 0), dual_c(tg->xi), 0, &G);
	*f = G.v;
	*slope = G.t;
	return st;
}

/* Primal (u, v) of one point of the unit square. */
static qaws_status scdf_solve(scdf_job const* job, double xu, double xv, double* out_u, double* out_v)
{
	scdf_target tg;
	double total = job->cumA[job->cells].v, target, lo, b1;
	unsigned int k = 0;
	qaws_dual1 B1;
	qaws_status st;
	if (xu < 0) xu = 0;
	if (xu > 1) xu = 1;
	if (xv < 0) xv = 0;
	if (xv > 1) xv = 1;
	target = xu * total;
	while (k + 1 < job->cells && job->cumA[k + 1].v <= target)
		k++;
	lo = job->u0 + job->hu * k;
	tg.xi = xu;
	tg.u = 0;
	st = scdf_newton(job, &tg, scdf_fu, lo, lo + job->hu, total, out_u);
	if (st != QAWS_STATUS_OK)
		return st;
	/* the v cell holding xi_v B(v1; u) */
	st = scdf_B(job, dual_c(*out_u), dual_c(job->v1), NULL, &B1);
	if (st != QAWS_STATUS_OK)
		return st;
	b1 = B1.v;
	k = 0;
	while (k + 1 < job->cells)
	{
		qaws_dual1 Bk;
		st = scdf_B(job, dual_c(*out_u), dual_c(job->v0 + job->hv * (k + 1)), NULL, &Bk);
		if (st != QAWS_STATUS_OK)
			return st;
		if (Bk.v > xv * b1)
			break;
		k++;
	}
	tg.xi = xv;
	tg.u = *out_u;
	lo = job->v0 + job->hv * k;
	return scdf_newton(job, &tg, scdf_gv, lo, lo + job->hv, b1, out_v);
}

/* First and second order rates of one sample along the job direction and
   the xi rates (xd may be NULL); out2 may be NULL. Returns
   QAWS_STATUS_NUMERICAL_FAILURE where the measure vanishes. */
static qaws_status scdf_rates(scdf_job const* job, double u, double v, double xu, double xv, double const* xd,
	qaws_surface_cdf_sample* out1, qaws_surface_cdf_sample* out2)
{
	qaws_dual1 F = dual_c(0), G = dual_c(0), ud, vd;
	qaws_dual1 xiu = qaws_dual1_make((qaws_scalar)xu, (qaws_scalar)(xd ? xd[0] : 0), 0);
	qaws_dual1 xiv = qaws_dual1_make((qaws_scalar)xv, (qaws_scalar)(xd ? xd[1] : 0), 0);
	double Fu, Gv, u1, u2 = 0, v1, v2 = 0;
	qaws_dual3 y[SCDF_Y];
	qaws_status st;
	/* slopes of the discrete equations */
	st = scdf_F(job, qaws_dual1_make((qaws_scalar)u, 1, 0), dual_c(xu), 0, &F);
	if (st == QAWS_STATUS_OK)
		st = scdf_G(job, dual_c(u), qaws_dual1_make((qaws_scalar)v, 1, 0), dual_c(xv), 0, &G);
	if (st != QAWS_STATUS_OK)
		return st;
	Fu = F.t;
	Gv = G.t;
	if (!(Fu > QAWS_EPSILON) || !(Gv > QAWS_EPSILON))
		return QAWS_STATUS_NUMERICAL_FAILURE;
	/* u: one dual Newton step per order */
	st = scdf_F(job, dual_c(u), xiu, 1, &F);
	if (st != QAWS_STATUS_OK)
		return st;
	u1 = -F.t / Fu;
	if (out2)
	{
		st = scdf_F(job, qaws_dual1_make((qaws_scalar)u, (qaws_scalar)u1, 0), xiu, 1, &F);
		if (st != QAWS_STATUS_OK)
			return st;
		u2 = -F.tt / Fu;
	}
	ud = qaws_dual1_make((qaws_scalar)u, (qaws_scalar)u1, (qaws_scalar)u2);
	/* v, with u moving */
	st = scdf_G(job, ud, dual_c(v), xiv, 1, &G);
	if (st != QAWS_STATUS_OK)
		return st;
	v1 = -G.t / Gv;
	if (out2)
	{
		st = scdf_G(job, ud, qaws_dual1_make((qaws_scalar)v, (qaws_scalar)v1, 0), xiv, 1, &G);
		if (st != QAWS_STATUS_OK)
			return st;
		v2 = -G.tt / Gv;
	}
	vd = qaws_dual1_make((qaws_scalar)v, (qaws_scalar)v1, (qaws_scalar)v2);
	st = scdf_jet(job, ud, vd, job->direction, y);
	if (st != QAWS_STATUS_OK)
		return st;
	out1->u = (qaws_scalar)u1;
	out1->v = (qaws_scalar)v1;
	out1->position = y[0].t;
	if (out2)
	{
		out2->u = (qaws_scalar)u2;
		out2->v = (qaws_scalar)v2;
		out2->position = y[0].tt;
	}
	return QAWS_STATUS_OK;
}

qaws_status qaws_surface_cdf_sample_tangent(
	qaws_diff_context const* ctx,
	qaws_surface const* surface,
	qaws_sample_measure_desc const* measure,
	qaws_scalar const* xi,
	qaws_scalar const* xi_tangent,
	unsigned int count,
	unsigned int cells,
	unsigned int quadrature,
	qaws_diff_views const* param_tangent,
	qaws_surface_cdf_sample* out_value,
	qaws_surface_cdf_sample* out_tangent,
	qaws_surface_cdf_sample* out_tangent2,
	qaws_scalar* out_total)
{
	scdf_job job;
	unsigned int i;
	qaws_status st;
	if (!xi && count)
		return QAWS_STATUS_INVALID_ARGUMENT;
	st = scdf_job_init(&job, ctx, surface, measure, cells, quadrature, param_tangent, NULL);
	if (st != QAWS_STATUS_OK)
		return st;
	if (out_total)
		*out_total = job.cumA[job.cells].v;
	for (i = 0; i < count && st == QAWS_STATUS_OK; i++)
	{
		double u, v, xd[2];
		qaws_surface_cdf_sample t1, t2;
		st = scdf_solve(&job, xi[2 * i], xi[2 * i + 1], &u, &v);
		if (st != QAWS_STATUS_OK)
			break;
		if (out_value)
		{
			qaws_surface_jet p;
			st = qaws_surface_eval_jet(surface, (qaws_scalar)u, (qaws_scalar)v, QAWS_SJET_P, &p);
			if (st != QAWS_STATUS_OK)
				break;
			out_value[i].u = (qaws_scalar)u;
			out_value[i].v = (qaws_scalar)v;
			out_value[i].position = p.d[0];
		}
		if (!out_tangent && !out_tangent2)
			continue;
		if (xi_tangent)
		{
			xd[0] = xi_tangent[2 * i];
			xd[1] = xi_tangent[2 * i + 1];
		}
		memset(&t1, 0, sizeof(t1));
		memset(&t2, 0, sizeof(t2));
		st = scdf_rates(&job, u, v, xi[2 * i], xi[2 * i + 1], xi_tangent ? xd : NULL, &t1, out_tangent2 ? &t2 : NULL);
		if (st == QAWS_STATUS_NUMERICAL_FAILURE)
		{
			qaws_internal_diff_report_note(ctx, QAWS_DIFF_SMOOTH, QAWS_DIFF_ILL_CONDITIONED, 0, i);
			memset(&t1, 0, sizeof(t1));
			memset(&t2, 0, sizeof(t2));
			st = QAWS_STATUS_OK;
		}
		else
			qaws_internal_diff_report_note(ctx, QAWS_DIFF_SMOOTH, QAWS_DIFF_VALID, 0, i);
		if (out_tangent) out_tangent[i] = t1;
		if (out_tangent2) out_tangent2[i] = t2;
	}
	scdf_job_free(&job);
	return st;
}

/*
 * Backward pass. Every coefficient is a dual number whose rate is its
 * derivative along the job direction; the plain adjoint reads the values
 * only. With `hvp` set, the pass is differentiated instead (forward over
 * reverse, for families linear in their fields): a node at the moving point
 * (u, v) pulling back C dw/dy adds
 *   J^T (C' dw/dy + C H_w ydot) + (J_u^T u' + J_v^T v') C dw/dy,
 * the last term being the same pullback shifted to the next jet channels.
 */
typedef struct scdf_pass
{
	scdf_job const* job;
	qaws_diff_views* sink;
	int hvp;
} scdf_pass;

/* jet index of the u and v derivatives of entry k */
static unsigned char const g_scdf_du[6] = { 1, 3, 4, 6, 7, 8 };
static unsigned char const g_scdf_dv[6] = { 2, 4, 5, 7, 8, 9 };

/* Channels of a pullback of the entries (S, S_u, S_v[, S_uu, S_uv, S_vv]),
   one order more when shifted for the HVP. */
static unsigned int scdf_channels(scdf_job const* job, int shifted)
{
	unsigned int order = (scdf_y_count(job) > 3 ? 2u : 1u) + (shifted ? 1u : 0u);
	return order >= 3 ? QAWS_SJET_ORDER3 : (order == 2 ? QAWS_SJET_ORDER2 : QAWS_SJET_ORDER1);
}

/* zdot^T H_w zdot along the jet values y (z may be NULL: zero). */
static double scdf_w_quad(scdf_job const* job, qaws_dual3 const* y, qaws_vec3 const* z, unsigned int i, unsigned int comp)
{
	qaws_dual3 s[SCDF_Y];
	unsigned int k, n = scdf_y_count(job);
	for (k = 0; k < n; k++)
		s[k] = qaws_dual3_make(y[k].v, z ? z[k] : qaws_v3_zero(), qaws_v3_zero());
	if (i < n)
		v3_set(&s[i].t, comp, v3_get(&s[i].t, comp) + QAWS_ONE);
	return scdf_w(job, s).tt;
}

/* Pullback of C dw/dy at the node (u, v); C, u, v carry their rates. */
static qaws_status scdf_node_pullback(scdf_pass const* ps, qaws_dual1 u, qaws_dual1 v, qaws_dual1 C)
{
	scdf_job const* job = ps->job;
	qaws_dual3 y[SCDF_Y], s[SCDF_Y];
	qaws_vec3 g[SCDF_Y], yd[SCDF_Y];
	qaws_surface_jet bar;
	unsigned int i, k, comp, n = scdf_y_count(job), ch;
	double q_d = 0;
	qaws_status st;
	if (C.v == 0 && (!ps->hvp || C.t == 0))
		return QAWS_STATUS_OK;
	if (!ps->hvp)
	{
		u = dual_c(u.v);
		v = dual_c(v.v);
	}
	st = scdf_jet(job, u, v, ps->hvp ? job->direction : NULL, y);
	if (st != QAWS_STATUS_OK)
		return st;
	for (k = 0; k < n; k++)
		yd[k] = y[k].t;
	if (ps->hvp)
		q_d = scdf_w_quad(job, y, yd, SCDF_Y, 0);
	memset(&bar, 0, sizeof(bar));
	for (i = 0; i < n; i++)
	{
		if (i == 0 && job->measure.kind != QAWS_MEASURE_DENSITY)
			continue;
		for (comp = 0; comp < 3; comp++)
		{
			double gi;
			for (k = 0; k < n; k++)
				s[k] = qaws_dual3_const(y[k].v);
			v3_set(&s[i].t, comp, QAWS_ONE);
			gi = scdf_w(job, s).t;
			v3_set(&g[i], comp, (qaws_scalar)gi);
			if (!ps->hvp)
				v3_set(&bar.d[i], comp, (qaws_scalar)(C.v * gi));
			else
			{
				/* (H_w ydot)_j by polarization */
				double hj = 0.5 * (scdf_w_quad(job, y, yd, i, comp) - q_d - scdf_w_quad(job, y, NULL, i, comp));
				v3_set(&bar.d[i], comp, (qaws_scalar)(C.t * gi + C.v * hj));
			}
		}
	}
	ch = scdf_channels(job, 0);
	if (ps->hvp)
	{
		for (i = 0; i < n; i++)
		{
			if (i == 0 && job->measure.kind != QAWS_MEASURE_DENSITY)
				continue;
			bar.d[g_scdf_du[i]] = v3_axpy(C.v * u.t, g[i], bar.d[g_scdf_du[i]]);
			bar.d[g_scdf_dv[i]] = v3_axpy(C.v * v.t, g[i], bar.d[g_scdf_dv[i]]);
		}
		ch = scdf_channels(job, 1);
	}
	bar.channels = ch;
	return qaws_surface_eval_adjoint(job->ctx, job->surface, (qaws_scalar)u.v, (qaws_scalar)v.v, ch, &bar, ps->sink, NULL, NULL);
}

/* Pullback c grad B(b; u): the full v cells, then the partial cell whose
   nodes and weights move with b. */
static qaws_status scdf_pullback_B(scdf_pass const* ps, qaws_dual1 u, qaws_dual1 b, qaws_dual1 c)
{
	scdf_job const* job = ps->job;
	unsigned int cc, q, k = (unsigned int)((b.v - job->v0) / job->hv);
	qaws_dual1 len;
	double vk;
	qaws_status st = QAWS_STATUS_OK;
	if (k >= job->cells)
		k = job->cells - 1;
	for (cc = 0; cc < k && st == QAWS_STATUS_OK; cc++)
		for (q = 0; q < job->n && st == QAWS_STATUS_OK; q++)
			st = scdf_node_pullback(ps, u, dual_c(job->v0 + job->hv * (cc + 0.5 + 0.5 * job->x[q])),
				dual_axpy(0.5 * job->hv * job->w[q], c, 0));
	vk = job->v0 + job->hv * k;
	len = dual_axpy(1, b, -vk);
	for (q = 0; q < job->n && st == QAWS_STATUS_OK; q++)
		st = scdf_node_pullback(ps, u, dual_axpy(0.5 * (1 + job->x[q]), len, vk),
			qaws_dual1_mul(c, dual_axpy(0.5 * job->w[q], len, 0)));
	return st;
}

/* Slope of an equation along one unknown (`which` 0: u, 1: v) and, for the
   HVP, its rate along (u', v', direction) by polarization of second order
   passes: (q(e + x) - q(x) - q(e)) / 2 + q-free terms, q the tt part. */
static qaws_status scdf_slope(scdf_job const* job, int eq, int which, double u, double v, double ud, double vd, double xi,
	int hvp, qaws_dual1* out)
{
	qaws_dual1 r = dual_c(0), a = dual_c(0), b = dual_c(0);
	double eu = which == 0 ? 1 : 0, ev = which == 1 ? 1 : 0;
	qaws_status st;
#define SCDF_EVAL(uu, vv, dir, res) \
	(eq == 0 ? scdf_F(job, qaws_dual1_make((qaws_scalar)u, (qaws_scalar)(uu), 0), dual_c(xi), dir, res) \
	         : scdf_G(job, qaws_dual1_make((qaws_scalar)u, (qaws_scalar)(uu), 0), qaws_dual1_make((qaws_scalar)v, (qaws_scalar)(vv), 0), \
	                  dual_c(xi), dir, res))
	st = SCDF_EVAL(eu, ev, 0, &r);
	if (st != QAWS_STATUS_OK || !hvp)
	{
		*out = dual_c(r.t);
		return st;
	}
	st = SCDF_EVAL(ud + eu, vd + ev, 1, &a);
	if (st == QAWS_STATUS_OK)
		st = SCDF_EVAL(ud, vd, 1, &b);
#undef SCDF_EVAL
	/* a.tt - b.tt = q(e) + 2 e^T H x, so the rate e^T H x = (a.tt - b.tt - r.tt) / 2 */
	*out = qaws_dual1_make(r.t, (qaws_scalar)(0.5 * ((double)a.tt - b.tt - r.tt)), 0);
	return st;
}

/*
 * Adjoint: the point pulls back J^T p_bar and adds p_bar . S_u, p_bar . S_v
 * to the (u, v) adjoints; the conditional equation G(v; u) = 0 then gives
 * mu = v_bar / G_v (theta: -mu grad G, u: -mu G_u, xi_v: mu B(v1; u)); the
 * marginal equation F(u) = 0 gives lambda = u_bar / F_u (theta: -lambda
 * grad F, xi_u: lambda A_total). The full u cells of grad A of the whole
 * batch are pulled back once with summed weights.
 */
static qaws_status scdf_backward(scdf_pass const* ps, qaws_scalar const* xi, unsigned int count,
	qaws_surface_cdf_sample const* adjoint, qaws_scalar* xi_adjoint)
{
	scdf_job const* job = ps->job;
	qaws_dual1* cellw;
	double total = job->cumA[job->cells].v;
	unsigned int i, c, q;
	qaws_status st = QAWS_STATUS_OK;
	cellw = (qaws_dual1*)qaws_internal_alloc(NULL, (unsigned long)(sizeof(qaws_dual1) * (job->cells + 1)));
	if (!cellw)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (c = 0; c <= job->cells; c++)
		cellw[c] = dual_c(0);
	for (i = 0; i < count && st == QAWS_STATUS_OK; i++)
	{
		double u, v, xu = xi[2 * i], xv = xi[2 * i + 1];
		qaws_dual1 ud, vd, ub, vb, Fu = dual_c(0), Gv = dual_c(0), Gu = dual_c(0), mu, lambda;
		qaws_surface_cdf_sample rate;
		unsigned int ku;
		st = scdf_solve(job, xu, xv, &u, &v);
		if (st != QAWS_STATUS_OK)
			break;
		memset(&rate, 0, sizeof(rate));
		if (ps->hvp)
		{
			st = scdf_rates(job, u, v, xu, xv, NULL, &rate, NULL);
			if (st == QAWS_STATUS_NUMERICAL_FAILURE)
			{
				st = QAWS_STATUS_OK;
				continue;
			}
			if (st != QAWS_STATUS_OK)
				break;
		}
		ud = qaws_dual1_make((qaws_scalar)u, rate.u, 0);
		vd = qaws_dual1_make((qaws_scalar)v, rate.v, 0);
		/* the position term */
		{
			qaws_surface_jet bar;
			qaws_scalar ua = QAWS_ZERO, va = QAWS_ZERO;
			memset(&bar, 0, sizeof(bar));
			if (!ps->hvp)
			{
				bar.d[0] = adjoint[i].position;
				bar.channels = QAWS_SJET_P;
				st = qaws_surface_eval_adjoint(job->ctx, job->surface, (qaws_scalar)u, (qaws_scalar)v, QAWS_SJET_P, &bar, ps->sink, &ua, &va);
				ub = dual_c((double)adjoint[i].u + ua);
				vb = dual_c((double)adjoint[i].v + va);
			}
			else
			{
				qaws_dual3 y[SCDF_Y];
				bar.d[1] = qaws_v3_scale(adjoint[i].position, rate.u);
				bar.d[2] = qaws_v3_scale(adjoint[i].position, rate.v);
				bar.channels = QAWS_SJET_ORDER1;
				st = qaws_surface_eval_adjoint(job->ctx, job->surface, (qaws_scalar)u, (qaws_scalar)v, QAWS_SJET_ORDER1, &bar, ps->sink,
					NULL, NULL);
				if (st == QAWS_STATUS_OK)
					st = scdf_jet(job, ud, vd, job->direction, y);
				ub = qaws_dual1_make((qaws_scalar)((double)adjoint[i].u + v3_dot(adjoint[i].position, y[1].v)),
					(qaws_scalar)v3_dot(adjoint[i].position, y[1].t), 0);
				vb = qaws_dual1_make((qaws_scalar)((double)adjoint[i].v + v3_dot(adjoint[i].position, y[2].v)),
					(qaws_scalar)v3_dot(adjoint[i].position, y[2].t), 0);
			}
			if (st != QAWS_STATUS_OK)
				break;
		}
		st = scdf_slope(job, 0, 0, u, v, rate.u, rate.v, xu, ps->hvp, &Fu);
		if (st == QAWS_STATUS_OK)
			st = scdf_slope(job, 1, 1, u, v, rate.u, rate.v, xv, ps->hvp, &Gv);
		if (st == QAWS_STATUS_OK)
			st = scdf_slope(job, 1, 0, u, v, rate.u, rate.v, xv, ps->hvp, &Gu);
		if (st != QAWS_STATUS_OK)
			break;
		if (!(Fu.v > QAWS_EPSILON) || !(Gv.v > QAWS_EPSILON))
		{
			qaws_internal_diff_report_note(job->ctx, QAWS_DIFF_SMOOTH, QAWS_DIFF_ILL_CONDITIONED, 0, i);
			continue;
		}
		qaws_internal_diff_report_note(job->ctx, QAWS_DIFF_SMOOTH, QAWS_DIFF_VALID, 0, i);
		/* conditional equation */
		mu = qaws_dual1_div(vb, Gv);
		ub = qaws_dual1_sub(ub, qaws_dual1_mul(mu, Gu));
		if (xi_adjoint)
		{
			qaws_dual1 B1;
			st = scdf_B(job, dual_c(u), dual_c(job->v1), NULL, &B1);
			if (st != QAWS_STATUS_OK)
				break;
			xi_adjoint[2 * i + 1] += (qaws_scalar)(mu.v * B1.v);
		}
		if (ps->sink)
		{
			st = scdf_pullback_B(ps, ud, vd, dual_axpy(-1, mu, 0));
			if (st == QAWS_STATUS_OK)
				st = scdf_pullback_B(ps, ud, dual_c(job->v1), dual_axpy(xv, mu, 0));
			if (st != QAWS_STATUS_OK)
				break;
		}
		/* marginal equation */
		lambda = qaws_dual1_div(ub, Fu);
		if (xi_adjoint)
			xi_adjoint[2 * i] += (qaws_scalar)(lambda.v * total);
		if (!ps->sink)
			continue;
		ku = (unsigned int)((u - job->u0) / job->hu);
		if (ku >= job->cells)
			ku = job->cells - 1;
		for (c = 0; c < job->cells; c++)
			cellw[c] = qaws_dual1_add(cellw[c], dual_axpy(xu - (c < ku ? 1.0 : 0.0), lambda, 0));
		{
			double uk = job->u0 + job->hu * ku;
			qaws_dual1 len = dual_axpy(1, ud, -uk);
			for (q = 0; q < job->n && st == QAWS_STATUS_OK; q++)
				st = scdf_pullback_B(ps, dual_axpy(0.5 * (1 + job->x[q]), len, uk), dual_c(job->v1),
					qaws_dual1_mul(dual_axpy(-1, lambda, 0), dual_axpy(0.5 * job->w[q], len, 0)));
		}
	}
	for (c = 0; c < job->cells && st == QAWS_STATUS_OK && ps->sink; c++)
		for (q = 0; q < job->n && st == QAWS_STATUS_OK; q++)
			st = scdf_pullback_B(ps, dual_c(job->u0 + job->hu * (c + 0.5 + 0.5 * job->x[q])), dual_c(job->v1),
				dual_axpy(0.5 * job->hu * job->w[q], cellw[c], 0));
	qaws_internal_dealloc(NULL, cellw);
	return st;
}

qaws_status qaws_surface_cdf_sample_adjoint(
	qaws_diff_context const* ctx,
	qaws_surface const* surface,
	qaws_sample_measure_desc const* measure,
	qaws_scalar const* xi,
	unsigned int count,
	unsigned int cells,
	unsigned int quadrature,
	qaws_surface_cdf_sample const* adjoint,
	qaws_diff_views* param_adjoint,
	qaws_scalar* xi_adjoint)
{
	scdf_job job;
	scdf_pass ps;
	qaws_diff_views rest;
	qaws_field_view kept[16];
	unsigned int f;
	qaws_status st;
	if ((!xi || !adjoint) && count)
		return QAWS_STATUS_INVALID_ARGUMENT;
	/* the reverse pass takes every field but the knots */
	if (param_adjoint)
	{
		rest = *param_adjoint;
		rest.fields = kept;
		rest.field_count = 0;
		for (f = 0; f < param_adjoint->field_count && rest.field_count < 16; f++)
			if (param_adjoint->fields[f].field != QAWS_FIELD_U_KNOTS && param_adjoint->fields[f].field != QAWS_FIELD_V_KNOTS)
				kept[rest.field_count++] = param_adjoint->fields[f];
	}
	st = scdf_job_init(&job, ctx, surface, measure, cells, quadrature, NULL, param_adjoint ? &rest : NULL);
	if (st != QAWS_STATUS_OK)
		return st;
	ps.job = &job;
	ps.sink = param_adjoint ? &rest : NULL;
	ps.hvp = 0;
	st = scdf_backward(&ps, xi, count, adjoint, xi_adjoint);
	scdf_job_free(&job);
	/*
	 * Knots move the cell grid and every quadrature node with the domain
	 * ends: their adjoints are exact forward passes, one per knot entry
	 * (sum adjoint . sample' along e_j), knots being few.
	 */
	for (f = 0; param_adjoint && f < param_adjoint->field_count && st == QAWS_STATUS_OK; f++)
	{
		qaws_field_view* kv = &param_adjoint->fields[f];
		qaws_surface_cdf_sample* t1;
		qaws_scalar* unit;
		qaws_field_view dv;
		qaws_diff_views dir;
		unsigned int e, i;
		if (kv->field != QAWS_FIELD_U_KNOTS && kv->field != QAWS_FIELD_V_KNOTS)
			continue;
		t1 = (qaws_surface_cdf_sample*)qaws_internal_alloc(NULL, (unsigned long)(sizeof(qaws_surface_cdf_sample) * (count + 1)));
		unit = (qaws_scalar*)qaws_internal_alloc(NULL, (unsigned long)(sizeof(qaws_scalar) * (kv->count + 1)));
		if (!t1 || !unit)
		{
			qaws_internal_dealloc(NULL, t1);
			qaws_internal_dealloc(NULL, unit);
			return QAWS_STATUS_ALLOCATION_FAILURE;
		}
		memset(unit, 0, sizeof(qaws_scalar) * kv->count);
		dv = qaws_field_view_make(kv->field, unit, kv->count, 1);
		dir.fields = &dv;
		dir.field_count = 1;
		dir.children = NULL;
		dir.child_count = 0;
		for (e = 0; e < kv->count && st == QAWS_STATUS_OK; e++)
		{
			double g = 0;
			qaws_scalar gs;
			if (!qaws_internal_view_element_active(kv, e))
				continue;
			unit[e] = QAWS_ONE;
			st = qaws_surface_cdf_sample_tangent(ctx, surface, measure, xi, NULL, count, cells, quadrature, &dir, NULL, t1, NULL, NULL);
			unit[e] = QAWS_ZERO;
			for (i = 0; i < count && st == QAWS_STATUS_OK; i++)
				g += (double)adjoint[i].u * t1[i].u + (double)adjoint[i].v * t1[i].v + v3_dot(adjoint[i].position, t1[i].position);
			gs = (qaws_scalar)g;
			if (st == QAWS_STATUS_OK)
				qaws_internal_view_add(kv, e, 1, &gs);
		}
		qaws_internal_dealloc(NULL, t1);
		qaws_internal_dealloc(NULL, unit);
	}
	return st;
}

/* The quadratic form q(x) = sum adjoint . sample''(x) of the surface warp. */
typedef struct scdf_quadratic
{
	qaws_diff_context const* ctx;
	qaws_surface const* surface;
	qaws_sample_measure_desc const* measure;
	qaws_scalar const* xi;
	unsigned int count, cells, quadrature;
	qaws_surface_cdf_sample const* adjoint;
	qaws_surface_cdf_sample* scratch;
} scdf_quadratic;

static qaws_status scdf_second(void const* user, qaws_diff_views const* dir, double* out)
{
	scdf_quadratic const* a = (scdf_quadratic const*)user;
	unsigned int i;
	qaws_status st = qaws_surface_cdf_sample_tangent(a->ctx, a->surface, a->measure, a->xi, NULL, a->count, a->cells, a->quadrature, dir,
		NULL, NULL, a->scratch, NULL);
	*out = 0;
	if (st != QAWS_STATUS_OK)
		return st;
	for (i = 0; i < a->count; i++)
		*out += (double)a->adjoint[i].u * a->scratch[i].u + (double)a->adjoint[i].v * a->scratch[i].v +
			v3_dot(a->adjoint[i].position, a->scratch[i].position);
	return QAWS_STATUS_OK;
}

/* HVP. Families linear in their fields differentiate the backward pass
   along the direction (forward over reverse, one pass). Others polarize the
   second order forward pass (hvp_polarize). */
qaws_status qaws_surface_cdf_sample_hvp(
	qaws_diff_context const* ctx,
	qaws_surface const* surface,
	qaws_sample_measure_desc const* measure,
	qaws_scalar const* xi,
	unsigned int count,
	unsigned int cells,
	unsigned int quadrature,
	qaws_surface_cdf_sample const* adjoint,
	qaws_diff_views const* direction,
	qaws_diff_views* out_hv)
{
	scdf_quadratic a;
	qaws_status st;
	if ((!xi || !adjoint) && count)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (!surface)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (!direction || !out_hv)
		return QAWS_STATUS_OK;
	if (out_hv->child_count || direction->child_count)
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	if ((qaws_surface_get_diff_capabilities(surface) & QAWS_CAP_LINEAR) &&
	    !qaws_diff_views_find(direction, QAWS_FIELD_U_KNOTS) && !qaws_diff_views_find(direction, QAWS_FIELD_V_KNOTS) &&
	    !qaws_diff_views_find(out_hv, QAWS_FIELD_U_KNOTS) && !qaws_diff_views_find(out_hv, QAWS_FIELD_V_KNOTS))
	{
		/* forward over reverse: one differentiated backward pass */
		scdf_job job;
		scdf_pass ps;
		st = scdf_job_init(&job, ctx, surface, measure, cells, quadrature, direction, out_hv);
		if (st != QAWS_STATUS_OK)
			return st;
		ps.job = &job;
		ps.sink = out_hv;
		ps.hvp = 1;
		st = scdf_backward(&ps, xi, count, adjoint, NULL);
		scdf_job_free(&job);
		return st;
	}
	a.ctx = ctx;
	a.surface = surface;
	a.measure = measure;
	a.xi = xi;
	a.count = count;
	a.cells = cells;
	a.quadrature = quadrature;
	a.adjoint = adjoint;
	a.scratch = (qaws_surface_cdf_sample*)qaws_internal_alloc(NULL, (unsigned long)(sizeof(qaws_surface_cdf_sample) * (count + 1)));
	if (!a.scratch)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	st = hvp_polarize(&a, scdf_second, direction, out_hv);
	qaws_internal_dealloc(NULL, a.scratch);
	return st;
}
