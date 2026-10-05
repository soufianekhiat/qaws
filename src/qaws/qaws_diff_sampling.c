#include "qaws_diff_sampling.h"
#include "qaws_diff.h"
#include "internal/qaws_internal_types.h"
#include "internal/qaws_internal_curve.h"
#include "internal/qaws_internal_diff.h"
#include "core/qaws_dual_core.h"
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

#define SAMPLING_CHANNELS (QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2)

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

/* Shared state of one call. */
typedef struct arc_job
{
	qaws_diff_context const* ctx;
	qaws_curve const* curve;
	qaws_diff_views const* direction;
	unsigned int n;
	double x[8], w[8];
	unsigned int spans;
	double* cum;      /* cum[k]: length of spans [0, k) */
	double* cum1;     /* first directional derivative of cum */
	double* cum2;     /* second directional derivative of cum */
} arc_job;

/* Composite rule: every interval is split into ARC_PIECES equal pieces of
   n Gauss points, so that the length (a square root of the speed) is
   integrated to near machine precision. */
#define ARC_PIECES 4

/* Node q (0 .. n * ARC_PIECES - 1) of [a, b] and its weight. */
static void arc_node(arc_job const* job, double a, double b, unsigned int q, double* t, double* W)
{
	unsigned int piece = q / job->n, k = q % job->n;
	double h = (b - a) / ARC_PIECES, lo = a + piece * h;
	*t = lo + 0.5 * h * (1 + job->x[k]);
	*W = 0.5 * h * job->w[k];
}

/* Jets (position, D1, D2) at t moving at rate t_dot, along the job
   direction; tangent2 may be NULL. */
static qaws_status arc_jets(arc_job const* job, qaws_scalar t, qaws_scalar t_dot, qaws_diff_views const* direction,
	qaws_curve_jet_3d* p, qaws_curve_jet_3d* tg, qaws_curve_jet_3d* tt)
{
	return qaws_internal_curve_tangent_any(job->ctx, job->curve, t, t_dot, SAMPLING_CHANNELS, direction, p, tg, tt);
}

/* Length of [a, b] and its first and second directional derivatives along
   the job direction (d1, d2 may be NULL when there is none). */
static qaws_status arc_integrate(arc_job const* job, double a, double b, double* v, double* d1, double* d2)
{
	unsigned int q;
	*v = 0;
	if (d1) *d1 = 0;
	if (d2) *d2 = 0;
	if (!(b > a))
		return QAWS_STATUS_OK;
	for (q = 0; q < job->n * ARC_PIECES; q++)
	{
		double t, W;
		arc_node(job, a, b, q, &t, &W);
		qaws_curve_jet_3d p, tg, tt;
		qaws_dual1 len;
		qaws_status st = arc_jets(job, (qaws_scalar)t, QAWS_ZERO, job->direction, &p, &tg, job->direction ? &tt : NULL);
		if (st != QAWS_STATUS_OK)
			return st;
		len = qaws_dual3_length(qaws_dual3_make(p.d[1], job->direction ? tg.d[1] : qaws_v3_zero(),
			job->direction ? tt.d[1] : qaws_v3_zero()));
		*v += W * len.v;
		if (d1) *d1 += W * len.t;
		if (d2) *d2 += W * len.tt;
	}
	return QAWS_STATUS_OK;
}

static void arc_job_free(arc_job* job)
{
	qaws_internal_dealloc(NULL, job->cum);
	job->cum = NULL;
}

/* Rejects knot views (they move the span boundaries) and builds the
   cumulative span lengths. */
static qaws_status arc_job_init(arc_job* job, qaws_diff_context const* ctx, qaws_curve const* curve,
	unsigned int quadrature, qaws_diff_views const* direction, qaws_diff_views const* other)
{
	unsigned int k;
	qaws_status st;
	memset(job, 0, sizeof(*job));
	if (!curve || curve->span_count == 0 || !curve->span_boundaries)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if ((direction && qaws_diff_views_find(direction, QAWS_FIELD_KNOTS)) ||
	    (other && qaws_diff_views_find(other, QAWS_FIELD_KNOTS)))
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	job->ctx = ctx;
	job->curve = curve;
	job->direction = direction;
	job->n = quadrature == 0 ? 8 : (quadrature < 2 ? 2 : (quadrature > 8 ? 8 : quadrature));
	sampling_gauss_rule(job->n, job->x, job->w);
	job->spans = curve->span_count;
	job->cum = (double*)qaws_internal_alloc(NULL, (unsigned long)(sizeof(double) * 3 * (job->spans + 1)));
	if (!job->cum)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	job->cum1 = job->cum + job->spans + 1;
	job->cum2 = job->cum1 + job->spans + 1;
	job->cum[0] = job->cum1[0] = job->cum2[0] = 0;
	for (k = 0; k < job->spans; k++)
	{
		double v, d1, d2;
		st = arc_integrate(job, curve->span_boundaries[k], curve->span_boundaries[k + 1], &v, &d1, &d2);
		if (st != QAWS_STATUS_OK)
		{
			arc_job_free(job);
			return st;
		}
		job->cum[k + 1] = job->cum[k] + v;
		job->cum1[k + 1] = job->cum1[k] + d1;
		job->cum2[k + 1] = job->cum2[k] + d2;
	}
	return QAWS_STATUS_OK;
}

static double arc_total(arc_job const* job)
{
	return job->cum[job->spans];
}

/* Parameter and span of the arc length sigma (clamped to the curve):
   safeguarded Newton on the span holding sigma. */
static qaws_status arc_solve(arc_job const* job, double sigma, double* out_t, unsigned int* out_span)
{
	double a, b, lo, hi, t, total = arc_total(job);
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
		double v, f, s;
		qaws_status st = arc_integrate(job, a, t, &v, NULL, NULL);
		if (st != QAWS_STATUS_OK)
			return st;
		f = job->cum[k] + v - sigma;
		if (f > 0) hi = t; else lo = t;
		if (fabs(f) <= 1e-15 * (total > 1 ? total : 1))
			break;
		st = arc_jets(job, (qaws_scalar)t, QAWS_ZERO, NULL, &p, &tg, NULL);
		if (st != QAWS_STATUS_OK)
			return st;
		s = sqrt(v3_dot(p.d[1], p.d[1]));
		t = s > 0 ? t - f / s : 0.5 * (lo + hi);
		if (!(t > lo && t < hi))
			t = 0.5 * (lo + hi);
		if (hi - lo <= 1e-15 * (fabs(b - a) + 1e-300))
			break;
	}
	*out_t = t;
	*out_span = k;
	return QAWS_STATUS_OK;
}

static double target_sigma(qaws_arc_length_target const* tg, double total)
{
	return (double)tg->distance + (double)tg->fraction * total;
}

/* ================================================================== */
/*  Forward                                                           */
/* ================================================================== */

qaws_status qaws_curve_arc_length_sample_tangent(
	qaws_diff_context const* ctx,
	qaws_curve const* curve,
	qaws_arc_length_target const* targets,
	qaws_scalar const* distance_tangent,
	unsigned int count,
	unsigned int quadrature,
	qaws_diff_views const* param_tangent,
	qaws_arc_length_sample* out_value,
	qaws_arc_length_sample* out_tangent,
	qaws_arc_length_sample* out_tangent2,
	qaws_scalar* out_total_length)
{
	arc_job job;
	unsigned int i;
	qaws_status st;
	if (!targets && count)
		return QAWS_STATUS_INVALID_ARGUMENT;
	st = arc_job_init(&job, ctx, curve, quadrature, param_tangent, NULL);
	if (st != QAWS_STATUS_OK)
		return st;
	if (out_total_length)
		*out_total_length = (qaws_scalar)arc_total(&job);
	for (i = 0; i < count && st == QAWS_STATUS_OK; i++)
	{
		double total = arc_total(&job), t, s, g, A, sig1, sig2, Ld1, Ld2, v, td, tdd;
		unsigned int k;
		qaws_curve_jet_3d p, tg, tt;
		st = arc_solve(&job, target_sigma(&targets[i], total), &t, &k);
		if (st != QAWS_STATUS_OK)
			break;
		st = arc_jets(&job, (qaws_scalar)t, QAWS_ZERO, param_tangent, &p, &tg, NULL);
		if (st != QAWS_STATUS_OK)
			break;
		if (out_value)
		{
			out_value[i].t = (qaws_scalar)t;
			out_value[i].position = p.d[0];
		}
		if (!out_tangent && !out_tangent2)
			continue;
		s = sqrt(v3_dot(p.d[1], p.d[1]));
		if (!(s > QAWS_EPSILON))
		{
			qaws_internal_diff_report_note(ctx, QAWS_DIFF_SMOOTH, QAWS_DIFF_ILL_CONDITIONED, 0, i);
			if (out_tangent) memset(&out_tangent[i], 0, sizeof(out_tangent[i]));
			if (out_tangent2) memset(&out_tangent2[i], 0, sizeof(out_tangent2[i]));
			continue;
		}
		qaws_internal_diff_report_note(ctx, QAWS_DIFF_SMOOTH, QAWS_DIFF_VALID, 0, i);
		st = arc_integrate(&job, curve->span_boundaries[k], t, &v, &Ld1, &Ld2);
		if (st != QAWS_STATUS_OK)
			break;
		Ld1 += job.cum1[k];
		Ld2 += job.cum2[k];
		g = v3_dot(p.d[1], tg.d[1]) / s;
		A = v3_dot(p.d[1], p.d[2]) / s;
		sig1 = (distance_tangent ? (double)distance_tangent[i] : 0) + (double)targets[i].fraction * job.cum1[job.spans];
		sig2 = (double)targets[i].fraction * job.cum2[job.spans];
		td = (sig1 - Ld1) / s;
		tdd = (sig2 - Ld2 - 2 * g * td - A * td * td) / s;
		/* the point along the curve moving at t' (straight path in (t, theta)),
		   plus C' t'' for the curvature of the path in t */
		st = arc_jets(&job, (qaws_scalar)t, (qaws_scalar)td, param_tangent, &p, &tg, out_tangent2 ? &tt : NULL);
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
	arc_job_free(&job);
	return st;
}

/* ================================================================== */
/*  Backward                                                          */
/* ================================================================== */

/* Per-sample state shared by the adjoint and the HVP. */
typedef struct arc_sample
{
	double t, s, lambda, lambda_dot, t_dot;
	unsigned int span;
	qaws_vec3 u;          /* unit tangent C' / |C'| */
} arc_sample;

/* Adds the D1-channel pullback bar_d1 at t. */
static qaws_status arc_add_d1(arc_job const* job, double t, qaws_vec3 bar_d1, qaws_diff_views* sink)
{
	qaws_curve_jet_3d bar;
	memset(&bar, 0, sizeof(bar));
	bar.d[1] = bar_d1;
	bar.channels = QAWS_EVAL_FLAG_D1;
	return qaws_internal_curve_adjoint_any(job->ctx, job->curve, (qaws_scalar)t, QAWS_EVAL_FLAG_D1, &bar, sink, NULL);
}

/*
 * Pulls back w1(tau) u(tau) + w2(tau) u_dot(tau) through the D1 channel
 * over the length integrals: every full span k gets
 *   w = F - sum of the samples whose span lies beyond k
 * and every sample subtracts its own weight over [a_k, t_i]. c1 / c2 are the
 * per-sample weights of the u and u_dot terms (c2 may be NULL: no u_dot).
 */
static qaws_status arc_pullback_lengths(arc_job const* job, arc_sample const* smp, unsigned int count,
	double const* c1, double const* c2, double F1, double F2, qaws_diff_views* sink)
{
	unsigned int S = job->spans, k, i, q;
	double* beyond = (double*)qaws_internal_alloc(NULL, (unsigned long)(sizeof(double) * 2 * (S + 1)));
	qaws_status st = QAWS_STATUS_OK;
	if (!beyond)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	memset(beyond, 0, sizeof(double) * 2 * (S + 1));
	/* beyond[k] = sum over samples with span > k */
	for (i = 0; i < count; i++)
		if (smp[i].s > QAWS_EPSILON && smp[i].span > 0)
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
		double w1 = F1 - beyond[k], w2 = F2 - beyond[S + 1 + k];
		double a = job->curve->span_boundaries[k], b = job->curve->span_boundaries[k + 1];
		if (w1 == 0 && (!c2 || w2 == 0))
			continue;
		for (q = 0; q < job->n * ARC_PIECES && st == QAWS_STATUS_OK; q++)
		{
			double t, W, s;
			arc_node(job, a, b, q, &t, &W);
			qaws_curve_jet_3d p, tg;
			qaws_vec3 u, bar;
			st = arc_jets(job, (qaws_scalar)t, QAWS_ZERO, c2 ? job->direction : NULL, &p, &tg, NULL);
			if (st != QAWS_STATUS_OK)
				break;
			s = sqrt(v3_dot(p.d[1], p.d[1]));
			if (!(s > QAWS_EPSILON))
				continue;
			u = v3_make(p.d[1].x / s, p.d[1].y / s, p.d[1].z / s);
			bar = v3_make(W * w1 * u.x, W * w1 * u.y, W * w1 * u.z);
			if (c2)
			{
				double ud = v3_dot(u, tg.d[1]);
				qaws_vec3 udot = v3_make((tg.d[1].x - ud * u.x) / s, (tg.d[1].y - ud * u.y) / s, (tg.d[1].z - ud * u.z) / s);
				bar = v3_axpy(W * w2, udot, bar);
			}
			st = arc_add_d1(job, t, bar, sink);
		}
	}
	/* partial spans of the samples: minus their weight over [a_k, t_i] */
	for (i = 0; i < count && st == QAWS_STATUS_OK; i++)
	{
		double a = job->curve->span_boundaries[smp[i].span], b = smp[i].t;
		if (!(smp[i].s > QAWS_EPSILON) || !(b > a))
			continue;
		for (q = 0; q < job->n * ARC_PIECES && st == QAWS_STATUS_OK; q++)
		{
			double t, W, s;
			arc_node(job, a, b, q, &t, &W);
			qaws_curve_jet_3d p, tg;
			qaws_vec3 u, bar;
			st = arc_jets(job, (qaws_scalar)t, QAWS_ZERO, c2 ? job->direction : NULL, &p, &tg, NULL);
			if (st != QAWS_STATUS_OK)
				break;
			s = sqrt(v3_dot(p.d[1], p.d[1]));
			if (!(s > QAWS_EPSILON))
				continue;
			u = v3_make(p.d[1].x / s, p.d[1].y / s, p.d[1].z / s);
			bar = v3_make(-W * c1[i] * u.x, -W * c1[i] * u.y, -W * c1[i] * u.z);
			if (c2)
			{
				double ud = v3_dot(u, tg.d[1]);
				qaws_vec3 udot = v3_make((tg.d[1].x - ud * u.x) / s, (tg.d[1].y - ud * u.y) / s, (tg.d[1].z - ud * u.z) / s);
				bar = v3_axpy(-W * c2[i], udot, bar);
			}
			st = arc_add_d1(job, t, bar, sink);
		}
	}
	qaws_internal_dealloc(NULL, beyond);
	return st;
}

/* Solves the samples and their lambda = (t_bar + p_bar . C') / |C'|; the
   position pullback J^T p_bar is added into sink when it is given. */
static qaws_status arc_samples(arc_job const* job, qaws_arc_length_target const* targets, unsigned int count,
	qaws_arc_length_sample const* adjoint, qaws_diff_views* sink, arc_sample* smp)
{
	unsigned int i;
	double total = arc_total(job);
	for (i = 0; i < count; i++)
	{
		qaws_curve_jet_3d p, tg;
		qaws_status st = arc_solve(job, target_sigma(&targets[i], total), &smp[i].t, &smp[i].span);
		qaws_scalar t_adj = QAWS_ZERO;
		if (st != QAWS_STATUS_OK)
			return st;
		st = arc_jets(job, (qaws_scalar)smp[i].t, QAWS_ZERO, NULL, &p, &tg, NULL);
		if (st != QAWS_STATUS_OK)
			return st;
		smp[i].s = sqrt(v3_dot(p.d[1], p.d[1]));
		smp[i].lambda = 0;
		smp[i].lambda_dot = 0;
		smp[i].t_dot = 0;
		smp[i].u = qaws_v3_zero();
		if (!(smp[i].s > QAWS_EPSILON))
		{
			qaws_internal_diff_report_note(job->ctx, QAWS_DIFF_SMOOTH, QAWS_DIFF_ILL_CONDITIONED, 0, i);
			continue;
		}
		qaws_internal_diff_report_note(job->ctx, QAWS_DIFF_SMOOTH, QAWS_DIFF_VALID, 0, i);
		smp[i].u = v3_make(p.d[1].x / smp[i].s, p.d[1].y / smp[i].s, p.d[1].z / smp[i].s);
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
		smp[i].lambda = ((double)adjoint[i].t + t_adj) / smp[i].s;
	}
	return QAWS_STATUS_OK;
}

qaws_status qaws_curve_arc_length_sample_adjoint(
	qaws_diff_context const* ctx,
	qaws_curve const* curve,
	qaws_arc_length_target const* targets,
	unsigned int count,
	unsigned int quadrature,
	qaws_arc_length_sample const* adjoint,
	qaws_diff_views* param_adjoint,
	qaws_scalar* distance_adjoint)
{
	arc_job job;
	arc_sample* smp;
	double* c1;
	double F = 0;
	unsigned int i;
	qaws_status st;
	if ((!targets || !adjoint) && count)
		return QAWS_STATUS_INVALID_ARGUMENT;
	st = arc_job_init(&job, ctx, curve, quadrature, NULL, param_adjoint);
	if (st != QAWS_STATUS_OK)
		return st;
	smp = (arc_sample*)qaws_internal_alloc(NULL, (unsigned long)((sizeof(arc_sample) + sizeof(double)) * (count + 1)));
	if (!smp)
	{
		arc_job_free(&job);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	c1 = (double*)(smp + count + 1);
	st = arc_samples(&job, targets, count, adjoint, param_adjoint, smp);
	for (i = 0; i < count && st == QAWS_STATUS_OK; i++)
	{
		c1[i] = smp[i].lambda;
		F += smp[i].lambda * (double)targets[i].fraction;
		if (distance_adjoint)
			distance_adjoint[i] += (qaws_scalar)smp[i].lambda;
	}
	if (st == QAWS_STATUS_OK && param_adjoint)
		st = arc_pullback_lengths(&job, smp, count, c1, NULL, F, 0, param_adjoint);
	qaws_internal_dealloc(NULL, smp);
	arc_job_free(&job);
	return st;
}

/*
 * HVP, forward over reverse. The gradient is
 *   G = sum_i J(t_i)^T p_bar_i + lambda_i (f_i grad L_total - grad L(t_i))
 * and, for a family linear in its fields (J = basis at t, grad L = integral
 * of B'^T u), its derivative along the direction is
 *   sum_i  B'(t_i)^T (p_bar_i - lambda_i u_i) t'_i
 *        + lambda'_i (f_i grad L_total - grad L(t_i))
 *        + lambda_i  (f_i grad' L_total - grad' L(t_i))
 * with grad' L the integral of B'^T u_dot, u_dot = (I - u u^T) dC' / |C'|,
 * lambda' = (p_bar . dC'(t, t') - lambda s') / s and s' = u . dC'(t, t').
 */
qaws_status qaws_curve_arc_length_sample_hvp(
	qaws_diff_context const* ctx,
	qaws_curve const* curve,
	qaws_arc_length_target const* targets,
	unsigned int count,
	unsigned int quadrature,
	qaws_arc_length_sample const* adjoint,
	qaws_diff_views const* direction,
	qaws_diff_views* out_hv)
{
	arc_job job;
	arc_sample* smp;
	double *c1, *c2;
	double F1 = 0, F2 = 0;
	unsigned int i;
	qaws_status st;
	if ((!targets || !adjoint) && count)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (!curve)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (!(qaws_curve_get_diff_capabilities(curve) & QAWS_CAP_LINEAR))
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	if (!direction || !out_hv)
		return QAWS_STATUS_OK;
	st = arc_job_init(&job, ctx, curve, quadrature, direction, out_hv);
	if (st != QAWS_STATUS_OK)
		return st;
	smp = (arc_sample*)qaws_internal_alloc(NULL, (unsigned long)((sizeof(arc_sample) + 2 * sizeof(double)) * (count + 1)));
	if (!smp)
	{
		arc_job_free(&job);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	c1 = (double*)(smp + count + 1);
	c2 = c1 + count + 1;
	st = arc_samples(&job, targets, count, adjoint, NULL, smp);
	for (i = 0; i < count && st == QAWS_STATUS_OK; i++)
	{
		qaws_curve_jet_3d p, tg;
		double v, Ld1, sdot;
		arc_sample* a = &smp[i];
		c1[i] = c2[i] = 0;
		if (!(a->s > QAWS_EPSILON))
			continue;
		st = arc_integrate(&job, curve->span_boundaries[a->span], a->t, &v, &Ld1, NULL);
		if (st != QAWS_STATUS_OK)
			break;
		Ld1 += job.cum1[a->span];
		a->t_dot = ((double)targets[i].fraction * job.cum1[job.spans] - Ld1) / a->s;
		st = arc_jets(&job, (qaws_scalar)a->t, (qaws_scalar)a->t_dot, direction, &p, &tg, NULL);
		if (st != QAWS_STATUS_OK)
			break;
		/* tg.d[1] = dC'(t, t'): the direction plus C'' t' */
		sdot = v3_dot(a->u, tg.d[1]);
		a->lambda_dot = (v3_dot(adjoint[i].position, tg.d[1]) - a->lambda * sdot) / a->s;
		c1[i] = a->lambda_dot;
		c2[i] = a->lambda;
		F1 += a->lambda_dot * (double)targets[i].fraction;
		F2 += a->lambda * (double)targets[i].fraction;
		/* moving evaluation point: B'(t)^T (p_bar - lambda u) t' */
		st = arc_add_d1(&job, a->t, v3_make(a->t_dot * (adjoint[i].position.x - a->lambda * a->u.x),
			a->t_dot * (adjoint[i].position.y - a->lambda * a->u.y),
			a->t_dot * (adjoint[i].position.z - a->lambda * a->u.z)), out_hv);
	}
	if (st == QAWS_STATUS_OK)
		st = arc_pullback_lengths(&job, smp, count, c1, c2, F1, F2, out_hv);
	qaws_internal_dealloc(NULL, smp);
	arc_job_free(&job);
	return st;
}
