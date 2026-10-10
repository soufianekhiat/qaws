#include "qaws_clothoid.h"
#include "qaws_curve.h"
#include "internal/qaws_internal_types.h"
#include "internal/qaws_internal_kinds.h"
#include "internal/qaws_internal_curve.h"
#include "internal/qaws_internal_diff.h"
#include "core/qaws_dual_core.h"
#include "qaws_diff.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ---------------------------------------------------------------------------
 * Impl struct
 * ------------------------------------------------------------------------- */


/* ---------------------------------------------------------------------------
 * Helpers
 * ------------------------------------------------------------------------- */

/* Compute theta(s) = start_angle + kappa_0 * s + (kappa_1 - kappa_0) * s^2 / (2 * L) */
static qaws_scalar clothoid_theta(qaws_clothoid_impl const *impl, qaws_scalar s)
{
	qaws_scalar kd = (impl->kappa_1 - impl->kappa_0) / impl->length;
	return impl->start_angle + impl->kappa_0 * s + kd * s * s * (qaws_scalar)0.5;
}

/* Compute position at arc-length s: composite 10-point Gauss-Legendre, with
   enough panels that the heading turns at most half a radian per panel. */
static void clothoid_position(qaws_clothoid_impl const *impl, qaws_scalar s,
	qaws_scalar *out_x, qaws_scalar *out_y)
{
	static double const gx[5] = { 0.1488743389816312, 0.4333953941292472, 0.6794095682990244,
		0.8650633666889845, 0.9739065285171717 };
	static double const gw[5] = { 0.2955242247147529, 0.2692667143361318, 0.2190863625159820,
		0.1494513491505806, 0.0666713443086881 };
	double kd = ((double)impl->kappa_1 - (double)impl->kappa_0) / (double)impl->length;
	double k0 = fabs((double)impl->kappa_0), k1 = fabs((double)impl->kappa_0 + kd * (double)s);
	double turn = fabs((double)s) * (k0 > k1 ? k0 : k1);
	unsigned int panels = 1 + (unsigned int)(turn * 2.0), p, i;
	double h = (double)s / (double)panels, sum_x = 0.0, sum_y = 0.0;

	if (panels > 4096)
		panels = 4096, h = (double)s / 4096.0;
	for (p = 0; p < panels; p++)
	{
		double mid = h * ((double)p + 0.5);
		for (i = 0; i < 10; i++)
		{
			double u = mid + 0.5 * h * (i < 5 ? -gx[i] : gx[i - 5]);
			double theta = (double)impl->start_angle + (double)impl->kappa_0 * u + kd * u * u * 0.5;
			double w = gw[i < 5 ? i : i - 5];
			sum_x += w * cos(theta);
			sum_y += w * sin(theta);
		}
	}
	*out_x = impl->origin_x + (qaws_scalar)(sum_x * h * 0.5);
	*out_y = impl->origin_y + (qaws_scalar)(sum_y * h * 0.5);
}

/* ---------------------------------------------------------------------------
 * Vtable: eval_span_2d
 * ------------------------------------------------------------------------- */

static qaws_status clothoid_eval_span_2d(
	qaws_curve const *curve,
	unsigned int span_index,
	qaws_scalar local_t,
	unsigned int eval_flags,
	qaws_eval_result_2d *out_result)
{
	qaws_clothoid_impl *impl = (qaws_clothoid_impl *)curve->impl;
	qaws_scalar L = impl->length;
	qaws_scalar s = local_t * L; /* local_t is [0,1]; map to arc-length */
	qaws_scalar theta;
	qaws_scalar cos_theta, sin_theta;
	qaws_scalar kappa, kappa_prime;

	(void)span_index;

	memset(out_result, 0, sizeof(*out_result));

	/* Position: numerical integration */
	if (eval_flags & QAWS_EVAL_FLAG_POSITION) {
		if (s <= (qaws_scalar)0) {
			out_result->position.x = impl->origin_x;
			out_result->position.y = impl->origin_y;
		} else {
			clothoid_position(impl, s, &out_result->position.x, &out_result->position.y);
		}
		out_result->valid_flags |= QAWS_EVAL_FLAG_POSITION;
	}

	/* Precompute theta and trig values for derivatives */
	theta = clothoid_theta(impl, s);
	cos_theta = (qaws_scalar)cos((double)theta);
	sin_theta = (qaws_scalar)sin((double)theta);

	/* Derivatives with respect to the curve parameter, the arc length s
	   (local_t only locates the point): D1 = (cos, sin) */
	if (eval_flags & (QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3))
	{
		out_result->d1.x = cos_theta;
		out_result->d1.y = sin_theta;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D1;
	}

	/* D2 = kappa (-sin, cos) */
	if (eval_flags & (QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3))
	{
		kappa = impl->kappa_0 + (impl->kappa_1 - impl->kappa_0) * s / L;
		out_result->d2.x = kappa * (-sin_theta);
		out_result->d2.y = kappa * cos_theta;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D2;
	}

	/* D3 = kappa' (-sin, cos) + kappa^2 (-cos, -sin) */
	if (eval_flags & QAWS_EVAL_FLAG_D3)
	{
		kappa = impl->kappa_0 + (impl->kappa_1 - impl->kappa_0) * s / L;
		kappa_prime = (impl->kappa_1 - impl->kappa_0) / L;
		out_result->d3.x = kappa_prime * (-sin_theta) + kappa * kappa * (-cos_theta);
		out_result->d3.y = kappa_prime * cos_theta + kappa * kappa * (-sin_theta);
		out_result->valid_flags |= QAWS_EVAL_FLAG_D3;
	}

	return QAWS_STATUS_OK;
}

/* ---------------------------------------------------------------------------
 * Vtable: eval_span_3d
 * ------------------------------------------------------------------------- */

static qaws_status clothoid_eval_span_3d(
	qaws_curve const *curve,
	unsigned int span_index,
	qaws_scalar local_t,
	unsigned int eval_flags,
	qaws_eval_result_3d *out_result)
{
	(void)curve;
	(void)span_index;
	(void)local_t;
	(void)eval_flags;
	(void)out_result;
	return QAWS_STATUS_INVALID_DIMENSION;
}

/* ---------------------------------------------------------------------------
 * Vtable: destroy
 * ------------------------------------------------------------------------- */

static void clothoid_destroy_impl(void *impl, qaws_allocator const* allocator)
{
	if (impl)
		qaws_internal_dealloc(allocator, impl);
}

/* ---------------------------------------------------------------------------
 * Vtable: property queries
 * ------------------------------------------------------------------------- */

static int clothoid_is_closed(qaws_curve const *curve)
{
	return qaws_internal_curve_ends_meet(curve);
}

static int clothoid_is_periodic(qaws_curve const *curve)
{
	(void)curve;
	return 0;
}

static int clothoid_is_rational(qaws_curve const *curve)
{
	(void)curve;
	return 0;
}

static qaws_continuity clothoid_get_continuity(qaws_curve const *curve)
{
	(void)curve;
	return QAWS_CONTINUITY_C3;
}

/* ---------------------------------------------------------------------------
 * Differentiation
 *
 * Fields: CENTER (the origin), ANGLE_START (theta_0), CURVATURE (kappa_0)
 * and CURVATURE_RATE (c = (kappa_1 - kappa_0) / L). With them the curve at
 * arc length s,
 *   C(s) = origin + int_0^s (cos, sin)(theta_0 + kappa_0 u + c u^2 / 2) du,
 * does not depend on L (the length only sets the domain). The jets are
 * evaluated in dual numbers, the position through the same Simpson rule as
 * the primal, so tangents and second tangents are exact for it.
 * ------------------------------------------------------------------------- */

#define CL_DIFF_CAPS (QAWS_CAP_TANGENT | QAWS_CAP_ADJOINT | QAWS_CAP_TANGENT2)

static qaws_dual1 cl_dual(qaws_scalar v, qaws_scalar t)
{
	return qaws_dual1_make(v, t, QAWS_ZERO);
}

static qaws_dual1 cl_sin(qaws_dual1 a)
{
	qaws_scalar s = (qaws_scalar)sin((double)a.v), c = (qaws_scalar)cos((double)a.v);
	return qaws_dual1_make(s, c * a.t, c * a.tt - s * a.t * a.t);
}

static qaws_dual1 cl_cos(qaws_dual1 a)
{
	qaws_scalar s = (qaws_scalar)sin((double)a.v), c = (qaws_scalar)cos((double)a.v);
	return qaws_dual1_make(c, -s * a.t, -s * a.tt - c * a.t * a.t);
}

/* Parameters as duals: origin, theta_0, kappa_0, rate. */
typedef struct cl_params
{
	qaws_dual1 ox, oy, th0, k0, c;
} cl_params;

static qaws_dual1 cl_theta(cl_params const *p, qaws_dual1 u)
{
	qaws_dual1 half_u2 = qaws_dual1_mul(qaws_dual1_make((qaws_scalar)0.5, 0, 0), qaws_dual1_mul(u, u));
	return qaws_dual1_add(p->th0, qaws_dual1_add(qaws_dual1_mul(p->k0, u), qaws_dual1_mul(p->c, half_u2)));
}

/* Jet (P, D1, D2, D3) as x / y duals at the dual arc length s. */
static void cl_jet(cl_params const *p, qaws_dual1 s, qaws_dual1 out[4][2])
{
	unsigned int n = 64, i;
	qaws_dual1 h = qaws_dual1_mul(s, qaws_dual1_make(QAWS_ONE / (qaws_scalar)n, 0, 0));
	qaws_dual1 sx = qaws_dual1_const(0), sy = qaws_dual1_const(0), th, kappa, ct, st;
	for (i = 0; i <= n; i++)
	{
		qaws_dual1 u = qaws_dual1_mul(h, qaws_dual1_const((qaws_scalar)i)), t = cl_theta(p, u);
		qaws_scalar w = (i == 0 || i == n) ? (qaws_scalar)1 : (i % 2 == 1) ? (qaws_scalar)4 : (qaws_scalar)2;
		sx = qaws_dual1_add(sx, qaws_dual1_mul(qaws_dual1_const(w), cl_cos(t)));
		sy = qaws_dual1_add(sy, qaws_dual1_mul(qaws_dual1_const(w), cl_sin(t)));
	}
	h = qaws_dual1_mul(h, qaws_dual1_const(QAWS_ONE / (qaws_scalar)3));
	out[0][0] = qaws_dual1_add(p->ox, qaws_dual1_mul(sx, h));
	out[0][1] = qaws_dual1_add(p->oy, qaws_dual1_mul(sy, h));
	th = cl_theta(p, s);
	ct = cl_cos(th);
	st = cl_sin(th);
	kappa = qaws_dual1_add(p->k0, qaws_dual1_mul(p->c, s));
	out[1][0] = ct;
	out[1][1] = st;
	out[2][0] = qaws_dual1_sub(qaws_dual1_const(0), qaws_dual1_mul(kappa, st));
	out[2][1] = qaws_dual1_mul(kappa, ct);
	{
		qaws_dual1 k2 = qaws_dual1_mul(kappa, kappa);
		out[3][0] = qaws_dual1_sub(qaws_dual1_sub(qaws_dual1_const(0), qaws_dual1_mul(p->c, st)), qaws_dual1_mul(k2, ct));
		out[3][1] = qaws_dual1_sub(qaws_dual1_mul(p->c, ct), qaws_dual1_mul(k2, st));
	}
}

static void cl_read(qaws_diff_views const *views, qaws_diff_field field, unsigned int comps, qaws_scalar *out)
{
	qaws_field_view const *v = views ? qaws_diff_views_find(views, field) : NULL;
	out[0] = out[1] = 0;
	if (v)
		qaws_internal_view_read(v, 0, comps, out);
}

static void cl_params_make(qaws_clothoid_impl const *impl, qaws_scalar const *d_origin, qaws_scalar d_th0,
	qaws_scalar d_k0, qaws_scalar d_c, cl_params *p)
{
	p->ox = cl_dual(impl->origin[0], d_origin[0]);
	p->oy = cl_dual(impl->origin[1], d_origin[1]);
	p->th0 = cl_dual(impl->start_angle, d_th0);
	p->k0 = cl_dual(impl->kappa_0, d_k0);
	p->c = cl_dual(impl->rate, d_c);
}

static qaws_status clothoid_tangent_span(
	qaws_diff_context const *ctx, qaws_curve const *curve, unsigned int span_index, qaws_scalar local_t,
	qaws_scalar t_dot, unsigned int channels, qaws_diff_views const *views,
	qaws_curve_jet_3d *primal, qaws_curve_jet_3d *tangent, qaws_curve_jet_3d *tangent2)
{
	qaws_clothoid_impl const *impl = (qaws_clothoid_impl const *)curve->impl;
	qaws_scalar d_origin[2], d_th0[2], d_k0[2], d_c[2];
	qaws_dual1 jet[4][2];
	cl_params p;
	unsigned int k;
	(void)ctx;
	(void)span_index;
	cl_read(views, QAWS_FIELD_CENTER, 2, d_origin);
	cl_read(views, QAWS_FIELD_ANGLE_START, 1, d_th0);
	cl_read(views, QAWS_FIELD_CURVATURE, 1, d_k0);
	cl_read(views, QAWS_FIELD_CURVATURE_RATE, 1, d_c);
	cl_params_make(impl, d_origin, d_th0[0], d_k0[0], d_c[0], &p);
	cl_jet(&p, cl_dual(local_t * impl->length, t_dot), jet);
	memset(primal, 0, sizeof(*primal));
	memset(tangent, 0, sizeof(*tangent));
	if (tangent2)
		memset(tangent2, 0, sizeof(*tangent2));
	for (k = 0; k <= QAWS_CURVE_JET_ORDER; k++)
	{
		if (!(channels & (1u << k)))
			continue;
		primal->d[k].x = jet[k][0].v;
		primal->d[k].y = jet[k][1].v;
		tangent->d[k].x = jet[k][0].t;
		tangent->d[k].y = jet[k][1].t;
		if (tangent2)
		{
			tangent2->d[k].x = jet[k][0].tt;
			tangent2->d[k].y = jet[k][1].tt;
		}
	}
	primal->channels = tangent->channels = channels;
	if (tangent2)
		tangent2->channels = channels;
	return QAWS_STATUS_OK;
}

/* <jet adjoint, jet tangent> for one seeded parameter (or the arc length). */
static qaws_scalar cl_pullback(qaws_clothoid_impl const *impl, qaws_scalar s, int which, unsigned int channels,
	qaws_curve_jet_3d const *bar)
{
	qaws_scalar zero[2] = { 0, 0 }, sum = 0;
	qaws_dual1 jet[4][2];
	cl_params p;
	unsigned int k;
	cl_params_make(impl, zero, which == 0 ? QAWS_ONE : 0, which == 1 ? QAWS_ONE : 0, which == 2 ? QAWS_ONE : 0, &p);
	cl_jet(&p, cl_dual(s, which == 3 ? QAWS_ONE : 0), jet);
	for (k = 0; k <= QAWS_CURVE_JET_ORDER; k++)
		if (channels & (1u << k))
			sum += bar->d[k].x * jet[k][0].t + bar->d[k].y * jet[k][1].t;
	return sum;
}

static qaws_status clothoid_adjoint_span(
	qaws_diff_context const *ctx, qaws_curve const *curve, unsigned int span_index, qaws_scalar local_t,
	unsigned int channels, qaws_curve_jet_3d const *jet_adjoint, qaws_diff_views *views, qaws_scalar *t_adjoint)
{
	qaws_clothoid_impl const *impl = (qaws_clothoid_impl const *)curve->impl;
	qaws_scalar s = local_t * impl->length;
	static qaws_diff_field const scalar_fields[3] = { QAWS_FIELD_ANGLE_START, QAWS_FIELD_CURVATURE, QAWS_FIELD_CURVATURE_RATE };
	unsigned int i;
	(void)ctx;
	(void)span_index;
	if (t_adjoint)
		*t_adjoint += cl_pullback(impl, s, 3, channels, jet_adjoint);
	if (!views)
		return QAWS_STATUS_OK;
	{
		/* the origin moves the position only */
		qaws_field_view *ov = qaws_diff_views_find(views, QAWS_FIELD_CENTER);
		if (ov && (channels & QAWS_EVAL_FLAG_POSITION))
		{
			qaws_scalar g[2];
			g[0] = jet_adjoint->d[0].x;
			g[1] = jet_adjoint->d[0].y;
			qaws_internal_view_add(ov, 0, 2, g);
		}
	}
	for (i = 0; i < 3; i++)
	{
		qaws_field_view *v = qaws_diff_views_find(views, scalar_fields[i]);
		if (v)
		{
			qaws_scalar g = cl_pullback(impl, s, (int)i, channels, jet_adjoint);
			qaws_internal_view_add(v, 0, 1, &g);
		}
	}
	return QAWS_STATUS_OK;
}

static unsigned int clothoid_describe_fields(qaws_curve const *curve, qaws_field_desc *out, unsigned int capacity)
{
	(void)curve;
	if (capacity >= 1)
		out[0] = qaws_internal_field_desc(QAWS_FIELD_CENTER, QAWS_VALUE_VEC2, 1, QAWS_DOMAIN_POSITION, QAWS_CONSTRAINT_NONE,
			QAWS_DIFF_SMOOTH, CL_DIFF_CAPS);
	if (capacity >= 2)
		out[1] = qaws_internal_field_desc(QAWS_FIELD_ANGLE_START, QAWS_VALUE_SCALAR, 1, QAWS_DOMAIN_ANGLE, QAWS_CONSTRAINT_NONE,
			QAWS_DIFF_SMOOTH, CL_DIFF_CAPS);
	if (capacity >= 3)
		out[2] = qaws_internal_field_desc(QAWS_FIELD_CURVATURE, QAWS_VALUE_SCALAR, 1, QAWS_DOMAIN_CURVATURE, QAWS_CONSTRAINT_NONE,
			QAWS_DIFF_SMOOTH, CL_DIFF_CAPS);
	if (capacity >= 4)
		out[3] = qaws_internal_field_desc(QAWS_FIELD_CURVATURE_RATE, QAWS_VALUE_SCALAR, 1, QAWS_DOMAIN_GENERIC, QAWS_CONSTRAINT_NONE,
			QAWS_DIFF_SMOOTH, CL_DIFF_CAPS);
	return 4;
}

static qaws_status clothoid_primal_field(qaws_curve const *curve, qaws_diff_field field,
	qaws_scalar const **out_data, unsigned int *out_count, unsigned int *out_components)
{
	qaws_clothoid_impl const *impl = (qaws_clothoid_impl const *)curve->impl;
	*out_count = 1;
	*out_components = 1;
	switch (field)
	{
	case QAWS_FIELD_CENTER: *out_data = impl->origin; *out_components = 2; return QAWS_STATUS_OK;
	case QAWS_FIELD_ANGLE_START: *out_data = &impl->start_angle; return QAWS_STATUS_OK;
	case QAWS_FIELD_CURVATURE: *out_data = &impl->kappa_0; return QAWS_STATUS_OK;
	case QAWS_FIELD_CURVATURE_RATE: *out_data = &impl->rate; return QAWS_STATUS_OK;
	default: return QAWS_STATUS_INVALID_ARGUMENT;
	}
}

static qaws_status clothoid_rebuild(qaws_curve const *curve, qaws_diff_views const *values, qaws_curve **out_curve)
{
	qaws_clothoid_impl const *impl = (qaws_clothoid_impl const *)curve->impl;
	qaws_clothoid_desc d;
	qaws_scalar origin[2], th0[2], k0[2], rate[2];
	qaws_field_view const *v;
	origin[0] = impl->origin[0];
	origin[1] = impl->origin[1];
	th0[0] = impl->start_angle;
	k0[0] = impl->kappa_0;
	rate[0] = impl->rate;
	if (values)
	{
		if ((v = qaws_diff_views_find(values, QAWS_FIELD_CENTER)) != NULL && v->data) qaws_internal_view_read(v, 0, 2, origin);
		if ((v = qaws_diff_views_find(values, QAWS_FIELD_ANGLE_START)) != NULL && v->data) qaws_internal_view_read(v, 0, 1, th0);
		if ((v = qaws_diff_views_find(values, QAWS_FIELD_CURVATURE)) != NULL && v->data) qaws_internal_view_read(v, 0, 1, k0);
		if ((v = qaws_diff_views_find(values, QAWS_FIELD_CURVATURE_RATE)) != NULL && v->data) qaws_internal_view_read(v, 0, 1, rate);
	}
	memset(&d, 0, sizeof(d));
	d.origin_x = origin[0];
	d.origin_y = origin[1];
	d.start_angle = th0[0];
	d.start_curvature = k0[0];
	d.end_curvature = k0[0] + rate[0] * impl->length;
	d.length = impl->length;
	return qaws_curve_create_clothoid(&d, out_curve);
}

static qaws_curve_diff_vtable const clothoid_diff_vtable = {
	CL_DIFF_CAPS,
	QAWS_DIFF_SMOOTH,
	clothoid_describe_fields,
	clothoid_primal_field,
	NULL,
	clothoid_tangent_span,
	clothoid_adjoint_span,
	NULL,
	clothoid_rebuild
};

/* ---------------------------------------------------------------------------
 * Vtable definition
 * ------------------------------------------------------------------------- */

static qaws_curve_vtable const clothoid_vtable = {
	clothoid_eval_span_2d,
	clothoid_eval_span_3d,
	clothoid_destroy_impl,
	clothoid_is_closed,
	clothoid_is_periodic,
	clothoid_is_rational,
	clothoid_get_continuity,
	&clothoid_diff_vtable
};

/* ---------------------------------------------------------------------------
 * Creation
 * ------------------------------------------------------------------------- */

qaws_status qaws_curve_create_clothoid(
	qaws_clothoid_desc const *desc,
	qaws_curve **out_curve)
{
	qaws_range parameter_range;
	qaws_curve *curve;
	qaws_clothoid_impl *impl;

	if (!desc)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (!out_curve)
		return QAWS_STATUS_INVALID_ARGUMENT;

	*out_curve = NULL;

	/* Length must be positive */
	if (desc->length <= (qaws_scalar)0)
		return QAWS_STATUS_INVALID_ARGUMENT;

	/* Parameter range: arc-length [0, L] */
	parameter_range.min_value = (qaws_scalar)0;
	parameter_range.max_value = desc->length;

	/* Allocate curve: single span, degree 3 convention, 2D only */
	curve = qaws_internal_curve_alloc(
		QAWS_CURVE_KIND_CLOTHOID,
		QAWS_DIMENSION_2D,
		3,
		1,
		parameter_range,
		&clothoid_vtable);

	if (!curve)
		return QAWS_STATUS_ALLOCATION_FAILURE;

	/* Set span boundaries: [0, L] */
	curve->span_boundaries[0] = (qaws_scalar)0;
	curve->span_boundaries[1] = desc->length;

	/* Allocate impl */
	impl = (qaws_clothoid_impl *)calloc(1, sizeof(qaws_clothoid_impl));
	if (!impl)
	{
		qaws_curve_destroy(curve);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	impl->origin_x = desc->origin_x;
	impl->origin_y = desc->origin_y;
	impl->start_angle = desc->start_angle;
	impl->kappa_0 = desc->start_curvature;
	impl->kappa_1 = desc->end_curvature;
	impl->length = desc->length;
	impl->origin[0] = desc->origin_x;
	impl->origin[1] = desc->origin_y;
	impl->rate = (desc->end_curvature - desc->start_curvature) / desc->length;

	curve->impl = impl;
	*out_curve = curve;

	return QAWS_STATUS_OK;
}
