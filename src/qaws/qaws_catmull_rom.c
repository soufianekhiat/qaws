#include "qaws_catmull_rom.h"
#include "qaws_curve.h"
#include "qaws_prepare.h"
#include "internal/qaws_internal_types.h"
#include "internal/qaws_internal_curve.h"
#include "internal/qaws_internal_basis.h"
#include "internal/qaws_internal_validation.h"
#include "core/qaws_cubic_poly_core.h"
#include "core/qaws_dual_core.h"
#include "internal/qaws_internal_diff.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "qaws_platform.h"

/* -------------------------------------------------------------------------- */
/*  Helpers                                                                   */
/* -------------------------------------------------------------------------- */

static qaws_scalar compute_alpha(qaws_parameterization param)
{
	switch (param) {
	case QAWS_PARAMETERIZATION_UNIFORM:
		return QAWS_ZERO;
	case QAWS_PARAMETERIZATION_CHORDAL:
		return QAWS_ONE;
	case QAWS_PARAMETERIZATION_CENTRIPETAL:
		return QAWS_LITERAL(0.5);
	case QAWS_PARAMETERIZATION_DEFAULT:
	default:
		return QAWS_LITERAL(0.5); /* default is centripetal */
	}
}

static void compute_knot_params(
	qaws_scalar const *control_points,
	unsigned int control_point_count,
	unsigned int dim_count,
	qaws_scalar alpha,
	qaws_scalar *knot_params)
{
	unsigned int i, d;

	knot_params[0] = QAWS_ZERO;

	for (i = 1; i < control_point_count; i++) {
		qaws_scalar dist_sq = QAWS_ZERO;
		for (d = 0; d < dim_count; d++) {
			qaws_scalar diff = control_points[i * dim_count + d]
			                 - control_points[(i - 1) * dim_count + d];
			dist_sq += diff * diff;
		}

		if (alpha == QAWS_ZERO) {
			/* Uniform: each interval = 1 */
			knot_params[i] = knot_params[i - 1] + QAWS_ONE;
		} else {
			qaws_scalar dist = QAWS_SQRT(dist_sq);
			qaws_scalar interval = QAWS_POW(dist, alpha);
			knot_params[i] = knot_params[i - 1] + interval;
		}
	}
}

/* -------------------------------------------------------------------------- */
/*  Vtable: eval_span (2D)                                                    */
/* -------------------------------------------------------------------------- */

static qaws_status catmull_rom_eval_span_2d(
	qaws_curve const *curve,
	unsigned int span_index,
	qaws_scalar local_t,
	unsigned int eval_flags,
	qaws_eval_result_2d *out_result)
{
	qaws_catmull_rom_impl const *impl =
		(qaws_catmull_rom_impl const *)curve->impl;
	unsigned int const dim_count = 2;
	qaws_scalar const *cx;
	qaws_scalar const *cy;
	qaws_scalar ax, bx, cx_coeff, dx_val;
	qaws_scalar ay, by, cy_coeff, dy_val;
	qaws_vec2 va, vb, vc, vd;
	qaws_eval_2d core;

	(void)curve;

	memset(out_result, 0, sizeof(*out_result));

	cx = &impl->segment_coeffs[span_index * dim_count * 4 + 0 * 4];
	cy = &impl->segment_coeffs[span_index * dim_count * 4 + 1 * 4];

	ax = cx[0]; bx = cx[1]; cx_coeff = cx[2]; dx_val = cx[3];
	ay = cy[0]; by = cy[1]; cy_coeff = cy[2]; dy_val = cy[3];

	va.x = ax; va.y = ay;
	vb.x = bx; vb.y = by;
	vc.x = cx_coeff; vc.y = cy_coeff;
	vd.x = dx_val; vd.y = dy_val;

	core = qaws_cubic_eval_2d(va, vb, vc, vd, local_t, eval_flags);

	if (eval_flags & QAWS_EVAL_FLAG_POSITION) {
		out_result->position = core.position;
		out_result->valid_flags |= QAWS_EVAL_FLAG_POSITION;
	}

	if (eval_flags & QAWS_EVAL_FLAG_D1) {
		out_result->d1 = core.d1;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D1;
	}

	if (eval_flags & QAWS_EVAL_FLAG_D2) {
		out_result->d2 = core.d2;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D2;
	}

	if (eval_flags & QAWS_EVAL_FLAG_D3) {
		out_result->d3.x = QAWS_LITERAL(6.0) * ax;
		out_result->d3.y = QAWS_LITERAL(6.0) * ay;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D3;
	}

	return QAWS_STATUS_OK;
}

/* -------------------------------------------------------------------------- */
/*  Vtable: eval_span (3D)                                                    */
/* -------------------------------------------------------------------------- */

static qaws_status catmull_rom_eval_span_3d(
	qaws_curve const *curve,
	unsigned int span_index,
	qaws_scalar local_t,
	unsigned int eval_flags,
	qaws_eval_result_3d *out_result)
{
	qaws_catmull_rom_impl const *impl =
		(qaws_catmull_rom_impl const *)curve->impl;
	unsigned int const dim_count = 3;
	qaws_scalar const *cx;
	qaws_scalar const *cy;
	qaws_scalar const *cz;
	qaws_scalar ax, bx, cx_coeff, dx_val;
	qaws_scalar ay, by, cy_coeff, dy_val;
	qaws_scalar az, bz, cz_coeff, dz_val;
	qaws_vec3 va, vb, vc, vd;
	qaws_eval_3d core;

	(void)curve;

	memset(out_result, 0, sizeof(*out_result));

	cx = &impl->segment_coeffs[span_index * dim_count * 4 + 0 * 4];
	cy = &impl->segment_coeffs[span_index * dim_count * 4 + 1 * 4];
	cz = &impl->segment_coeffs[span_index * dim_count * 4 + 2 * 4];

	ax = cx[0]; bx = cx[1]; cx_coeff = cx[2]; dx_val = cx[3];
	ay = cy[0]; by = cy[1]; cy_coeff = cy[2]; dy_val = cy[3];
	az = cz[0]; bz = cz[1]; cz_coeff = cz[2]; dz_val = cz[3];

	va.x = ax; va.y = ay; va.z = az;
	vb.x = bx; vb.y = by; vb.z = bz;
	vc.x = cx_coeff; vc.y = cy_coeff; vc.z = cz_coeff;
	vd.x = dx_val; vd.y = dy_val; vd.z = dz_val;

	core = qaws_cubic_eval_3d(va, vb, vc, vd, local_t, eval_flags);

	if (eval_flags & QAWS_EVAL_FLAG_POSITION) {
		out_result->position = core.position;
		out_result->valid_flags |= QAWS_EVAL_FLAG_POSITION;
	}

	if (eval_flags & QAWS_EVAL_FLAG_D1) {
		out_result->d1 = core.d1;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D1;
	}

	if (eval_flags & QAWS_EVAL_FLAG_D2) {
		out_result->d2 = core.d2;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D2;
	}

	if (eval_flags & QAWS_EVAL_FLAG_D3) {
		out_result->d3.x = QAWS_LITERAL(6.0) * ax;
		out_result->d3.y = QAWS_LITERAL(6.0) * ay;
		out_result->d3.z = QAWS_LITERAL(6.0) * az;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D3;
	}

	return QAWS_STATUS_OK;
}

/* -------------------------------------------------------------------------- */
/*  Vtable: lifecycle and query functions                                     */
/* -------------------------------------------------------------------------- */

static void catmull_rom_destroy_impl(void *impl, qaws_allocator const* allocator)
{
	qaws_catmull_rom_impl *cr = (qaws_catmull_rom_impl *)impl;
	if (cr) {
		qaws_internal_dealloc(allocator, cr->control_points);
		qaws_internal_dealloc(allocator, cr->knot_params);
		qaws_internal_dealloc(allocator, cr->segment_coeffs);
		qaws_internal_dealloc(allocator, cr);
	}
}

static int catmull_rom_is_closed(qaws_curve const *curve)
{
	qaws_catmull_rom_impl const *impl =
		(qaws_catmull_rom_impl const *)curve->impl;
	return impl->closed;
}

static int catmull_rom_is_periodic(qaws_curve const *curve)
{
	qaws_catmull_rom_impl const *impl =
		(qaws_catmull_rom_impl const *)curve->impl;
	return impl->closed;
}

static int catmull_rom_is_rational(qaws_curve const *curve)
{
	(void)curve;
	return 0;
}

static qaws_continuity catmull_rom_get_continuity(qaws_curve const *curve)
{
	(void)curve;
	return QAWS_CONTINUITY_C1;
}

/* -------------------------------------------------------------------------- */
/*  Vtable                                                                    */
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
/*  Differential rules                                                        */
/*                                                                            */
/*  Span coefficients come from the four surrounding points and knot          */
/*  intervals |P_{i+1} - P_i|^alpha, so centripetal and chordal curves are    */
/*  not linear in their points. The coefficient construction is replayed in   */
/*  second-order dual numbers along the point tangent; adjoints sweep each    */
/*  point coordinate of the span.                                             */
/* -------------------------------------------------------------------------- */

#define CR_DIFF_CAPS (QAWS_CAP_TANGENT | QAWS_CAP_ADJOINT | QAWS_CAP_TANGENT2)

static void cr_span_points(qaws_catmull_rom_impl const *impl, unsigned int span, unsigned int *idx)
{
	unsigned int N = impl->control_point_count;
	if (impl->closed) {
		idx[0] = (span == 0) ? N - 1 : span - 1;
		idx[1] = span;
		idx[2] = (span + 1) % N;
		idx[3] = (span + 2) % N;
	} else {
		idx[0] = span; idx[1] = span + 1; idx[2] = span + 2; idx[3] = span + 3;
	}
}

/* |a - b|^alpha in dual numbers (alpha 0, 1/2 or 1). */
static qaws_dual1 cr_interval(qaws_dual1 const *a, qaws_dual1 const *b, unsigned int dim, qaws_scalar alpha)
{
	qaws_dual1 s = qaws_dual1_const(QAWS_ZERO), d;
	unsigned int c;
	if (alpha == QAWS_ZERO)
		return qaws_dual1_const(QAWS_ONE);
	for (c = 0; c < dim; c++) {
		qaws_dual1 diff = qaws_dual1_sub(a[c], b[c]);
		s = qaws_dual1_add(s, qaws_dual1_mul(diff, diff));
	}
	d = qaws_dual1_sqrt(s);
	return alpha == QAWS_ONE ? d : qaws_dual1_sqrt(d);
}

static int cr_nonzero(qaws_dual1 x)
{
	return !(x.v < QAWS_LITERAL(1e-10) && x.v > -QAWS_LITERAL(1e-10));
}

/* coeff[c][0..3] = a, b, c, d of the span polynomial a l^3 + b l^2 + c l + d,
   with point tangents p_dot[4][3] (row i = point idx[i]). */
static void cr_span_coeffs_dual(qaws_catmull_rom_impl const *impl, unsigned int dim, unsigned int span,
	qaws_scalar const p_dot[4][3], qaws_dual1 coeff[3][4])
{
	qaws_dual1 P[4][3], t0, t2, t3, dt10, dt21, dt20, dt32, dt31;
	qaws_scalar alpha = compute_alpha(impl->parameterization);
	unsigned int idx[4], i, c;

	cr_span_points(impl, span, idx);
	for (i = 0; i < 4; i++)
		for (c = 0; c < dim; c++)
			P[i][c] = qaws_dual1_make(impl->control_points[idx[i] * dim + c], p_dot[i][c], QAWS_ZERO);

	/* t1 = 0, t0 = -|P1 - P0|^a, t2 = |P2 - P1|^a, t3 = t2 + |P3 - P2|^a */
	t0 = qaws_dual1_sub(qaws_dual1_const(QAWS_ZERO), cr_interval(P[1], P[0], dim, alpha));
	t2 = cr_interval(P[2], P[1], dim, alpha);
	t3 = qaws_dual1_add(t2, cr_interval(P[3], P[2], dim, alpha));
	dt10 = qaws_dual1_sub(qaws_dual1_const(QAWS_ZERO), t0);
	dt21 = t2;
	dt20 = qaws_dual1_sub(t2, t0);
	dt32 = qaws_dual1_sub(t3, t2);
	dt31 = t3;

	for (c = 0; c < dim; c++) {
		qaws_dual1 m1 = qaws_dual1_const(QAWS_ZERO), m2 = qaws_dual1_const(QAWS_ZERO), M1, M2, d21;
		if (cr_nonzero(dt21) && cr_nonzero(dt20) && cr_nonzero(dt10))
			m1 = qaws_dual1_add(qaws_dual1_sub(
				qaws_dual1_div(qaws_dual1_sub(P[2][c], P[1][c]), dt21),
				qaws_dual1_div(qaws_dual1_sub(P[2][c], P[0][c]), dt20)),
				qaws_dual1_div(qaws_dual1_sub(P[1][c], P[0][c]), dt10));
		if (cr_nonzero(dt32) && cr_nonzero(dt31) && cr_nonzero(dt21))
			m2 = qaws_dual1_add(qaws_dual1_sub(
				qaws_dual1_div(qaws_dual1_sub(P[3][c], P[2][c]), dt32),
				qaws_dual1_div(qaws_dual1_sub(P[3][c], P[1][c]), dt31)),
				qaws_dual1_div(qaws_dual1_sub(P[2][c], P[1][c]), dt21));
		M1 = qaws_dual1_mul(m1, dt21);
		M2 = qaws_dual1_mul(m2, dt21);
		d21 = qaws_dual1_sub(P[2][c], P[1][c]);
		coeff[c][3] = P[1][c];
		coeff[c][2] = M1;
		/* 3 (p2 - p1) - 2 M1 - M2 */
		coeff[c][1] = qaws_dual1_sub(qaws_dual1_sub(qaws_dual1_make(3 * d21.v, 3 * d21.t, 3 * d21.tt),
			qaws_dual1_make(2 * M1.v, 2 * M1.t, 2 * M1.tt)), M2);
		/* 2 (p1 - p2) + M1 + M2 */
		coeff[c][0] = qaws_dual1_add(qaws_dual1_add(qaws_dual1_make(-2 * d21.v, -2 * d21.t, -2 * d21.tt), M1), M2);
	}
	for (c = dim; c < 3; c++)
		for (i = 0; i < 4; i++)
			coeff[c][i] = qaws_dual1_const(QAWS_ZERO);
}

/* k-th derivative of a l^3 + b l^2 + c l + d at l (k = 0..4). */
static qaws_dual1 cr_poly_derivative(qaws_dual1 const *q, qaws_scalar l, unsigned int k)
{
	qaws_scalar w[4];
	unsigned int i;
	qaws_dual1 r = qaws_dual1_const(QAWS_ZERO);
	switch (k) {
	case 0: w[0] = l * l * l; w[1] = l * l; w[2] = l; w[3] = QAWS_ONE; break;
	case 1: w[0] = 3 * l * l; w[1] = 2 * l; w[2] = QAWS_ONE; w[3] = QAWS_ZERO; break;
	case 2: w[0] = 6 * l; w[1] = 2; w[2] = QAWS_ZERO; w[3] = QAWS_ZERO; break;
	case 3: w[0] = 6; w[1] = QAWS_ZERO; w[2] = QAWS_ZERO; w[3] = QAWS_ZERO; break;
	default: return r;
	}
	for (i = 0; i < 4; i++)
		r = qaws_dual1_add(r, qaws_dual1_make(w[i] * q[i].v, w[i] * q[i].t, w[i] * q[i].tt));
	return r;
}

static void cr_set(qaws_vec3 *v, unsigned int c, qaws_scalar x)
{
	if (c == 0) v->x = x; else if (c == 1) v->y = x; else v->z = x;
}

static qaws_scalar cr_get(qaws_vec3 const *v, unsigned int c)
{
	return c == 0 ? v->x : (c == 1 ? v->y : v->z);
}

static void cr_read_tangents(qaws_catmull_rom_impl const *impl, unsigned int dim, unsigned int span,
	qaws_diff_views const *views, qaws_scalar p_dot[4][3])
{
	qaws_field_view const *pv = views ? qaws_diff_views_find(views, QAWS_FIELD_POINTS) : NULL;
	unsigned int idx[4], i;
	memset(p_dot, 0, sizeof(qaws_scalar) * 12);
	if (!pv)
		return;
	cr_span_points(impl, span, idx);
	for (i = 0; i < 4; i++)
		qaws_internal_view_read(pv, idx[i], dim, p_dot[i]);
}

static qaws_status catmull_rom_tangent_span(
	qaws_diff_context const *ctx, qaws_curve const *curve, unsigned int span_index, qaws_scalar local_t,
	qaws_scalar t_dot, unsigned int channels, qaws_diff_views const *views,
	qaws_curve_jet_3d *primal, qaws_curve_jet_3d *tangent, qaws_curve_jet_3d *tangent2)
{
	qaws_catmull_rom_impl const *impl = (qaws_catmull_rom_impl const *)curve->impl;
	unsigned int dim = (unsigned int)curve->dimension, c, k;
	qaws_scalar p_dot[4][3];
	qaws_dual1 coeff[3][4];
	(void)ctx;

	cr_read_tangents(impl, dim, span_index, views, p_dot);
	cr_span_coeffs_dual(impl, dim, span_index, (qaws_scalar const (*)[3])p_dot, coeff);
	memset(primal, 0, sizeof(*primal));
	memset(tangent, 0, sizeof(*tangent));
	if (tangent2)
		memset(tangent2, 0, sizeof(*tangent2));
	for (k = 0; k <= QAWS_CURVE_JET_ORDER; k++) {
		if (!(channels & (1u << k)))
			continue;
		for (c = 0; c < 3; c++) {
			qaws_dual1 q0 = cr_poly_derivative(coeff[c], local_t, k);
			qaws_dual1 q1 = cr_poly_derivative(coeff[c], local_t, k + 1);
			qaws_dual1 q2 = cr_poly_derivative(coeff[c], local_t, k + 2);
			cr_set(&primal->d[k], c, q0.v);
			cr_set(&tangent->d[k], c, q0.t + t_dot * q1.v);
			if (tangent2)
				cr_set(&tangent2->d[k], c, q0.tt + QAWS_LITERAL(2.0) * t_dot * q1.t + t_dot * t_dot * q2.v);
		}
	}
	primal->channels = tangent->channels = channels;
	if (tangent2)
		tangent2->channels = channels;
	return QAWS_STATUS_OK;
}

static qaws_status catmull_rom_adjoint_span(
	qaws_diff_context const *ctx, qaws_curve const *curve, unsigned int span_index, qaws_scalar local_t,
	unsigned int channels, qaws_curve_jet_3d const *jet_adjoint, qaws_diff_views *views, qaws_scalar *t_adjoint)
{
	qaws_catmull_rom_impl const *impl = (qaws_catmull_rom_impl const *)curve->impl;
	qaws_field_view *pv = views ? qaws_diff_views_find(views, QAWS_FIELD_POINTS) : NULL;
	unsigned int dim = (unsigned int)curve->dimension, idx[4], i, c, k, cc;
	qaws_scalar p_dot[4][3];
	qaws_dual1 coeff[3][4];
	(void)ctx;

	memset(p_dot, 0, sizeof(p_dot));
	if (t_adjoint) {
		cr_span_coeffs_dual(impl, dim, span_index, (qaws_scalar const (*)[3])p_dot, coeff);
		for (k = 0; k <= QAWS_CURVE_JET_ORDER; k++)
			if (channels & (1u << k))
				for (cc = 0; cc < dim; cc++)
					*t_adjoint += cr_get(&jet_adjoint->d[k], cc) * cr_poly_derivative(coeff[cc], local_t, k + 1).v;
	}
	if (!pv)
		return QAWS_STATUS_OK;
	cr_span_points(impl, span_index, idx);
	for (i = 0; i < 4; i++)
		for (c = 0; c < dim; c++) {
			qaws_scalar g[3] = { 0, 0, 0 };
			p_dot[i][c] = QAWS_ONE;
			cr_span_coeffs_dual(impl, dim, span_index, (qaws_scalar const (*)[3])p_dot, coeff);
			p_dot[i][c] = QAWS_ZERO;
			for (k = 0; k <= QAWS_CURVE_JET_ORDER; k++)
				if (channels & (1u << k))
					for (cc = 0; cc < dim; cc++)
						g[c] += cr_get(&jet_adjoint->d[k], cc) * cr_poly_derivative(coeff[cc], local_t, k).t;
			/* add only this coordinate (closed curves may repeat a point) */
			{
				qaws_scalar one[3] = { 0, 0, 0 };
				one[c] = g[c];
				qaws_internal_view_add(pv, idx[i], dim, one);
			}
		}
	return QAWS_STATUS_OK;
}

static unsigned int catmull_rom_describe_fields(qaws_curve const *curve, qaws_field_desc *out, unsigned int capacity)
{
	qaws_catmull_rom_impl const *impl = (qaws_catmull_rom_impl const *)curve->impl;
	if (capacity >= 1)
		out[0] = qaws_internal_field_desc(QAWS_FIELD_POINTS, (qaws_value_type)curve->dimension,
			impl->control_point_count, QAWS_DOMAIN_POSITION, QAWS_CONSTRAINT_NONE, QAWS_DIFF_SMOOTH, CR_DIFF_CAPS);
	return 1;
}

static qaws_status catmull_rom_primal_field(qaws_curve const *curve, qaws_diff_field field,
	qaws_scalar const **out_data, unsigned int *out_count, unsigned int *out_components)
{
	qaws_catmull_rom_impl const *impl = (qaws_catmull_rom_impl const *)curve->impl;
	if (field != QAWS_FIELD_POINTS)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_data = impl->control_points;
	*out_count = impl->control_point_count;
	*out_components = (unsigned int)curve->dimension;
	return QAWS_STATUS_OK;
}

static qaws_status catmull_rom_rebuild(qaws_curve const *curve, qaws_diff_views const *values, qaws_curve **out_curve)
{
	qaws_catmull_rom_impl const *impl = (qaws_catmull_rom_impl const *)curve->impl;
	qaws_catmull_rom_desc d;
	qaws_scalar const *pts = NULL;
	qaws_scalar *owned = NULL;
	qaws_status st = qaws_internal_field_override(values, QAWS_FIELD_POINTS, impl->control_points,
		impl->control_point_count, (unsigned int)curve->dimension, &pts, &owned);
	if (st != QAWS_STATUS_OK)
		return st;
	memset(&d, 0, sizeof(d));
	d.dimension = curve->dimension;
	d.control_points = pts;
	d.control_point_count = impl->control_point_count;
	d.parameterization = impl->parameterization;
	d.closed = impl->closed;
	st = qaws_curve_create_catmull_rom(&d, out_curve);
	qaws_internal_dealloc(NULL, owned);
	return st;
}

static qaws_curve_diff_vtable const catmull_rom_diff_vtable = {
	CR_DIFF_CAPS,
	QAWS_DIFF_PIECEWISE_SMOOTH,
	catmull_rom_describe_fields,
	catmull_rom_primal_field,
	NULL,
	catmull_rom_tangent_span,
	catmull_rom_adjoint_span,
	NULL,
	catmull_rom_rebuild
};

static qaws_curve_vtable const catmull_rom_vtable = {
	catmull_rom_eval_span_2d,
	catmull_rom_eval_span_3d,
	catmull_rom_destroy_impl,
	catmull_rom_is_closed,
	catmull_rom_is_periodic,
	catmull_rom_is_rational,
	catmull_rom_get_continuity,
	&catmull_rom_diff_vtable
};

/* -------------------------------------------------------------------------- */
/*  Creation                                                                  */
/* -------------------------------------------------------------------------- */

qaws_status qaws_curve_create_catmull_rom_ex(
	qaws_catmull_rom_desc const *desc,
	qaws_allocator const *allocator,
	qaws_curve **out_curve)
{
	qaws_catmull_rom_impl *impl = NULL;
	qaws_curve *curve = NULL;
	unsigned int dim_count;
	unsigned int span_count;
	unsigned int N;
	unsigned int i, s;
	qaws_scalar alpha;
	qaws_range range;

	/* --- Validation ---------------------------------------------------- */

	if (!desc || !out_curve)
		return QAWS_STATUS_INVALID_ARGUMENT;

	*out_curve = NULL;

	if (desc->dimension != QAWS_DIMENSION_2D &&
	    desc->dimension != QAWS_DIMENSION_3D)
		return QAWS_STATUS_INVALID_DIMENSION;

	if (!desc->control_points)
		return QAWS_STATUS_INVALID_ARGUMENT;

	dim_count = (desc->dimension == QAWS_DIMENSION_2D) ? 2u : 3u;
	N = desc->control_point_count;

	if (desc->closed) {
		if (N < 3)
			return QAWS_STATUS_INVALID_CONTROL_POINT_COUNT;
		span_count = N;
	} else {
		if (N < 4)
			return QAWS_STATUS_INVALID_CONTROL_POINT_COUNT;
		span_count = N - 3;
	}

	/* --- Allocate implementation --------------------------------------- */

	impl = (qaws_catmull_rom_impl *)qaws_internal_alloc(allocator, (unsigned long)sizeof(qaws_catmull_rom_impl));
	if (!impl)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	memset(impl, 0, sizeof(qaws_catmull_rom_impl));

	impl->control_point_count = N;
	impl->parameterization = desc->parameterization;
	impl->closed = desc->closed ? 1 : 0;

	/* Copy control points */
	impl->control_points = (qaws_scalar *)qaws_internal_alloc(allocator,
		(unsigned long)(sizeof(qaws_scalar) * (size_t)(N * dim_count)));
	if (!impl->control_points) {
		qaws_internal_dealloc(allocator, impl);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	memcpy(impl->control_points, desc->control_points,
	       sizeof(qaws_scalar) * (size_t)(N * dim_count));

	/* Compute knot parameters */
	impl->knot_params = (qaws_scalar *)qaws_internal_alloc(allocator,
		(unsigned long)(sizeof(qaws_scalar) * (size_t)N));
	if (!impl->knot_params) {
		qaws_internal_dealloc(allocator, impl->control_points);
		qaws_internal_dealloc(allocator, impl);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	alpha = compute_alpha(impl->parameterization);
	compute_knot_params(impl->control_points, N, dim_count,
	                    alpha, impl->knot_params);

	/* Compute segment coefficients via prepare functions */
	{
		size_t coeffs_size = (size_t)(span_count * dim_count * 4) * sizeof(qaws_scalar);
		impl->segment_coeffs = (qaws_scalar *)qaws_internal_alloc(allocator,
			(unsigned long)coeffs_size);
		if (!impl->segment_coeffs) {
			qaws_internal_dealloc(allocator, impl->knot_params);
			qaws_internal_dealloc(allocator, impl->control_points);
			qaws_internal_dealloc(allocator, impl);
			return QAWS_STATUS_ALLOCATION_FAILURE;
		}

		if (dim_count == 2) {
			qaws_vec2 *tmp_a, *tmp_b, *tmp_c, *tmp_d;
			unsigned int prep_span_count = 0;
			size_t vec_size = (size_t)span_count * sizeof(qaws_vec2);
			tmp_a = (qaws_vec2 *)qaws_internal_alloc(allocator, (unsigned long)(vec_size * 4));
			if (!tmp_a) {
				qaws_internal_dealloc(allocator, impl->segment_coeffs);
				qaws_internal_dealloc(allocator, impl->knot_params);
				qaws_internal_dealloc(allocator, impl->control_points);
				qaws_internal_dealloc(allocator, impl);
				return QAWS_STATUS_ALLOCATION_FAILURE;
			}
			tmp_b = tmp_a + span_count;
			tmp_c = tmp_b + span_count;
			tmp_d = tmp_c + span_count;

			qaws_catmull_rom_prepare_2d(
				(const qaws_vec2 *)impl->control_points, N,
				impl->parameterization, impl->closed,
				tmp_a, tmp_b, tmp_c, tmp_d, &prep_span_count);

			for (s = 0; s < span_count; s++) {
				impl->segment_coeffs[s * 8 + 0] = tmp_a[s].x;
				impl->segment_coeffs[s * 8 + 1] = tmp_b[s].x;
				impl->segment_coeffs[s * 8 + 2] = tmp_c[s].x;
				impl->segment_coeffs[s * 8 + 3] = tmp_d[s].x;
				impl->segment_coeffs[s * 8 + 4] = tmp_a[s].y;
				impl->segment_coeffs[s * 8 + 5] = tmp_b[s].y;
				impl->segment_coeffs[s * 8 + 6] = tmp_c[s].y;
				impl->segment_coeffs[s * 8 + 7] = tmp_d[s].y;
			}
			qaws_internal_dealloc(allocator, tmp_a);
		} else {
			qaws_vec3 *tmp_a, *tmp_b, *tmp_c, *tmp_d;
			unsigned int prep_span_count = 0;
			size_t vec_size = (size_t)span_count * sizeof(qaws_vec3);
			tmp_a = (qaws_vec3 *)qaws_internal_alloc(allocator, (unsigned long)(vec_size * 4));
			if (!tmp_a) {
				qaws_internal_dealloc(allocator, impl->segment_coeffs);
				qaws_internal_dealloc(allocator, impl->knot_params);
				qaws_internal_dealloc(allocator, impl->control_points);
				qaws_internal_dealloc(allocator, impl);
				return QAWS_STATUS_ALLOCATION_FAILURE;
			}
			tmp_b = tmp_a + span_count;
			tmp_c = tmp_b + span_count;
			tmp_d = tmp_c + span_count;

			qaws_catmull_rom_prepare_3d(
				(const qaws_vec3 *)impl->control_points, N,
				impl->parameterization, impl->closed,
				tmp_a, tmp_b, tmp_c, tmp_d, &prep_span_count);

			for (s = 0; s < span_count; s++) {
				impl->segment_coeffs[s * 12 + 0]  = tmp_a[s].x;
				impl->segment_coeffs[s * 12 + 1]  = tmp_b[s].x;
				impl->segment_coeffs[s * 12 + 2]  = tmp_c[s].x;
				impl->segment_coeffs[s * 12 + 3]  = tmp_d[s].x;
				impl->segment_coeffs[s * 12 + 4]  = tmp_a[s].y;
				impl->segment_coeffs[s * 12 + 5]  = tmp_b[s].y;
				impl->segment_coeffs[s * 12 + 6]  = tmp_c[s].y;
				impl->segment_coeffs[s * 12 + 7]  = tmp_d[s].y;
				impl->segment_coeffs[s * 12 + 8]  = tmp_a[s].z;
				impl->segment_coeffs[s * 12 + 9]  = tmp_b[s].z;
				impl->segment_coeffs[s * 12 + 10] = tmp_c[s].z;
				impl->segment_coeffs[s * 12 + 11] = tmp_d[s].z;
			}
			qaws_internal_dealloc(allocator, tmp_a);
		}
	}

	/* --- Allocate curve ------------------------------------------------ */

	range.min_value = QAWS_ZERO;
	range.max_value = (qaws_scalar)span_count;

	curve = qaws_internal_curve_alloc_ex(
		QAWS_CURVE_KIND_CATMULL_ROM,
		desc->dimension,
		3,            /* degree */
		span_count,
		range,
		&catmull_rom_vtable,
		allocator);
	if (!curve) {
		qaws_internal_dealloc(allocator, impl->segment_coeffs);
		qaws_internal_dealloc(allocator, impl->knot_params);
		qaws_internal_dealloc(allocator, impl->control_points);
		qaws_internal_dealloc(allocator, impl);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	/* Set span boundaries: uniform [0..span_count] */
	for (i = 0; i <= span_count; i++) {
		curve->span_boundaries[i] = (qaws_scalar)i;
	}

	curve->impl = impl;
	*out_curve = curve;

	return QAWS_STATUS_OK;
}

qaws_status qaws_curve_create_catmull_rom(
	qaws_catmull_rom_desc const *desc,
	qaws_curve **out_curve)
{
	return qaws_curve_create_catmull_rom_ex(desc, NULL, out_curve);
}
