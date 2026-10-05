#include "qaws_surface_revolution.h"
#include "qaws_eval.h"
#include "qaws_inspect.h"
#include "internal/qaws_internal_surface.h"
#include "internal/qaws_internal_curve.h"
#include "internal/qaws_internal_diff.h"
#include "core/qaws_dual_core.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "qaws_platform.h"

typedef struct qaws_surface_revolution_impl
{
	qaws_curve const* profile;
	qaws_range prof_range;
	qaws_vec3 origin;
	qaws_vec3 X_axis;   /* perpendicular to rotation axis */
	qaws_vec3 Y_axis;   /* perpendicular to rotation axis and X_axis */
	qaws_vec3 Z_axis;   /* normalized rotation axis direction */
	qaws_scalar angle;
} qaws_surface_revolution_impl;

/* S(u,v) = origin + r(v)*cos(theta)*X + r(v)*sin(theta)*Y + h(v)*Z
   where theta = u * angle, r = profile_x(v), h = profile_y(v).

   dS/du = angle * r * (-sin(theta)*X + cos(theta)*Y)
   dS/dv = (dr/dv*cos(theta)*X + dr/dv*sin(theta)*Y + dh/dv*Z) * s_prof

   d2S/du2 = -angle^2 * r * (cos(theta)*X + sin(theta)*Y)
   d2S/dv2 = (d2r/dv2*cos(theta)*X + d2r/dv2*sin(theta)*Y + d2h/dv2*Z) * s_prof^2
   d2S/dudv = angle * dr/dv * (-sin(theta)*X + cos(theta)*Y) * s_prof */
static qaws_status revolution_surface_eval(
	qaws_surface const* surface,
	qaws_scalar u,
	qaws_scalar v,
	unsigned int eval_flags,
	qaws_surface_eval_result* out_result)
{
	qaws_surface_revolution_impl const* impl =
		(qaws_surface_revolution_impl const*)surface->impl;
	qaws_scalar s_prof = impl->prof_range.max_value - impl->prof_range.min_value;
	qaws_scalar t_prof = impl->prof_range.min_value + v * s_prof;
	qaws_scalar theta = u * impl->angle;
	qaws_scalar cos_t = QAWS_COS(theta);
	qaws_scalar sin_t = QAWS_SIN(theta);

	unsigned int pf = QAWS_EVAL_FLAG_POSITION;
	qaws_eval_result_2d pr;
	qaws_status s;
	qaws_scalar r, h;

	if (eval_flags & (QAWS_SURFACE_EVAL_DV | QAWS_SURFACE_EVAL_DUV))
		pf |= QAWS_EVAL_FLAG_D1;
	if (eval_flags & QAWS_SURFACE_EVAL_DVV)
		pf |= QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2;

	memset(&pr, 0, sizeof(pr));
	s = qaws_curve_evaluate_2d(impl->profile, t_prof, pf, &pr);
	if (s != QAWS_STATUS_OK) return s;

	r = pr.position.x;
	h = pr.position.y;

	/* Position: origin + r*cos(theta)*X + r*sin(theta)*Y + h*Z */
	if (eval_flags & QAWS_SURFACE_EVAL_POSITION)
	{
		out_result->position.x = impl->origin.x
			+ r * cos_t * impl->X_axis.x
			+ r * sin_t * impl->Y_axis.x
			+ h * impl->Z_axis.x;
		out_result->position.y = impl->origin.y
			+ r * cos_t * impl->X_axis.y
			+ r * sin_t * impl->Y_axis.y
			+ h * impl->Z_axis.y;
		out_result->position.z = impl->origin.z
			+ r * cos_t * impl->X_axis.z
			+ r * sin_t * impl->Y_axis.z
			+ h * impl->Z_axis.z;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_POSITION;
	}

	/* dS/du = angle * r * (-sin(theta)*X + cos(theta)*Y) */
	if (eval_flags & QAWS_SURFACE_EVAL_DU)
	{
		qaws_scalar ar = impl->angle * r;
		out_result->du.x = ar * (-sin_t * impl->X_axis.x + cos_t * impl->Y_axis.x);
		out_result->du.y = ar * (-sin_t * impl->X_axis.y + cos_t * impl->Y_axis.y);
		out_result->du.z = ar * (-sin_t * impl->X_axis.z + cos_t * impl->Y_axis.z);
		out_result->valid_flags |= QAWS_SURFACE_EVAL_DU;
	}

	/* dS/dv = (dr/dv*cos(theta)*X + dr/dv*sin(theta)*Y + dh/dv*Z) * s_prof */
	if (eval_flags & QAWS_SURFACE_EVAL_DV)
	{
		qaws_scalar dr = pr.d1.x * s_prof;
		qaws_scalar dh = pr.d1.y * s_prof;
		out_result->dv.x = dr * cos_t * impl->X_axis.x
			+ dr * sin_t * impl->Y_axis.x
			+ dh * impl->Z_axis.x;
		out_result->dv.y = dr * cos_t * impl->X_axis.y
			+ dr * sin_t * impl->Y_axis.y
			+ dh * impl->Z_axis.y;
		out_result->dv.z = dr * cos_t * impl->X_axis.z
			+ dr * sin_t * impl->Y_axis.z
			+ dh * impl->Z_axis.z;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_DV;
	}

	/* d2S/du2 = -angle^2 * r * (cos(theta)*X + sin(theta)*Y) */
	if (eval_flags & QAWS_SURFACE_EVAL_DUU)
	{
		qaws_scalar a2r = -(impl->angle * impl->angle) * r;
		out_result->duu.x = a2r * (cos_t * impl->X_axis.x + sin_t * impl->Y_axis.x);
		out_result->duu.y = a2r * (cos_t * impl->X_axis.y + sin_t * impl->Y_axis.y);
		out_result->duu.z = a2r * (cos_t * impl->X_axis.z + sin_t * impl->Y_axis.z);
		out_result->valid_flags |= QAWS_SURFACE_EVAL_DUU;
	}

	/* d2S/dv2 = (d2r/dv2*cos(theta)*X + d2r/dv2*sin(theta)*Y + d2h/dv2*Z) * s_prof^2 */
	if (eval_flags & QAWS_SURFACE_EVAL_DVV)
	{
		qaws_scalar sp2 = s_prof * s_prof;
		qaws_scalar ddr = pr.d2.x * sp2;
		qaws_scalar ddh = pr.d2.y * sp2;
		out_result->dvv.x = ddr * cos_t * impl->X_axis.x
			+ ddr * sin_t * impl->Y_axis.x
			+ ddh * impl->Z_axis.x;
		out_result->dvv.y = ddr * cos_t * impl->X_axis.y
			+ ddr * sin_t * impl->Y_axis.y
			+ ddh * impl->Z_axis.y;
		out_result->dvv.z = ddr * cos_t * impl->X_axis.z
			+ ddr * sin_t * impl->Y_axis.z
			+ ddh * impl->Z_axis.z;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_DVV;
	}

	/* d2S/dudv = angle * dr/dv * (-sin(theta)*X + cos(theta)*Y) * s_prof */
	if (eval_flags & QAWS_SURFACE_EVAL_DUV)
	{
		qaws_scalar a_dr = impl->angle * pr.d1.x * s_prof;
		out_result->duv.x = a_dr * (-sin_t * impl->X_axis.x + cos_t * impl->Y_axis.x);
		out_result->duv.y = a_dr * (-sin_t * impl->X_axis.y + cos_t * impl->Y_axis.y);
		out_result->duv.z = a_dr * (-sin_t * impl->X_axis.z + cos_t * impl->Y_axis.z);
		out_result->valid_flags |= QAWS_SURFACE_EVAL_DUV;
	}

	/* Normal */
	if (eval_flags & QAWS_SURFACE_EVAL_NORMAL)
	{
		qaws_internal_surface_normal(out_result->du, out_result->dv, &out_result->normal);
		out_result->valid_flags |= QAWS_SURFACE_EVAL_NORMAL;
	}

	return QAWS_STATUS_OK;
}

static void revolution_surface_destroy(void* impl, qaws_allocator const* allocator)
{
	qaws_internal_dealloc(allocator, impl);
}

static int revolution_surface_is_rational(qaws_surface const* s)
{
	(void)s;
	return 0;
}

/* -------------------------------------------------------------------------- */
/*  Differential rules                                                        */
/*                                                                            */
/*  S(u,v) = O + r(t) e(theta) + h(t) Z,  theta = u alpha,  t = t0 + v s,     */
/*  e = cos X + sin Y.  Own fields: center (O) and angle_end (alpha).         */
/*  Child 0: the 2D profile (r, h). The frame X, Y, Z is fixed.               */
/* -------------------------------------------------------------------------- */

#define REVOLUTION_DIFF_CAPS (QAWS_CAP_TANGENT | QAWS_CAP_ADJOINT | QAWS_CAP_TANGENT2)

static unsigned char const g_revolution_jet_a[QAWS_SURFACE_JET_COUNT] = { 0, 1, 0, 2, 1, 0, 3, 2, 1, 0 };
static unsigned char const g_revolution_jet_b[QAWS_SURFACE_JET_COUNT] = { 0, 0, 1, 0, 1, 2, 0, 1, 2, 3 };

static unsigned int revolution_describe_fields(qaws_surface const* surface,
	qaws_field_desc* out, unsigned int capacity)
{
	(void)surface;
	if (capacity >= 1)
		out[0] = qaws_internal_field_desc(QAWS_FIELD_CENTER, QAWS_VALUE_VEC3, 1,
			QAWS_DOMAIN_POSITION, QAWS_CONSTRAINT_NONE, QAWS_DIFF_SMOOTH, REVOLUTION_DIFF_CAPS);
	if (capacity >= 2)
		out[1] = qaws_internal_field_desc(QAWS_FIELD_ANGLE_END, QAWS_VALUE_SCALAR, 1,
			QAWS_DOMAIN_ANGLE, QAWS_CONSTRAINT_NONE, QAWS_DIFF_SMOOTH, REVOLUTION_DIFF_CAPS);
	return 2;
}

static qaws_status revolution_primal_field(qaws_surface const* surface, qaws_diff_field field,
	qaws_scalar const** out_data, unsigned int* out_count, unsigned int* out_components)
{
	qaws_surface_revolution_impl const* impl = (qaws_surface_revolution_impl const*)surface->impl;
	*out_count = 1;
	if (field == QAWS_FIELD_CENTER)
	{
		*out_data = &impl->origin.x;
		*out_components = 3;
		return QAWS_STATUS_OK;
	}
	if (field == QAWS_FIELD_ANGLE_END)
	{
		*out_data = &impl->angle;
		*out_components = 1;
		return QAWS_STATUS_OK;
	}
	return QAWS_STATUS_INVALID_ARGUMENT;
}

static unsigned int revolution_children(qaws_surface const* surface, qaws_diff_child* out, unsigned int capacity)
{
	qaws_surface_revolution_impl const* impl = (qaws_surface_revolution_impl const*)surface->impl;
	if (capacity >= 1)
	{
		out[0].curve = impl->profile;
		out[0].surface = NULL;
	}
	return 1;
}

static qaws_dual1 revolution_cos(qaws_dual1 x)
{
	qaws_scalar c = QAWS_COS(x.v), s = QAWS_SIN(x.v);
	return qaws_dual1_make(c, -s * x.t, -c * x.t * x.t - s * x.tt);
}

static qaws_dual1 revolution_sin(qaws_dual1 x)
{
	qaws_scalar c = QAWS_COS(x.v), s = QAWS_SIN(x.v);
	return qaws_dual1_make(s, c * x.t, -s * x.t * x.t + c * x.tt);
}

/* All ten partials from dual inputs; r[b], h[b] are d^b/dt^b of the profile. */
static void revolution_jets(
	qaws_surface_revolution_impl const* impl,
	qaws_dual1 u, qaws_dual1 alpha, qaws_dual3 origin,
	qaws_dual1 const* r, qaws_dual1 const* h, qaws_scalar s,
	qaws_dual3* S)
{
	qaws_dual1 theta = qaws_dual1_mul(u, alpha);
	qaws_dual1 c = revolution_cos(theta), sn = revolution_sin(theta);
	qaws_dual3 X = qaws_dual3_const(impl->X_axis), Y = qaws_dual3_const(impl->Y_axis), Z = qaws_dual3_const(impl->Z_axis);
	qaws_dual3 e[4];
	qaws_dual1 ap[4];
	qaws_scalar sp[4];
	unsigned int k;

	e[0] = qaws_dual3_add(qaws_dual3_scale(X, c), qaws_dual3_scale(Y, sn));
	e[1] = qaws_dual3_sub(qaws_dual3_scale(Y, c), qaws_dual3_scale(X, sn));
	e[2] = qaws_dual3_mul_const(e[0], -QAWS_ONE);
	e[3] = qaws_dual3_mul_const(e[1], -QAWS_ONE);
	ap[0] = qaws_dual1_const(QAWS_ONE);
	sp[0] = QAWS_ONE;
	for (k = 1; k < 4; k++)
	{
		ap[k] = qaws_dual1_mul(ap[k - 1], alpha);
		sp[k] = sp[k - 1] * s;
	}

	for (k = 0; k < QAWS_SURFACE_JET_COUNT; k++)
	{
		unsigned int a = g_revolution_jet_a[k], b = g_revolution_jet_b[k];
		qaws_dual1 coef = qaws_dual1_mul(ap[a], r[b]);
		coef = qaws_dual1_make(coef.v * sp[b], coef.t * sp[b], coef.tt * sp[b]);
		S[k] = qaws_dual3_scale(e[a], coef);
		if (a == 0)
			S[k] = qaws_dual3_add(S[k], qaws_dual3_scale(Z, qaws_dual1_make(h[b].v * sp[b], h[b].t * sp[b], h[b].tt * sp[b])));
		if (a == 0 && b == 0)
			S[k] = qaws_dual3_add(S[k], origin);
	}
}

static qaws_scalar revolution_view_scalar(qaws_diff_views const* views, qaws_diff_field field)
{
	qaws_field_view const* fv = views ? qaws_diff_views_find(views, field) : NULL;
	qaws_scalar x = QAWS_ZERO;
	if (fv)
		qaws_internal_view_read(fv, 0, 1, &x);
	return x;
}

static qaws_status revolution_tangent(
	qaws_diff_context const* ctx, qaws_surface const* surface,
	qaws_scalar u, qaws_scalar v, qaws_scalar u_dot, qaws_scalar v_dot,
	unsigned int channels, qaws_diff_views const* views,
	qaws_surface_jet* primal, qaws_surface_jet* tangent, qaws_surface_jet* tangent2)
{
	qaws_surface_revolution_impl const* impl = (qaws_surface_revolution_impl const*)surface->impl;
	qaws_scalar s = impl->prof_range.max_value - impl->prof_range.min_value;
	qaws_curve_jet_3d cp, ct, ctt;
	qaws_dual1 r[4], h[4];
	qaws_dual3 S[QAWS_SURFACE_JET_COUNT], origin;
	qaws_scalar o_dot[3] = { 0, 0, 0 };
	qaws_field_view const* ov = views ? qaws_diff_views_find(views, QAWS_FIELD_CENTER) : NULL;
	unsigned int k;
	qaws_status st;

	st = qaws_internal_curve_tangent_any(ctx, impl->profile, impl->prof_range.min_value + v * s, s * v_dot,
		0xFu, qaws_internal_child_views(views, 0), &cp, &ct, tangent2 ? &ctt : NULL);
	if (st != QAWS_STATUS_OK)
		return st;
	for (k = 0; k < 4; k++)
	{
		r[k] = qaws_dual1_make(cp.d[k].x, ct.d[k].x, tangent2 ? ctt.d[k].x : QAWS_ZERO);
		h[k] = qaws_dual1_make(cp.d[k].y, ct.d[k].y, tangent2 ? ctt.d[k].y : QAWS_ZERO);
	}
	if (ov)
		qaws_internal_view_read(ov, 0, 3, o_dot);
	origin = qaws_dual3_make(impl->origin, qaws_v3(o_dot[0], o_dot[1], o_dot[2]), qaws_v3_zero());

	revolution_jets(impl, qaws_dual1_make(u, u_dot, QAWS_ZERO),
		qaws_dual1_make(impl->angle, revolution_view_scalar(views, QAWS_FIELD_ANGLE_END), QAWS_ZERO),
		origin, r, h, s, S);

	memset(primal, 0, sizeof(*primal));
	memset(tangent, 0, sizeof(*tangent));
	if (tangent2)
		memset(tangent2, 0, sizeof(*tangent2));
	for (k = 0; k < QAWS_SURFACE_JET_COUNT; k++)
	{
		if (!(channels & (1u << k)))
			continue;
		primal->d[k] = S[k].v;
		tangent->d[k] = S[k].t;
		if (tangent2)
			tangent2->d[k] = S[k].tt;
	}
	primal->channels = tangent->channels = channels;
	if (tangent2)
		tangent2->channels = channels;
	return QAWS_STATUS_OK;
}

/* Exact pullback by forward seeds: 8 profile jet scalars, the angle, the
   three origin components and u. */
static qaws_status revolution_adjoint(
	qaws_diff_context const* ctx, qaws_surface const* surface,
	qaws_scalar u, qaws_scalar v, unsigned int channels,
	qaws_surface_jet const* ybar, qaws_diff_views* views,
	qaws_scalar* u_adjoint, qaws_scalar* v_adjoint)
{
	qaws_surface_revolution_impl const* impl = (qaws_surface_revolution_impl const*)surface->impl;
	qaws_scalar s = impl->prof_range.max_value - impl->prof_range.min_value;
	qaws_scalar t = impl->prof_range.min_value + v * s;
	qaws_field_view* ov = views ? qaws_diff_views_find(views, QAWS_FIELD_CENTER) : NULL;
	qaws_field_view* av = views ? qaws_diff_views_find(views, QAWS_FIELD_ANGLE_END) : NULL;
	qaws_curve_jet_3d cp, ct, cbar;
	qaws_scalar bars[13];
	qaws_scalar tbar = QAWS_ZERO;
	unsigned int seed, k;
	qaws_status st;

	st = qaws_internal_curve_tangent_any(ctx, impl->profile, t, QAWS_ZERO, 0xFu, NULL, &cp, &ct, NULL);
	if (st != QAWS_STATUS_OK)
		return st;

	for (seed = 0; seed < 13; seed++)
	{
		qaws_dual1 r[4], h[4];
		qaws_dual3 S[QAWS_SURFACE_JET_COUNT];
		qaws_vec3 od = qaws_v3_zero();
		qaws_scalar acc = QAWS_ZERO;
		for (k = 0; k < 4; k++)
		{
			r[k] = qaws_dual1_make(cp.d[k].x, seed == k ? QAWS_ONE : QAWS_ZERO, QAWS_ZERO);
			h[k] = qaws_dual1_make(cp.d[k].y, seed == 4 + k ? QAWS_ONE : QAWS_ZERO, QAWS_ZERO);
		}
		if (seed == 9) od.x = QAWS_ONE;
		if (seed == 10) od.y = QAWS_ONE;
		if (seed == 11) od.z = QAWS_ONE;
		revolution_jets(impl,
			qaws_dual1_make(u, seed == 12 ? QAWS_ONE : QAWS_ZERO, QAWS_ZERO),
			qaws_dual1_make(impl->angle, seed == 8 ? QAWS_ONE : QAWS_ZERO, QAWS_ZERO),
			qaws_dual3_make(impl->origin, od, qaws_v3_zero()), r, h, s, S);
		for (k = 0; k < QAWS_SURFACE_JET_COUNT; k++)
			if (channels & (1u << k))
				acc += qaws_v3_dot(ybar->d[k], S[k].t);
		bars[seed] = acc;
	}

	memset(&cbar, 0, sizeof(cbar));
	for (k = 0; k < 4; k++)
	{
		cbar.d[k].x = bars[k];
		cbar.d[k].y = bars[4 + k];
	}
	cbar.channels = 0xFu;
	st = qaws_internal_curve_adjoint_any(ctx, impl->profile, t, 0xFu, &cbar,
		(qaws_diff_views*)qaws_internal_child_views(views, 0), v_adjoint ? &tbar : NULL);
	if (st != QAWS_STATUS_OK)
		return st;
	if (v_adjoint)
		*v_adjoint += s * tbar;
	if (u_adjoint)
		*u_adjoint += bars[12];
	if (av)
		qaws_internal_view_add(av, 0, 1, &bars[8]);
	if (ov)
		qaws_internal_view_add(ov, 0, 3, &bars[9]);
	return QAWS_STATUS_OK;
}

static qaws_surface_diff_vtable const revolution_surface_diff_vtable = {
	REVOLUTION_DIFF_CAPS,
	QAWS_DIFF_PIECEWISE_SMOOTH,
	revolution_describe_fields,
	revolution_primal_field,
	NULL,
	NULL,
	revolution_tangent,
	revolution_adjoint,
	revolution_children
};

static qaws_surface_vtable const revolution_surface_vtable = {
	revolution_surface_eval,
	revolution_surface_destroy,
	revolution_surface_is_rational,
	&revolution_surface_diff_vtable
};

qaws_status qaws_surface_create_revolution(
	qaws_surface_revolution_desc const* desc,
	qaws_surface** out_surface)
{
	qaws_surface* surface;
	qaws_surface_revolution_impl* impl;
	qaws_range u_range, v_range;
	qaws_scalar axis_len;
	qaws_scalar angle;
	qaws_vec3 Z, X, Y;

	if (!desc || !out_surface) return QAWS_STATUS_INVALID_ARGUMENT;
	if (!desc->profile) return QAWS_STATUS_INVALID_ARGUMENT;
	if (qaws_curve_get_dimension(desc->profile) != QAWS_DIMENSION_2D)
		return QAWS_STATUS_INVALID_DIMENSION;

	axis_len = QAWS_SQRT(desc->axis_direction.x * desc->axis_direction.x
		+ desc->axis_direction.y * desc->axis_direction.y
		+ desc->axis_direction.z * desc->axis_direction.z);
	if (axis_len < QAWS_LITERAL(1e-12))
		return QAWS_STATUS_INVALID_ARGUMENT;

	/* Normalize axis direction -> Z_axis */
	Z.x = desc->axis_direction.x / axis_len;
	Z.y = desc->axis_direction.y / axis_len;
	Z.z = desc->axis_direction.z / axis_len;

	/* Build perpendicular X_axis: pick the coordinate axis least aligned with Z */
	{
		qaws_scalar ax = QAWS_FABS(Z.x);
		qaws_scalar ay = QAWS_FABS(Z.y);
		qaws_scalar az = QAWS_FABS(Z.z);
		qaws_vec3 up;
		qaws_scalar dot, xlen;

		if (ax <= ay && ax <= az)
		{
			up.x = QAWS_ONE; up.y = QAWS_ZERO; up.z = QAWS_ZERO;
		}
		else if (ay <= az)
		{
			up.x = QAWS_ZERO; up.y = QAWS_ONE; up.z = QAWS_ZERO;
		}
		else
		{
			up.x = QAWS_ZERO; up.y = QAWS_ZERO; up.z = QAWS_ONE;
		}

		/* X = up - dot(up, Z) * Z, then normalize */
		dot = up.x * Z.x + up.y * Z.y + up.z * Z.z;
		X.x = up.x - dot * Z.x;
		X.y = up.y - dot * Z.y;
		X.z = up.z - dot * Z.z;
		xlen = QAWS_SQRT(X.x * X.x + X.y * X.y + X.z * X.z);
		X.x /= xlen; X.y /= xlen; X.z /= xlen;
	}

	/* Y = cross(Z, X) */
	Y.x = Z.y * X.z - Z.z * X.y;
	Y.y = Z.z * X.x - Z.x * X.z;
	Y.z = Z.x * X.y - Z.y * X.x;

	/* Default angle: 2*PI for full revolution */
	angle = desc->angle;
	if (QAWS_FABS(angle) < QAWS_LITERAL(1e-12))
		angle = QAWS_LITERAL(2.0) * QAWS_LITERAL(3.14159265358979323846);

	u_range.min_value = 0; u_range.max_value = 1;
	v_range.min_value = 0; v_range.max_value = 1;

	surface = qaws_internal_surface_alloc(
		QAWS_SURFACE_KIND_REVOLUTION,
		0, 0, u_range, v_range,
		&revolution_surface_vtable);
	if (!surface) return QAWS_STATUS_ALLOCATION_FAILURE;

	impl = (qaws_surface_revolution_impl*)malloc(sizeof(qaws_surface_revolution_impl));
	if (!impl)
	{
		qaws_internal_surface_free(surface);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	impl->profile = desc->profile;
	impl->prof_range = qaws_curve_get_parameter_range(desc->profile);
	impl->origin = desc->axis_origin;
	impl->X_axis = X;
	impl->Y_axis = Y;
	impl->Z_axis = Z;
	impl->angle = angle;

	surface->impl = impl;
	*out_surface = surface;
	return QAWS_STATUS_OK;
}
