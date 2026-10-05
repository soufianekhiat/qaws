#include "qaws_surface_revolution.h"
#include "qaws_eval.h"
#include "qaws_inspect.h"
#include "internal/qaws_internal_surface.h"
#include "internal/qaws_internal_curve.h"
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

static qaws_surface_vtable const revolution_surface_vtable = {
	revolution_surface_eval,
	revolution_surface_destroy,
	revolution_surface_is_rational,
	NULL /* diff */
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
