#include "qaws_surface_extrusion.h"
#include "qaws_eval.h"
#include "qaws_inspect.h"
#include "internal/qaws_internal_surface.h"
#include "internal/qaws_internal_curve.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "qaws_platform.h"

typedef struct qaws_surface_extrusion_impl
{
	qaws_curve const* profile;
	qaws_range prof_range;
	qaws_vec3 dir;          /* normalized direction */
	qaws_scalar length;
	int profile_is_2d;
} qaws_surface_extrusion_impl;

/* S(u,v) = profile_3d(u) + v * dir * length
   dS/du = profile'(u) * s_prof  (lifted to 3D if 2D)
   dS/dv = dir * length
   dS/duu = profile''(u) * s_prof^2
   dS/dvv = 0
   dS/duv = 0 */
static qaws_status extrusion_surface_eval(
	qaws_surface const* surface,
	qaws_scalar u,
	qaws_scalar v,
	unsigned int eval_flags,
	qaws_surface_eval_result* out_result)
{
	qaws_surface_extrusion_impl const* impl =
		(qaws_surface_extrusion_impl const*)surface->impl;
	qaws_scalar s_prof = impl->prof_range.max_value - impl->prof_range.min_value;
	qaws_scalar t_prof = impl->prof_range.min_value + u * s_prof;
	qaws_vec3 prof_pos;
	qaws_vec3 prof_d1;
	qaws_vec3 prof_d2;

	prof_pos.x = 0; prof_pos.y = 0; prof_pos.z = 0;
	prof_d1.x = 0; prof_d1.y = 0; prof_d1.z = 0;
	prof_d2.x = 0; prof_d2.y = 0; prof_d2.z = 0;

	if (impl->profile_is_2d)
	{
		unsigned int pf = QAWS_EVAL_FLAG_POSITION;
		qaws_eval_result_2d pr;
		qaws_status s;

		if (eval_flags & (QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_DUV
			| QAWS_SURFACE_EVAL_NORMAL))
			pf |= QAWS_EVAL_FLAG_D1;
		if (eval_flags & QAWS_SURFACE_EVAL_DUU)
			pf |= QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2;

		memset(&pr, 0, sizeof(pr));
		s = qaws_curve_evaluate_2d(impl->profile, t_prof, pf, &pr);
		if (s != QAWS_STATUS_OK) return s;

		/* Lift 2D to 3D: (x,y) -> (x, y, 0) */
		prof_pos.x = pr.position.x;
		prof_pos.y = pr.position.y;
		prof_d1.x = pr.d1.x;
		prof_d1.y = pr.d1.y;
		prof_d2.x = pr.d2.x;
		prof_d2.y = pr.d2.y;
	}
	else
	{
		unsigned int pf = QAWS_EVAL_FLAG_POSITION;
		qaws_eval_result_3d pr;
		qaws_status s;

		if (eval_flags & (QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_DUV
			| QAWS_SURFACE_EVAL_NORMAL))
			pf |= QAWS_EVAL_FLAG_D1;
		if (eval_flags & QAWS_SURFACE_EVAL_DUU)
			pf |= QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2;

		memset(&pr, 0, sizeof(pr));
		s = qaws_curve_evaluate_3d(impl->profile, t_prof, pf, &pr);
		if (s != QAWS_STATUS_OK) return s;

		prof_pos = pr.position;
		prof_d1 = pr.d1;
		prof_d2 = pr.d2;
	}

	/* Position: profile_3d(u) + v * dir * length */
	if (eval_flags & QAWS_SURFACE_EVAL_POSITION)
	{
		out_result->position.x = prof_pos.x + v * impl->dir.x * impl->length;
		out_result->position.y = prof_pos.y + v * impl->dir.y * impl->length;
		out_result->position.z = prof_pos.z + v * impl->dir.z * impl->length;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_POSITION;
	}

	/* dS/du = profile'(u) * s_prof */
	if (eval_flags & QAWS_SURFACE_EVAL_DU)
	{
		out_result->du.x = prof_d1.x * s_prof;
		out_result->du.y = prof_d1.y * s_prof;
		out_result->du.z = prof_d1.z * s_prof;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_DU;
	}

	/* dS/dv = dir * length */
	if (eval_flags & QAWS_SURFACE_EVAL_DV)
	{
		out_result->dv.x = impl->dir.x * impl->length;
		out_result->dv.y = impl->dir.y * impl->length;
		out_result->dv.z = impl->dir.z * impl->length;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_DV;
	}

	/* dS/duu = profile''(u) * s_prof^2 */
	if (eval_flags & QAWS_SURFACE_EVAL_DUU)
	{
		qaws_scalar sp2 = s_prof * s_prof;
		out_result->duu.x = prof_d2.x * sp2;
		out_result->duu.y = prof_d2.y * sp2;
		out_result->duu.z = prof_d2.z * sp2;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_DUU;
	}

	/* dS/dvv = 0 */
	if (eval_flags & QAWS_SURFACE_EVAL_DVV)
	{
		out_result->dvv.x = 0; out_result->dvv.y = 0; out_result->dvv.z = 0;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_DVV;
	}

	/* dS/duv = 0 */
	if (eval_flags & QAWS_SURFACE_EVAL_DUV)
	{
		out_result->duv.x = 0; out_result->duv.y = 0; out_result->duv.z = 0;
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

static void extrusion_surface_destroy(void* impl, qaws_allocator const* allocator)
{
	qaws_internal_dealloc(allocator, impl);
}

static int extrusion_surface_is_rational(qaws_surface const* s)
{
	(void)s;
	return 0;
}

static qaws_surface_vtable const extrusion_surface_vtable = {
	extrusion_surface_eval,
	extrusion_surface_destroy,
	extrusion_surface_is_rational
};

qaws_status qaws_surface_create_extrusion(
	qaws_surface_extrusion_desc const* desc,
	qaws_surface** out_surface)
{
	qaws_surface* surface;
	qaws_surface_extrusion_impl* impl;
	qaws_range u_range, v_range;
	qaws_scalar dir_len;
	qaws_dimension dim;

	if (!desc || !out_surface) return QAWS_STATUS_INVALID_ARGUMENT;
	if (!desc->profile) return QAWS_STATUS_INVALID_ARGUMENT;

	dim = qaws_curve_get_dimension(desc->profile);
	if (dim != QAWS_DIMENSION_2D && dim != QAWS_DIMENSION_3D)
		return QAWS_STATUS_INVALID_DIMENSION;

	dir_len = QAWS_SQRT(desc->direction.x * desc->direction.x
		+ desc->direction.y * desc->direction.y
		+ desc->direction.z * desc->direction.z);
	if (dir_len < QAWS_LITERAL(1e-12))
		return QAWS_STATUS_INVALID_ARGUMENT;

	u_range.min_value = 0; u_range.max_value = 1;
	v_range.min_value = 0; v_range.max_value = 1;

	surface = qaws_internal_surface_alloc(
		QAWS_SURFACE_KIND_EXTRUSION,
		0, 0, u_range, v_range,
		&extrusion_surface_vtable);
	if (!surface) return QAWS_STATUS_ALLOCATION_FAILURE;

	impl = (qaws_surface_extrusion_impl*)malloc(sizeof(qaws_surface_extrusion_impl));
	if (!impl)
	{
		qaws_internal_surface_free(surface);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	impl->profile = desc->profile;
	impl->prof_range = qaws_curve_get_parameter_range(desc->profile);
	impl->profile_is_2d = (dim == QAWS_DIMENSION_2D) ? 1 : 0;

	/* Normalize direction */
	impl->dir.x = desc->direction.x / dir_len;
	impl->dir.y = desc->direction.y / dir_len;
	impl->dir.z = desc->direction.z / dir_len;

	/* Length: if desc->length > 0 use it, otherwise use the magnitude of direction */
	impl->length = (desc->length > QAWS_ZERO) ? desc->length : dir_len;

	surface->impl = impl;
	*out_surface = surface;
	return QAWS_STATUS_OK;
}
