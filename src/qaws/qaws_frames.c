/*
 * Rotation-minimizing frames along curves (qaws_inspect.h): the double-reflection step of core/qaws_frame_core.h from sample to
 * sample, the normal kept perpendicular to the tangent and of unit length.
 */

#include "qaws_inspect.h"
#include "qaws_eval.h"
#include "internal/qaws_internal_types.h"
#include "core/qaws_frame_core.h"
#include <math.h>

static int fr_unit(qaws_vec3* v)
{
	qaws_scalar l = QAWS_SQRT(qaws_frame_dot(*v, *v));
	if (!(l > QAWS_ZERO))
		return 0;
	v->x /= l;
	v->y /= l;
	v->z /= l;
	return 1;
}

/* r made perpendicular to the unit t and of unit length; 0 when it lies along t */
static int fr_perpendicular(qaws_vec3* r, qaws_vec3 t)
{
	qaws_scalar d = qaws_frame_dot(*r, t);
	qaws_scalar l0 = QAWS_SQRT(qaws_frame_dot(*r, *r));
	r->x -= d * t.x;
	r->y -= d * t.y;
	r->z -= d * t.z;
	return QAWS_SQRT(qaws_frame_dot(*r, *r)) > QAWS_LITERAL(1e-9) * l0 && fr_unit(r);
}

qaws_status qaws_curve_compute_rotation_minimizing_frames_3d(qaws_curve const* curve, qaws_scalar const* parameters, unsigned int count,
	qaws_vec3 const* initial_normal, qaws_scalar* out_frames)
{
	qaws_vec3 x0 = { 0, 0, 0 }, t0 = { 0, 0, 0 }, r0 = { 0, 0, 0 };
	unsigned int i;
	if (!curve || !parameters || !out_frames || count == 0)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (curve->dimension != QAWS_DIMENSION_3D)
		return QAWS_STATUS_INVALID_DIMENSION;
	for (i = 1; i < count; i++)
		if (!(parameters[i] >= parameters[i - 1]))
			return QAWS_STATUS_INVALID_ARGUMENT;
	for (i = 0; i < count; i++)
	{
		qaws_eval_result_3d e;
		qaws_vec3 t, r, b;
		qaws_status st = qaws_curve_evaluate_3d(&e, curve, parameters[i], QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1);
		if (st != QAWS_STATUS_OK)
			return st;
		t = e.d1;
		if (!fr_unit(&t))
			return QAWS_STATUS_DEGENERATE_CURVE;
		if (i == 0)
		{
			r.x = r.y = r.z = QAWS_ZERO;
			if (initial_normal)
				r = *initial_normal;
			if (!initial_normal || !fr_perpendicular(&r, t))
			{
				/* the axis least along the tangent, made perpendicular */
				qaws_scalar ax = QAWS_FABS(t.x), ay = QAWS_FABS(t.y), az = QAWS_FABS(t.z);
				r.x = ax <= ay && ax <= az ? QAWS_ONE : QAWS_ZERO;
				r.y = r.x == QAWS_ZERO && ay <= az ? QAWS_ONE : QAWS_ZERO;
				r.z = r.x == QAWS_ZERO && r.y == QAWS_ZERO ? QAWS_ONE : QAWS_ZERO;
				fr_perpendicular(&r, t);
			}
		}
		else
		{
			r = qaws_double_reflection(x0, t0, r0, e.position, t);
			if (!fr_perpendicular(&r, t))
				return QAWS_STATUS_NUMERICAL_FAILURE;
		}
		b.x = t.y * r.z - t.z * r.y;
		b.y = t.z * r.x - t.x * r.z;
		b.z = t.x * r.y - t.y * r.x;
		out_frames[9 * i + 0] = t.x; out_frames[9 * i + 1] = t.y; out_frames[9 * i + 2] = t.z;
		out_frames[9 * i + 3] = r.x; out_frames[9 * i + 4] = r.y; out_frames[9 * i + 5] = r.z;
		out_frames[9 * i + 6] = b.x; out_frames[9 * i + 7] = b.y; out_frames[9 * i + 8] = b.z;
		x0 = e.position;
		t0 = t;
		r0 = r;
	}
	return QAWS_STATUS_OK;
}
