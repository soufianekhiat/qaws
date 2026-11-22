/*
 * qaws_transform: points mapped through a placement (qaws_transform.h).
 */
#include "qaws_platform.h"
#include "qaws_transform.h"
#include "internal/qaws_internal_transform.h"

static qaws_status tf_map(qaws_scalar* out, qaws_transform const* x, qaws_scalar const* p, unsigned int n, int inverse)
{
	qaws_tf f;
	unsigned int i;
	if (!x || (n && (!p || !out)))
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (!qaws_tf_init(&f, x))
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (i = 0; i < n; i++)
	{
		double a[3] = { p[3 * i], p[3 * i + 1], p[3 * i + 2] }, b[3];
		if (!isfinite(a[0]) || !isfinite(a[1]) || !isfinite(a[2]))
			return QAWS_STATUS_INVALID_ARGUMENT;
		if (inverse)
			qaws_tf_to_frame(&f, a, b);
		else
			qaws_tf_to_world(&f, a, b);
		out[3 * i] = (qaws_scalar)b[0], out[3 * i + 1] = (qaws_scalar)b[1], out[3 * i + 2] = (qaws_scalar)b[2];
	}
	return QAWS_STATUS_OK;
}

qaws_status qaws_transform_points(qaws_scalar* out_points, qaws_transform const* transform, qaws_scalar const* points, unsigned int count)
{
	return tf_map(out_points, transform, points, count, 0);
}

qaws_status qaws_transform_points_inverse(qaws_scalar* out_points, qaws_transform const* transform, qaws_scalar const* points,
	unsigned int count)
{
	return tf_map(out_points, transform, points, count, 1);
}
