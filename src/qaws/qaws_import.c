#include "qaws_import.h"
#include "internal/qaws_internal_memory.h"
#include "qaws_catmull_rom.h"
#include "qaws_trajectory.h"
#include <stdlib.h>
#include <string.h>

/* -------------------------------------------------------------------------- */
/*  Polyline import                                                           */
/* -------------------------------------------------------------------------- */

qaws_status qaws_polyline_import_catmull_rom(
	qaws_dimension dimension,
	void const* points,
	unsigned int point_count,
	qaws_parameterization parameterization,
	int closed,
	qaws_curve** out_curve)
{
	qaws_catmull_rom_desc desc;

	if (!points || !out_curve)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (point_count < 2)
		return QAWS_STATUS_INVALID_ARGUMENT;

	desc.dimension = dimension;
	desc.control_points = points;
	desc.control_point_count = point_count;
	desc.parameterization = parameterization;
	desc.closed = closed;

	return qaws_curve_create_catmull_rom(&desc, out_curve);
}

qaws_status qaws_polyline_import_trajectory(
	qaws_dimension dimension,
	void const* points,
	unsigned int point_count,
	int closed,
	qaws_curve** out_curve)
{
	qaws_trajectory_desc desc;
	qaws_scalar times_stack[64];
	qaws_scalar* times = NULL;
	qaws_status s;
	unsigned int i;

	if (!points || !out_curve)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (point_count < 2)
		return QAWS_STATUS_INVALID_ARGUMENT;

	/* Generate uniform key times [0, 1, 2, ..., n-1] */
	if (point_count <= 64)
		times = times_stack;
	else
	{
		times = (qaws_scalar*)qaws_internal_malloc(qaws_internal_mul2(point_count, sizeof(qaws_scalar)));
		if (!times)
			return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (i = 0; i < point_count; i++)
		times[i] = (qaws_scalar)i;

	memset(&desc, 0, sizeof(desc));
	desc.dimension = dimension;
	desc.degree = 3;
	desc.key_positions = points;
	desc.key_count = point_count;
	desc.key_times = times;
	desc.key_time_count = point_count;
	desc.closed = closed;

	s = qaws_curve_create_trajectory(&desc, out_curve);

	if (times != times_stack)
		qaws_internal_free(times);

	return s;
}
