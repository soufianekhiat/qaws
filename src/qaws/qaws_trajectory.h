#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_TRAJECTORY_H_F32
#define QAWS_TRAJECTORY_H_F32
#define QAWS_TRAJECTORY_H_BODY
#endif
#else
#ifndef QAWS_TRAJECTORY_H_F64
#define QAWS_TRAJECTORY_H_F64
#define QAWS_TRAJECTORY_H_BODY
#endif
#endif
#ifdef QAWS_TRAJECTORY_H_BODY
#undef QAWS_TRAJECTORY_H_BODY
#define QAWS_TRAJECTORY_H

#include "qaws_types.h"
#include "qaws_status.h"

typedef struct qaws_trajectory_desc
{
	qaws_dimension dimension;
	unsigned int degree;
	void const* key_positions;
	unsigned int key_count;
	qaws_scalar const* key_times;
	unsigned int key_time_count;
	void const* key_velocities;
	unsigned int key_velocity_count;
	void const* key_accelerations;
	unsigned int key_acceleration_count;
	int closed;
} qaws_trajectory_desc;

qaws_status qaws_curve_create_trajectory(
	qaws_curve** out_curve,
	qaws_trajectory_desc const* desc);

#endif /* QAWS_TRAJECTORY_H */
