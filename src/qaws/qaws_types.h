#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_TYPES_H_F32
#define QAWS_TYPES_H_F32
#define QAWS_TYPES_H_BODY
#endif
#else
#ifndef QAWS_TYPES_H_F64
#define QAWS_TYPES_H_F64
#define QAWS_TYPES_H_BODY
#endif
#endif
#ifdef QAWS_TYPES_H_BODY
#undef QAWS_TYPES_H_BODY
#define QAWS_TYPES_H

/*
 * Thread safety guarantees:
 *
 * - qaws_curve objects are immutable after creation. Multiple threads
 *   may concurrently evaluate, sample, inspect, or traverse the same
 *   curve without synchronization.
 *
 * - qaws_traversal objects are safe for concurrent reads (evaluate,
 *   map queries) but qaws_traversal_advance_* and qaws_traversal_reset
 *   mutate internal state and require external synchronization.
 *
 * - Creation and destruction functions (create_*, destroy) are not
 *   thread-safe on the same object. Do not destroy a curve while
 *   another thread is using it.
 */

#include "qaws_platform.h"   /* (qaws_scalar: the current precision's) */

typedef struct qaws_vec2
{
	qaws_scalar x;
	qaws_scalar y;
} qaws_vec2;

typedef struct qaws_vec3
{
	qaws_scalar x;
	qaws_scalar y;
	qaws_scalar z;
} qaws_vec3;

#ifndef QAWS_ENUM_QAWS_DIMENSION
#define QAWS_ENUM_QAWS_DIMENSION
typedef enum qaws_dimension
{
	QAWS_DIMENSION_2D = 2,
	QAWS_DIMENSION_3D = 3
} qaws_dimension;
#endif

#ifndef QAWS_ENUM_QAWS_CURVE_KIND
#define QAWS_ENUM_QAWS_CURVE_KIND
typedef enum qaws_curve_kind
{
	QAWS_CURVE_KIND_INVALID = 0,
	QAWS_CURVE_KIND_BEZIER,
	QAWS_CURVE_KIND_HERMITE,
	QAWS_CURVE_KIND_CATMULL_ROM,
	QAWS_CURVE_KIND_BSPLINE,
	QAWS_CURVE_KIND_NURBS,
	QAWS_CURVE_KIND_TRAJECTORY,
	QAWS_CURVE_KIND_YUKSEL,
	QAWS_CURVE_KIND_RATIONAL_BEZIER,
	QAWS_CURVE_KIND_COMPOSITE,
	QAWS_CURVE_KIND_ARC,
	QAWS_CURVE_KIND_POLYNOMIAL,
	QAWS_CURVE_KIND_CLOTHOID,
	QAWS_CURVE_KIND_SUBDIVISION,
	QAWS_CURVE_KIND_REPARAMETERIZED
} qaws_curve_kind;
#endif

#ifndef QAWS_ENUM_QAWS_PARAMETERIZATION
#define QAWS_ENUM_QAWS_PARAMETERIZATION
typedef enum qaws_parameterization
{
	QAWS_PARAMETERIZATION_DEFAULT = 0,
	QAWS_PARAMETERIZATION_UNIFORM,
	QAWS_PARAMETERIZATION_CHORDAL,
	QAWS_PARAMETERIZATION_CENTRIPETAL
} qaws_parameterization;
#endif

#ifndef QAWS_ENUM_QAWS_CONTINUITY
#define QAWS_ENUM_QAWS_CONTINUITY
typedef enum qaws_continuity
{
	QAWS_CONTINUITY_NONE = 0,
	QAWS_CONTINUITY_C0,
	QAWS_CONTINUITY_C1,
	QAWS_CONTINUITY_C2,
	QAWS_CONTINUITY_C3
} qaws_continuity;
#endif

#ifndef QAWS_ENUM_QAWS_TRAVERSAL_MODE
#define QAWS_ENUM_QAWS_TRAVERSAL_MODE
typedef enum qaws_traversal_mode
{
	QAWS_TRAVERSAL_MODE_PARAMETER = 0,
	QAWS_TRAVERSAL_MODE_ARC_LENGTH,
	QAWS_TRAVERSAL_MODE_TIME
} qaws_traversal_mode;
#endif

#ifndef QAWS_ENUM_QAWS_MOTION_PROFILE
#define QAWS_ENUM_QAWS_MOTION_PROFILE
typedef enum qaws_motion_profile
{
	QAWS_MOTION_PROFILE_NONE = 0,
	QAWS_MOTION_PROFILE_CONSTANT_SPEED,
	QAWS_MOTION_PROFILE_CONSTANT_ACCELERATION,
	QAWS_MOTION_PROFILE_TRAPEZOIDAL_SPEED,
	QAWS_MOTION_PROFILE_S_CURVE,
	QAWS_MOTION_PROFILE_CUSTOM
} qaws_motion_profile;
#endif

/* Custom speed function: returns speed given distance along the curve. */
typedef qaws_scalar (*qaws_speed_fn)(qaws_scalar distance, void* user_data);

#ifndef QAWS_ENUM_QAWS_EASING
#define QAWS_ENUM_QAWS_EASING
typedef enum qaws_easing
{
	QAWS_EASING_LINEAR = 0,
	QAWS_EASING_QUAD_IN,
	QAWS_EASING_QUAD_OUT,
	QAWS_EASING_QUAD_IN_OUT,
	QAWS_EASING_CUBIC_IN,
	QAWS_EASING_CUBIC_OUT,
	QAWS_EASING_CUBIC_IN_OUT,
	QAWS_EASING_SINE_IN,
	QAWS_EASING_SINE_OUT,
	QAWS_EASING_SINE_IN_OUT
} qaws_easing;
#endif

#ifndef QAWS_ENUM_QAWS_WRAP_MODE
#define QAWS_ENUM_QAWS_WRAP_MODE
typedef enum qaws_wrap_mode
{
	QAWS_WRAP_CLAMP = 0,
	QAWS_WRAP_LOOP,
	QAWS_WRAP_PING_PONG
} qaws_wrap_mode;
#endif

#ifndef QAWS_ENUM_QAWS_EVAL_FLAGS
#define QAWS_ENUM_QAWS_EVAL_FLAGS
typedef enum qaws_eval_flags
{
	QAWS_EVAL_FLAG_NONE     = 0,
	QAWS_EVAL_FLAG_POSITION = 1 << 0,
	QAWS_EVAL_FLAG_D1       = 1 << 1,
	QAWS_EVAL_FLAG_D2       = 1 << 2,
	QAWS_EVAL_FLAG_D3       = 1 << 3
} qaws_eval_flags;
#endif

typedef struct qaws_range
{
	qaws_scalar min_value;
	qaws_scalar max_value;
} qaws_range;

typedef struct qaws_eval_result_2d
{
	qaws_vec2 position;
	qaws_vec2 d1;
	qaws_vec2 d2;
	qaws_vec2 d3;
	unsigned int valid_flags;
} qaws_eval_result_2d;

typedef struct qaws_eval_result_3d
{
	qaws_vec3 position;
	qaws_vec3 d1;
	qaws_vec3 d2;
	qaws_vec3 d3;
	unsigned int valid_flags;
} qaws_eval_result_3d;

typedef struct qaws_sampling_desc
{
	qaws_traversal_mode traversal_mode;
	unsigned int sample_count;
	qaws_scalar error_tolerance;
	int use_adaptive_sampling;
} qaws_sampling_desc;

typedef struct qaws_traversal_desc
{
	qaws_traversal_mode traversal_mode;
	qaws_motion_profile motion_profile;
	qaws_scalar speed;
	qaws_scalar acceleration;
	qaws_scalar max_speed;
	qaws_scalar start_time;
	qaws_scalar end_time;
	int clamp_to_domain;
	qaws_easing easing;
	qaws_wrap_mode wrap_mode;
	qaws_scalar jerk;                      /* max jerk for S-curve profile */
	qaws_speed_fn custom_speed_fn;         /* custom speed callback */
	void* custom_speed_fn_user_data;       /* user data for custom speed callback */
} qaws_traversal_desc;

typedef struct qaws_curve qaws_curve;
typedef struct qaws_traversal qaws_traversal;

/*
 * Custom allocator. Pass to _ex creation functions to override
 * malloc/free. If alloc returns NULL, creation fails with
 * QAWS_STATUS_ALLOCATION_FAILURE.
 *
 * alloc(size, user_data) — allocate size bytes (like malloc)
 * dealloc(ptr, user_data) — free pointer (like free)
 */
typedef struct qaws_allocator
{
	void* (*alloc)(unsigned long size, void* user_data);
	void  (*dealloc)(void* ptr, void* user_data);
	void* user_data;
} qaws_allocator;

/*
 * Stack-allocated curve for single-span curves (up to degree 7, 3D).
 * No heap allocation. Use qaws_curve_init_*_inline() creation functions.
 *
 * The curve pointer is retrieved via qaws_curve_inline_get_curve().
 * Do NOT call qaws_curve_destroy() on inline curves.
 *
 * Maximum storage: 8 control points * 3 components = 24 scalars for CPs,
 * plus impl struct overhead. Total fixed budget: 512 bytes.
 */
#define QAWS_INLINE_STORAGE_BYTES 512

typedef struct qaws_curve_inline
{
	/* Opaque storage — do not access directly.
	   Layout is set by the init functions. */
	char storage[QAWS_INLINE_STORAGE_BYTES];
} qaws_curve_inline;

/* Get the usable qaws_curve pointer from an initialized inline curve. */
qaws_curve* qaws_curve_inline_get_curve(qaws_curve_inline* ic);

#endif /* QAWS_TYPES_H */
