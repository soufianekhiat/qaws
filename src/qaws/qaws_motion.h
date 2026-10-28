#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_MOTION_H_F32
#define QAWS_MOTION_H_F32
#define QAWS_MOTION_H_BODY
#endif
#else
#ifndef QAWS_MOTION_H_F64
#define QAWS_MOTION_H_F64
#define QAWS_MOTION_H_BODY
#endif
#endif
#ifdef QAWS_MOTION_H_BODY
#undef QAWS_MOTION_H_BODY
#define QAWS_MOTION_H

#include "qaws_types.h"
#include "qaws_status.h"
#include "qaws_curve.h"

/*
 * Motions over t in [0, 1]: the input of sweeps (a body moving, swept).
 *
 * Rigid motions move points of the body's frame: x(t) = R(t) p + T(t).
 *   - From two poses (each a position, a unit quaternion, its linear and
 *     angular velocities, and with QAWS_MOTION_QUINTIC its accelerations):
 *     T(t) is the Hermite polynomial through them (cubic from positions and
 *     velocities, quintic with accelerations too), so a cubic path given by
 *     its ends' values and derivatives is reproduced exactly. R(t) is
 *     R0 Cay(c(t)): c(t) a Hermite polynomial (of the same degree) in the
 *     Gibbs vector (axis tan(angle / 2)) of the rotation from pose 0, and
 *     Cay(c) = I + 2 ([c]x + [c]x^2) / (1 + c . c) the Cayley map, so the
 *     rotation is a rational function of t (its entries polynomials over
 *     1 + c . c): exact, no trigonometry. c(0) = 0, c(1) the relative
 *     rotation's, the ends' derivatives from the angular velocities (world
 *     frame) by the Gibbs vector's kinematics,
 *         omega = R0 (2 / (1 + c . c)) (c' + c x c'),
 *     and its derivative for the accelerations.
 *     Turning a half turn or more: a quaternion says where the body ends, not
 *     how many turns it took; the angular velocities say it. The total angle
 *     is the relative rotation's (about its axis) plus the whole turns that
 *     bring it nearest the angular velocities' mean along that axis (a spin
 *     of three turns ends where it started, its velocities 6 pi). Past a
 *     half turn the motion is cut into equal pieces of an eighth of a turn at
 *     most, each as above: the angle along the axis a Hermite polynomial
 *     (its rates the angular velocities' along the axis), their other parts
 *     blended linearly, the translation's polynomial cut at the joints
 *     (exact), so positions, orientations and their rates match where the
 *     pieces meet. A half turn the velocities do not orient (either way
 *     round as near) is refused.
 *   - Along a curve (3D): T(t) the curve at the parameter mapped from t
 *     (its range's start at 0, end at 1), R(t) a given constant rotation;
 *     or, rotation-minimizing (qaws_motion_create_curve_frames), the given
 *     rotation at 0 carried along the curve with no turn about its tangent
 *     (a pipe's cross-section, a tool following its path): R(t) = F(t)
 *     F(0)^T R0, F the frame of the unit tangent and a normal transported by
 *     double reflection (Wang et al. 2008) from 512 samples, one reflection
 *     on from the sample before t; the angular velocity c' x c'' / |c'|^2. A
 *     planar curve turns the body about the plane's normal only.
 *
 * Vertex motions move each of count vertices on its own Hermite polynomial
 * (positions and velocities at both ends, accelerations with
 * QAWS_MOTION_QUINTIC): a deforming mesh, each trajectory a Bezier curve in
 * time.
 *
 * Thread-safe on an immutable motion.
 */

typedef struct qaws_motion qaws_motion;

#define QAWS_MOTION_QUINTIC 1u      /* the accelerations matched too (quintic Hermite) */

typedef struct qaws_motion_pose
{
	qaws_vec3 position;
	qaws_scalar rotation[4];        /* a quaternion (w, x, y, z), normalized when made */
	qaws_vec3 velocity;
	qaws_vec3 angular_velocity;     /* in the world frame */
	qaws_vec3 acceleration;         /* with QAWS_MOTION_QUINTIC */
	qaws_vec3 angular_acceleration; /* with QAWS_MOTION_QUINTIC, in the world frame */
} qaws_motion_pose;

/* A rigid motion at t: x = rotation p + position (rotation by rows), its velocity and angular velocity (world frame) */
typedef struct qaws_motion_state
{
	qaws_vec3 position;
	qaws_scalar rotation[9];
	qaws_vec3 velocity;
	qaws_vec3 angular_velocity;
} qaws_motion_state;

/* INVALID_ARGUMENT for a value not finite, a zero quaternion, a half turn its angular velocities do not orient */
qaws_status qaws_motion_create_poses(qaws_motion** out_motion, qaws_motion_pose const* start, qaws_motion_pose const* end, unsigned int flags);

/* rotation: a quaternion (w, x, y, z), NULL the identity. The curve is borrowed: keep it until the motion is destroyed.
   INVALID_ARGUMENT for a curve not 3D */
qaws_status qaws_motion_create_curve(qaws_motion** out_motion, qaws_curve const* curve, qaws_scalar const* rotation);

#ifndef QAWS_ENUM_QAWS_MOTION_FRAME
#define QAWS_ENUM_QAWS_MOTION_FRAME
typedef enum qaws_motion_frame
{
	QAWS_MOTION_FRAME_FIXED = 0,             /* the rotation held (qaws_motion_create_curve) */
	QAWS_MOTION_FRAME_ROTATION_MINIMIZING    /* the rotation carried along the curve, no turn about its tangent */
} qaws_motion_frame;
#endif

/* As qaws_motion_create_curve, the rotation (NULL: the identity) the body's at t = 0, turned along the curve by the frame.
   INVALID_ARGUMENT as qaws_motion_create_curve, or an unknown frame */
qaws_status qaws_motion_create_curve_frames(qaws_motion** out_motion, qaws_curve const* curve, qaws_scalar const* rotation,
	qaws_motion_frame frame);

/* count vertices (3 each): positions and velocities at t = 0 and 1; accelerations (may be NULL) with QAWS_MOTION_QUINTIC */
qaws_status qaws_motion_create_vertices(qaws_motion** out_motion, qaws_scalar const* positions0, qaws_scalar const* positions1,
	qaws_scalar const* velocities0, qaws_scalar const* velocities1, qaws_scalar const* accelerations0, qaws_scalar const* accelerations1,
	unsigned int count, unsigned int flags);

void qaws_motion_destroy(qaws_motion* motion);

/* 1 for a vertex motion, 0 for a rigid one */
int qaws_motion_is_per_vertex(qaws_motion const* motion);
unsigned int qaws_motion_get_vertex_count(qaws_motion const* motion);

/* a rigid motion at t (INVALID_ARGUMENT for a vertex motion) */
qaws_status qaws_motion_evaluate(qaws_motion_state* out_state, qaws_motion const* motion, qaws_scalar t);

/* count points (3 each) moved to t: a rigid motion's R(t) p + T(t); a vertex motion's own vertices (points NULL, count its
   vertex count). out_velocities (may be NULL): their velocities */
qaws_status qaws_motion_transform_points(qaws_scalar* out_points, qaws_scalar* out_velocities, qaws_motion const* motion, qaws_scalar t,
	qaws_scalar const* points, unsigned int count);

#endif /* QAWS_MOTION_H */
