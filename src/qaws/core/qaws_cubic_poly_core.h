#ifndef QAWS_CUBIC_POLY_CORE_H
#define QAWS_CUBIC_POLY_CORE_H

#include "../qaws_platform.h"
#include "../qaws_core_types.h"

/* ===================================================================
 * Cubic polynomial: P(t) = a*t^3 + b*t^2 + c*t + d
 *
 * Shared eval code for Hermite, Catmull-Rom, and Trajectory.
 * All three pre-compute a,b,c,d per span per dimension at creation,
 * then this function is the shared eval hot path.
 * =================================================================== */

QAWS_INLINE qaws_eval_2d qaws_cubic_eval_2d(
    qaws_vec2 a, qaws_vec2 b,
    qaws_vec2 c, qaws_vec2 d,
    qaws_scalar t, int eval_flags)
{
    qaws_scalar t2 = t * t;
    qaws_scalar t3 = t2 * t;
    qaws_eval_2d result;

    result.position.x = a.x * t3 + b.x * t2 + c.x * t + d.x;
    result.position.y = a.y * t3 + b.y * t2 + c.y * t + d.y;

    if ((eval_flags & 0x2) != 0) {
        result.d1.x = QAWS_LITERAL(3.0) * a.x * t2
                    + QAWS_LITERAL(2.0) * b.x * t + c.x;
        result.d1.y = QAWS_LITERAL(3.0) * a.y * t2
                    + QAWS_LITERAL(2.0) * b.y * t + c.y;
    } else {
        result.d1.x = QAWS_ZERO;
        result.d1.y = QAWS_ZERO;
    }

    if ((eval_flags & 0x4) != 0) {
        result.d2.x = QAWS_LITERAL(6.0) * a.x * t + QAWS_LITERAL(2.0) * b.x;
        result.d2.y = QAWS_LITERAL(6.0) * a.y * t + QAWS_LITERAL(2.0) * b.y;
    } else {
        result.d2.x = QAWS_ZERO;
        result.d2.y = QAWS_ZERO;
    }

    return result;
}

/* ===================================================================
 * Cubic polynomial evaluation (3D)
 * =================================================================== */

QAWS_INLINE qaws_eval_3d qaws_cubic_eval_3d(
    qaws_vec3 a, qaws_vec3 b,
    qaws_vec3 c, qaws_vec3 d,
    qaws_scalar t, int eval_flags)
{
    qaws_scalar t2 = t * t;
    qaws_scalar t3 = t2 * t;
    qaws_eval_3d result;

    result.position.x = a.x * t3 + b.x * t2 + c.x * t + d.x;
    result.position.y = a.y * t3 + b.y * t2 + c.y * t + d.y;
    result.position.z = a.z * t3 + b.z * t2 + c.z * t + d.z;

    if ((eval_flags & 0x2) != 0) {
        result.d1.x = QAWS_LITERAL(3.0) * a.x * t2
                    + QAWS_LITERAL(2.0) * b.x * t + c.x;
        result.d1.y = QAWS_LITERAL(3.0) * a.y * t2
                    + QAWS_LITERAL(2.0) * b.y * t + c.y;
        result.d1.z = QAWS_LITERAL(3.0) * a.z * t2
                    + QAWS_LITERAL(2.0) * b.z * t + c.z;
    } else {
        result.d1.x = QAWS_ZERO;
        result.d1.y = QAWS_ZERO;
        result.d1.z = QAWS_ZERO;
    }

    if ((eval_flags & 0x4) != 0) {
        result.d2.x = QAWS_LITERAL(6.0) * a.x * t + QAWS_LITERAL(2.0) * b.x;
        result.d2.y = QAWS_LITERAL(6.0) * a.y * t + QAWS_LITERAL(2.0) * b.y;
        result.d2.z = QAWS_LITERAL(6.0) * a.z * t + QAWS_LITERAL(2.0) * b.z;
    } else {
        result.d2.x = QAWS_ZERO;
        result.d2.y = QAWS_ZERO;
        result.d2.z = QAWS_ZERO;
    }

    return result;
}

#endif /* QAWS_CUBIC_POLY_CORE_H */
