#ifndef QAWS_ARC_CORE_H
#define QAWS_ARC_CORE_H

#include "../qaws_platform.h"
#include "../qaws_core_types.h"

/* ===================================================================
 * Arc evaluation (2D)
 *
 * center:      arc center
 * axis_u:      cos direction vector (scaled by radius_x)
 * axis_v:      sin direction vector (scaled by radius_y)
 * theta_start: starting angle
 * theta_range: angular sweep
 * t:           parameter in [0,1]
 * =================================================================== */

QAWS_INLINE qaws_vec2 qaws_arc_eval_2d(
    qaws_vec2 center,
    qaws_vec2 axis_u,
    qaws_vec2 axis_v,
    qaws_scalar theta_start,
    qaws_scalar theta_range,
    qaws_scalar t)
{
    qaws_scalar theta = theta_start + t * theta_range;
    qaws_scalar ct = QAWS_COS(theta);
    qaws_scalar st = QAWS_SIN(theta);
    qaws_vec2 result;

    result.x = center.x + ct * axis_u.x + st * axis_v.x;
    result.y = center.y + ct * axis_u.y + st * axis_v.y;
    return result;
}

/* ===================================================================
 * Arc evaluation (3D)
 * =================================================================== */

QAWS_INLINE qaws_vec3 qaws_arc_eval_3d(
    qaws_vec3 center,
    qaws_vec3 axis_u,
    qaws_vec3 axis_v,
    qaws_scalar theta_start,
    qaws_scalar theta_range,
    qaws_scalar t)
{
    qaws_scalar theta = theta_start + t * theta_range;
    qaws_scalar ct = QAWS_COS(theta);
    qaws_scalar st = QAWS_SIN(theta);
    qaws_vec3 result;

    result.x = center.x + ct * axis_u.x + st * axis_v.x;
    result.y = center.y + ct * axis_u.y + st * axis_v.y;
    result.z = center.z + ct * axis_u.z + st * axis_v.z;
    return result;
}

/* ===================================================================
 * Full arc evaluation with derivatives (2D)
 *
 * D1 = (-sin*u + cos*v) * theta_range
 * D2 = (-cos*u - sin*v) * theta_range^2
 * =================================================================== */

QAWS_INLINE qaws_eval_2d qaws_arc_eval_full_2d(
    qaws_vec2 center,
    qaws_vec2 axis_u,
    qaws_vec2 axis_v,
    qaws_scalar theta_start,
    qaws_scalar theta_range,
    qaws_scalar t,
    int eval_flags)
{
    qaws_scalar theta = theta_start + t * theta_range;
    qaws_scalar ct = QAWS_COS(theta);
    qaws_scalar st = QAWS_SIN(theta);
    qaws_eval_2d result;

    result.position.x = center.x + ct * axis_u.x + st * axis_v.x;
    result.position.y = center.y + ct * axis_u.y + st * axis_v.y;

    if ((eval_flags & 0x2) != 0) {
        result.d1.x = (-st * axis_u.x + ct * axis_v.x) * theta_range;
        result.d1.y = (-st * axis_u.y + ct * axis_v.y) * theta_range;
    } else {
        result.d1.x = QAWS_ZERO;
        result.d1.y = QAWS_ZERO;
    }

    if ((eval_flags & 0x4) != 0) {
        qaws_scalar tr2 = theta_range * theta_range;
        result.d2.x = (-ct * axis_u.x - st * axis_v.x) * tr2;
        result.d2.y = (-ct * axis_u.y - st * axis_v.y) * tr2;
    } else {
        result.d2.x = QAWS_ZERO;
        result.d2.y = QAWS_ZERO;
    }

    return result;
}

/* ===================================================================
 * Full arc evaluation with derivatives (3D)
 * =================================================================== */

QAWS_INLINE qaws_eval_3d qaws_arc_eval_full_3d(
    qaws_vec3 center,
    qaws_vec3 axis_u,
    qaws_vec3 axis_v,
    qaws_scalar theta_start,
    qaws_scalar theta_range,
    qaws_scalar t,
    int eval_flags)
{
    qaws_scalar theta = theta_start + t * theta_range;
    qaws_scalar ct = QAWS_COS(theta);
    qaws_scalar st = QAWS_SIN(theta);
    qaws_eval_3d result;

    result.position.x = center.x + ct * axis_u.x + st * axis_v.x;
    result.position.y = center.y + ct * axis_u.y + st * axis_v.y;
    result.position.z = center.z + ct * axis_u.z + st * axis_v.z;

    if ((eval_flags & 0x2) != 0) {
        result.d1.x = (-st * axis_u.x + ct * axis_v.x) * theta_range;
        result.d1.y = (-st * axis_u.y + ct * axis_v.y) * theta_range;
        result.d1.z = (-st * axis_u.z + ct * axis_v.z) * theta_range;
    } else {
        result.d1.x = QAWS_ZERO;
        result.d1.y = QAWS_ZERO;
        result.d1.z = QAWS_ZERO;
    }

    if ((eval_flags & 0x4) != 0) {
        qaws_scalar tr2 = theta_range * theta_range;
        result.d2.x = (-ct * axis_u.x - st * axis_v.x) * tr2;
        result.d2.y = (-ct * axis_u.y - st * axis_v.y) * tr2;
        result.d2.z = (-ct * axis_u.z - st * axis_v.z) * tr2;
    } else {
        result.d2.x = QAWS_ZERO;
        result.d2.y = QAWS_ZERO;
        result.d2.z = QAWS_ZERO;
    }

    return result;
}

#endif /* QAWS_ARC_CORE_H */
