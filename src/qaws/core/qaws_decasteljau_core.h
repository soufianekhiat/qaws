#ifndef QAWS_DECASTELJAU_CORE_H
#define QAWS_DECASTELJAU_CORE_H

#include "../qaws_platform.h"
#include "../qaws_core_types.h"

/* ===================================================================
 * Cubic Bezier: position only (2D)
 *
 * Four CPs passed by value — no arrays needed.
 * =================================================================== */

QAWS_INLINE qaws_vec2 qaws_bezier3_eval_2d(
    qaws_vec2 p0, qaws_vec2 p1,
    qaws_vec2 p2, qaws_vec2 p3,
    qaws_scalar t)
{
    qaws_scalar u = QAWS_ONE - t;
    qaws_vec2 q0, q1, q2, r0, r1, result;

    q0.x = u * p0.x + t * p1.x;  q0.y = u * p0.y + t * p1.y;
    q1.x = u * p1.x + t * p2.x;  q1.y = u * p1.y + t * p2.y;
    q2.x = u * p2.x + t * p3.x;  q2.y = u * p2.y + t * p3.y;

    r0.x = u * q0.x + t * q1.x;  r0.y = u * q0.y + t * q1.y;
    r1.x = u * q1.x + t * q2.x;  r1.y = u * q1.y + t * q2.y;

    result.x = u * r0.x + t * r1.x;
    result.y = u * r0.y + t * r1.y;
    return result;
}

/* ===================================================================
 * Cubic Bezier: position only (3D)
 * =================================================================== */

QAWS_INLINE qaws_vec3 qaws_bezier3_eval_3d(
    qaws_vec3 p0, qaws_vec3 p1,
    qaws_vec3 p2, qaws_vec3 p3,
    qaws_scalar t)
{
    qaws_scalar u = QAWS_ONE - t;
    qaws_vec3 q0, q1, q2, r0, r1, result;

    q0.x = u * p0.x + t * p1.x;  q0.y = u * p0.y + t * p1.y;  q0.z = u * p0.z + t * p1.z;
    q1.x = u * p1.x + t * p2.x;  q1.y = u * p1.y + t * p2.y;  q1.z = u * p1.z + t * p2.z;
    q2.x = u * p2.x + t * p3.x;  q2.y = u * p2.y + t * p3.y;  q2.z = u * p2.z + t * p3.z;

    r0.x = u * q0.x + t * q1.x;  r0.y = u * q0.y + t * q1.y;  r0.z = u * q0.z + t * q1.z;
    r1.x = u * q1.x + t * q2.x;  r1.y = u * q1.y + t * q2.y;  r1.z = u * q1.z + t * q2.z;

    result.x = u * r0.x + t * r1.x;
    result.y = u * r0.y + t * r1.y;
    result.z = u * r0.z + t * r1.z;
    return result;
}

/* ===================================================================
 * Arbitrary-degree De Casteljau (2D, position only)
 *
 * Control points in flat array: [x0,y0, x1,y1, ..., xn,yn]
 * degree must be <= QAWS_CORE_MAX_DEGREE.
 * =================================================================== */

QAWS_INLINE qaws_vec2 qaws_decasteljau_2d(
    qaws_scalar cp[QAWS_CORE_MAX_POINTS * 2],
    int degree,
    qaws_scalar t)
{
    qaws_scalar work[QAWS_CORE_MAX_POINTS * 2];
    qaws_scalar u = QAWS_ONE - t;
    qaws_vec2 result;
    int r, i;

    work[0] = work[1] = QAWS_ZERO;   /* (degree < 0: no garbage read) */
    for (i = 0; i <= degree; i++) {
        work[i * 2]     = cp[i * 2];
        work[i * 2 + 1] = cp[i * 2 + 1];
    }

    for (r = 1; r <= degree; r++) {
        for (i = 0; i <= degree - r; i++) {
            work[i * 2]     = u * work[i * 2]     + t * work[(i + 1) * 2];
            work[i * 2 + 1] = u * work[i * 2 + 1] + t * work[(i + 1) * 2 + 1];
        }
    }

    result.x = work[0];
    result.y = work[1];
    return result;
}

/* ===================================================================
 * Arbitrary-degree De Casteljau (3D, position only)
 * =================================================================== */

QAWS_INLINE qaws_vec3 qaws_decasteljau_3d(
    qaws_scalar cp[QAWS_CORE_MAX_POINTS * 3],
    int degree,
    qaws_scalar t)
{
    qaws_scalar work[QAWS_CORE_MAX_POINTS * 3];
    qaws_scalar u = QAWS_ONE - t;
    qaws_vec3 result;
    int r, i;

    for (i = 0; i <= degree; i++) {
        work[i * 3]     = cp[i * 3];
        work[i * 3 + 1] = cp[i * 3 + 1];
        work[i * 3 + 2] = cp[i * 3 + 2];
    }

    for (r = 1; r <= degree; r++) {
        for (i = 0; i <= degree - r; i++) {
            work[i * 3]     = u * work[i * 3]     + t * work[(i + 1) * 3];
            work[i * 3 + 1] = u * work[i * 3 + 1] + t * work[(i + 1) * 3 + 1];
            work[i * 3 + 2] = u * work[i * 3 + 2] + t * work[(i + 1) * 3 + 2];
        }
    }

    result.x = work[0];
    result.y = work[1];
    result.z = work[2];
    return result;
}

/* ===================================================================
 * Bezier derivative control points (2D)
 *
 * Given degree-N CPs, writes degree-(N-1) differenced CPs.
 * =================================================================== */

QAWS_INLINE void qaws_bezier_deriv_points_2d(
    QAWS_OUT qaws_scalar out_dp[QAWS_CORE_MAX_POINTS * 2],
    qaws_scalar cp[QAWS_CORE_MAX_POINTS * 2],
    int degree)
{
    qaws_scalar deg = QAWS_ITOF(degree);
    int i;
    for (i = 0; i < degree; i++) {
        out_dp[i * 2]     = deg * (cp[(i + 1) * 2]     - cp[i * 2]);
        out_dp[i * 2 + 1] = deg * (cp[(i + 1) * 2 + 1] - cp[i * 2 + 1]);
    }
}

/* ===================================================================
 * Bezier derivative control points (3D)
 * =================================================================== */

QAWS_INLINE void qaws_bezier_deriv_points_3d(
    QAWS_OUT qaws_scalar out_dp[QAWS_CORE_MAX_POINTS * 3],
    qaws_scalar cp[QAWS_CORE_MAX_POINTS * 3],
    int degree)
{
    qaws_scalar deg = QAWS_ITOF(degree);
    int i;
    for (i = 0; i < degree; i++) {
        out_dp[i * 3]     = deg * (cp[(i + 1) * 3]     - cp[i * 3]);
        out_dp[i * 3 + 1] = deg * (cp[(i + 1) * 3 + 1] - cp[i * 3 + 1]);
        out_dp[i * 3 + 2] = deg * (cp[(i + 1) * 3 + 2] - cp[i * 3 + 2]);
    }
}

/* ===================================================================
 * Full Bezier evaluation with derivatives (2D)
 *
 * Returns position + D1 + D2 via qaws_eval_2d.
 * eval_flags uses QAWS_EVAL_FLAG_* bits (or just pass non-zero for all).
 * =================================================================== */

QAWS_INLINE qaws_eval_2d qaws_decasteljau_eval_2d(
    qaws_scalar cp[QAWS_CORE_MAX_POINTS * 2],
    int degree,
    qaws_scalar t,
    int eval_flags)
{
    qaws_eval_2d result;
    qaws_vec2 pos;
    qaws_scalar d1cp[QAWS_CORE_MAX_POINTS * 2];
    qaws_scalar d2cp[QAWS_CORE_MAX_POINTS * 2];
    qaws_vec2 d1val, d2val;

    pos = qaws_decasteljau_2d(cp, degree, t);
    result.position = pos;
    result.d1.x = QAWS_ZERO; result.d1.y = QAWS_ZERO;
    result.d2.x = QAWS_ZERO; result.d2.y = QAWS_ZERO;

    if (((eval_flags & 0x2) != 0) && degree >= 1) {
        qaws_bezier_deriv_points_2d(d1cp, cp, degree);
        d1val = qaws_decasteljau_2d(d1cp, degree - 1, t);
        result.d1 = d1val;
    }

    if (((eval_flags & 0x4) != 0) && degree >= 2) {
        if ((eval_flags & 0x2) == 0) {
            qaws_bezier_deriv_points_2d(d1cp, cp, degree);
        }
        qaws_bezier_deriv_points_2d(d2cp, d1cp, degree - 1);
        d2val = qaws_decasteljau_2d(d2cp, degree - 2, t);
        result.d2 = d2val;
    }

    return result;
}

/* ===================================================================
 * Full Bezier evaluation with derivatives (3D)
 * =================================================================== */

QAWS_INLINE qaws_eval_3d qaws_decasteljau_eval_3d(
    qaws_scalar cp[QAWS_CORE_MAX_POINTS * 3],
    int degree,
    qaws_scalar t,
    int eval_flags)
{
    qaws_eval_3d result;
    qaws_vec3 pos;
    qaws_scalar d1cp[QAWS_CORE_MAX_POINTS * 3];
    qaws_scalar d2cp[QAWS_CORE_MAX_POINTS * 3];
    qaws_vec3 d1val, d2val;

    pos = qaws_decasteljau_3d(cp, degree, t);
    result.position = pos;
    result.d1.x = QAWS_ZERO; result.d1.y = QAWS_ZERO; result.d1.z = QAWS_ZERO;
    result.d2.x = QAWS_ZERO; result.d2.y = QAWS_ZERO; result.d2.z = QAWS_ZERO;

    if (((eval_flags & 0x2) != 0) && degree >= 1) {
        qaws_bezier_deriv_points_3d(d1cp, cp, degree);
        d1val = qaws_decasteljau_3d(d1cp, degree - 1, t);
        result.d1 = d1val;
    }

    if (((eval_flags & 0x4) != 0) && degree >= 2) {
        if ((eval_flags & 0x2) == 0) {
            qaws_bezier_deriv_points_3d(d1cp, cp, degree);
        }
        qaws_bezier_deriv_points_3d(d2cp, d1cp, degree - 1);
        d2val = qaws_decasteljau_3d(d2cp, degree - 2, t);
        result.d2 = d2val;
    }

    return result;
}

#endif /* QAWS_DECASTELJAU_CORE_H */
